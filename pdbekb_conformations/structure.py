"""
structure.py -- the (optional) superposition-dependent layer.

PDBe-KB deliberately does NOT ship an RMSD matrix: the GLOCON score is
superposition-agnostic (it compares internal Ca-Ca distance patterns) and is
therefore not convertible to RMSD. If you want Ca-RMSD between the
representative structures of two clusters you must superpose them yourself.
Per the PDBe team, GESAMT (their internal tool) is currently unreliable, so we
superpose locally here.

Two things this module gives you for a pair of cluster representatives:
  * ca_rmsd  -- Kabsch Ca-RMSD over the residues both structures share, mapped
                onto a common UniProt numbering via SIFTS (a fair comparison).
  * per-residue Ca deviation -- localises WHERE the conformational change is,
                which is exactly what you want to steer binder sampling.

A thin wrapper for US-align / TM-align is included as well; the PDBe team
suggested TM-score or lDDT as more robust than RMSD for cross-structure
comparison, and you already use US-align elsewhere.

NOTE: this layer needs to download mmCIF files from PDBe/RCSB, so it runs on a
machine with access to those servers (not required for the clustering tables).
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

try:
    import gemmi
    _HAVE_GEMMI = True
except Exception:                       # pragma: no cover
    _HAVE_GEMMI = False


# --------------------------------------------------------------------------- #
# SIFTS: author residue number -> UniProt residue number, per chain
# --------------------------------------------------------------------------- #
def auth_to_uniprot(sifts: dict, pdb: str, acc: str, chain: str) -> dict[int, int]:
    """Build {label_seq_id: uniprot_number} for one pdb chain.

    Uses the PDBe `mappings/uniprot_segments` payload. Within a SIFTS segment
    the mapping is 1:1 and linear, so we walk each segment via its
    `residue_number` (mmCIF label_seq_id) endpoints -- NOT
    `author_residue_number`: confirmed on live data that the author number is
    null for a substantial fraction of mappings (disordered/construct
    termini with no PDB author numbering), while `residue_number` is always
    populated. Keying on label_seq_id instead, matched against gemmi's
    `Residue.label_seq` (see `ca_from_chain`), recovers those chains instead
    of silently dropping them."""
    out: dict[int, int] = {}
    try:
        mappings = sifts[pdb.lower()]["UniProt"][acc]["mappings"]
    except (KeyError, TypeError):
        return out
    for mp in mappings:
        if str(mp.get("chain_id")) != str(chain):
            continue
        unp0 = int(mp["unp_start"])
        lbl0_raw = mp["start"]["residue_number"]
        lbl1_raw = mp["end"]["residue_number"]
        if lbl0_raw is None or lbl1_raw is None:
            log.debug("%s chain %s: SIFTS mapping has null residue_number, "
                     "skipping (unp %s-%s)", pdb, chain, mp.get("unp_start"), mp.get("unp_end"))
            continue
        lbl0, lbl1 = int(lbl0_raw), int(lbl1_raw)
        for lbl in range(lbl0, lbl1 + 1):
            out[lbl] = unp0 + (lbl - lbl0)
    return out


# --------------------------------------------------------------------------- #
# load a chain + extract Ca coordinates keyed by UniProt residue number
# --------------------------------------------------------------------------- #
def load_chain(mmcif_bytes: bytes, chain: str):
    """Return (gemmi.Structure, gemmi.Chain) for `chain`; (st, None) if absent.

    Keep the returned Structure alive for as long as you use the Chain."""
    if not _HAVE_GEMMI:
        raise RuntimeError("gemmi is required for coordinate extraction")
    with tempfile.NamedTemporaryFile(suffix=".cif") as tmp:
        tmp.write(mmcif_bytes)
        tmp.flush()
        st = gemmi.read_structure(tmp.name)
    st.setup_entities()
    model = st[0]
    for ch in model:
        if ch.name == chain:
            return st, ch
    return st, None


def ca_from_chain(chain, a2u: dict[int, int], lo: int, hi: int) -> dict[int, np.ndarray]:
    """{uniprot_resnum: xyz} for CA atoms of an already-loaded chain in [lo, hi].

    `a2u` is keyed by label_seq_id (see `auth_to_uniprot`), matched here via
    gemmi's `Residue.label_seq` -- not `seqid.num` (author numbering), which
    SIFTS does not always report."""
    coords: dict[int, np.ndarray] = {}
    if chain is None:
        return coords
    for res in chain:
        unp = a2u.get(res.label_seq)
        if unp is None or not (lo <= unp <= hi):
            continue
        ca = res.find_atom("CA", "*")
        if ca is not None:
            coords[unp] = np.array([ca.pos.x, ca.pos.y, ca.pos.z])
    return coords


def ca_by_uniprot(mmcif_bytes: bytes, chain: str, a2u: dict[int, int],
                  lo: int, hi: int) -> dict[int, np.ndarray]:
    """Convenience wrapper: load a chain and return its CA coords."""
    _, ch = load_chain(mmcif_bytes, chain)
    return ca_from_chain(ch, a2u, lo, hi)


# --------------------------------------------------------------------------- #
# Kabsch superposition (returns a re-usable affine transform)
# --------------------------------------------------------------------------- #
IDENTITY = (np.eye(3), np.zeros(3))


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Optimal rigid transform mapping P onto Q.

    Returns (R, t, rmsd) so that  X_aligned = (R @ X.T).T + t  for any point
    set X in P's frame. Returning the affine (R, t) -- not just the fitted
    points -- is what lets us re-apply the fit to whole structures when saving
    aligned coordinates."""
    Pmean, Qmean = P.mean(0), Q.mean(0)
    H = (P - Pmean).T @ (Q - Qmean)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Qmean - R @ Pmean
    P_aligned = (R @ P.T).T + t
    rmsd = float(np.sqrt(((P_aligned - Q) ** 2).sum(1).mean()))
    return R, t, rmsd


def compare_reps(ca_a: dict[int, np.ndarray],
                 ca_b: dict[int, np.ndarray]) -> Optional[dict]:
    """Ca-RMSD + per-residue deviation between two representatives.

    Superposes on the common UniProt residues, then reports overall RMSD and a
    per-residue deviation profile (deviation after global superposition)."""
    common = sorted(set(ca_a) & set(ca_b))
    if len(common) < 3:
        return None
    A = np.array([ca_a[i] for i in common])
    B = np.array([ca_b[i] for i in common])
    R, t, rmsd = kabsch(A, B)
    A_aln = (R @ A.T).T + t
    dev = np.sqrt(((A_aln - B) ** 2).sum(1))
    profile = pd.DataFrame({"uniprot_resnum": common, "ca_deviation_A": dev})
    peak = profile.loc[profile.ca_deviation_A.idxmax()]
    return {
        "n_common_ca": len(common),
        "ca_rmsd": rmsd,
        "max_ca_deviation_A": float(dev.max()),
        "peak_residue": int(peak.uniprot_resnum),
        "moving_region": _contiguous_hotspots(profile),
        "profile": profile,
    }


def _contiguous_hotspots(profile: pd.DataFrame, z: float = 1.5) -> str:
    """Residue ranges whose deviation exceeds mean + z*std -> 'the moving region'."""
    thr = profile.ca_deviation_A.mean() + z * profile.ca_deviation_A.std()
    hot = profile[profile.ca_deviation_A > thr].uniprot_resnum.tolist()
    if not hot:
        return ""
    ranges, s, p = [], hot[0], hot[0]
    for r in hot[1:]:
        if r == p + 1:
            p = r
        else:
            ranges.append((s, p)); s = p = r
    ranges.append((s, p))
    return ";".join(f"{a}-{b}" if a != b else f"{a}" for a, b in ranges)


# --------------------------------------------------------------------------- #
# write out (optionally transformed) segment coordinates
# --------------------------------------------------------------------------- #
def write_segment_pdb(chain, a2u: dict[int, int], lo: int, hi: int, path: Path,
                      R: np.ndarray = None, t: np.ndarray = None,
                      full_chain: bool = False, ca_only: bool = False) -> int:
    """Write the segment residues of `chain` to `path` (PDB), applying the rigid
    transform (R, t) to every atom if given. Returns the residue count written.

    full_chain=True keeps the whole chain (useful to see flanking motion);
    otherwise only residues mapping into [lo, hi] are written."""
    if not _HAVE_GEMMI:
        raise RuntimeError("gemmi is required to write structures")
    R = np.eye(3) if R is None else R
    t = np.zeros(3) if t is None else t
    st = gemmi.Structure()
    model = gemmi.Model("1")
    ch = gemmi.Chain(chain.name)
    n = 0
    for res in chain:
        unp = a2u.get(res.label_seq)
        if not full_chain and (unp is None or not (lo <= unp <= hi)):
            continue
        nr = gemmi.Residue()
        nr.name, nr.seqid, nr.subchain = res.name, res.seqid, res.subchain
        for atom in res:
            if ca_only and atom.name != "CA":
                continue
            na = gemmi.Atom()
            na.name, na.element = atom.name, atom.element
            na.b_iso, na.occ, na.altloc = atom.b_iso, atom.occ, atom.altloc
            xyz = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
            xt = R @ xyz + t
            na.pos = gemmi.Position(*xt)
            nr.add_atom(na)
        if len(nr):
            ch.add_residue(nr); n += 1
    model.add_chain(ch)
    st.add_model(model)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(path))
    return n


def align_representatives(reps: dict, ref_id, out_dir: Path,
                          lo: int, hi: int, full_chain: bool = False) -> list[dict]:
    """Superpose every representative onto the reference cluster's rep and save
    each as a PDB, all in the reference frame so they overlay directly.

    `reps[cluster_id]` = {'chain', 'a2u', 'ca', 'pdb', 'chain_name'}.
    Returns a manifest (one dict per saved structure)."""
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    ref = reps[ref_id]
    manifest = []
    for cid, r in reps.items():
        if cid == ref_id:
            R, t = IDENTITY
            fit_rmsd, n_fit = 0.0, len(r["ca"])
        else:
            common = sorted(set(r["ca"]) & set(ref["ca"]))
            if len(common) < 3:
                log.warning("cluster %s: <3 common Ca with reference, skipping", cid)
                continue
            A = np.array([r["ca"][i] for i in common])
            B = np.array([ref["ca"][i] for i in common])
            R, t, fit_rmsd = kabsch(A, B)
            n_fit = len(common)
        fname = out_dir / f"cluster{cid}_{r['pdb']}_{r['chain_name']}.pdb"
        n_res = write_segment_pdb(r["chain"], r["a2u"], lo, hi, fname, R, t,
                                  full_chain=full_chain)
        manifest.append({"cluster_id": cid, "pdb_id": r["pdb"],
                         "auth_chain": r["chain_name"],
                         "is_reference": cid == ref_id,
                         "n_residues_written": n_res, "n_fit_ca": n_fit,
                         "fit_rmsd_to_ref": round(float(fit_rmsd), 3),
                         "path": str(fname)})
    return manifest


# --------------------------------------------------------------------------- #
# optional external structural aligner (US-align / TM-align) -> TM-score
# --------------------------------------------------------------------------- #
def usalign(pdb_a: Path, pdb_b: Path, exe: str = "USalign",
            save_superposition: Path = None) -> Optional[dict]:
    """Run US-align on two (single-chain) PDBs; parse TM-score, RMSD, aligned
    length. TM-score is superposition-independent, so pre-aligned inputs give
    identical results. `save_superposition` (a path prefix) additionally writes
    US-align's own superposed coordinates via `-o`."""
    if shutil.which(exe) is None:
        log.info("%s not on PATH; skipping TM-score", exe)
        return None
    cmd = [exe, str(pdb_a), str(pdb_b)]
    if save_superposition is not None:
        cmd += ["-o", str(save_superposition)]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=600, check=True)
    except (subprocess.SubprocessError, OSError) as exc:
        log.warning("US-align failed on %s vs %s: %s", pdb_a, pdb_b, exc)
        return None
    out: dict = {}
    for line in res.stdout.splitlines():
        s = line.strip()
        if s.startswith("Aligned length"):
            # Real output pads after "=" (e.g. "Aligned length=  96, RMSD=
            # 0.54, ..."), so a plain whitespace .split() separates "RMSD="
            # from its value into two tokens -- match with regex instead.
            m = re.search(r"RMSD=\s*([\d.]+)", s)
            if m:
                out["tm_rmsd"] = float(m.group(1))
            m = re.search(r"Aligned length=\s*(\d+)", s)
            if m:
                out["aligned_length"] = int(m.group(1))
        elif s.startswith("TM-score=") and "Structure_1" in s:
            out["tm_score_norm1"] = float(s.split("=")[1].split()[0])
        elif s.startswith("TM-score=") and "Structure_2" in s:
            out["tm_score_norm2"] = float(s.split("=")[1].split()[0])
    if "tm_score_norm1" in out and "tm_score_norm2" in out:
        out["tm_score"] = min(out["tm_score_norm1"], out["tm_score_norm2"])
    return out or None

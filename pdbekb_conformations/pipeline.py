"""
pipeline.py -- orchestration.

build_dataset(accessions) walks the API + FTP for a list of UniProt accessions
and writes a tidy, analysis-ready dataset:

    <out>/segments.csv        one row per UniProt segment (+ GLOCON magnitude)
    <out>/clusters.csv        one row per conformational cluster (+ representative)
    <out>/members.csv         one row per PDB chain, labelled by cluster
    <out>/protein_summary.csv one row per protein -> ranking of conformational spread
    <out>/matrices/<acc>_<s>_<e>.npz   normalised score+linkage+labels per segment

add_rmsd(...) then (optionally) superposes cluster representatives locally and
appends Ca-RMSD + a per-residue deviation profile.

flag_bindome(...) joins the per-protein summary back onto a Bindome accession
list to answer: "which Bindome domains have reported alternative conformations?"
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from .fetch import Fetcher
from . import parse, structure

log = logging.getLogger(__name__)


def build_dataset(accessions: Iterable[str],
                  out_dir: str | Path = "conformation_data",
                  cache_dir: str | Path = "cache",
                  pause_s: float = 0.2) -> dict[str, pd.DataFrame]:
    out = Path(out_dir)
    (out / "matrices").mkdir(parents=True, exist_ok=True)
    fx = Fetcher(cache_dir=cache_dir, pause_s=pause_s)

    seg_all, clu_all, mem_all = [], [], []

    for acc in accessions:
        payload = fx.superposition(acc)
        if not payload:
            log.info("%s: no superposition data", acc)
            continue
        tabs = parse.parse_superposition(acc, payload)
        segs = tabs["segments"]
        if segs.empty:
            continue

        # enrich each segment with GLOCON magnitude from the score matrix
        extra = []
        for s in segs.itertuples():
            score_npz = fx.score_npz(acc, s.segment_start, s.segment_end)
            link_npz = fx.linkage_npz(acc, s.segment_start, s.segment_end)
            score = labels = Z = None
            if score_npz:
                Z_from_score = None
                score = parse.find_score_matrix(score_npz)
                labels = parse.find_labels(score_npz, score.shape[0]) if score is not None else None
            if link_npz:
                Z = parse.find_linkage(link_npz)
                if labels is None and Z is not None:
                    labels = parse.find_labels(link_npz, Z.shape[0] + 1)
            m = parse.summarise_segment_matrix(score, labels, tabs["members"],
                                               s.segment_start, s.segment_end)
            extra.append({**{"segment_start": s.segment_start,
                             "segment_end": s.segment_end}, **m,
                          "has_score_matrix": score is not None,
                          "has_linkage": Z is not None})
            # persist normalised matrices for downstream re-use
            np.savez_compressed(
                out / "matrices" / f"{acc}_{s.segment_start}_{s.segment_end}.npz",
                score=score if score is not None else np.array([]),
                linkage=Z if Z is not None else np.array([]),
                labels=labels if labels is not None else np.array([]),
            )
        segs = segs.merge(pd.DataFrame(extra), on=["segment_start", "segment_end"])
        seg_all.append(segs)
        clu_all.append(tabs["clusters"])
        mem_all.append(tabs["members"])

    segments = _concat(seg_all)
    clusters = _concat(clu_all)
    members = _concat(mem_all)
    summary = _protein_summary(segments)

    segments.to_csv(out / "segments.csv", index=False)
    clusters.to_csv(out / "clusters.csv", index=False)
    members.to_csv(out / "members.csv", index=False)
    summary.to_csv(out / "protein_summary.csv", index=False)
    log.info("wrote dataset to %s (%d proteins, %d segments, %d clusters)",
             out, summary.shape[0], len(segments), len(clusters))
    return {"segments": segments, "clusters": clusters,
            "members": members, "protein_summary": summary}


def _concat(frames):
    frames = [f for f in frames if f is not None and not f.empty]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _protein_summary(segments: pd.DataFrame) -> pd.DataFrame:
    """One row per protein: how much conformational spread does it show?"""
    if segments.empty:
        return pd.DataFrame()
    g = segments.groupby("uniprot")
    out = g.agg(
        n_segments=("segment_start", "count"),
        max_clusters=("n_clusters", "max"),
        total_chains=("n_chains", "sum"),
        max_glocon=("glocon_max", "max"),
        max_between_cluster_glocon=("glocon_between_cluster_max", "max"),
    ).reset_index()
    out["multi_conformation"] = out.max_clusters > 1
    return out.sort_values(["max_between_cluster_glocon", "max_glocon"],
                           ascending=False, na_position="last").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# optional RMSD-between-representatives layer
# --------------------------------------------------------------------------- #
def add_rmsd(dataset: dict[str, pd.DataFrame],
             out_dir: str | Path = "conformation_data",
             cache_dir: str | Path = "cache",
             tm_score: bool = False,
             save_aligned: bool = True,
             full_chain: bool = False,
             usalign_exe: str = "USalign") -> pd.DataFrame:
    """For every multi-cluster segment: superpose the cluster representatives,
    save the aligned structures, and record Ca-RMSD (+ optional TM-score) and a
    per-residue deviation profile for each representative pair.

    Parameters
    ----------
    tm_score      also run US-align to record TM-score / TM-RMSD (needs the
                  `USalign` binary on PATH).
    save_aligned  write each segment's representatives, superposed onto the
                  reference cluster, as PDBs under aligned_structures/ so they
                  overlay directly in PyMOL/ChimeraX. A manifest.csv indexes them.
    full_chain    save the whole chain rather than just the segment residues.
    """
    out = Path(out_dir)
    (out / "residue_deviation").mkdir(exist_ok=True, parents=True)
    aln_root = out / "aligned_structures"
    fx = Fetcher(cache_dir=cache_dir)
    clusters = dataset["clusters"]
    rows, manifest_rows = [], []

    for (acc, s, e), grp in clusters.groupby(["uniprot", "segment_start", "segment_end"]):
        grp = grp.dropna(subset=["representative_pdb"])
        if len(grp) < 2:
            continue

        # load each representative once (structure + Ca coords in UniProt numbering)
        reps: dict = {}
        for r in grp.itertuples():
            cif = fx.mmcif(r.representative_pdb)
            sifts = fx.sifts_uniprot_segments(r.representative_pdb)
            if cif is None or sifts is None:
                continue
            a2u = structure.auth_to_uniprot(sifts, r.representative_pdb, acc,
                                            r.representative_auth_chain)
            st, chain = structure.load_chain(cif, r.representative_auth_chain)
            ca = structure.ca_from_chain(chain, a2u, s, e)
            if len(ca) < 3:
                continue
            reps[r.cluster_id] = {"chain": chain, "a2u": a2u, "ca": ca, "_st": st,
                                  "pdb": r.representative_pdb,
                                  "chain_name": r.representative_auth_chain,
                                  "n_members": int(r.n_members)}
        if len(reps) < 2:
            continue

        # reference = most populated cluster (tie -> lowest cluster_id)
        ref_id = max(reps, key=lambda c: (reps[c]["n_members"], -c))

        # save aligned structures (all in the reference frame)
        seg_dir = aln_root / f"{acc}_{s}_{e}"
        if save_aligned:
            man = structure.align_representatives(reps, ref_id, seg_dir, s, e,
                                                  full_chain=full_chain)
            for m in man:
                manifest_rows.append({"uniprot": acc, "segment_start": s,
                                      "segment_end": e, "reference_cluster": ref_id, **m})

        # pairwise metrics (RMSD is frame-invariant; TM-score via US-align)
        ids = sorted(reps)
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                ci, cj = ids[i], ids[j]
                cmp = structure.compare_reps(reps[ci]["ca"], reps[cj]["ca"])
                if cmp is None:
                    continue
                cmp["profile"].to_csv(
                    out / "residue_deviation" / f"{acc}_{s}_{e}_c{ci}_c{cj}.csv",
                    index=False)
                row = {"uniprot": acc, "segment_start": s, "segment_end": e,
                       "cluster_i": ci, "cluster_j": cj,
                       "n_common_ca": cmp["n_common_ca"],
                       "ca_rmsd": cmp["ca_rmsd"],
                       "max_ca_deviation_A": cmp["max_ca_deviation_A"],
                       "peak_residue": cmp["peak_residue"],
                       "moving_region": cmp["moving_region"]}
                if tm_score and save_aligned:
                    pa = seg_dir / f"cluster{ci}_{reps[ci]['pdb']}_{reps[ci]['chain_name']}.pdb"
                    pb = seg_dir / f"cluster{cj}_{reps[cj]['pdb']}_{reps[cj]['chain_name']}.pdb"
                    tm = structure.usalign(pa, pb, exe=usalign_exe)
                    if tm:
                        row.update({"tm_score": tm.get("tm_score"),
                                    "tm_rmsd": tm.get("tm_rmsd"),
                                    "aligned_length": tm.get("aligned_length")})
                rows.append(row)

    rmsd = pd.DataFrame(rows)
    rmsd.to_csv(out / "rmsd_pairs.csv", index=False)
    if save_aligned and manifest_rows:
        pd.DataFrame(manifest_rows).to_csv(aln_root / "manifest.csv", index=False)
    log.info("RMSD layer: %d representative pairs, aligned structures under %s",
             len(rmsd), aln_root)
    return rmsd


# --------------------------------------------------------------------------- #
# Bindome join
# --------------------------------------------------------------------------- #
def flag_bindome(bindome_accessions: Iterable[str],
                 dataset: dict[str, pd.DataFrame],
                 rmsd: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Answer: which Bindome domains have reported alternative conformations?"""
    summary = dataset["protein_summary"].set_index("uniprot")
    rows = []
    rmsd_by_acc = rmsd.groupby("uniprot").ca_rmsd.max() if rmsd is not None and not rmsd.empty else None
    for acc in bindome_accessions:
        if acc in summary.index:
            s = summary.loc[acc]
            rows.append({"uniprot": acc,
                         "in_pdbekb": True,
                         "multi_conformation": bool(s.multi_conformation),
                         "max_clusters": int(s.max_clusters),
                         "max_between_cluster_glocon": s.max_between_cluster_glocon,
                         "max_ca_rmsd": (float(rmsd_by_acc[acc])
                                         if rmsd_by_acc is not None and acc in rmsd_by_acc
                                         else np.nan)})
        else:
            rows.append({"uniprot": acc, "in_pdbekb": False,
                         "multi_conformation": False, "max_clusters": 0,
                         "max_between_cluster_glocon": np.nan, "max_ca_rmsd": np.nan})
    return pd.DataFrame(rows).sort_values(
        ["multi_conformation", "max_between_cluster_glocon"],
        ascending=False, na_position="last").reset_index(drop=True)

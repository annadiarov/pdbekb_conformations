"""
parse.py -- turn the raw PDBe-KB payloads into tidy tables.

Two jobs:

A) Parse the superposition API JSON into three long-format frames:
       segments  (one row per UniProt segment)
       clusters  (one row per conformational cluster)
       members   (one row per PDB chain instance)

B) Read the FTP `.npz` matrices *introspectively*. The internal array names in
   those files are not officially documented, so instead of hard-coding key
   names we detect arrays by shape/dtype:
       - a linkage matrix is the (N-1, 4) float array whose 4th column is an
         integer cluster-size column ending at N (scipy.cluster.hierarchy format)
       - a score/dissimilarity matrix is either a square symmetric (N, N) array
         or a condensed length-N(N-1)/2 vector
       - labels are a 1-D array of length N (chain identifiers)
   This makes the loader robust if PDBe rename keys or ship extra arrays.
"""

from __future__ import annotations

import logging
from typing import Any, Optional

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform

log = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# A. API JSON -> tables
# --------------------------------------------------------------------------- #
def _member_fields(m: Any) -> dict:
    """Normalise one cluster-member record. Handles minor key-name variation."""
    if not isinstance(m, dict):  # some payloads use a bare "pdbid_chain" string
        s = str(m)
        pdb, _, chain = s.partition("_")
        return {"pdb_id": pdb.lower(), "auth_asym_id": chain or None,
                "struct_asym_id": None, "is_representative": False}
    get = lambda *ks: next((m[k] for k in ks if k in m and m[k] is not None), None)
    return {
        "pdb_id": (get("pdb_id", "pdbId", "pdb") or "").lower() or None,
        "auth_asym_id": get("auth_asym_id", "chain_id", "auth_chain_id", "chain"),
        "struct_asym_id": get("struct_asym_id", "struct_asym", "label_asym_id"),
        "is_representative": bool(get("is_representative", "representative", "rep") or False),
    }


def parse_superposition(acc: str, payload: dict) -> dict[str, pd.DataFrame]:
    """API JSON -> {'segments','clusters','members'} DataFrames for one accession."""
    seg_rows, clu_rows, mem_rows = [], [], []

    segments = payload.get(acc, payload) if isinstance(payload, dict) else payload
    if isinstance(segments, dict):                      # tolerate {acc: {...}} nesting
        segments = segments.get(acc, list(segments.values()))
    if not isinstance(segments, list):
        segments = [segments]

    for seg in segments:
        start = int(seg.get("segment_start", seg.get("start")))
        end = int(seg.get("segment_end", seg.get("end")))
        clusters = seg.get("clusters", seg.get("cluster", []))
        n_chains = 0
        for ci, cluster in enumerate(clusters):
            members = cluster if isinstance(cluster, list) else cluster.get("members", [])
            rep_pdb = rep_chain = None
            for m in members:
                f = _member_fields(m)
                if f["is_representative"]:
                    rep_pdb, rep_chain = f["pdb_id"], f["auth_asym_id"]
                mem_rows.append({"uniprot": acc, "segment_start": start,
                                 "segment_end": end, "cluster_id": ci, **f})
                n_chains += 1
            # if the API did not flag a representative, fall back to first member
            if rep_pdb is None and members:
                f0 = _member_fields(members[0])
                rep_pdb, rep_chain = f0["pdb_id"], f0["auth_asym_id"]
            clu_rows.append({"uniprot": acc, "segment_start": start,
                             "segment_end": end, "cluster_id": ci,
                             "n_members": len(members),
                             "representative_pdb": rep_pdb,
                             "representative_auth_chain": rep_chain})
        seg_rows.append({"uniprot": acc, "segment_start": start, "segment_end": end,
                         "segment_length": end - start + 1,
                         "n_clusters": len(clusters), "n_chains": n_chains})

    return {
        "segments": pd.DataFrame(seg_rows),
        "clusters": pd.DataFrame(clu_rows),
        "members": pd.DataFrame(mem_rows),
    }


# --------------------------------------------------------------------------- #
# B. npz introspection
# --------------------------------------------------------------------------- #
def _looks_like_linkage(a: np.ndarray) -> bool:
    if a.ndim != 2 or a.shape[1] != 4 or a.shape[0] < 1:
        return False
    col3 = a[:, 3]
    n = a.shape[0] + 1
    return bool(np.allclose(col3, np.round(col3))
                and col3.min() >= 1 and col3.max() == n
                and a[:, 2].min() >= 0)


def find_linkage(npz: dict) -> Optional[np.ndarray]:
    for v in npz.values():
        v = np.asarray(v)
        if _looks_like_linkage(v):
            return v.astype(float)
    return None


def find_labels(npz: dict, n: int) -> Optional[np.ndarray]:
    """A 1-D array of length n that is NOT purely numeric = chain labels."""
    best = None
    for v in npz.values():
        v = np.asarray(v)
        if v.ndim == 1 and v.shape[0] == n:
            if v.dtype.kind in ("U", "S", "O"):
                return v.astype(str)
            best = best or v            # numeric fallback (indices)
    return best.astype(str) if best is not None else None


def find_score_matrix(npz: dict, n: Optional[int] = None) -> Optional[np.ndarray]:
    """Return a square symmetric dissimilarity matrix (expanding condensed form).

    Confirmed against a live GLOCON `*_score_data.npz`: PDBe ships `scores` as
    a numeric (N, N) array with only the **upper triangle** populated (lower
    triangle and diagonal are 0), alongside an unrelated (N, N) *string*
    `labels` array (pairwise "pdb_A_to_pdb_B" annotations, not a label
    vector -- see `find_labels`). So this must (a) skip non-numeric arrays
    before any numeric comparison, and (b) accept a triangular matrix by
    mirroring it, not just an already-symmetric one."""
    # square first
    for v in npz.values():
        v = np.asarray(v)
        if v.dtype.kind not in "fiu":
            continue
        if v.ndim == 2 and v.shape[0] == v.shape[1] >= 2:
            if n is None or v.shape[0] == n:
                v = v.astype(float)
                if np.allclose(v, v.T, atol=1e-6):
                    return v
                upper = np.triu(v, k=1)
                lower = np.tril(v, k=-1)
                if np.allclose(lower, 0, atol=1e-6) and not np.allclose(upper, 0, atol=1e-6):
                    return upper + upper.T
                if np.allclose(upper, 0, atol=1e-6) and not np.allclose(lower, 0, atol=1e-6):
                    return lower + lower.T
    # condensed
    for v in npz.values():
        v = np.asarray(v)
        if v.ndim == 1 and v.dtype.kind in "fiu":
            m = v.shape[0]
            k = int((1 + np.sqrt(1 + 8 * m)) / 2)   # solve k(k-1)/2 = m
            if k * (k - 1) // 2 == m and (n is None or k == n):
                return squareform(v.astype(float))
    return None


def clusters_from_linkage(Z: np.ndarray, n_clusters: int) -> np.ndarray:
    """Reproduce discrete clusters from the linkage matrix (sanity / re-threshold).

    Not scipy's ``fcluster(..., criterion="maxclust")``: that codepath is
    broken for very small trees (verified: silently returns a single cluster
    for any t when there are exactly 3 leaves, scipy 1.10.1). Instead, cut
    the dendrogram directly -- perform the ``n - n_clusters`` smallest-height
    merges (Z's rows are already sorted by increasing distance) via
    union-find, which is exactly what "maxclust" is supposed to compute."""
    n = Z.shape[0] + 1
    n_clusters = max(1, min(n_clusters, n))
    n_merges = n - n_clusters

    id_to_node = {k: {k} for k in range(n)}   # node id -> set of leaf indices
    next_id = n
    for i, j, _dist, _size in Z[:n_merges, :4]:
        i, j = int(i), int(j)
        merged = id_to_node.pop(i) | id_to_node.pop(j)
        id_to_node[next_id] = merged
        next_id += 1

    labels = np.empty(n, dtype=int)
    for cid, (_node, leaves) in enumerate(id_to_node.items(), start=1):
        for leaf in leaves:
            labels[leaf] = cid
    return labels


def normalise_label(lbl: str) -> str:
    """'1ABC_A' / '1abc.A' / '1abc:A' -> '1abc_a' so npz labels match API members."""
    s = str(lbl).strip().lower().replace(".", "_").replace(":", "_").replace("/", "_")
    return s


def summarise_segment_matrix(score: Optional[np.ndarray],
                             labels: Optional[np.ndarray],
                             members: pd.DataFrame,
                             start: int, end: int) -> dict:
    """GLOCON magnitude metrics for one segment.

    Returns overall max/mean and, if labels can be matched to clusters, the
    max BETWEEN-cluster GLOCON (the biologically interesting number: how far
    apart the distinct conformations are)."""
    out = {"glocon_max": np.nan, "glocon_mean": np.nan,
           "glocon_between_cluster_max": np.nan, "matrix_n": 0}
    if score is None:
        return out
    iu = np.triu_indices_from(score, k=1)
    out["glocon_max"] = float(score[iu].max())
    out["glocon_mean"] = float(score[iu].mean())
    out["matrix_n"] = int(score.shape[0])

    if labels is None:
        return out
    # map each matrix row -> cluster_id via the members table.
    # Confirmed against live FTP npz files: the score/linkage `labels` arrays
    # key chains by struct_asym_id (label_asym_id), not auth_asym_id -- these
    # differ whenever a PDB entry's author and internal chain IDs diverge
    # (common for multi-chain entries; verified matches drop from 100% to as
    # low as ~40% using auth_asym_id alone). Index by both so either
    # convention resolves.
    seg_mem = members[(members.segment_start == start) & (members.segment_end == end)]
    lut = {}
    for r in seg_mem.itertuples():
        lut[normalise_label(f"{r.pdb_id}_{r.struct_asym_id}")] = r.cluster_id
        lut.setdefault(normalise_label(f"{r.pdb_id}_{r.auth_asym_id}"), r.cluster_id)
    row_cluster = np.array([lut.get(normalise_label(l), -1) for l in labels])
    if (row_cluster >= 0).sum() >= 2 and len(set(row_cluster[row_cluster >= 0])) >= 2:
        between = [score[i, j] for i, j in zip(*iu)
                   if row_cluster[i] != row_cluster[j]
                   and row_cluster[i] >= 0 and row_cluster[j] >= 0]
        if between:
            out["glocon_between_cluster_max"] = float(max(between))
    return out

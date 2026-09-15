"""
Self-contained sanity test. No network: fabricates API JSON + npz matrices +
coordinates and checks every piece of logic the real pipeline relies on.
Run:  python tests/test_parsing_synthetic.py
"""
import io
import sys
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pdbekb_conformations import parse, structure   # noqa: E402


def test_api_parsing():
    acc = "P0TEST"
    payload = {acc: [
        {"segment_start": 19, "segment_end": 138, "clusters": [
            [{"pdb_id": "1AAA", "auth_asym_id": "A", "is_representative": True},
             {"pdb_id": "1AAB", "auth_asym_id": "B", "is_representative": False}],
            [{"pdb_id": "2BBB", "auth_asym_id": "A", "is_representative": True}],
        ]},
        {"segment_start": 200, "segment_end": 260, "clusters": [
            [{"pdb_id": "3CCC", "auth_asym_id": "A", "is_representative": True}],
        ]},
    ]}
    t = parse.parse_superposition(acc, payload)
    assert list(t["segments"].n_clusters) == [2, 1]
    assert t["clusters"].shape[0] == 3
    assert t["members"].shape[0] == 4
    rep = t["clusters"].query("segment_start==19 and cluster_id==0").iloc[0]
    assert rep.representative_pdb == "1aaa" and rep.representative_auth_chain == "A"
    print("[ok] API parsing -> segments/clusters/members")
    return acc, t


def _make_npz(labels, points):
    """Fake a PDBe segment npz: score matrix + scipy linkage + labels."""
    D = squareform(pdist(points))
    Z = linkage(pdist(points), method="average")   # UPGMA, like PDBe
    buf = io.BytesIO()
    np.savez(buf, glocon=D, linkage_matrix=Z, chain_labels=np.array(labels))
    buf.seek(0)
    with np.load(buf, allow_pickle=True) as z:
        return {k: z[k] for k in z.files}


def test_npz_introspection_and_glocon(acc, tables):
    # two tight clusters far apart in "GLOCON space"
    rng = np.random.default_rng(0)
    cl0 = rng.normal(0, 0.1, (2, 5))
    cl1 = rng.normal(10, 0.1, (1, 5))
    pts = np.vstack([cl0, cl1])
    labels = ["1aaa_A", "1aab_B", "2bbb_A"]      # match the API members
    npz = _make_npz(labels, pts)

    score = parse.find_score_matrix(npz)
    Z = parse.find_linkage(npz)
    lab = parse.find_labels(npz, score.shape[0])
    assert score.shape == (3, 3)
    assert Z is not None and Z.shape == (2, 4)
    assert list(lab) == labels
    # re-deriving 2 clusters from the linkage must recover the true grouping
    assignments = parse.clusters_from_linkage(Z, 2)
    assert len(set(assignments)) == 2

    summ = parse.summarise_segment_matrix(score, lab, tables["members"], 19, 138)
    assert summ["matrix_n"] == 3
    # between-cluster GLOCON (~10) must exceed overall mean, and be finite
    assert summ["glocon_between_cluster_max"] > summ["glocon_mean"]
    print(f"[ok] npz introspection + GLOCON summary "
          f"(between-cluster max={summ['glocon_between_cluster_max']:.2f})")


def test_kabsch_recovers_rmsd():
    rng = np.random.default_rng(1)
    Q = rng.normal(0, 5, (40, 3))
    # rigid transform of Q -> should give ~0 RMSD after superposition
    theta = 0.7
    Rtrue = np.array([[np.cos(theta), -np.sin(theta), 0],
                      [np.sin(theta),  np.cos(theta), 0],
                      [0, 0, 1]])
    P = (Rtrue @ Q.T).T + np.array([12, -4, 7])
    R, t, rmsd = structure.kabsch(P, Q)
    assert rmsd < 1e-6, rmsd
    # the returned affine must actually map P back onto Q
    P_aln = (R @ P.T).T + t
    assert np.allclose(P_aln, Q, atol=1e-6)
    print(f"[ok] Kabsch returns re-usable affine (rmsd={rmsd:.2e} A)")


def test_per_residue_localisation():
    """Two 'conformers' identical except a hinge region that swings out."""
    rng = np.random.default_rng(2)
    n = 60
    base = np.cumsum(rng.normal(0, 1.2, (n, 3)), axis=0)   # a fake backbone
    ca_a = {19 + i: base[i].copy() for i in range(n)}
    moved = base.copy()
    moved[30:40] += np.array([8.0, 0, 0])                  # residues 49-58 move
    ca_b = {19 + i: moved[i].copy() for i in range(n)}
    cmp = structure.compare_reps(ca_a, ca_b)
    assert cmp is not None
    assert 49 <= cmp["peak_residue"] <= 58
    reg = cmp["moving_region"]
    assert reg, "expected a moving region to be reported"
    print(f"[ok] per-residue deviation localises hinge "
          f"(peak={cmp['peak_residue']}, region={reg}, rmsd={cmp['ca_rmsd']:.2f} A)")


def test_align_and_save_roundtrip(tmp):
    """Build two toy gemmi chains, align+save, and confirm the written
    coordinates are actually superposed (RMSD ~ 0 after saving)."""
    import gemmi

    def make_chain(name, coords):
        st = gemmi.Structure(); m = gemmi.Model("1"); ch = gemmi.Chain(name)
        for k, xyz in enumerate(coords, start=19):
            r = gemmi.Residue(); r.name = "ALA"; r.seqid = gemmi.SeqId(k, " ")
            r.label_seq = k   # a2u is keyed by label_seq_id, like real SIFTS data
            a = gemmi.Atom(); a.name = "CA"; a.element = gemmi.Element("C")
            a.pos = gemmi.Position(*xyz); r.add_atom(a); ch.add_residue(r)
        m.add_chain(ch); st.add_model(m); st.setup_entities()
        return st, st[0][0]

    rng = np.random.default_rng(3)
    base = rng.normal(0, 4, (20, 3))
    theta = 1.1
    Rt = np.array([[np.cos(theta), -np.sin(theta), 0],
                   [np.sin(theta),  np.cos(theta), 0], [0, 0, 1]])
    moved = (Rt @ base.T).T + np.array([20, 5, -3])   # same shape, displaced

    st0, ch0 = make_chain("A", base)
    st1, ch1 = make_chain("A", moved)
    a2u = {19 + i: 19 + i for i in range(20)}
    reps = {
        0: {"chain": ch0, "a2u": a2u, "ca": {19 + i: base[i] for i in range(20)},
            "pdb": "1ref", "chain_name": "A", "_st": st0},
        1: {"chain": ch1, "a2u": a2u, "ca": {19 + i: moved[i] for i in range(20)},
            "pdb": "2mov", "chain_name": "A", "_st": st1},
    }
    man = structure.align_representatives(reps, 0, Path(tmp), 19, 38)
    assert len(man) == 2
    # reload the two written PDBs and check they now overlay
    import gemmi as g
    def read_ca(p):
        s = g.read_structure(str(p))
        return np.array([r.find_atom("CA", "*").pos.tolist() for r in s[0][0]])
    ref = next(m for m in man if m["is_reference"])
    mov = next(m for m in man if not m["is_reference"])
    d = np.sqrt(((read_ca(ref["path"]) - read_ca(mov["path"])) ** 2).sum(1).mean())
    assert d < 1e-4, d
    assert mov["fit_rmsd_to_ref"] < 1e-3
    print(f"[ok] aligned PDBs written and overlay correctly (rmsd={d:.2e} A)")


if __name__ == "__main__":
    import tempfile
    acc, tables = test_api_parsing()
    test_npz_introspection_and_glocon(acc, tables)
    test_kabsch_recovers_rmsd()
    test_per_residue_localisation()
    with tempfile.TemporaryDirectory() as td:
        test_align_and_save_roundtrip(td)
    print("\nAll synthetic checks passed.")

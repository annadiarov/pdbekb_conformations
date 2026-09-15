#!/usr/bin/env python3
"""
A small, end-to-end example of the Python API.

Runs the two stages on two accessions and prints, for each one, the question
the pipeline is meant to answer:

  1. does this protein have more than one reported conformation?
  2. how far apart are the conformations (GLOCON, Ca-RMSD)?
  3. which residues move?

Usage
-----
    conda activate pdbekb
    python examples/quickstart.py                 # writes examples/output/

Fetches everything live, into a cache of its own under examples/output/cache,
so it is a genuine cold-start run that neither touches nor depends on the repo's
`cache/`. Needs outbound internet; takes a couple of minutes, most of it
downloading the four mmCIFs stage 2 needs. Delete examples/output/ to rerun cold.

TM-score needs the US-align binary (`conda install -c bioconda usalign`); if it
isn't found the example still completes, just without the TM-score columns.

The equivalent command line is:

    python run_pipeline.py --accession-file examples/example_accessions.txt \
        --out examples/output --cache examples/output/cache --rmsd --tmscore \
        --bindome-flag examples/example_accessions.txt -v
"""
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))       # the package is used from the repo, not installed

from pdbekb_conformations import pipeline  # noqa: E402

OUT = REPO / "examples" / "output"
ACCESSIONS = ["O00214", "O43488"]

# A cache of this example's own, so the run is a real cold start rather than a
# replay of whatever the repo's `cache/` happens to hold.
CACHE = OUT / "cache"


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    # ------------------------------------------------------------------ stage 1
    # Clustering tables + GLOCON. Network-bound but cheap: no structures are
    # downloaded. Writes segments/clusters/members/protein_summary.csv and
    # matrices/*.npz under `out_dir`, and returns the same tables as DataFrames.
    ds = pipeline.build_dataset(
        ACCESSIONS,
        out_dir=OUT,
        cache_dir=CACHE,
        pause_s=0.2,          # be polite to EBI; 0 is fine once CACHE is warm
    )

    summary = ds["protein_summary"]
    print("\n=== 1. Which proteins have more than one conformation? ===")
    print(summary[["uniprot", "n_segments", "max_clusters",
                   "max_between_cluster_glocon", "multi_conformation"]]
          .to_string(index=False))

    # `multi_conformation` is per protein; the segment is the real unit, because
    # a protein can be rigid in one domain and flexible in another.
    print("\n=== 2. How different are they? (per segment, GLOCON) ===")
    print(ds["segments"][["uniprot", "segment_start", "segment_end", "n_clusters",
                          "n_chains", "glocon_between_cluster_max"]]
          .to_string(index=False))

    # Each cluster has a representative PDB chain -- that is what stage 2 superposes.
    print("\n=== ...and the cluster representatives ===")
    print(ds["clusters"].to_string(index=False))

    # ------------------------------------------------------------------ stage 2
    # Local Kabsch superposition of the representatives -> true Ca-RMSD and a
    # per-residue deviation profile. Only multi-cluster segments do any work, so
    # O43488 (1 cluster) is silently skipped.
    rmsd = pipeline.add_rmsd(
        ds,
        out_dir=OUT,
        cache_dir=CACHE,
        tm_score=True,        # needs USalign; skipped with a log note if absent
        save_aligned=True,    # writes examples/output/aligned_structures/
    )

    print("\n=== 3. Which residues actually move? ===")
    cols = ["uniprot", "segment_start", "segment_end", "cluster_i", "cluster_j",
            "n_common_ca", "ca_rmsd", "max_ca_deviation_A", "peak_residue",
            "moving_region"]
    # tm_score/tm_rmsd/aligned_length only exist if US-align actually ran.
    tm_cols = [c for c in ("tm_score", "tm_rmsd", "aligned_length") if c in rmsd.columns]
    if tm_cols:
        cols += tm_cols
    else:
        print("(no TM-score columns -- USalign was not found; see requirements.txt)")
    print(rmsd[cols].to_string(index=False))

    # `moving_region` is the pre-summarised answer. The full profile behind it --
    # one row per residue -- is on disk; here is the top of it for the segment
    # with the largest local displacement.
    worst = rmsd.loc[rmsd.max_ca_deviation_A.idxmax()]
    prof_path = (OUT / "residue_deviation" /
                 f"{worst.uniprot}_{worst.segment_start}_{worst.segment_end}"
                 f"_c{worst.cluster_i}_c{worst.cluster_j}.csv")
    prof = pd.read_csv(prof_path)
    print(f"\nPer-residue profile for {worst.uniprot} "
          f"{worst.segment_start}-{worst.segment_end} "
          f"(cluster {worst.cluster_i} vs {worst.cluster_j}), "
          f"{len(prof)} residues, 10 most displaced:")
    print(prof.nlargest(10, "ca_deviation_A").to_string(index=False))
    print(f"\n(residue numbers are UniProt numbering; full file: {prof_path})")

    # ------------------------------------------------------- the one-row answer
    # A yes/no row per accession, including accessions PDBe-KB has nothing for.
    flags = pipeline.flag_bindome(ACCESSIONS, ds, rmsd)
    flags.to_csv(OUT / "bindome_conformation_flags.csv", index=False)
    print("\n=== Per-accession summary ===")
    print(flags.to_string(index=False))

    print(f"\nAll outputs under {OUT}/")
    print("Superposed representatives (drop a folder into PyMOL to see the motion):")
    for p in sorted((OUT / "aligned_structures").rglob("*.pdb")):
        print(f"  {p.relative_to(REPO)}")


if __name__ == "__main__":
    main()

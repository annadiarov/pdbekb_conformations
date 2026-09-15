#!/usr/bin/env python3
"""
Command-line driver.

Examples
--------
# 1. Build the clustering tables for a few accessions
python run_pipeline.py --accessions Q14676 Q96Y14 P38398 --out conformation_data

# 2. Same, reading accessions (one per line) from a Bindome list
python run_pipeline.py --accession-file bindome_accessions.txt --out conformation_data

# 3. Also compute local Ca-RMSD between cluster representatives (needs network
#    access to PDBe/RCSB for mmCIF; add --tmscore to also run US-align)
python run_pipeline.py --accession-file bindome_accessions.txt --rmsd --tmscore

# 4. Flag which Bindome accessions have alternative conformations
python run_pipeline.py --accession-file bindome_accessions.txt --rmsd \
       --bindome-flag bindome_accessions.txt
"""
import argparse
import logging
from pathlib import Path

from pdbekb_conformations import pipeline


def read_list(path):
    return [l.strip() for l in Path(path).read_text().splitlines()
            if l.strip() and not l.startswith("#")]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--accessions", nargs="+", help="UniProt accessions")
    src.add_argument("--accession-file", help="file with one accession per line")
    ap.add_argument("--out", default="conformation_data")
    ap.add_argument("--cache", default="cache")
    ap.add_argument("--pause", type=float, default=0.2, help="seconds between requests")
    ap.add_argument("--rmsd", action="store_true", help="compute Ca-RMSD between reps")
    ap.add_argument("--tmscore", action="store_true",
                    help="also record TM-score via US-align (needs USalign on PATH)")
    ap.add_argument("--no-save-aligned", action="store_true",
                    help="do not write the superposed representative PDBs")
    ap.add_argument("--full-chain", action="store_true",
                    help="save whole chains rather than just the segment residues")
    ap.add_argument("--usalign-exe", default="USalign",
                    help="US-align binary; a name resolved on PATH (the default, "
                         "which is what `conda install -c bioconda usalign` "
                         "provides) or an explicit path to the binary")
    ap.add_argument("--bindome-flag", help="accession list to flag for alt conformations")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")

    accs = args.accessions or read_list(args.accession_file)
    ds = pipeline.build_dataset(accs, out_dir=args.out, cache_dir=args.cache,
                                pause_s=args.pause)
    print(f"\nProteins with >1 conformational cluster in any segment: "
          f"{int(ds['protein_summary'].multi_conformation.sum())}/"
          f"{len(ds['protein_summary'])}")

    rmsd = None
    if args.rmsd:
        rmsd = pipeline.add_rmsd(ds, out_dir=args.out, cache_dir=args.cache,
                                 tm_score=args.tmscore,
                                 save_aligned=not args.no_save_aligned,
                                 full_chain=args.full_chain,
                                 usalign_exe=args.usalign_exe)
        # Report what was actually written: --tmscore is a request, but US-align
        # is optional and add_rmsd skips it (with a log note) if the binary
        # isn't there, so key off the columns rather than the flag.
        if args.tmscore:
            extra = (" (+TM-score)" if "tm_score" in rmsd.columns else
                     f" (no TM-score: '{args.usalign_exe}' not found; see --usalign-exe)")
        else:
            extra = ""
        print(f"Computed {len(rmsd)} representative-pair RMSDs{extra} "
              f"-> {args.out}/rmsd_pairs.csv")
        if not args.no_save_aligned:
            print(f"Aligned representative PDBs -> {args.out}/aligned_structures/")

    if args.bindome_flag:
        flags = pipeline.flag_bindome(read_list(args.bindome_flag), ds, rmsd)
        outp = Path(args.out) / "bindome_conformation_flags.csv"
        flags.to_csv(outp, index=False)
        print(f"Bindome flags -> {outp} "
              f"({int(flags.multi_conformation.sum())} with alt conformations)")


if __name__ == "__main__":
    main()

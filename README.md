# pdbekb_conformations

Given a list of UniProt accessions, harvest **PDBe-KB's conformational
clustering data** and answer:

1. which protein domains have more than one reported conformation,
2. how different are those conformations (GLOCON magnitude, and optionally a
   real Cα-RMSD / TM-score), and
3. **which residues actually move** — the part you'd use to steer binder
   sampling toward the flexible region.

---

## 1. Setup

```bash
conda create -n pdbekb python=3.11
conda activate pdbekb
pip install -r requirements.txt
```

**US-align** is optional — only `--tmscore` uses it. It is an external binary,
not a pip package, so it is not in `requirements.txt`:

```bash
conda install -c bioconda usalign        # installs `USalign` on PATH
```

That puts it where the code looks by default, so `--tmscore` then works with no
further configuration. If you'd rather build it yourself, or you're on a machine
where it's already installed somewhere unusual, point `--usalign-exe` at the
binary (§5). Without it, everything except `--tmscore` still runs; the TM-score
columns are simply skipped, with a note in the log.

## 2. Input

A list of UniProt accessions — nothing else. Either:

```bash
--accessions Q14676 Q96Y14 P38398      # a few, on the command line
--accession-file my_accessions.txt     # one accession per line; blank lines
                                        # and lines starting with '#' are skipped
```

Everything else is fetched live over the network (PDBe-KB API + FTP, and —
only if you pass `--rmsd` — PDBe/RCSB mmCIF files and the PDBe SIFTS API), so
the machine you run this on needs outbound internet. See [§6](#6-where-to-run-it-hpc-note)
if that's not true of your compute nodes.

## 3. Quick start

```bash
# Stage 1 only: clustering + GLOCON tables for your whole list (cheap, no structures)
python run_pipeline.py --accession-file list_bioemu.txt --out conformation_data -v

# Stage 2 as well, for a shortlist: local Ca-RMSD + TM-score + aligned PDBs
python run_pipeline.py --accession-file top_movers.txt --out conformation_data \
    --rmsd --tmscore -v

# Add a per-accession yes/no summary on top of either of the above
python run_pipeline.py --accession-file list_bioemu.txt --out conformation_data \
    --bindome-flag list_bioemu.txt -v
```

`cache/` is reused across runs (incremental and offline once populated), so
re-running any of the above is cheap.

**Worked example.** [`examples/`](examples/) runs both stages on two contrasting
accessions, fetching everything live into a cache of its own, and walks through
how to read each output file:

```bash
python examples/quickstart.py     # or the equivalent CLI, see examples/README.md
```

## 4. Outputs

Everything lands under `--out` (default `conformation_data/`). Two stages:

* **Stage 1** (always runs): the clustering tables + GLOCON magnitude.
* **Stage 2** (only with `--rmsd`): local superposition, Cα-RMSD, and
  per-residue deviation.

```
conformation_data/
├── segments.csv                                     stage 1
├── clusters.csv                                      "
├── members.csv                                       "
├── protein_summary.csv                                "
├── matrices/{acc}_{start}_{end}.npz                   "
├── rmsd_pairs.csv                                   stage 2
├── residue_deviation/{acc}_{s}_{e}_c{i}_c{j}.csv       "
├── aligned_structures/{acc}_{s}_{e}/cluster{id}_{pdb}_{chain}.pdb
├── aligned_structures/manifest.csv                     "
└── bindome_conformation_flags.csv                   --bindome-flag
```

Long-format and relational — one row per finest unit, easy to filter and
join on `uniprot` / `segment_start` / `segment_end` / `cluster_id`.

| File | One row = | Key columns | Needs |
|---|---|---|---|
| `segments.csv` | one UniProt segment | `uniprot, segment_start, segment_end, n_clusters, n_chains, glocon_max, glocon_mean, glocon_between_cluster_max, matrix_n, has_score_matrix, has_linkage` | — |
| `clusters.csv` | one conformational cluster | `uniprot, segment_start, segment_end, cluster_id, n_members, representative_pdb, representative_auth_chain` | — |
| `members.csv` | one PDB chain in a cluster | `uniprot, segment_start, segment_end, cluster_id, pdb_id, auth_asym_id, struct_asym_id, is_representative` | — |
| `protein_summary.csv` | one protein, ranked | `uniprot, n_segments, max_clusters, total_chains, max_glocon, max_between_cluster_glocon, multi_conformation` | — |
| `matrices/{acc}_{s}_{e}.npz` | one segment | arrays `score` (GLOCON matrix), `linkage` (scipy linkage tree), `labels` | — |
| `rmsd_pairs.csv` | one cluster-pair in a segment | `uniprot, segment_start, segment_end, cluster_i, cluster_j, n_common_ca, ca_rmsd, max_ca_deviation_A, peak_residue, moving_region` (+`tm_score, tm_rmsd, aligned_length` with `--tmscore`) | `--rmsd` |
| `residue_deviation/{acc}_{s}_{e}_c{i}_c{j}.csv` | one residue | `uniprot_resnum, ca_deviation_A` | `--rmsd` |
| `aligned_structures/{acc}_{s}_{e}/cluster{id}_{pdb}_{chain}.pdb` | one representative, superposed | coordinates only | `--rmsd` |
| `aligned_structures/manifest.csv` | one saved structure | `uniprot, segment_start, segment_end, reference_cluster, cluster_id, pdb_id, auth_chain, is_reference, n_residues_written, n_fit_ca, fit_rmsd_to_ref, path` | `--rmsd` |
| `bindome_conformation_flags.csv` | one accession from your `--bindome-flag` list | `uniprot, in_pdbekb, multi_conformation, max_clusters, max_between_cluster_glocon, max_ca_rmsd` | `--bindome-flag` |

**Things that aren't obvious from the column names:**

* **`glocon_between_cluster_max`** (in `segments.csv`/`protein_summary.csv`)
  is the number that actually matters for ranking: the largest GLOCON
  distance *between* two different clusters, i.e. how far apart the
  distinct conformations are — not just the spread within one cluster.
  `protein_summary.csv` is pre-sorted by this column.
* **`moving_region`** / **`peak_residue`** (in `rmsd_pairs.csv`) are the
  pre-summarized answer to "which residues differ": `moving_region` lists
  the UniProt residue ranges whose Cα moved more than
  mean + 1.5·std after superposing the two representatives (e.g.
  `"48-49;56;103-104"`); `peak_residue` is the single worst one. For the
  full per-residue profile behind that summary, see `residue_deviation/*.csv`.
* **All residue numbers everywhere in the stage-2 output** (`uniprot_resnum`,
  `peak_residue`, `moving_region`, `segment_start`/`segment_end`) are
  **UniProt sequence numbering** (via SIFTS), not PDB author numbering — so
  they're directly comparable across different PDB entries of the same
  protein and map straight onto the UniProt sequence.
* **`aligned_structures/`**: every cluster representative in a segment is
  superposed onto that segment's *reference* cluster (its most populated
  one, `reference_cluster` in `manifest.csv`) and written out in that shared
  frame — load a segment's folder straight into PyMOL/ChimeraX and the
  conformers overlay with no further fitting. Only the segment's own
  residues are written by default; pass `--full-chain` to keep whole chains
  (useful to see flanking motion). `ca_rmsd` in `rmsd_pairs.csv` is
  frame-invariant, so it doesn't depend on this saved frame.
* **`--bindome-flag FILE`** doesn't fetch anything new — it looks `FILE`'s
  accessions up in the dataset *already built in the same command* and
  writes one row per accession, including `in_pdbekb: False` rows for ones
  with no data (unlike `protein_summary.csv`, which just omits them). So
  `FILE` should be the same list (or a subset of it) you passed to
  `--accessions`/`--accession-file` — an accession never fetched in this run
  will read as `in_pdbekb: False` regardless of whether PDBe-KB actually has
  it.

## 5. All flags

| Flag | Default | Meaning |
|---|---|---|
| `--accessions ACC [ACC ...]` | — | Accessions on the command line. Mutually exclusive with `--accession-file`. |
| `--accession-file FILE` | — | One accession per line (`#` comments and blank lines skipped). |
| `--out DIR` | `conformation_data` | Output directory (§4). |
| `--cache DIR` | `cache` | HTTP response cache; reused across runs. |
| `--pause SECONDS` | `0.2` | Delay between requests to EBI. Raise it for big runs; `--pause 0` is fine when reading purely from a warm cache (e.g. on an offline compute node). |
| `--rmsd` | off | Run stage 2: local Cα-RMSD + per-residue deviation + aligned PDBs. Needs network access to PDBe/RCSB. |
| `--tmscore` | off | Adds `tm_score`/`tm_rmsd`/`aligned_length` to `rmsd_pairs.csv` via US-align. Only takes effect together with `--rmsd` — passed alone it does nothing. If the binary can't be found the run continues without the TM-score columns. |
| `--usalign-exe PATH` | `USalign` | The US-align binary: a name looked up on `PATH` (the default works after `conda install -c bioconda usalign`) or an explicit path. Only used with `--tmscore`. |
| `--no-save-aligned` | off | With `--rmsd`, skip writing the superposed PDB files (`rmsd_pairs.csv`/`residue_deviation/` are still written). |
| `--full-chain` | off | With `--rmsd`, save whole chains in `aligned_structures/` rather than just the segment's residues. |
| `--bindome-flag FILE` | — | Write `bindome_conformation_flags.csv` for this accession list (§4 note above). |
| `-v`, `--verbose` | off | INFO-level logging (recommended — otherwise only warnings/errors print). |

## 6. Where to run it (HPC note)

Two stages, two different profiles:

* **Stage 1** (clustering tables): network-bound, tiny CPU — a laptop or a
  login node is fine even for a large list. Budget disk for `cache/`.
* **Stage 2** (`--rmsd`): downloads mmCIFs (bandwidth + potentially many GB
  of cache) and runs many small local superpositions — a cluster helps here,
  mainly for disk and parallel local alignment, not for hitting EBI harder.

The usual HPC gotcha: **compute nodes often have no outbound internet**
(true on many EPFL/SCITAS partitions). So decouple the two stages:

1. Run **stage 1 on a login node / laptop** (has internet), caching to
   shared scratch.
2. Run **stage 2 as an array job on compute nodes** (no internet needed —
   everything reads from the warm cache) with `--pause 0`, sharding the
   accession list across tasks.

Practical pattern: run stage 1 on everything (cheap), rank on
`protein_summary.csv`'s `max_between_cluster_glocon`, then run stage 2 only
on the top movers.

## 7. Method notes

* **GLOCON, not RMSD.** PDBe-KB clusters on GLOCON — a global backbone
  dissimilarity from pairwise Cα *distances* (superposition-agnostic), via
  UPGMA. It is **not convertible to RMSD**, and PDBe does not ship an RMSD
  matrix. `--rmsd` is a separate, optional local step (Kabsch superposition
  of cluster representatives) that gives you a true Cα-RMSD, and
  `--tmscore` adds TM-score via US-align on top of that (both flagged by the
  PDBe team as more robust than RMSD for cross-structure comparison).
* **Never enumerate the FTP tree.** The FTP paths for the precomputed
  matrices contain a `segment_<start>_<end>` fragment you can't know in
  advance — the fix is to always call the superposition Graph-API first
  (`GET .../graph-api/uniprot/superposition/{accession}`), which returns
  every segment's `(start, end)` *and* its clusters/members; those values
  then determine the FTP paths deterministically.
* **Format quirks confirmed against live data**, baked into the parser so
  you don't need to re-derive them: the FTP `*_score_data.npz` score matrix
  is upper-triangular (not symmetric); its chain labels are keyed by
  `struct_asym_id`, not `auth_asym_id`; and PDBe's SIFTS mapping leaves
  `author_residue_number` null for a large fraction of chains, so residue
  correspondence is built from `residue_number` (label_seq_id) instead,
  matched against gemmi's `Residue.label_seq`.
* A benchmarking set of distinct monomer conformers exists at
  `ftp.ebi.ac.uk/pub/databases/pdbe-kb/benchmarking/distinct-monomer-conformers/`
  if you want a curated positive-control set.

## 8. Tests

```bash
python tests/test_parsing_synthetic.py
```

A self-contained, no-network sanity check (fabricated API JSON + npz
matrices + coordinates) covering the parsing, clustering, Kabsch, and
per-residue localisation logic. Run as a plain script, not via `pytest`
(the functions deliberately chain outputs between checks rather than using
fixtures).

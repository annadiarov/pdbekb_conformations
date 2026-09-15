# A small worked example

Two accessions, both stages, from a cold start.

```bash
conda activate pdbekb
python examples/quickstart.py            # Python API, writes examples/output/
```

or the same thing from the command line:

```bash
python run_pipeline.py --accession-file examples/example_accessions.txt \
    --out examples/output --cache examples/output/cache --rmsd --tmscore \
    --bindome-flag examples/example_accessions.txt -v
```

It caches into `examples/output/cache` rather than the repo's `cache/`, so this
is a real first run: it fetches from PDBe-KB and downloads the four mmCIFs that
stage 2 needs (~7 MB), and so needs outbound internet. `rm -rf examples/output`
resets it to cold.

TM-score needs the US-align binary, which is not a pip dependency:

```bash
conda install -c bioconda usalign        # installs `USalign` on PATH
```

Without it the example still finishes; the `tm_score` / `tm_rmsd` /
`aligned_length` columns are simply absent, and the log says why.

The list is deliberately a contrast:

| Accession | Protein | Why it's here |
|---|---|---|
| `O00214` | galectin-8 (LGALS8) | 2 segments, 2 conformational clusters each — the "yes, it moves" case |
| `O43488` | aflatoxin B1 aldehyde reductase (AKR7A2) | 1 segment, 1 cluster — the "nothing reported" case |

Any other accession works the same way.

## What comes out, and how to read it

**1. Does it move?** — `protein_summary.csv`

```
uniprot  n_segments  max_clusters  max_between_cluster_glocon  multi_conformation
 O00214           2             2               118451.296875                True
 O43488           1             1                         NaN               False
```

`O43488` has one cluster, so there is no *between*-cluster distance to report
and `max_between_cluster_glocon` is empty — an empty value here means "only one
conformation", not "missing data".

**2. Where, and how far?** — `segments.csv`, then `rmsd_pairs.csv`

The segment, not the protein, is the unit that matters: a protein can be rigid
in one domain and flexible in another. Galectin-8's two segments both split
into 2 clusters, but by very different amounts:

```
uniprot  segment_start  segment_end  n_clusters  n_chains  glocon_between_cluster_max
 O00214              1          317           2        59               118451.296875
 O00214            171          317           2         7                 2627.399902
```

GLOCON is superposition-agnostic and **not** convertible to Å (§7 of the main
README), so it ranks but doesn't measure. `--rmsd` adds the measurement, by
superposing the two clusters' representative chains:

```
uniprot  seg    cluster_i cluster_j  n_common_ca  ca_rmsd  max_ca_deviation_A  peak_residue  moving_region  tm_score  tm_rmsd  aligned_length
 O00214  1-317          0         1          148 0.956638            5.668759             8   8-11;154-155   0.50314     0.74             147
 O00214  171-317        0         1          136 3.183565           24.416922           177        176-177   0.89098     1.12             134
```

Note the metrics disagree in an informative way. The 1–317 segment ranks far
higher on GLOCON but has the *lower* Cα-RMSD — GLOCON is computed over all 59
chains in the segment, while the RMSD compares only the two representatives.
Use GLOCON to rank candidates and RMSD to characterise the pair you picked.

The last three columns come from `--tmscore` (US-align) and disagree again, for
a different reason: 1–317 has the *worse* TM-score (0.50 vs 0.89) despite the
better RMSD, because TM-score is normalised by length and penalises the parts
that don't superpose, while `ca_rmsd` averages over only the 148 Cα the two
structures share. Three numbers, three questions: how far apart PDBe-KB's
clustering put them, how far the shared atoms move, and whether they're still
the same fold.

**3. Which residues?** — `moving_region`, and `residue_deviation/*.csv`

`moving_region` is the summarised answer (residues more than mean + 1.5·std
from the mean deviation); the per-residue profile behind it is on disk:

```
 uniprot_resnum  ca_deviation_A
            177       24.416922
            176       23.810185
            253        3.152309
            193        2.843178
```

All numbers are **UniProt** numbering, so they map straight onto the sequence
and are comparable across PDB entries.

This particular hit is also a good lesson in reading the output: residues
176–177 sit at the very start of the 171–317 segment, and a 24 Å swing over two
terminal residues with a ~2.5 Å background is a floppy chain end, not a domain
rearrangement. Before handing a `moving_region` to binder sampling, check
whether it is terminal (like this one) or interior and contiguous (like
`8-11;154-155` in the other segment). The `ca_rmsd` and `n_common_ca` columns
give the context needed for that call.

**4. The one-row-per-accession answer** — `bindome_conformation_flags.csv`

```
uniprot  in_pdbekb  multi_conformation  max_clusters  max_between_cluster_glocon  max_ca_rmsd
 O00214       True                True             2               118451.296875     3.183565
 O43488       True               False             1                         NaN          NaN
```

This is the table to join back onto a Bindome list. Accessions PDBe-KB knows
nothing about appear here as `in_pdbekb: False` rather than being dropped.

**5. See it** — `aligned_structures/`

```
examples/output/aligned_structures/O00214_171_317/cluster0_8hl9_A.pdb
examples/output/aligned_structures/O00214_171_317/cluster1_2yro_A.pdb
```

Both are already in the same frame, so `pymol examples/output/aligned_structures/O00214_171_317/*.pdb`
overlays them with no further fitting.

## Scaling this up

For a real list, don't run `--rmsd` on everything: stage 1 is cheap and stage 2
downloads mmCIFs. Run stage 1 over the whole list, rank on
`protein_summary.csv`'s `max_between_cluster_glocon`, and run stage 2 on the top
of that ranking only. §6 of the main README has the HPC pattern for splitting
the two stages across a login node and compute nodes.

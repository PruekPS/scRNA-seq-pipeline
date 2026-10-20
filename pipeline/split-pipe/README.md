# split-pipe (Parse Biosciences) — not included here

**split-pipe** is Parse Biosciences' own command-line pipeline for Evercode
data: barcode correction, alignment, and counting, the same logic as
Parse's cloud tool Trailmaker. It is licensed software (Parse's intellectual
property), available to Evercode customers through Parse's
[Pipeline Download Request Form](https://www.parsebiosciences.com/pipeline-download/).
**No split-pipe code, installer, or run script is distributed in this
repository.** Install and run it by following Parse's own documentation
(Support Suite: "Pipeline Installation", "Running the Pipeline", "Pipeline
Options"), on an HPC compute node rather than a login node.

This pipeline starts where split-pipe stops. What it expects from
split-pipe's output (checked against split-pipe v1.9.0 and real Parse output):

## Whole transcriptome (WT)

- Each sublibrary is processed separately, then merged with split-pipe's
  combine mode. Cells are named by their barcode wells, `NN_NN_NN`; after
  combining, each name gets a sublibrary tag `__s1`, `__s2`, ... from the
  sublibrary's position in the combine list.
- `scripts/00_load_parse_output.py` reads a `DGE_filtered/` folder
  (`count_matrix.mtx`, `all_genes.csv`, `cell_metadata.csv`). Point it at
  either `<combined>/all-sample/DGE_filtered` or the per-sample folders,
  never both (all-sample holds the same cells again).
- Study design (patient, time point, tissue, ...) is not in split-pipe's
  output: pass it with `--sample-sheet` (a CSV with a `sample` column whose
  names match the split-pipe sample names exactly).
- The `__sN` tag becomes `obs["sublibrary"]`, the technical batch to give
  Harmony (`02b_harmony_integration.py --batch-key sublibrary`).

## Evercode TCR

- Each TCR sublibrary is run against its matching WT sublibrary (split-pipe's
  parent-directory option), then combined. **Combine TCR sublibraries in the
  same order as the WT sublibraries**, because the `__sN` tags come from list
  position; `scripts/12_tcr_attach.py` checks every sublibrary and stops if
  one does not line up.
- `scripts/12_tcr_attach.py` reads the per-sample AIRR tables
  (`<sample>/TCR_filtered/tcr_annotation_airr.tsv`, or `.tsv.gz`), skipping
  the duplicate `all-sample` folder.
- TCR libraries need reads long enough to cover the CDR3 (read 1) and the cell
  barcodes (read 2); see Parse's Evercode TCR documentation for the required
  sequencing configuration.

## Testing without real data

`tests/synthetic_splitpipe.py` writes synthetic data in this output layout and
`tests/integration_smoke.sh` runs the whole chain on it, checking the results
against the planted truth.

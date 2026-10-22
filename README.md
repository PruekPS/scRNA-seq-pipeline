# Parse Evercode scRNA-seq: from split-pipe output to biology, step by step

This repository takes **Parse Biosciences Evercode single-cell data** — whole
transcriptome, plus Evercode TCR — from Parse's demultiplexed output to
quality-controlled, annotated cells and the analyses built on them:
cell-cell communication, trajectories, differential expression between
conditions, and T-cell clonotypes. The heavy lifting runs on an HPC cluster
with SLURM (built on UF's [HiPerGator](https://www.rc.ufl.edu/)); small tests
run on a laptop.

It is written so you can run it **without a coding assistant**, and so you
understand *why* each step exists, not just which command to type.

Work through the pages in order.

---

## What you will end up with

```
results/my_run/
  04b_annotated.h5ad        every cell QC-filtered, clustered, and labelled (or honestly "Ambiguous")
  cellchat_results/         ligand-receptor signalling between cell types
  pseudobulk_de/            per-cell-type differential expression, volcano plots, heatmap
  tcr/13_tcr_results/       clonotypes, repertoire diversity, shared and expanded clones
```

plus a plot at every step to check the result is sane.

---

## The pages, in order

| # | Page | What it covers |
|---|------|----------------|
| 1 | [Getting started](docs/01-getting-started.md) | accounts, logging in, placeholders, the one HPC rule |
| 2 | [Installing the software](docs/02-installation.md) | the conda environments, CellChat's extras, verifying the install |
| 3 | [Bringing in your data](docs/03-your-data.md) | Trailmaker or split-pipe output, the sample sheet, sublibraries, TCR data |
| 4 | [Test small before you test big](docs/04-test-small-first.md) | public and synthetic test data, the end-to-end test |
| 5 | [The core pipeline](docs/05-core-pipeline.md) | load → QC → normalize → cluster → annotate, and what to check at each step |
| 6 | [Batch correction](docs/06-batch-correction.md) | Harmony, and choosing the right batch column |
| 7 | [Cell-cell communication](docs/07-cell-communication.md) | CellChat |
| 8 | [Trajectories and fate mapping](docs/08-trajectory.md) | PAGA, CytoTRACE, CellRank |
| 9 | [Differential expression](docs/09-differential-expression.md) | pseudobulk DESeq2, matching the design to your experiment; optional TF activity |
| 10 | [T-cell receptors](docs/10-tcr.md) | clonotypes, repertoire diversity, shared and expanded clones |
| 11 | [Running at scale with SLURM](docs/11-running-on-slurm.md) | the job scripts, chaining, watching jobs, memory |
| 12 | [When it goes wrong](docs/12-troubleshooting.md) | real errors, what they mean, what to do |

**Learning the field?** [`teaching/`](teaching/README.md) explains the
concepts behind each stage and has exercises.
**Want the reasoning and the evidence?** [`docs/reference/`](docs/reference/README.md)
has the design decisions and a validation report per step, each checked
against an alternative method on real or public data.

---

## The short version

Once installed (page 2), a split-pipe run end to end is four jobs:

```bash
cd /blue/<GROUP>/<USER>/scRNA-seq-pipeline
ACC="--account=<GROUP> --qos=<your-qos>"; OUT=results/my_run

SAMPLE_SHEET=samples.csv sbatch $ACC pipeline/run_pipeline.slurm \
    <combined>/all-sample/DGE_filtered $OUT all-sample           # load, QC, normalize
sbatch $ACC pipeline/run_downstream.slurm $OUT sublibrary        # Harmony, annotate, CellChat
COVARIATES=patient sbatch $ACC pipeline/run_pseudobulk_de.slurm \
    $OUT timepoint 6h pre                                         # DE, paired by patient
sbatch $ACC pipeline/run_tcr.slurm $OUT/tcr $OUT/04b_annotated.h5ad \
    <tcr_combined> sample patient                                 # T-cell receptors
```

Each job needs the previous one's output: wait for it to finish, or chain them
(page 11). Anything in `<angle brackets>` is a placeholder for your own value;
page 1 explains each one, page 11 every argument.

---

## What has been validated, and how

| Steps | Status |
|---|---|
| load, QC, normalize, Harmony, cluster, annotate, CellChat | **validated on real Parse data at full scale**: Parse's public 1M PBMC dataset, 987,585 cells → 907,696 after QC → 694,803 confidently annotated → 291 significant interactions across 19 pathways ([report](docs/reference/real-data-validation.md)) |
| differential expression (incl. paired designs), T-cell receptors | **tested on synthetic data in split-pipe's exact output format**, with a planted truth: 19/19 checks pass ([page 4](docs/04-test-small-first.md)); not yet on real data |
| trajectories | built and checked against CellRank's documentation; not yet run ([design](docs/reference/trajectory-design.md)) |

Each validated step is also compared against an alternative method — fixed vs.
adaptive QC thresholds, Leiden vs. Louvain, Wilcoxon vs. t-test, CellChat's
permutation test vs. naive ranking — with the numbers in
[`docs/reference/`](docs/reference/README.md).

---

## What is in this repository

```
scripts/          the pipeline (Python; R for CellChat, DESeq2 and Trailmaker conversion)
  00_load_parse_output.py      split-pipe/Trailmaker output -> .h5ad, sample sheet, sublibrary
  00a_check_cell_calls.py      optional: are split-pipe's per-sample cell calls reliable? (page 3)
  00b_seurat_to_dge_format.R   Trailmaker Seurat export -> flat files
  01_qc_filter.py              adaptive QC thresholds, doublets per sample
  02_normalize_pca.py          normalize, highly variable genes, PCA
  02b_harmony_integration.py   batch correction
  03_cluster.py                Leiden clustering at several resolutions, UMAP
  04b_better_annotation.py     confidence-gated cell-type labels (human, mouse, dog)
  04c_celltypist_annotation.py optional: reference-based labels (CellTypist, human immune)
  05_refine_annotation.py, 05b_apply_label_map.py        optional: per-lineage refinement, doublet flags, curated label maps
  05c_score_gene_programs.py   optional: published gene programs per cell / per sample (T-cell set built in)
  06_export_for_cellchat.py, 07_cellchat_analysis.R      cell-cell communication
  08_trajectory_paga.py, 09_cellrank_fate_mapping.py     trajectories, fate probabilities
  10_pseudobulk_prep.py, 11_pseudobulk_deseq2.R          differential expression
  14a_fetch_collectri.py, 14_tf_activity_ulm.py          optional: transcription-factor activity, human or mouse (page 9)
  15_gsea_preranked.py         optional: preranked GSEA on the DE results (page 9)
  12_tcr_attach.py, 13_tcr_analysis.py                   T-cell receptor clonotypes
pipeline/         SLURM job scripts that chain the steps (page 11)
  split-pipe/     what this pipeline expects from split-pipe's output (no split-pipe code)
  download/       fetch and verify a sequencing core's data release (page 11)
configs/          example label maps for 05b
envs/             conda environments and their installers (page 2)
tests/            synthetic test data and the end-to-end test (page 4)
docs/             the twelve pages above; docs/reference/ for design and validation
teaching/         concepts behind each stage, and exercises
```

---

## About split-pipe

split-pipe is Parse Biosciences' licensed software, available to Evercode
customers from Parse
([Pipeline Download Request Form](https://www.parsebiosciences.com/pipeline-download/)).
**No split-pipe code, installer or run script is distributed here.** Install
and run it following Parse's documentation;
[`pipeline/split-pipe/README.md`](pipeline/split-pipe/README.md) describes what
this pipeline expects from its output.

## Method credits

The analysis is built on **scanpy** (Wolf et al., *Genome Biology* 2018),
**Scrublet** (Wolock et al., *Cell Systems* 2019), **Harmony** (Korsunsky et
al., *Nature Methods* 2019), **CellChat** (Jin et al., *Nature Communications*
2021), **DESeq2** (Love et al., *Genome Biology* 2014), **PAGA** (Wolf et al.,
*Genome Biology* 2019), **CellRank 2** (Weiler et al., *Nature Methods* 2024)
and **scirpy** (Sturm et al., *Bioinformatics* 2020). Parse Evercode chemistry
is based on SPLiT-seq (Rosenberg et al., *Science* 2018). Full references with
DOIs: [design decisions](docs/reference/design-decisions.md#references). This
repository is the wrapper that runs them on Parse data at scale and ties the
pieces together.

## Scope

**Covered:** Parse Evercode WT and TCR data, from Parse's demultiplexed output
through annotation, cell-cell communication, trajectories, differential
expression and clonotypes. **Not covered:** demultiplexing itself (Parse's
software), RNA velocity (standard Parse output lacks spliced/unspliced counts),
reference-based annotation, and the biological interpretation of your
results — a scientific call only you can make.

## License

MIT — see [LICENSE](LICENSE). This covers the code and documentation in this
repository only; split-pipe and the other tools it runs keep their own licenses.

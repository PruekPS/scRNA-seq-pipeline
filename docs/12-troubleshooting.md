# 12. When it goes wrong

Every problem below actually happened while building this pipeline. Find the
symptom, read what it means, do what it says.

---

## First: read the log

```bash
tail -50 logs/<name>-<jobid>.out
sacct -j <jobid> -X --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

A job that dies **within seconds with almost no output** is almost always a
script problem — wrong path, missing argument, wrong folder — not a compute
problem. The log nearly always has the Python or R error spelled out.

---

## Jobs

| Symptom | What it means | What to do |
|---|---|---|
| `python: can't open file '.../scripts/...'` | you ran `sbatch` from the wrong folder | `cd` into the repository, resubmit (page 11) |
| state `OUT_OF_MEMORY`, or "Killed" in the log | the job used more memory than requested | check `MaxRSS`; then see the next table before asking for more |
| `PENDING` for a long time | waiting for a free machine | normal; do not resubmit. Smaller memory/time requests start sooner |

## Running out of memory

More memory is not always the fix. Two real cases where it was not:

| Case | Real cause | Real fix |
|---|---|---|
| QC died at 400 GB on 987,585 cells | doublet detection ran on all 24 samples at once | per sample (`--batch-key sample`, the default) — correct *and* far less memory |
| CellChat died with 8 workers | each worker is a separate R process with its own copy of the data, so 8 workers ≈ 8× the memory | `--workers=1` (the default); it finished in ~1.5 h |

If a parallel job runs out of memory where a smaller one did not, try fewer
workers before asking for more memory.

---

## Installs

| Symptom | What is going on | What to do |
|---|---|---|
| `dependencies X are not available` installing an R package from GitHub | the installer looks only on CRAN; X is a Bioconductor package | install X with conda first (already done in `envs/install_r_extras.sh`) |
| a C/C++ compile fails with a linker error | the cluster's compilers lack a library the source build needs | use the precompiled conda-forge / bioconda build instead |
| `namespace X ... >= Y is required` | conda's newest build of X is older than needed | check `conda search X`; install that one from CRAN instead (done for NMF) |
| `No module named 'pkg_resources'` reading TCR tables | setuptools 81+ removed it; scirpy's `airr` still needs it | `setuptools<81` (pinned in `environment.yml`, or run `envs/install_tcr.slurm`) |
| `cellrank` will not install into `scrna` | CellRank needs Python ≥3.12, `scrna` is 3.11 | use the separate `scrna-traj` environment (page 2) |
| an install "succeeded" but the script fails to import | a clean exit is not proof it worked | always run the verify step on page 2 |

---

## Pipeline steps

| Symptom | Meaning | What to do |
|---|---|---|
| loader stops: `Samples with cells but no --sample-sheet row` | sheet names differ from split-pipe's sample names | copy the names exactly (case, underscores) |
| some samples have only a few hundred cells, each with tens of thousands of transcripts | split-pipe's per-sample cell cutoff landed above the main population | run `00a_check_cell_calls.py`; if it flags samples, load `DGE_unfiltered` with `--min-counts`/`--min-genes` (page 3) |
| mouse data: QC removes almost nothing for mitochondria, and `0 genes matched mt-prefix 'MT-'` | species left at the human default | `--species mouse` (step 01) or `SPECIES=mouse` for `run_pipeline.slurm` |
| tissue cells (lung, tumour, stroma) labelled "Platelet" or "HSPC" by 04b | 04b's marker sets are blood-focused | read `top_markers.csv` and label those clusters yourself, or extend the marker dictionary |
| mouse TF activity: few TFs scored, or results odd | human CollecTRI used on mouse genes | fetch the mouse network: `14a_fetch_collectri.py ... --organism mouse` (page 9) |
| loader warns `DGE_unfiltered holds every barcode ... no --min-counts` | unfiltered folder loaded without a cutoff, so empty barcodes come in too | add `--min-counts N --min-genes M` (page 3) |
| about a third of genes have IDs instead of names | those genes have a blank `gene_name` in Parse's output | expected; the loader falls back to `gene_id` |
| Harmony step: `Value passed for key 'X_pca_harmony' is of incorrect shape` | harmonypy 2.0.0 is installed | the environment pins 0.0.10; reinstall it |
| many `Ambiguous (A / B)` labels | the marker sets cannot separate those clusters confidently | look at `top_markers.csv` for those clusters; for non-blood tissue, extend the marker dictionary in `04b_better_annotation.py` |
| step 12 stops: `Only N% of TCR cell ids match` | WT and TCR cell names differ (usually the sublibrary tag) | follow the option the message names (`--tcr-suffix-format`, `--strip-tcr-suffix`) |
| step 12 stops: one sublibrary `matched only N%` | TCR was combined in a different sublibrary order than WT | re-combine TCR in the WT order |
| step 12 warns few TCR cells are in T-cell clusters | wrong cell matching, doublets, or annotation without T-cell labels | do not use the TCR results until explained |
| step 11: `SKIP '<cell type>': design ... is not estimable` | the covariate lines up exactly with the condition for that cell type | that cell type cannot be tested with this design; drop or change `--covariates` |
| step 11: `SKIP ... need >= 2 each` | fewer than 2 samples in a condition for that cell type | expected with few samples; not an error |

---

## Something in the results looks wrong

Stop and check before building on it. A step that finishes but produces
nonsense is worse than one that crashes, because it looks like a result. Good
first questions: does the cell count still make sense? Do the plots for that
step look like the ones in [docs/reference/](reference/README.md)? Does the
result hold in every sample, or only when pooled?

---

That is the end of the guide. For the concepts behind each step, and
exercises, see [teaching/](../teaching/README.md). For the design decisions and
validation reports, see [docs/reference/](reference/README.md).

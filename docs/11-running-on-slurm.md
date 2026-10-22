# 11. Running at scale with SLURM

Once every step works on small data (pages 4–5), run the real dataset as
batch jobs. The scripts in `pipeline/` chain several steps into one job, all
writing into one output folder.

---

## Every time you log in

```bash
cd /blue/<GROUP>/<USER>/scRNA-seq-pipeline
mkdir -p logs
```

> **Always `cd` into the repository before `sbatch`.** The job scripts find
> `scripts/` relative to where you ran `sbatch` from (`$SLURM_SUBMIT_DIR`), not
> where the job file lives. Submit from your home folder and the job fails in
> seconds with "no such file" — this happened while building the pipeline.

---

## The job scripts, in run order

| Job | Runs | Arguments | Resources |
|---|---|---|---|
| `run_pipeline.slurm` | 00 load → 01 QC → 02 PCA → [02b Harmony] → 03 cluster → 04b annotation → [04c CellTypist] | `DGE_DIR OUT_DIR [sample_name] [mad] [sample_sheet.csv]` + optional settings below | 16 CPU, 128 GB, 24 h — ~150k cells in ~30 min; use `--mem=400gb` for ~1M cells |
| `run_downstream.slurm` | 02b Harmony → 03 cluster → 04b annotation → 06/07 CellChat | `OUT_DIR [batch_key] [resolution] [skip_harmony 0/1] [species]` | 16 CPU, 400 GB, 24 h — proven on ~1M cells |
| `run_pseudobulk_de.slurm` | 10 → 11 differential expression | `OUT_DIR CONDITION_COL A B [annotation_col] [sample_col]` + optional `COVARIATES=` | 4 CPU, 64 GB, 4 h |
| `run_trajectory.slurm` | 08 → 09 trajectories (uses `scrna-traj`) | `OUT_DIR [annotation_col] ROOT_CLUSTER [n_states]` | 8 CPU, 200 GB, 12 h — not yet sized on real data |
| `run_tcr.slurm` | 12 → 13 T-cell receptors | `OUT_DIR WT_H5AD AIRR_ROOT [sample_col] [group_col] [annotation_col]` | 4 CPU, 64 GB, 4 h — not yet sized on real data |

The resources are set inside each file; ~1M cells needed 400 GB.
A run a tenth that size needs much less — see "Memory" below. 16 CPUs can
exceed a group's regular QoS limit; submit to the burst QoS (`<group>-b`) if
the job is rejected with `QOSGrpCpuLimit`.

**`run_pipeline.slurm` settings** (environment variables at submit time, e.g.
`sbatch --export=ALL,SPECIES=mouse,HARMONY_KEY=sublibrary ...`):

| Setting | Effect |
|---|---|
| `SPECIES=mouse` | mouse mitochondrial prefix (`mt-`) in QC and mouse markers in 04b. **Set it for mouse data**, otherwise the mitochondrial filter finds no genes |
| `SAMPLE_SHEET=…`, `SAMPLE_SHEET_KEY=…` | join the study design (page 3); the key column defaults to `sample` |
| `MIN_COUNTS=500 MIN_GENES=300` | one cell cutoff for every sample, with a `DGE_unfiltered` folder (page 3) |
| `EXCLUDE_SUBLIBS="s1"` | drop a sublibrary with a library-level artifact |
| `HARMONY_KEY=sublibrary` | Harmony (02b) before clustering |
| `RESOLUTION=0.6` | the Leiden resolution 04b and 04c label |
| `CELLTYPIST_MODEL=…/Immune_All_Low.pkl` (+ `CELLTYPIST_MODEL2=…`) | run 04c CellTypist after 04b (human immune models) |

`run_downstream.slurm` starts again from `02_pca.h5ad` (Harmony → cluster → 04b
→ CellChat) and replaces `04b_annotated.h5ad`; skip `HARMONY_KEY` in the first
job if you plan to run it.

---

## A full run, end to end

```bash
ACC="--account=<GROUP> --qos=<your-qos>"
OUT=/blue/<GROUP>/<USER>/results/my_run

# 1. load (with sample sheet), QC, normalize, first clustering
SAMPLE_SHEET=samples.csv sbatch $ACC --output=logs/pipeline-%j.out \
  pipeline/run_pipeline.slurm <combined>/all-sample/DGE_filtered $OUT all-sample

# 2. Harmony by sublibrary, re-cluster, confident annotation, CellChat
sbatch $ACC --output=logs/downstream-%j.out \
  pipeline/run_downstream.slurm $OUT sublibrary 0.6 0 human

# 3. differential expression: 6h vs pre, paired by patient
COVARIATES=patient sbatch $ACC --output=logs/de-%j.out \
  pipeline/run_pseudobulk_de.slurm $OUT timepoint 6h pre

# 4. T-cell receptors
sbatch $ACC --output=logs/tcr-%j.out \
  pipeline/run_tcr.slurm $OUT/tcr $OUT/04b_annotated.h5ad <tcr_combined> sample patient
```

Steps 3 and 4 need step 2's `04b_annotated.h5ad`. Extra options for step 12
go in `TCR_ATTACH_EXTRA` (for example `TCR_ATTACH_EXTRA="--strip-tcr-suffix"`).

---

## Downloading data from a sequencing core

`pipeline/download/fetch_core_data.slurm` pulls a core's data release (an HTTPS
folder behind a password, e.g. the UMN Genomics Center) and checks it. It
compares every file with the core's checksum list (skipping the line where the
list checksums itself), tests each FASTQ for truncation, and writes its own
`MANIFEST.md5` for later checks:

```bash
sbatch --account=<GROUP> --qos=<GROUP> --output=<DEST>/logs/fetch-%j.out \
  pipeline/download/fetch_core_data.slurm <URL> <DEST>
```

The password never goes on the command line or into the repository. Before
submitting, put it in a private `~/.core_download.wgetrc` (`chmod 600`), as the
script's header shows, and delete that file once the download is verified. The
job ends with `all checks passed`. Re-running it is safe: it fetches only
missing or changed files.

---

## Run the next job only if the previous one succeeded

```bash
J1=$(SAMPLE_SHEET=samples.csv sbatch --parsable $ACC --output=logs/pipeline-%j.out \
     pipeline/run_pipeline.slurm <combined>/all-sample/DGE_filtered $OUT all-sample)
sbatch $ACC --dependency=afterok:$J1 --output=logs/downstream-%j.out \
     pipeline/run_downstream.slurm $OUT sublibrary
```

The second job waits as `PENDING (Dependency)` and never starts if the first
one fails.

---

## Watching jobs

```bash
squeue -u $USER                                             # running and waiting jobs
tail -f logs/pipeline-<jobid>.out                           # live output
sacct -j <jobid> -X --format=JobID,State,ExitCode,Elapsed,MaxRSS   # after it ends
```

| `squeue` state | Means |
|---|---|
| `PD` | pending — waiting for a free machine; normal, do not resubmit |
| `R` | running |
| (gone) | finished one way or the other — check `sacct` |

---

## Memory

`MaxRSS` in `sacct` is the peak memory the job actually used. Compare it with
what you asked for:

- far below → ask for less next time; smaller requests start sooner
- the job died with `OUT_OF_MEMORY` → see page 12 before simply asking for more

Check what nodes actually have before guessing a bigger number:

```bash
sinfo -N -o "%N %P %m" | sort -u -k3 -n | tail
```

---

Next: [12. When it goes wrong](12-troubleshooting.md)

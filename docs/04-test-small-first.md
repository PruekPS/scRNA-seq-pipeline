# 4. Test small before you test big

Before your full dataset, run the chain on a few thousand cells. This is not
optional caution; it is the most useful habit in this guide.

> **Why it is not optional.** While building this pipeline, the QC step was first
> run directly on ~1 million cells, with no small test. It ran for 24 hours of
> allocated compute before dying of out-of-memory, a bug a 5-minute test on 50,000
> cells would have caught. Scale up only after the small run succeeds *and* its
> output looks sane.

> **Before you start**
> - The `scrna` environment is installed ([page 2](02-installation.md)).
> - You are in the code folder (`cd /blue/<GROUP>/<USER>/scRNA-seq-pipeline`).
> - `$PROJECT` and `$OUT` are set (page 1).

---

## Step 1: work on a compute node

Small tests are quick, but they are still computation, so not for the login node.
Ask SLURM for an interactive session:

```bash
srun --account=<GROUP> --qos=<your-qos> --cpus-per-task=4 --mem=16gb \
     --time=02:00:00 --pty bash -i
```

**You should see:** after a short wait, the prompt shows a compute node's name
instead of `login...`. A new session starts fresh, so set everything again:

```bash
cd /blue/<GROUP>/<USER>/scRNA-seq-pipeline
PROJECT=/blue/<GROUP>/<USER>/my_project
OUT=$PROJECT/results
module load conda
conda activate scrna
```

**You should see:** `(scrna)` at the start of the prompt. Type `exit` when you are
done, to give the machine back.

---

## Three kinds of test data

| Data | Size | What it tests | What it cannot test |
|---|---|---|---|
| Synthetic split-pipe output (`tests/synthetic_splitpipe.py`) | 2,000 cells, 4 samples, 3 sublibraries, WT + TCR | the whole chain on split-pipe's exact output layout, against a **planted truth** | real biology; real file quirks |
| PBMC3k (public, 10x Genomics) | 2,700 cells, 1 sample | every step runs; results can be compared with published ones | anything Parse-specific |
| A subset of your own data | 20–50k cells | memory, runtime, real quirks | |

---

## Step 2: the end-to-end test (do this one first)

From the code folder, in your interactive session. The test makes its own folders:

```bash
bash tests/integration_smoke.sh tests/data/integration
```

It builds synthetic data in split-pipe's output layout, with known cell types, known
T-cell clones and known differentially expressed genes. It runs the real scripts, then
**checks every result against the planted truth**. It takes a few minutes and ends with a
list like this:

```text
PASS  00a flags the planted unstable cutoff: flagged [('s1', 'P2_6h'), ('s2', 'P2_6h'), ('s3', 'P2_6h')]
PASS  00c flags the planted duplicate library: ...
PASS  04b annotation matches planted types: 99.7% of 933 confidently labeled cells
PASS  12 mis-ordered sublibraries refused: no output written
PASS  11 significant genes are all planted 6h genes: false hits: []
...
32/32 checks passed
```

**It worked if:** the last line says every check passed. Run it after installing, and
again after changing any script. Any `FAIL` line says what went wrong. About half of its
cells end up labelled *Ambiguous*: its fake CD4 and CD8 T cells differ in only two genes,
and the annotator rightly declines to guess. That is a limit of the synthetic data.

To see what the figures look like, open `tests/data/integration/run/qc_plots/qc_violin.png`
(download it, or use your cluster's file browser).

---

## Step 3: a real public dataset (PBMC3k, optional)

Everything in `tests/data/` is ignored by git, so it can never be published by accident.
Create the folder, and download the data (this needs internet, and also creates a `data/`
folder in the code folder):

```bash
mkdir -p tests/data
python -c "import scanpy as sc; sc.datasets.pbmc3k().write_h5ad('tests/data/pbmc3k_raw.h5ad')"
python scripts/01_qc_filter.py --input tests/data/pbmc3k_raw.h5ad \
  --output tests/data/pbmc3k_qc.h5ad --plot-dir tests/data/pbmc3k_qc_plots
```

**You should see:** `Final: ... cells`, about 2,600 of the 2,700. Continue through the steps
on page 5, using `tests/data/` for the files.

---

## When is a small run "sane"?

Before scaling up, look at the plots each step writes ([page 5](05-core-pipeline.md) says
what to look for at each step). A run that finishes but produces nonsense is worse than one
that crashes, because it looks like a result.

---

Next: [5. The core pipeline, step by step](05-core-pipeline.md)

# 4. Test small before you test big

Before your full dataset, run the chain on a few thousand cells. This is not
optional caution; it is the most useful habit in this guide.

> **Why it is not optional.** While building this pipeline, the QC step was
> first run directly on ~1 million cells, with no small test. It ran for 24
> hours of allocated compute before dying of out-of-memory — a bug a 5-minute
> test on 50,000 cells would have caught. Scale up only after the small run
> succeeds *and* its output looks sane.

---

## Work in an interactive session on a compute node

Small tests are quick, but they are still computation — not for the login
node. Ask SLURM for an interactive shell:

```bash
srun --account=<GROUP> --qos=<your-qos> --cpus-per-task=4 --mem=16gb \
     --time=02:00:00 --pty bash -i
cd /blue/<GROUP>/<USER>/scRNA-seq-pipeline
module load conda && conda activate scrna
```

Your prompt now shows a compute node's name instead of the login node's.

---

## Three kinds of test data

| Data | Size | What it tests | What it cannot test |
|---|---|---|---|
| PBMC3k (public, 10x Genomics) | 2,700 cells, 1 sample | every step runs; the standard tutorial dataset, so results can be compared with published ones | anything Parse-specific (split-pool layout, many samples) |
| Synthetic split-pipe output (`tests/synthetic_splitpipe.py`) | 2,000 cells, 4 samples, 3 sublibraries, WT + TCR | the whole chain on split-pipe's exact output layout, against a **planted truth** | real biology; real file quirks |
| A subset of your own data | 20–50k cells | memory, runtime, real quirks | — |

### PBMC3k

```bash
python -c "import scanpy as sc; sc.datasets.pbmc3k().write_h5ad('tests/data/pbmc3k_raw.h5ad')"
python scripts/01_qc_filter.py --input tests/data/pbmc3k_raw.h5ad \
  --output tests/data/pbmc3k_qc.h5ad --plot-dir tests/data/pbmc3k_qc_plots
```

and continue through the steps on page 5. Everything in `tests/data/` is
ignored by git.

### The end-to-end test

```bash
bash tests/integration_smoke.sh tests/data/integration
```

It builds synthetic data in split-pipe's output layout, with known cell types,
known T-cell clones and known differentially expressed genes. It then runs the
real scripts — load with a sample sheet, QC, Harmony by sublibrary, clustering,
annotation, TCR, paired DE — and **checks every result against the planted
truth**, ending with a list like:

```
PASS  00a flags the planted unstable cutoff: flagged [('s1', 'P2_6h'), ('s2', 'P2_6h'), ('s3', 'P2_6h')]
PASS  04b annotation matches planted types: 99.7% of 933 confidently labeled cells
PASS  12 combined: match rate: 98.1% of 903 TCR cells (QC removed 38 WT cells)
PASS  12 mis-ordered sublibraries refused: no output written
PASS  13 three expanded P1 clones found once each and shared pre/6h
PASS  11 significant genes are all planted 6h genes: false hits: []
...
17/17 checks passed
```

It takes about two minutes. Run it after installing, and again after changing
any script. (About half of its cells end up labelled *Ambiguous*: its fake CD4
and CD8 T cells differ in only two genes, and the annotator rightly declines to
guess. That is a limit of the synthetic data, not of the annotator.)

---

## When is a small run "sane"?

Before scaling up, look at the plots each step writes (page 5 says what to
look for at each step). A run that finishes but produces nonsense is worse than
one that crashes, because it looks like a result.

---

Next: [5. The core pipeline, step by step](05-core-pipeline.md)

# 2. Installing the software

One-time setup: get the code, then build the software environments. Every
install below runs as a **job on a compute node**, never on the login node
(page 1).

---

## Get the code

```bash
ssh <your-username>@<cluster-login-host>
cd /blue/<GROUP>/<USER>/
git clone https://github.com/PruekPS/scRNA-seq-pipeline.git
cd scRNA-seq-pipeline
```

Check you are in the right folder — this should print `OK`:

```bash
test -f scripts/01_qc_filter.py && echo "OK, right folder" || echo "WRONG folder"
```

### What is in the repository

| Path | What it holds |
|---|---|
| `scripts/` | the pipeline steps, numbered in run order `00_` to `13_`; each is a standalone command with its usage written at the top of the file |
| `pipeline/` | SLURM job scripts that chain several steps into one job (page 11) |
| `envs/` | the software environments and their installers (this page) |
| `tests/` | synthetic test data generators and the end-to-end test (page 4) |
| `docs/` | this guide; `docs/reference/` holds design decisions and validation reports |
| `teaching/` | the concepts behind each stage, and exercises |

> Read a script's top docstring (the text in triple quotes) before running it.
> Each one documents its inputs, outputs, and a real usage example.

---

## Environment 1: `scrna` — the main pipeline

Python 3.11 with scanpy, Harmony, scirpy (TCR), plus R with Seurat and DESeq2.
Used by every step except 08–09.

```bash
sbatch --account=<GROUP> --qos=<your-qos> --job-name=env-scrna \
  --cpus-per-task=4 --mem=32gb --time=02:00:00 \
  --output=env-scrna-%j.out \
  --wrap="module load conda && mamba env create -f envs/environment.yml"
```

**Why versions are pinned.** An earlier, unpinned version of this file
installed a working `scanpy` on one machine and a five-year-old, broken one on
the cluster — same file, different result. The install finished with no
error; it only failed when a script tried to import scanpy. Two pins in
`environment.yml` exist because of real failures:

| Pin | Why |
|---|---|
| `harmonypy=0.0.10` | version 2.0.0 returns a malformed array on large data (failed on a 907,696-cell run) |
| `setuptools<81` | scirpy's dependency `airr` still needs `pkg_resources`, which setuptools 81 removed |

**An install that exits cleanly is not proof it worked** — always run the
verify step below.

Two optional steps need pip packages, listed at the end of `environment.yml`:
CellTypist for 04c and gseapy for 15. They were added to an existing
environment without changing anything else
(`python -m pip install celltypist==1.6.3`;
`python -m pip install --no-deps gseapy==1.1.3`). Step 14 (TF activity) needs
nothing extra: it re-implements decoupler's method in numpy, because installing
decoupler would downgrade this environment's numpy.

## CellChat's extra pieces

CellChat (page 7) needs packages conda cannot provide: Bioconductor packages,
one package whose conda build is too old, and two that exist only on GitHub.
`envs/install_r_extras.sh` installs them; every line is commented with the
real failure it fixes.

```bash
sbatch --account=<GROUP> --qos=<your-qos> --job-name=env-r-extras \
  --cpus-per-task=4 --mem=16gb --time=01:30:00 \
  --wrap="bash envs/install_r_extras.sh"
```

## Environment 2: `scrna-traj` — only for trajectories (page 8)

CellRank 2 needs Python 3.12 or newer; the validated `scrna` environment is
Python 3.11. Rather than risk the working environment, the trajectory steps
get their own:

```bash
sbatch --account=<GROUP> --qos=<your-qos> --job-name=env-traj \
  --cpus-per-task=4 --mem=32gb --time=02:00:00 \
  --wrap="module load conda && mamba env create -f envs/environment_trajectory.yml"
```

Skip this if you will not run steps 08–09.

---

## Verify

Run as a tiny job (or in an interactive session on a compute node):

```bash
sbatch --account=<GROUP> --qos=<your-qos> --cpus-per-task=1 --mem=8gb --time=00:20:00 \
  --output=env-verify-%j.out --wrap='
module load conda
conda run -n scrna python -c "import scanpy, scirpy, harmonypy; print(\"scanpy\", scanpy.__version__, \"scirpy\", scirpy.__version__)"
conda run -n scrna Rscript -e "suppressMessages({library(CellChat); library(DESeq2)}); cat(\"CellChat + DESeq2 load OK\\n\")"
conda run -n scrna-traj python -c "import cellrank; print(\"cellrank\", cellrank.__version__)"'
```

Then read `env-verify-<jobid>.out`. If `scrna-traj` was skipped, its line
fails — that is expected.

---

## Already have an older `scrna` environment?

Two installers add the newer pieces to an **existing** environment without
risking it. Each records a rollback point, solves the install as a dry run
first, and **refuses to install anything** if that would replace a package the
validated pipeline depends on (Python, R, harmonypy, scanpy, numpy, …). After
installing, it runs small real tests (a tiny DESeq2 run, a scirpy clonotype
call) and checks CellChat still loads.

```bash
# DESeq2 into scrna, and builds scrna-traj (CellRank)
sbatch --account=<GROUP> --qos=<your-qos> envs/install_trajectory_de.slurm /path/to/envs/scrna

# scirpy (+ setuptools<81) into scrna, for the TCR steps
sbatch --account=<GROUP> --qos=<your-qos> envs/install_tcr.slurm /path/to/envs/scrna
```

Find the environment path with `conda env list`.

---

Next: [3. Bringing in your data](03-your-data.md)

# 2. Installing the software

One-time setup: get the code, make your project folders, then build the software
environments. Every install runs as a **job on a compute node**, never on the
login node ([page 1](01-getting-started.md)). Allow an afternoon, mostly waiting.

> **Before you start**
> - You are logged in to the cluster.
> - You know your `<GROUP>` and `<your-qos>` (page 1).

---

## Step 1: get the code

`git clone` downloads the code into a new folder called `scRNA-seq-pipeline`, so
first move to the place you want it:

```bash
cd /blue/<GROUP>/<USER>/
git clone https://github.com/PruekPS/scRNA-seq-pipeline.git
cd scRNA-seq-pipeline
```

Check you are in the right folder:

```bash
test -f scripts/01_qc_filter.py && echo "OK, right folder" || echo "WRONG folder"
```

**You should see:** `OK, right folder`.

To get updates later: `cd` into this folder and run `git pull`.

### What is in the repository

| Path | What it holds |
|---|---|
| `scripts/` | the pipeline steps, numbered in run order; each is a standalone command with its usage written at the top of the file |
| `pipeline/` | SLURM job files that chain several steps into one job (page 11), and small helpers |
| `envs/` | the software environments and their installers (this page) |
| `tests/` | synthetic test data and the end-to-end test (page 4) |
| `docs/` | this guide; `docs/reference/` holds design decisions and validation reports |
| `teaching/` | the concepts behind each stage, and exercises |

---

## Step 2: make your project folders

This creates `logs/`, `results/`, `references/`, `figures/` and a `samples.csv` to
fill in. It takes a second and is safe to run again, because it never overwrites
anything. Give the **full path**, starting with `/`:

```bash
bash pipeline/new_project.sh /blue/<GROUP>/<USER>/my_project
```

**You should see:** `Created:` and the list of folders.

> **Why `logs/` must exist first.** SLURM writes a job's output to the file named
> in `--output`, but it will not create the folder. If the folder is missing, the
> job fails at once and leaves no log at all, which is confusing to debug.

Now set your saved values (page 1):

```bash
PROJECT=/blue/<GROUP>/<USER>/my_project
OUT=$PROJECT/results
```

---

## Step 3: tell conda to install on your large storage

Conda environments take many gigabytes. By default they go in your home folder,
which is usually small and can fill up partway through an install. Skip this if
your cluster gives a large home folder.

Create the folders:

```bash
mkdir -p /blue/<GROUP>/<USER>/.conda/envs /blue/<GROUP>/<USER>/.conda/pkgs
```

Open conda's settings file:

```bash
nano ~/.condarc
```

Type these lines, with your own values, then press Ctrl+O, Enter, Ctrl+X:

```yaml
envs_dirs:
  - /blue/<GROUP>/<USER>/.conda/envs
pkgs_dirs:
  - /blue/<GROUP>/<USER>/.conda/pkgs
```

---

## Step 4: build the main environment, `scrna`

An *environment* is a self-contained set of software at known versions. `scrna`
holds Python 3.11 with scanpy, Harmony and scirpy, plus R with Seurat and DESeq2.
Every step except 08–09 uses it.

Make sure you are in the code folder, then submit:

```bash
cd /blue/<GROUP>/<USER>/scRNA-seq-pipeline
sbatch --account=<GROUP> --qos=<your-qos> --job-name=env-scrna \
  --cpus-per-task=4 --mem=32gb --time=04:00:00 \
  --output=$PROJECT/logs/env-scrna-%j.out \
  --wrap="module load conda && mamba env create -f envs/environment.yml"
```

**You should see:** `Submitted batch job 12345`. Then watch it:

```bash
squeue -u $USER                                # is it still running?
tail -50 $PROJECT/logs/env-scrna-12345.out     # use your own job number
```

**It worked if:** the log ends with lines about `conda activate scrna` and no
error. Step 7 confirms it properly. If your cluster has no `mamba`, use
`conda env create` in the command, which is slower. If `module load conda` fails,
ask your cluster's support how they provide conda.

**Why versions are pinned.** An earlier, unpinned version of this file installed a
working `scanpy` on one machine and a five-year-old, broken one on the cluster:
same file, different result. The install finished with no error; it only failed
when a script tried to import scanpy. Two pins exist because of real failures:

| Pin | Why |
|---|---|
| `harmonypy=0.0.10` | version 2.0.0 returns a malformed array on large data (failed on a 907,696-cell run) |
| `setuptools<81` | scirpy's dependency `airr` still needs `pkg_resources`, which setuptools 81 removed |

**An install that exits cleanly is not proof it worked.** Always do step 7.

Step 14 (TF activity) needs nothing extra: it re-implements decoupler's method in
numpy, because installing decoupler would downgrade this environment's numpy.

---

## Step 5: CellChat's extra pieces (only for page 7)

CellChat needs packages conda cannot provide. `envs/install_r_extras.sh` installs
them; its comments explain the real failure each line fixes. Wait until step 4 has
finished, then, from the code folder:

```bash
sbatch --account=<GROUP> --qos=<your-qos> --job-name=env-r-extras \
  --cpus-per-task=4 --mem=16gb --time=01:30:00 \
  --output=$PROJECT/logs/env-r-extras-%j.out \
  --wrap="bash envs/install_r_extras.sh"
```

**It worked if:** the log ends with `CELLCHAT_AND_PRESTO_OK`.

---

## Step 6: the trajectory environment, `scrna-traj` (only for page 8)

CellRank 2 needs Python 3.12 or newer, while `scrna` is Python 3.11. Rather than
risk the working environment, the trajectory steps get their own:

```bash
sbatch --account=<GROUP> --qos=<your-qos> --job-name=env-traj \
  --cpus-per-task=4 --mem=32gb --time=02:00:00 \
  --output=$PROJECT/logs/env-traj-%j.out \
  --wrap="module load conda && mamba env create -f envs/environment_trajectory.yml"
```

---

## Step 7: check everything

This job checks the environments, the packages and your project folders, and
prints what to do about anything missing. It changes nothing. From the code folder:

```bash
sbatch --account=<GROUP> --qos=<your-qos> --mem=8gb --time=00:30:00 \
  --output=$PROJECT/logs/check-%j.out pipeline/check_setup.sh $PROJECT
```

Read the result once the job has finished, with your job number:

```bash
cat $PROJECT/logs/check-12345.out
```

**You should see:** lines starting `OK`, ending `All good: no problems found.`
Lines starting `note` are fine: they mean an optional piece is not installed. Each
`FIX` line says what to do. It will say `samples.csv still has the example rows`
until you do [page 3](03-your-data.md).

---

## Already have an older `scrna` environment?

Two installers add the newer pieces to an **existing** environment without risking
it. Each records a rollback point, solves the install as a dry run first, and
**refuses to install anything** if that would replace a package the validated
pipeline depends on. Find the environment's path with `conda env list`.

```bash
# DESeq2 into scrna, and builds scrna-traj (CellRank)
sbatch --account=<GROUP> --qos=<your-qos> --output=$PROJECT/logs/env-de-%j.out \
  envs/install_trajectory_de.slurm /path/to/envs/scrna

# scirpy (+ setuptools<81) into scrna, for the TCR steps
sbatch --account=<GROUP> --qos=<your-qos> --output=$PROJECT/logs/env-tcr-%j.out \
  envs/install_tcr.slurm /path/to/envs/scrna
```

Two optional steps need pip packages, listed at the end of `environment.yml`. To
add them to an existing environment without changing anything else, activate it
inside a job and run `python -m pip install celltypist==1.6.3` (for 04c) and
`python -m pip install --no-deps gseapy==1.1.3` (for 15).

---

Next: [3. Bringing in your data](03-your-data.md)

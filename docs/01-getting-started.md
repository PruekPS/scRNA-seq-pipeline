# 1. Getting started: accounts, folders, and the one HPC rule

This pipeline runs its heavy steps on an HPC cluster with SLURM (it was built
and tested on UF's [HiPerGator](https://www.rc.ufl.edu/)). This page gets you
onto the cluster and sets up the folders every later page uses. New to the
terminal? Read [page 0](00-command-line-basics.md) first.

---

## What you need before starting

- **A cluster account** and a **group allocation** from your PI or lab. Together
  they give you the two values every job needs: `--account` and `--qos`. Nothing
  in this guide can substitute for them; ask your PI or HPC admin.
- **SSH access** that works, including two-factor approval (Duo on HiPerGator).
- **Data**: your own Parse Evercode output (Trailmaker export or split-pipe
  output, [page 3](03-your-data.md)), or none yet: [page 4](04-test-small-first.md)
  shows how to practise on test data first.

---

## Placeholders

Anything in `<angle brackets>` is a value you replace with your own, removing the
brackets:

| Placeholder | What it is | How to find it |
|---|---|---|
| `<your-username>` | your cluster login name | from your account email |
| `<cluster-login-host>` | the login address | `hpg.rc.ufl.edu` on HiPerGator |
| `<GROUP>` | your SLURM account, usually your lab's name | `sacctmgr show assoc user=$USER format=account,qos` |
| `<your-qos>` | your SLURM QoS, often `<GROUP>` or `<GROUP>-b` | same command |
| `/blue/<GROUP>/<USER>` | your project storage | HiPerGator's `/blue` area; on other clusters, ask where large files go. Your home folder is usually too small |

---

## Step 1: log in

Run this in the terminal on your own computer:

```bash
ssh <your-username>@<cluster-login-host>
```

Approve the two-factor prompt. You land on a **login node**. The prompt now shows
a name like `login7`.

---

## The one rule: never compute on the login node

Login nodes are shared by everyone on the cluster at once. They are for light
work only: editing files, submitting jobs, checking the queue. Anything that uses
real CPU or memory for more than a few seconds (running a pipeline step,
installing software) goes into a job on a **compute node**, submitted with `sbatch`.

> **A real example of getting this wrong.** While building this pipeline, a
> package install with several minutes of compilation was run on the login node
> in the background instead of through `sbatch`. That slows the machine for
> everyone logged in. The fix was the same install, submitted as a small job.

---

## Step 2: choose where things go

The guide uses two folders, kept apart so that updating the code can never touch
your data:

```text
/blue/<GROUP>/<USER>/
  scRNA-seq-pipeline/      the code. You download this once on page 2.
  my_project/              one folder per project.
    samples.csv            your sample sheet (page 3)
    logs/                  what each job printed
    results/               everything the pipeline writes
    references/            downloaded networks and gene sets
    figures/               publication figures
```

You will create these on [page 2](02-installation.md), once the code is there:
it includes a small script that makes the project folders for you.

---

## Step 3: know the three values you set at each login

Saved values disappear when you log out, so **paste these each time you log in**,
with your own paths:

```bash
cd /blue/<GROUP>/<USER>/scRNA-seq-pipeline     # always run commands from here
PROJECT=/blue/<GROUP>/<USER>/my_project
OUT=$PROJECT/results
```

Every later page assumes you have done this. Each page begins with a "Before you
start" box saying which folder you must be in and what must already exist.

---

## How the rest of the guide is organised

| Pages | What they do |
|---|---|
| 2–4 | install, bring in your data, test small first |
| 5 | the core chain: load → QC → normalize → cluster → annotate |
| 6–10 | optional stages: batch correction, cell-cell communication, trajectories, differential expression, T-cell receptors |
| 11–12 | running at scale with SLURM, and what to do when things break |
| 13 | publication figures in Nature style |

The *why* behind each stage is in [`teaching/`](../teaching/README.md). The design
decisions and validation reports are in [`docs/reference/`](reference/README.md).

---

Next: [2. Installing the software](02-installation.md)

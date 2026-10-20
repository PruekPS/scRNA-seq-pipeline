# 1. Getting started — accounts, logging in, and the one HPC rule

This pipeline runs its heavy steps on an HPC cluster with SLURM (it was built
and tested on UF's [HiPerGator](https://www.rc.ufl.edu/)). Small tests also run
on a laptop. This page gets you onto the cluster and explains the placeholders
used everywhere in this guide.

---

## What you need before starting

- **A cluster account** and a **group allocation** from your PI or lab — that
  gives you the two values every job needs: `--account` and `--qos`.
  Nothing in this guide can substitute for them; ask your PI or HPC admin.
- **SSH access** that works, including two-factor approval (Duo on HiPerGator).
- **Basic command-line comfort**: `cd`, `ls`, `cat`, editing a text file.
  Page 2 covers the rest.
- **Data**: your own Parse Evercode output (Trailmaker export or split-pipe
  output, page 3) — or none yet: page 4 shows how to practise on public and
  synthetic data first.

---

## Placeholders

Anything in `<angle brackets>` is a value you replace with your own:

| Placeholder | What it is | How to find it |
|---|---|---|
| `<your-username>` | your cluster login name | from your account email |
| `<cluster-login-host>` | the login address | `hpg.rc.ufl.edu` on HiPerGator |
| `<GROUP>` / `<your-account>` | your SLURM account | `sacctmgr show assoc user=$USER format=account,qos` |
| `<your-qos>` | your SLURM QoS (often `<GROUP>` or `<GROUP>-b`) | same command |
| `/blue/<GROUP>/<USER>` | your project storage | HiPerGator's `/blue` area; other clusters differ |

---

## Logging in

```bash
ssh <your-username>@<cluster-login-host>
```

Approve the two-factor prompt. You land on a **login node** — a shared
machine, used by everyone on the cluster at once.

---

## The one rule: never compute on the login node

Login nodes are for light work only: editing files, submitting jobs, checking
the queue. Anything that uses real CPU or memory for more than a few seconds —
running a pipeline step, installing or compiling software — goes into a job on
a **compute node**, submitted with `sbatch`.

> **A real example of getting this wrong.** While building this pipeline, a
> package install with several minutes of C++/Fortran compilation was run on
> the login node in the background (`nohup`) instead of through `sbatch`. That
> slows the machine for everyone logged in. The fix was the same install,
> submitted as a small job — only *where* it ran was the mistake.

Rule of thumb: if a command takes more than a couple of seconds or noticeable
memory, it goes in an `sbatch` job — even a one-off setup step.

---

## How the rest of the guide is organised

| Pages | What they do |
|---|---|
| 2–4 | install, bring in your data, test small first |
| 5 | the core chain: load → QC → normalize → cluster → annotate |
| 6–10 | optional stages: batch correction, cell-cell communication, trajectories, differential expression, T-cell receptors |
| 11–12 | running at scale with SLURM, and what to do when things break |

The *why* behind each stage — the concepts, at more depth than these pages —
is in [`teaching/`](../teaching/README.md). The design decisions and validation
reports are in [`docs/reference/`](reference/README.md).

---

Next: [2. Installing the software](02-installation.md)

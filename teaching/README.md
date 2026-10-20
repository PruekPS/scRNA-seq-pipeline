# Teaching: single-cell RNA-seq analysis, from the reasoning up

This folder is for students learning single-cell analysis with this pipeline.
The [guide in `docs/`](../docs/01-getting-started.md) says *what to run*; this
folder explains *why each step exists*, so you can make the judgement calls
yourself — and gives exercises to check that you can.

You do not need a coding assistant or anyone driving for you. You do need a
cluster account (see [guide page 1](../docs/01-getting-started.md)).

---

## What you will be able to explain

By the end you should be able to explain, not just run:

- why QC uses **adaptive thresholds** instead of fixed cutoffs, and what that
  trades off;
- what a **doublet** is, and why it is detected per sample;
- why PCA needs **highly variable gene selection** and scaling first;
- what a **batch effect** is, when to correct it, and when correcting would
  destroy your result;
- why clustering has a **resolution** and no single "correct" number of
  clusters;
- how marker-based **annotation** computes a label, and why an honest "not
  sure" beats a confident wrong answer;
- what CellChat's **communication probability** is — and is not;
- why **pseudobulk** is the right way to compare conditions, and how to match
  the statistical **design** to the experiment;
- what a T-cell **clonotype** is and how to compare repertoires fairly;
- enough **SLURM** to submit, watch and debug your own jobs.

---

## How to use this folder

Work through it alongside the guide, one stage at a time:

| Step | Do this in `docs/` | Read this here |
|---|---|---|
| 1 | [set up](../docs/01-getting-started.md), [install](../docs/02-installation.md), [bring in data](../docs/03-your-data.md) | — |
| 2 | [test small](../docs/04-test-small-first.md) on PBMC3k | [Concepts: why test small](concepts.md#testing-small) |
| 3 | [core pipeline](../docs/05-core-pipeline.md), one step at a time | [Concepts: QC](concepts.md#quality-control) → [normalization and PCA](concepts.md#normalization-highly-variable-genes-and-pca) → [clustering](concepts.md#clustering) → [annotation](concepts.md#annotation) |
| 4 | [batch correction](../docs/06-batch-correction.md) | [Concepts: batch effects](concepts.md#batch-effects) |
| 5 | the optional stages you need ([7](../docs/07-cell-communication.md)–[10](../docs/10-tcr.md)) | [cell-cell communication](concepts.md#cell-cell-communication), [trajectories](concepts.md#trajectories), [differential expression](concepts.md#differential-expression), [T-cell receptors](concepts.md#t-cell-receptors) |
| 6 | [run at scale](../docs/11-running-on-slurm.md) | [Concepts: computing on a cluster](concepts.md#computing-on-a-cluster) |
| 7 | — | [Exercises](exercises.md) |

Look at the plots every step writes. That is where the learning happens.

---

## Further reading in this repository

The [validation reports](../docs/reference/README.md) show each step checked
against an alternative method on real or public data, with the numbers — a
good model of how to justify an analysis choice. The
[design decisions](../docs/reference/design-decisions.md) page shows what was
planned, what changed once real data arrived, and why.

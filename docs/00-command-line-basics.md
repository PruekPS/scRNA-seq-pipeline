# 0. Command-line basics: start here if you have never used a terminal

You do not need to know how to program to run this pipeline. You need to type
(or paste) commands into a terminal, change a few values to your own, and read
what comes back. This page teaches exactly that, in about 20 minutes.

Already comfortable with `cd`, `mkdir`, `nano` and `sbatch`? Skip to
[page 1](01-getting-started.md).

---

## What a terminal is

A terminal is a window where you type a command and press **Enter**. The
computer runs it and prints the result. You will use it to connect to the
cluster, a large shared computer that does the heavy calculation.

| Your computer | Open this |
|---|---|
| Mac | **Terminal** (press Cmd+Space, type `Terminal`) |
| Windows | **PowerShell** (press the Windows key, type `PowerShell`) |
| Linux | **Terminal** |

---

## The commands you will use

Try these once you are logged in to the cluster ([page 1](01-getting-started.md)).

| Command | What it does | Example |
|---|---|---|
| `pwd` | shows the folder you are in ("print working directory") | `pwd` |
| `ls` | lists what is in this folder | `ls` |
| `ls -lh` | the same, with file sizes | `ls -lh results/` |
| `cd FOLDER` | moves into a folder | `cd results` |
| `cd ..` | moves up one folder | `cd ..` |
| `mkdir -p FOLDER` | creates a folder (and any missing parent folders) | `mkdir -p my_project/logs` |
| `nano FILE` | opens a simple text editor | `nano samples.csv` |
| `cat FILE` | prints a whole small file | `cat samples.csv` |
| `head FILE` | prints the first 10 lines | `head samples.csv` |
| `tail -50 FILE` | prints the last 50 lines, the place to look in a log | `tail -50 logs/pipeline-123.out` |
| `cp A B` | copies file A to B | `cp samples.csv samples_old.csv` |

**Using nano:** type your text and move around with the arrow keys. To save,
press **Ctrl+O** then **Enter**. To leave, press **Ctrl+X**. The mouse does not
move the cursor.

**Pasting:** on a Mac use Cmd+V. In most other terminals, use right-click or
Ctrl+Shift+V.

**If something seems stuck:** press **Ctrl+C** to stop the running command.

---

## How to read the commands in this guide

### Replace anything in angle brackets

`<GROUP>` means "your own value goes here". Delete the angle brackets too.

```text
Guide:     cd /blue/<GROUP>/<USER>/
You type:  cd /blue/smithlab/jdoe/
```

### Long commands continue with a backslash

A `\` at the end of a line means "this command continues on the next line".
Paste the whole block, then press Enter once. The `\` must be the very last
character on its line, with no space after it.

```bash
python scripts/01_qc_filter.py --input in.h5ad \
  --output out.h5ad
```

### Words starting with `$` are saved values

You can save a long path under a short name, and the guide does this to avoid
retyping. Set it once:

```bash
PROJECT=/blue/smithlab/jdoe/my_project
```

There are no spaces around the `=`. Then `$PROJECT/logs` means
`/blue/smithlab/jdoe/my_project/logs`. Check a saved value with
`echo $PROJECT`.

**These are forgotten when you log out.** If a command complains about a path
such as `/logs/...` that starts oddly, the value is probably empty: set it
again.

### Lines starting with `#` are notes

Your terminal ignores everything after a `#`. You can paste them.

---

## Paths: where a file is

- A **full path** starts with `/` and works from anywhere:
  `/blue/smithlab/jdoe/my_project/samples.csv`.
- A **relative path** is measured from the folder you are in. `results/x.h5ad`
  means "inside the `results` folder here". It only works from the right folder,
  which is the most common cause of `No such file or directory`.

Avoid spaces in file and folder names: use `my_project`, not `my project`. Capital
letters matter on the cluster: `Samples.csv` and `samples.csv` are different files.

---

## Moving files between your computer and the cluster

The easiest way is often a web page: on HiPerGator, use
[Open OnDemand](https://ondemand.rc.ufl.edu/) and its Files menu to upload and
download by dragging. Ask your cluster's support whether it has something similar.

From a terminal **on your own computer** (not on the cluster):

```bash
# copy a sample sheet up to the cluster
scp samples.csv <your-username>@<cluster-login-host>:/blue/<GROUP>/<USER>/my_project/

# copy a results folder down to the current folder on your computer
scp -r <your-username>@<cluster-login-host>:/blue/<GROUP>/<USER>/my_project/results/qc_plots .
```

---

## Jobs: how work runs on a cluster

The machine you log in to is shared and is not for calculations. You describe
the work in a **job**, submit it with `sbatch`, and it waits in a queue until a
machine is free. You can log out while it runs.

| Command | What it does |
|---|---|
| `sbatch JOBFILE ...` | submits a job, and prints its job number |
| `squeue -u $USER` | lists your running and waiting jobs |
| `scancel 12345` | cancels job 12345 |
| `tail -50 logs/x-12345.out` | reads what the job has printed so far |

`PD` in the `squeue` list means waiting, `R` means running. When a job leaves the
list it has finished, one way or the other: read its log to see how it went.

---

## What the files in this pipeline are

| File | What it is |
|---|---|
| `.h5ad` | your cells and genes, plus everything added so far. Each step reads one and writes the next. Large, and not meant to be opened by hand |
| `.csv` | a plain table; opens in Excel |
| `.png`, `.tiff` | images to look at |
| `.pdf`, `.svg` | the same figures with editable text, for Illustrator or Inkscape |
| `.py`, `.R` | the scripts that do the work |
| `.slurm` | a job file: several steps chained together, for `sbatch` |

---

## When an error appears

Errors look alarming, but most of the text is detail. **Read the last few lines
first:** that is usually where the actual message is. Then look it up on
[page 12](12-troubleshooting.md).

---

Next: [1. Getting started](01-getting-started.md)

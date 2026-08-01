# split-pipe — command-line alternative to Trailmaker

Trailmaker (cloud) is the primary demultiplexing path for this project. This
folder holds the setup for **split-pipe**, Parse Biosciences' command-line
pipeline, kept ready as a HiPerGator-native alternative — same barcode
correction / alignment / quantification logic as Trailmaker's Pipeline
module, run as a batch SLURM job instead of through the web UI.

## Why keep both

- Trailmaker: no-code, good for quick looks and sharing with collaborators.
- split-pipe: data and compute never leave HiPerGator, fully scriptable,
  batchable across many sublibraries, fits this repo's step-numbered style.

## Getting the software

split-pipe is licensed, customer-gated software — it is **not** vendored in
this repo (see `.gitignore`) and cannot be fetched by an agent without your
Parse account credentials. To obtain it:

1. Request access via Parse's [Pipeline Download Request
   Form](https://www.parsebiosciences.com/pipeline-download/) (requires a
   Parse Biosciences customer account tied to your Evercode kit purchase).
2. The most recent version publicly documented in third-party HPC deployment
   notes at the time this was written is **1.5.0**; confirm the true current
   version from the download portal itself, since Parse does not publish a
   public changelog.
3. Download the tarball to your own machine or directly to HiPerGator (e.g.
   via the [OnDemand file
   browser](https://ondemand.rc.ufl.edu/pun/sys/dashboard/files/fs//blue/esayour/puttawatsuphapru/))
   into `/blue/esayour/puttawatsuphapru/software/split-pipe/`.

## Install on HiPerGator

Everything below stays under `/blue/esayour/puttawatsuphapru/` — no
system-wide installs, no shared/other-user paths.

```bash
ssh hipergator
cd /blue/esayour/puttawatsuphapru
bash scRNA-pipeline/pipeline/split-pipe/install_hipergator.sh \
     /blue/esayour/puttawatsuphapru/software/split-pipe/ParseBiosciences-Pipeline.<version>.zip
```

See [`install_hipergator.sh`](install_hipergator.sh) — it creates a dedicated
mamba env, unzips the pipeline, and pip-installs it exactly as Parse's own
install script does, all scoped to your `/blue` directory.

## Running it

[`run_split_pipe_template.slurm`](run_split_pipe_template.slurm) is a SLURM
batch template: modest CPU/RAM allocation (no GPU — this step doesn't use
one), reads raw FASTQ from a `raw/` input path and writes to an `out/` path,
both meant to sit under your `/blue/esayour/puttawatsuphapru/` project space.
Fill in the placeholders (kit chemistry, sample names, genome reference)
before submitting with `sbatch`.

## Generating input FASTQ from a raw Illumina run (bcl2fastq)

If starting from raw BCL files rather than FASTQ, Parse's guidance (see
[Generating fastq files with
bcl2fastq](https://support.parsebiosciences.com/hc/en-us/articles/360058912372-Generating-fastq-files-with-bcl2fastq),
vendor documentation) is:

- Run `bcl2fastq` with **`--no-lane-splitting`** so reads aren't split across
  lanes before split-pipe sees them.
- Do **not** set an adapter-trimming sequence — Parse explicitly warns this
  can trim into the cell-barcode region of Read 2, corrupting demultiplexing
  downstream.
- Use a `SampleSheet.csv` with the Evercode kit's documented index sequences
  (from the kit's user manual), one row per sublibrary.

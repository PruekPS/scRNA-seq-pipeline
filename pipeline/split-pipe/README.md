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
2. Confirmed current version as of this writing: **1.9.0** (zip is
   `ParseBiosciences-Pipeline.<version>.zip`; Parse does not publish a public
   changelog, so re-check the download portal for anything newer before
   installing).
3. Download the zip to your own machine or directly to HiPerGator (e.g. via
   UF Research Computing's OnDemand file browser) into your project's
   `software/split-pipe/` directory under your HiPerGator group/user storage
   path (`/blue/<group>/<user>/` — see your own local notes for the actual
   value; not published in this repo). Never commit the zip itself (already
   covered by `.gitignore`'s generic `*.zip` rule) — it's licensed software.

## Install on HiPerGator

Everything below stays under your HiPerGator group/user project root — no
system-wide installs, no shared/other-user paths.

```bash
ssh <your-hipergator-alias>
cd /blue/<group>/<user>
bash scRNA-pipeline/pipeline/split-pipe/install_hipergator.sh \
     /blue/<group>/<user>/software/split-pipe/ParseBiosciences-Pipeline.<version>.zip \
     immune   # optional: only for Evercode TCR/BCR data
```

See [`install_hipergator.sh`](install_hipergator.sh) — it creates a dedicated
conda env (`python>=3.12.8`, per split-pipe's own `setup.py`), unzips the
pipeline, runs Parse's `install_dependencies_conda.sh` into that env
(installs samtools/STAR/gcc/numba etc. — there is no `env_conda.yml` in this
version), then `pip install .`, all scoped to your project root (set via the
`SPLITPIPE_PROJECT_ROOT` environment variable — see the script). If it's the
first time this HiPerGator account has used conda's default/conda-forge
channels non-interactively, `install_dependencies_conda.sh` can fail with
`CondaToSNonInteractiveError` — the script's own comments show the one-time
`conda tos accept` fix. For BCR/TCR work (Evercode TCR/BCR kits), pass
`immune` as a second argument: the script then runs Parse's
`install_immune_dbs.sh` **before** `pip install .`, since the IgBlast
databases it writes into `splitpipe/immune_dir/` are packaged at pip-install
time — running it afterwards (as older notes here suggested) leaves the
installed package without them.

## Running it

[`run_split_pipe_template.slurm`](run_split_pipe_template.slurm) is a SLURM
batch template: modest CPU/RAM allocation (no GPU — this step doesn't use
one), reads raw FASTQ from a `raw/` input path and writes to an `out/` path,
both meant to sit under your HiPerGator group/user project space.
Fill in the placeholders (account, QoS, paths, kit chemistry, sample names,
genome reference) before submitting with `sbatch`.

**Chemistry, kit, and sample wells are all required together** (confirmed
from split-pipe 1.9.0's own docs) — `--chemistry` and `--kit` alone are not
enough to get a per-sample split; `--sample <name> <wells>` maps each of
your biological samples to its plate wells (e.g. `--sample donor1 A1-A6
--sample donor2 A7-A12`). As of 1.9.0: valid `--chemistry` values are `v1`,
`v2`, `v3`, `v4`; valid `--kit` values are `WT_mini`, `WT`, `WT_mega`,
`WT_mega_384`, `WT_penta`, `WT_penta_384`, `MAX_WT_S`, `MAX_WT_M`,
`MAX_WT_L`, `MAX_WT_L_384`. Confirm rather than trust this list blindly —
Parse updates it between releases:

```bash
conda run -n splitpipe split-pipe --chem_list
conda run -n splitpipe split-pipe --kit_list
conda run -n splitpipe split-pipe --explain   # full mode/flag reference, incl. well-range syntax
```

A real single-sublibrary call, from Parse's own worked example (multiple
FASTQ files per read = multiple sequencing lanes, run through together
rather than pre-concatenated; omit `--fq2` entirely if its filenames differ
from `--fq1` only by `R1`/`R2`):

```bash
conda run -n splitpipe split-pipe \
  --mode all --chemistry v3 --kit WT \
  --fq1 S1_L001_R1_001.fastq.gz S1_L002_R1_001.fastq.gz \
  --fq2 S1_L001_R2_001.fastq.gz S1_L002_R2_001.fastq.gz \
  --genome_dir /path/to/genomes/hg38 \
  --output_dir /path/to/analysis/S1-out \
  --sample donor1 A1-A6 --sample donor2 A7-A12
```

**Species-mixed vs. species-separate samples matters for `mkref`.** Only
build one combined genome (`--genome_name hg38 GRCm39 --fasta ... --genes
...`, both lists ordered identically) if cells from two species are
literally mixed in the same well. If different species occupy separate
wells as different samples (not this project's usual case, but relevant if
mouse and canine samples ever share a sublibrary), run each species through
its own `mkref` genome and its own `--mode all` call instead — a combined
genome increases multi-mapping loss for species-separate designs. Parse
recommends also lowering the multi-mapping cutoff for that case:
`echo "post_min_map_frac 0.01" > parfile.txt` then add `--parfile
parfile.txt` to the `--mode all` call.

**Version compatibility:** reference genomes (and other outputs) built with
a split-pipe version older than 1.9.0 are not compatible with 1.9.0 and
need to be regenerated with `--mode mkref` again before reuse.

## Building a genome reference (`--mode mkref`)

Before `--mode all` can run, `--genome_dir` needs a split-pipe-built index:

```bash
conda run -n splitpipe split-pipe --mode mkref \
  --genome_name <NAME> --fasta <genome.fa> --genes <annotation.gtf.gz> \
  --output_dir <PROJECT_ROOT>/references/<NAME>
```

Requirements worth checking before a build fails partway through (from
Parse's own "Adding Custom Sequences and Gene Annotation File Formatting"
support article — full detail in the licensed PDF kept alongside this repo,
not committed):

- `--genes` accepts GTF or GFF3 (GTF ≡ GFF2), **not** UCSC-formatted GTF —
  use Ensembl or NCBI annotations instead.
- The GTF's seqname (column 1) must exactly match the FASTA's SeqID —
  `>CP049705.1` in the FASTA requires `CP049705.1`, not a longer description,
  in the GTF.
- Column 3 needs both a `gene` and an `exon` line per feature (STAR requires
  both); column 9 needs `gene_id`. Must be tab-delimited — four spaces is a
  common accidental substitution and produces a `Problem parsing gtf
  attributes` error.
- This matters for the mouse/canine annotation this pipeline already
  supports ([`scripts/04b_better_annotation.py`](../../scripts/04b_better_annotation.py))
  — a non-Ensembl/NCBI mouse or canine GTF is the likely failure mode here,
  not anything in our own scripts.
- Adding a custom/reporter sequence (e.g. GFP) means appending it to both
  the FASTA (as a new record) and the GTF (one `gene` + one `exon` line,
  start/end matching the appended sequence's actual length), then re-running
  `mkref`. See the support article for the worked example if this comes up.

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

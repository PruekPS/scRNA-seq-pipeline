#!/bin/bash
#SBATCH --job-name=split-pipe-demux
#SBATCH --account=<GROUP>   # TODO: set to your HiPerGator group; confirm with `sacctmgr show assoc user=$USER format=account,qos`
#SBATCH --qos=<GROUP>       # TODO: set to your HiPerGator QoS — often <group> or <group>-b (burst); check the command above
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8   # TODO: pick a tier below, matching --mem
#SBATCH --mem=64gb          # TODO: pick a tier below, matching --cpus-per-task
#SBATCH --time=08:00:00
#SBATCH --output=<PROJECT_ROOT>/logs/split-pipe-%j.out   # TODO: set to your /blue/<group>/<user> path
#SBATCH --mail-type=END,FAIL
#
# Non-GPU allocation — split-pipe's alignment/quantification step is CPU +
# RAM bound, not GPU bound. Sizing per Parse's own hardware guidance for a
# SINGLE sublibrary (per-sublibrary total reads, not per-sample):
#   <=100M reads   -> 8 cpus,  64gb
#   100M-500M reads -> 16 cpus, 128gb
#   500M-1B reads   -> 24-32 cpus, 256gb
# Defaults below are the smallest tier -- raise both together to match your
# actual sublibrary's read count (`zcat *_R1*.fastq.gz | wc -l`, divide by 4).
#
# Fill in the placeholders below, then: sbatch run_split_pipe_template.slurm

set -euo pipefail

PROJECT_ROOT="<PROJECT_ROOT>"   # TODO: set to your /blue/<group>/<user> path (matches --output above)
ENV_NAME="splitpipe"

FASTQ_DIR="${PROJECT_ROOT}/raw/<RUN_NAME>"          # TODO: set input FASTQ path
OUT_DIR="${PROJECT_ROOT}/analysis/<RUN_NAME>/split-pipe-out"  # TODO: set output path
GENOME_DIR="${PROJECT_ROOT}/references/<GENOME_INDEX>"        # TODO: split-pipe genome index
SUBLIB_NAME="<SUBLIB_NAME>"                          # TODO: sublibrary name (this FASTQ pair)
CHEMISTRY="<CHEMISTRY>"                              # TODO: v1, v2, v3, or v4 as of split-pipe 1.9.0 -- confirm with `split-pipe --chem_list`
KIT="<KIT>"                                          # TODO: one of WT_mini, WT, WT_mega, WT_mega_384, WT_penta,
                                                      #       WT_penta_384, MAX_WT_S, MAX_WT_M, MAX_WT_L, MAX_WT_L_384
                                                      #       as of split-pipe 1.9.0 -- confirm with `split-pipe --kit_list`

mkdir -p "${OUT_DIR}" "${PROJECT_ROOT}/logs"

module load conda 2>/dev/null || module load mamba 2>/dev/null || true

# --sample maps each biological sample to its plate wells within this
# sublibrary -- repeat --sample for each one, e.g.:
#   --sample donor1 A1-A6 --sample donor2 A7-A12
# Well-range syntax and the full mode/flag reference:
#   conda run -n splitpipe split-pipe --explain
conda run -n "${ENV_NAME}" split-pipe \
  --mode all \
  --chemistry "${CHEMISTRY}" \
  --kit "${KIT}" \
  --genome_dir "${GENOME_DIR}" \
  --fq1 "${FASTQ_DIR}/${SUBLIB_NAME}_R1_001.fastq.gz" \
  --fq2 "${FASTQ_DIR}/${SUBLIB_NAME}_R2_001.fastq.gz" \
  --output_dir "${OUT_DIR}" \
  --sample "<SAMPLE1_NAME>" "<SAMPLE1_WELLS>" \
  --sample "<SAMPLE2_NAME>" "<SAMPLE2_WELLS>" \
  --nthreads "${SLURM_CPUS_PER_TASK}"

#!/bin/bash
#SBATCH --job-name=split-pipe-demux
#SBATCH --account=esayour   # TODO: verify — inferred from /blue/esayour/ path; confirm with `sacctmgr show assoc user=$USER format=account,qos`
#SBATCH --qos=esayour       # TODO: verify — HiPerGator QoS is often <group> or <group>-b (burst); check the command above
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48gb
#SBATCH --time=08:00:00
#SBATCH --output=/blue/esayour/puttawatsuphapru/logs/split-pipe-%j.out
#SBATCH --mail-type=END,FAIL
#
# Modest, non-GPU allocation — split-pipe's alignment/quantification step is
# CPU + RAM bound, not GPU bound. Resize --cpus-per-task / --mem to match
# your actual sublibrary size rather than requesting more by default.
#
# Fill in the placeholders below, then: sbatch run_split_pipe_template.slurm

set -euo pipefail

PROJECT_ROOT="/blue/esayour/puttawatsuphapru"
ENV_NAME="splitpipe"

FASTQ_DIR="${PROJECT_ROOT}/raw/<RUN_NAME>"          # TODO: set input FASTQ path
OUT_DIR="${PROJECT_ROOT}/analysis/<RUN_NAME>/split-pipe-out"  # TODO: set output path
GENOME_DIR="${PROJECT_ROOT}/references/<GENOME_INDEX>"        # TODO: split-pipe genome index
CHEMISTRY="<v2|v3|WT_mega>"                          # TODO: Evercode kit chemistry
SAMPLE_NAME="<SAMPLE_NAME>"                          # TODO: sublibrary/sample name

mkdir -p "${OUT_DIR}" "${PROJECT_ROOT}/logs"

module load conda 2>/dev/null || module load mamba 2>/dev/null || true

conda run -n "${ENV_NAME}" split-pipe \
  --mode all \
  --chemistry "${CHEMISTRY}" \
  --genome_dir "${GENOME_DIR}" \
  --fq1 "${FASTQ_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz" \
  --fq2 "${FASTQ_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz" \
  --output_dir "${OUT_DIR}" \
  --nthreads "${SLURM_CPUS_PER_TASK}"

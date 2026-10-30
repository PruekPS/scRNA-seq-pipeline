#!/bin/bash
# Creates the folders for a new project and a sample sheet to fill in. Safe to re-run:
# nothing that already exists is changed. Runs in seconds, fine on a login node.
#
# Usage: bash pipeline/new_project.sh /full/path/to/my_project
set -euo pipefail
PROJECT="${1:?Usage: bash pipeline/new_project.sh /full/path/to/my_project}"
case "${PROJECT}" in
  /*) ;;
  *) echo "ERROR: give the full path, starting with / (you gave: ${PROJECT})" >&2; exit 1 ;;
esac

mkdir -p "${PROJECT}"/{logs,results,references,figures}
if [[ ! -e "${PROJECT}/samples.csv" ]]; then
  cat > "${PROJECT}/samples.csv" <<'CSV'
sample,patient,condition
REPLACE_ME_sample1,patient1,control
REPLACE_ME_sample2,patient1,treated
CSV
fi
cat <<MSG
Created:
  ${PROJECT}/logs         job logs (SLURM cannot create this folder itself, so it must exist first)
  ${PROJECT}/results      pipeline output
  ${PROJECT}/references   downloaded networks and gene sets
  ${PROJECT}/figures      publication figures
  ${PROJECT}/samples.csv  sample sheet: open it with   nano ${PROJECT}/samples.csv

Copy these two lines into your terminal each time you log in:
  PROJECT=${PROJECT}
  OUT=\$PROJECT/results
MSG

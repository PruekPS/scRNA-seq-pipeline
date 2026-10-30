#!/bin/bash
# Checks that everything is in place, and says what to do about anything missing. It changes
# nothing. It imports the Python packages, which takes a minute or so, so run it as a job:
#   sbatch --account=<GROUP> --qos=<QOS> --mem=8gb --time=00:30:00 --output=check-%j.out pipeline/check_setup.sh
# then read check-<jobid>.out. Optional: pass your project folder as the first argument.
#SBATCH --job-name=scrna-check
set -uo pipefail
PROJECT="${1:-}"
problems=0
ok()   { echo "  OK    $*"; }
fix()  { echo "  FIX   $*"; problems=$((problems + 1)); }
note() { echo "  note  $*"; }

echo "== Folder"
if [[ -f scripts/01_qc_filter.py ]]; then ok "you are in the pipeline folder: $(pwd)"
else fix "not in the pipeline folder. Run: cd /path/to/scRNA-seq-pipeline"; echo "$problems problem(s)."; exit 1; fi

echo "== Conda"
module load conda 2>/dev/null || module load mamba 2>/dev/null || true
if command -v conda >/dev/null 2>&1; then ok "conda: $(conda --version)"
else fix "conda not found. See docs/02-installation.md"; echo "$problems problem(s)."; exit 1; fi
eval "$(conda shell.bash hook)"

check_env() { # name, python imports, what it is for
  if ! conda activate "$1" >/dev/null 2>&1; then
    if [[ "$3" == required ]]; then fix "environment '$1' not found. Build it: docs/02-installation.md"
    else note "environment '$1' not installed. Only needed for $3"; fi
    return
  fi
  if python -W ignore -c "import $2" >/dev/null 2>&1; then ok "environment '$1': packages load"
  else fix "environment '$1' exists but packages do not load. Run: python -c 'import $2' to see the error"; fi
}
check_env scrna "scanpy, anndata, harmonypy, scirpy, celltypist, gseapy" required
if conda activate scrna >/dev/null 2>&1; then
  version=$(python -c "import importlib.metadata as m; print(m.version('harmonypy'))" 2>/dev/null || echo missing)
  [[ "$version" == 0.0.* ]] && ok "harmonypy $version" || fix "harmonypy is $version, needs 0.0.10: conda install -n scrna harmonypy=0.0.10"
  font=$(python -W ignore -c "import sys; sys.path.insert(0,'scripts'); import figstyle; print(figstyle.apply_style())" 2>/dev/null || echo unknown)
  [[ "$font" == Arial || "$font" == Helvetica ]] && ok "figure font: $font" \
    || note "figure font: $font (Arial not installed here; see scripts/figstyle.py. Fine for checking results)"
  if command -v Rscript >/dev/null 2>&1; then
    Rscript -e 'suppressMessages(library(DESeq2)); suppressMessages(library(ComplexHeatmap))' >/dev/null 2>&1 \
      && ok "R: DESeq2 and ComplexHeatmap load" || fix "R packages for step 11 do not load: docs/02-installation.md"
    Rscript -e 'suppressMessages(library(CellChat))' >/dev/null 2>&1 \
      && ok "R: CellChat loads" || note "CellChat not installed. Only needed for page 7: envs/install_r_extras.sh"
  else
    fix "Rscript not found in the scrna environment"
  fi
  conda deactivate
fi
check_env scrna-traj "cellrank, scanpy" "trajectories (page 8)"

if [[ -n "$PROJECT" ]]; then
  echo "== Project: $PROJECT"
  for d in logs results references figures; do
    [[ -d "$PROJECT/$d" ]] && ok "folder $d exists" || fix "missing folder. Run: bash pipeline/new_project.sh $PROJECT"
  done
  if [[ -f "$PROJECT/samples.csv" ]]; then
    grep -q REPLACE_ME "$PROJECT/samples.csv" && fix "samples.csv still has the example rows. Edit it: nano $PROJECT/samples.csv" \
      || ok "samples.csv filled in"
  else
    fix "no samples.csv. Run: bash pipeline/new_project.sh $PROJECT"
  fi
fi

echo
if (( problems == 0 )); then echo "All good: no problems found."; else echo "$problems problem(s) to fix, listed above as FIX."; fi

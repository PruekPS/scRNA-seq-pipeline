#!/bin/bash
# End-to-end integration smoke test on synthetic data in split-pipe's output layout
# (tests/synthetic_splitpipe.py): WT load -> QC -> PCA -> Harmony (by sublibrary) ->
# cluster -> annotate -> TCR attach/analysis -> paired pseudobulk DE, then checks the
# results against the planted truth (tests/check_integration.py).
#
# Run from the repo root, inside the `scrna` env, ON A COMPUTE NODE, e.g.:
#   sbatch --account=<g> --qos=<q> --cpus-per-task=4 --mem=32gb --time=01:00:00 \
#     --wrap "conda activate scrna && bash tests/integration_smoke.sh tests/data/integration"
# Takes a few minutes; ~2,000 cells.
set -euo pipefail
OUT="${1:?Usage: bash tests/integration_smoke.sh OUT_DIR}"
F="${OUT}/fixture"; R="${OUT}/run"
rm -rf "${OUT}"; mkdir -p "${R}"

step() { echo; echo "######## $*"; }

step "fixture (split-pipe v1.9.0 output layout)"
python tests/synthetic_splitpipe.py "${F}"

step "00 load all-sample/DGE_filtered + sample sheet"
python scripts/00_load_parse_output.py --sample-dir "${F}/wt_comb/all-sample/DGE_filtered" all-sample \
  --sample-sheet "${F}/sample_sheet.csv" --output "${R}/00_raw.h5ad"

step "00a check split-pipe cell calls (planted unstable cutoff must be flagged)"
python scripts/00a_check_cell_calls.py --run-dir "${F}"/wt_persub/s* --output-dir "${R}/cell_call_check"

step "00 load all-sample/DGE_unfiltered with one cutoff (optional path)"
python scripts/00_load_parse_output.py --sample-dir "${F}/wt_comb/all-sample/DGE_unfiltered" all-sample \
  --min-counts 200 --min-genes 100 --sample-sheet "${F}/sample_sheet.csv" --output "${R}/00_raw_unfiltered.h5ad"

step "01 QC (per-sample Scrublet)"
python scripts/01_qc_filter.py --input "${R}/00_raw.h5ad" --output "${R}/01_filtered.h5ad" --plot-dir "${R}/qc_plots"

step "02 normalize / PCA"
python scripts/02_normalize_pca.py --input "${R}/01_filtered.h5ad" --output "${R}/02_pca.h5ad" --plot-dir "${R}/pca_plots"

step "02b Harmony, batch = sublibrary"
python scripts/02b_harmony_integration.py --input "${R}/02_pca.h5ad" --output "${R}/02b_harmony.h5ad" \
  --plot-dir "${R}/harmony_plots" --batch-key sublibrary

step "03 cluster on Harmony embedding"
python scripts/03_cluster.py --input "${R}/02b_harmony.h5ad" --output "${R}/03_clustered.h5ad" \
  --plot-dir "${R}/cluster_plots" --use-rep X_pca_harmony

step "04b annotate"
python scripts/04b_better_annotation.py --input "${R}/03_clustered.h5ad" --output "${R}/04b_annotated.h5ad" \
  --plot-dir "${R}/marker_plots" --resolution 0.6 --species human

step "12 TCR attach -- combined TCR output"
python scripts/12_tcr_attach.py --wt-input "${R}/04b_annotated.h5ad" --airr-root "${F}/tcr_comb" \
  --output "${R}/12_tcr.h5ad" --plot-dir "${R}/tcr_plots"

step "12 TCR attach -- per-sublibrary TCR outputs + --tcr-suffix-format"
AIRR_ARGS=()
for d in "${F}"/tcr_persub/s*; do
  n="${d##*/s}"
  AIRR_ARGS+=(--airr "${d}/all-sample/TCR_filtered/tcr_annotation_airr.tsv" "${n}")
done
python scripts/12_tcr_attach.py --wt-input "${R}/04b_annotated.h5ad" "${AIRR_ARGS[@]}" \
  --tcr-suffix-format '__s{label}' --output "${R}/12_tcr_persub.h5ad" --plot-dir "${R}/tcr_plots_persub"

step "12 TCR attach -- sublibraries combined in the WRONG order (must be refused)"
if python scripts/12_tcr_attach.py --wt-input "${R}/04b_annotated.h5ad" --airr-root "${F}/tcr_comb_misordered" \
     --output "${R}/12_should_not_exist.h5ad" --plot-dir "${R}/tcr_plots_bad"; then
  echo "FAIL: mis-ordered TCR was accepted"; exit 1
fi
[ ! -e "${R}/12_should_not_exist.h5ad" ] && echo "OK: mis-ordered TCR refused, nothing written"

step "13 TCR analysis"
python scripts/13_tcr_analysis.py --input "${R}/12_tcr.h5ad" --output-dir "${R}/13_tcr" --group-col patient

step "10/11 pseudobulk DE, paired by patient"
python scripts/10_pseudobulk_prep.py --input "${R}/04b_annotated.h5ad" --output-dir "${R}/pseudobulk" \
  --condition-col timepoint --covariate-cols patient
Rscript scripts/11_pseudobulk_deseq2.R --input-dir="${R}/pseudobulk" --output-dir="${R}/de" \
  --condition-a=6h --condition-b=pre --covariates=patient

step "check against planted truth"
python tests/check_integration.py "${OUT}"

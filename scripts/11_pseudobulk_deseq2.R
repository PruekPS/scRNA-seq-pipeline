#!/usr/bin/env Rscript
# Step 11: Pseudobulk differential expression with DESeq2, per cell type,
# plus volcano plots and a top-DE-gene heatmap.
#
# Reads the flat MTX/CSV files written by scripts/10_pseudobulk_prep.py
# (pseudobulk counts, one column per (sample, cell_type, condition) group)
# and runs DESeq2 (Love, Huber & Anders, Genome Biology 2014,
# doi:10.1186/s13059-014-0550-8) separately for each cell type, comparing
# two levels of whatever --condition-col scripts/10 was given. Dataset-
# agnostic by design -- --condition-a/--condition-b take any two level
# names present in pseudobulk_metadata.csv's condition column, nothing
# here assumes a specific comparison.
#
# Why per-cell-type, not one combined model: pooling all cell types into
# a single DESeq2 design would confound cell-type composition differences
# with the condition effect being tested -- the standard practice for
# scRNA-seq pseudobulk DE (see Murphy & Skene, Nat Commun 2022, cited in
# scripts/10_pseudobulk_prep.py) is one model per cell type.
#
# Cell types with fewer than 2 pseudo-samples in EITHER condition level
# are skipped with a warning (DESeq2 cannot estimate dispersion from a
# single replicate) rather than crashing the whole run.
#
# Design is per-experiment, chosen at run time with --covariates (columns
# that scripts/10 carried into pseudobulk_metadata.csv via --covariate-cols):
#   (none)                 ~ condition               unpaired groups
#   --covariates=patient   ~ patient + condition     paired before/after per patient
#   --covariates=batch     ~ batch + condition       known processing batch
# The condition term is always last, so the reported log2 fold change is the
# condition effect adjusted for the covariates. Before fitting each cell type
# the design matrix is checked for full rank; if a covariate is confounded
# with the condition (e.g. every patient appears in only one condition) that
# cell type is skipped with an explanation instead of DESeq2 erroring out.
#
# Smoke-tested 2026-10-19 on HiPerGator with tests/synthetic_tcr.py's planted
# DE (G0-G9 up 4x at "6h", G10-G39 up 3x in one patient): unpaired found 38
# significant hits, ~patient + condition found 62, all of them planted
# condition genes (0 patient genes, 0 others); a covariate confounded with the
# condition was refused per cell type. NOT yet run on real data.
#
# Usage (local):
#   Rscript 11_pseudobulk_deseq2.R --input-dir=pseudobulk/ --output-dir=de_results/ \
#     --condition-a=disease --condition-b=healthy [--covariates=patient]
#   paired shortcut (after scripts/10 --pair-col): ... --pair-col=pair
#
# --pair-col=pair fits ~ pair + [covariates +] condition like --covariates=pair, and in
# addition drops pseudo-samples with no pair and keeps, per cell type, only pairs seen in
# both conditions (a pair seen once is absorbed by its own term and adds nothing to the
# within-pair estimate); --min-replicates then counts complete pairs.
#
# Usage (HPC): see pipeline/run_pseudobulk_de.slurm.

suppressMessages({
  library(DESeq2)
  library(Matrix)
  library(ggplot2)
  library(ComplexHeatmap)
})

parse_args <- function(args) {
  out <- list(padj.threshold = 0.05, lfc.threshold = 1.0, top.n.heatmap = 30, min.replicates = 2, covariates = "")
  for (a in args) {
    if (!startsWith(a, "--") || !grepl("=", a)) next
    kv <- sub("^--", "", a)
    parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
    key <- gsub("-", ".", parts[1])
    val <- paste(parts[-1], collapse = "=")
    out[[key]] <- val
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$input.dir) || is.null(args$output.dir) || is.null(args$condition.a) || is.null(args$condition.b)) {
  stop("Usage: Rscript 11_pseudobulk_deseq2.R --input-dir=DIR --output-dir=DIR ",
       "--condition-a=LEVEL --condition-b=LEVEL [--padj-threshold=0.05] [--lfc-threshold=1.0] ",
       "[--top-n-heatmap=30] [--min-replicates=2] [--covariates=col1,col2] [--pair-col=pair]")
}
input_dir <- args$input.dir
output_dir <- args$output.dir
condition_a <- args$condition.a
condition_b <- args$condition.b
padj_threshold <- as.numeric(args$padj.threshold)
lfc_threshold <- as.numeric(args$lfc.threshold)
top_n_heatmap <- as.integer(args$top.n.heatmap)
min_replicates <- as.integer(args$min.replicates)
covariates <- trimws(strsplit(args$covariates, ",", fixed = TRUE)[[1]])
covariates <- covariates[nzchar(covariates)]
pair_col <- args$pair.col  # NULL = no pairing
if (!is.null(pair_col)) covariates <- unique(c(pair_col, covariates))
design_formula <- as.formula(paste("~", paste(c(covariates, "condition"), collapse = " + ")))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "volcano_plots"), recursive = TRUE, showWarnings = FALSE)

cat("=== Loading pseudobulk data from", input_dir, "===\n")
counts <- as.matrix(Matrix::readMM(file.path(input_dir, "pseudobulk_counts.mtx")))
genes <- read.csv(file.path(input_dir, "genes.csv"))$gene
meta <- read.csv(file.path(input_dir, "pseudobulk_metadata.csv"))

rownames(counts) <- genes
colnames(counts) <- meta$pseudo_sample_id
mode(counts) <- "integer"

meta <- meta[meta$condition %in% c(condition_a, condition_b), ]
if (!is.null(pair_col)) {
  if (!pair_col %in% names(meta)) {
    stop("--pair-col '", pair_col, "' not in pseudobulk_metadata.csv -- re-run scripts/10 with --pair-col")
  }
  unpaired <- is.na(meta[[pair_col]]) | meta[[pair_col]] == ""
  if (any(unpaired)) cat(sprintf("Paired mode: dropping %d pseudo-samples with no '%s'\n", sum(unpaired), pair_col))
  meta <- meta[!unpaired, ]
  meta[[pair_col]] <- as.character(meta[[pair_col]])
}
counts <- counts[, meta$pseudo_sample_id, drop = FALSE]
meta$condition <- factor(meta$condition, levels = c(condition_b, condition_a))  # condition_b = reference
missing_cov <- setdiff(covariates, names(meta))
if (length(missing_cov) > 0) {
  stop("--covariates not found in pseudobulk_metadata.csv: ", paste(missing_cov, collapse = ", "),
       ". Re-run scripts/10 with --covariate-cols=", paste(covariates, collapse = ","))
}
for (cv in covariates) if (is.character(meta[[cv]]) || is.logical(meta[[cv]])) meta[[cv]] <- factor(meta[[cv]])
rownames(meta) <- meta$pseudo_sample_id
cat("Design:", deparse(design_formula), "\n")

cat(sprintf("Loaded %d genes x %d pseudo-samples across %d cell types (comparing '%s' vs '%s')\n",
            nrow(counts), ncol(counts), length(unique(meta$cell_type)), condition_a, condition_b))

all_results <- list()
vst_by_celltype <- list()
n_skip_design <- 0

for (ct in sort(unique(meta$cell_type))) {
  ct_meta <- meta[meta$cell_type == ct, ]
  if (!is.null(pair_col)) {
    pr <- as.character(ct_meta[[pair_col]])
    complete <- intersect(pr[ct_meta$condition == condition_a], pr[ct_meta$condition == condition_b])
    dropped <- setdiff(unique(pr), complete)
    if (length(dropped) > 0) cat(sprintf("  '%s': dropping incomplete pair(s) %s\n", ct, paste(dropped, collapse = ", ")))
    ct_meta <- ct_meta[pr %in% complete, ]
    if (length(complete) < min_replicates) {
      cat(sprintf("SKIP '%s': only %d complete pair(s) (need >= %d)\n", ct, length(complete), min_replicates))
      next
    }
  }
  n_a <- sum(ct_meta$condition == condition_a)
  n_b <- sum(ct_meta$condition == condition_b)
  if (n_a < min_replicates || n_b < min_replicates) {
    cat(sprintf("SKIP '%s': only %d '%s' / %d '%s' pseudo-samples (need >= %d each)\n",
                ct, n_a, condition_a, n_b, condition_b, min_replicates))
    next
  }

  cat(sprintf("=== DESeq2: %s (%d vs %d pseudo-samples) ===\n", ct, n_a, n_b))
  ct_counts <- counts[, ct_meta$pseudo_sample_id, drop = FALSE]
  keep_genes <- rowSums(ct_counts) > 0
  ct_counts <- ct_counts[keep_genes, , drop = FALSE]

  ct_meta[] <- lapply(ct_meta, function(x) if (is.factor(x)) droplevels(x) else x)
  mm <- model.matrix(design_formula, data = ct_meta)
  if (qr(mm)$rank < ncol(mm)) {
    cat(sprintf("SKIP '%s': design %s is not estimable here -- a covariate is confounded with the condition ",
                ct, deparse(design_formula)),
        "(or has a level present in only one pseudo-sample). Check which samples survived --min-cells-per-group.\n")
    n_skip_design <- n_skip_design + 1
    next
  }
  if (nrow(mm) <= ncol(mm)) {
    cat(sprintf("SKIP '%s': %d pseudo-samples but %d model coefficients -- no residual degrees of freedom.\n",
                ct, nrow(mm), ncol(mm)))
    n_skip_design <- n_skip_design + 1
    next
  }

  dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = design_formula)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", condition_a, condition_b), alpha = padj_threshold)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$cell_type <- ct
  res_df$design <- deparse(design_formula)
  res_df <- res_df[order(res_df$padj), ]

  out_csv <- file.path(output_dir, paste0("de_", gsub("[^A-Za-z0-9]", "_", ct), ".csv"))
  write.csv(res_df, out_csv, row.names = FALSE)
  all_results[[ct]] <- res_df

  res_df$significant <- !is.na(res_df$padj) & res_df$padj < padj_threshold & abs(res_df$log2FoldChange) > lfc_threshold
  n_sig <- sum(res_df$significant, na.rm = TRUE)
  cat(sprintf("  %d significant genes (padj<%.3g, |log2FC|>%.2g) -> %s\n",
              n_sig, padj_threshold, lfc_threshold, out_csv))

  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pmax(padj, 1e-300)), color = significant)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#1f6f6b")) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey40") +
    labs(title = paste0(ct, ": ", condition_a, " vs ", condition_b),
         x = "log2 fold change", y = "-log10(adjusted p-value)") +
    theme_minimal() + theme(legend.position = "none")
  ggsave(file.path(output_dir, "volcano_plots", paste0("volcano_", gsub("[^A-Za-z0-9]", "_", ct), ".png")),
         p, width = 6, height = 5, dpi = 150)

  vsd <- tryCatch(vst(dds, blind = FALSE), error = function(e) varianceStabilizingTransformation(dds, blind = FALSE))
  vst_by_celltype[[ct]] <- assay(vsd)
}

if (length(all_results) == 0) {
  if (n_skip_design > 0) {
    stop("No cell type could be fitted: ", n_skip_design, " were skipped because the design ",
         deparse(design_formula), " is not estimable (see the SKIP lines above). Drop or change --covariates.")
  }
  stop("No cell type had enough replicates in both conditions (need >= ", min_replicates,
       " pseudo-samples per condition) -- nothing to report. Check pseudobulk_metadata.csv.")
}

cat("=== Building combined heatmap of top DE genes across all cell types ===\n")
combined_df <- do.call(rbind, all_results)
combined_df <- combined_df[order(combined_df$padj), ]
top_overall <- head(unique(combined_df$gene[!is.na(combined_df$padj)]), top_n_heatmap)

# Build a genes x pseudo-samples matrix from each cell type's VST values,
# restricted to the shared top-gene set (genes not tested in a given cell
# type -- e.g. filtered out for zero counts there -- show as NA, not zero).
heatmap_mat <- matrix(NA_real_, nrow = length(top_overall), ncol = ncol(counts),
                       dimnames = list(top_overall, colnames(counts)))
for (ct in names(vst_by_celltype)) {
  ct_vst <- vst_by_celltype[[ct]]
  ct_cols <- colnames(ct_vst)  # paired mode may have dropped incomplete pairs
  shared_genes <- intersect(top_overall, rownames(ct_vst))
  heatmap_mat[shared_genes, ct_cols] <- ct_vst[shared_genes, ct_cols, drop = FALSE]
}
# Only pseudo-samples of cell types that were actually fitted: skipped cell types have no
# VST values, and their all-NA columns used to crash the column clustering (hclust on NA).
fitted_cols <- unlist(lapply(vst_by_celltype, colnames), use.names = FALSE)
heatmap_mat <- heatmap_mat[, colnames(heatmap_mat) %in% fitted_cols, drop = FALSE]
heatmap_mat <- heatmap_mat[rowSums(!is.na(heatmap_mat)) > 0, , drop = FALSE]
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))  # z-score per gene across pseudo-samples
heatmap_mat_scaled[is.nan(heatmap_mat_scaled)] <- 0  # a gene constant across samples has no z-score; show as 0

col_ct <- meta$cell_type[match(colnames(heatmap_mat_scaled), meta$pseudo_sample_id)]
col_cond <- meta$condition[match(colnames(heatmap_mat_scaled), meta$pseudo_sample_id)]
col_anno <- HeatmapAnnotation(condition = col_cond, cell_type = col_ct)

# Genes filtered out in some cell types leave NA gaps; hierarchical clustering can't use
# them, so in that case keep columns grouped by cell type (split) and ordered by condition.
has_na <- anyNA(heatmap_mat_scaled)
if (has_na) cat("  heatmap: some genes untested in some cell types (grey); columns grouped, not clustered\n")
png(file.path(output_dir, "top_de_genes_heatmap.png"), width = 2400, height = max(1200, 40 * nrow(heatmap_mat_scaled)), res = 200)
draw(Heatmap(heatmap_mat_scaled, name = "z-score (VST)",
             top_annotation = col_anno,
             column_split = col_ct,
             column_order = order(col_ct, col_cond),
             cluster_columns = !has_na,
             cluster_rows = !has_na,
             show_column_names = FALSE,
             row_names_gp = grid::gpar(fontsize = 8),
             column_title_gp = grid::gpar(fontsize = 8),
             na_col = "grey90"))
dev.off()

write.csv(combined_df, file.path(output_dir, "de_all_celltypes_combined.csv"), row.names = FALSE)

cat("=== Done. Outputs in", output_dir, "===\n")
cat("  de_<cell_type>.csv          -- per-cell-type DESeq2 results\n")
cat("  de_all_celltypes_combined.csv -- all cell types concatenated\n")
cat("  volcano_plots/volcano_<cell_type>.png\n")
cat("  top_de_genes_heatmap.png\n")

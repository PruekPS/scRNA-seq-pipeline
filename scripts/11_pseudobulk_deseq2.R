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
# Paired design (--pair-col=pair, after scripts/10 was run with --pair-col):
# fits ~ pair + condition, so the condition effect is estimated within each
# pair (e.g. the same patient's pre and post draws) and between-pair
# baseline differences are absorbed by the pair term -- the standard DESeq2
# approach for paired samples (DESeq2 vignette, "Multi-factor designs").
# Per cell type, only complete pairs (>= 1 pseudo-sample in each condition)
# are kept; pseudo-samples with no pair are dropped; --min-replicates then
# counts complete pairs.
#
# NOT YET TEST-RUN as of this commit -- dataset-agnostic by design, so
# there is nothing to validate against until a real condition comparison
# exists. Validate against tests/data/ (with a synthetic condition column
# if needed) before trusting output on real data.
#
# Usage (local):
#   Rscript 11_pseudobulk_deseq2.R --input-dir pseudobulk/ --output-dir de_results/ \
#     --condition-a disease --condition-b healthy
#
# Usage (HPC): see pipeline/run_pseudobulk_de.slurm.

suppressMessages({
  library(DESeq2)
  library(Matrix)
  library(ggplot2)
  library(ComplexHeatmap)
})

parse_args <- function(args) {
  out <- list(padj.threshold = 0.05, lfc.threshold = 1.0, top.n.heatmap = 30, min.replicates = 2)
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
       "[--top-n-heatmap=30] [--min-replicates=2] [--pair-col=pair]")
}
pair_col <- args$pair.col  # NULL = unpaired ~condition
input_dir <- args$input.dir
output_dir <- args$output.dir
condition_a <- args$condition.a
condition_b <- args$condition.b
padj_threshold <- as.numeric(args$padj.threshold)
lfc_threshold <- as.numeric(args$lfc.threshold)
top_n_heatmap <- as.integer(args$top.n.heatmap)
min_replicates <- as.integer(args$min.replicates)

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
  if (!pair_col %in% colnames(meta)) {
    stop("--pair-col '", pair_col, "' not in pseudobulk_metadata.csv (columns: ",
         paste(colnames(meta), collapse = ", "), ") -- rerun scripts/10 with --pair-col")
  }
  n_unpaired <- sum(is.na(meta[[pair_col]]) | meta[[pair_col]] == "")
  if (n_unpaired > 0) cat(sprintf("Paired mode: dropping %d pseudo-samples with no '%s'\n", n_unpaired, pair_col))
  meta <- meta[!(is.na(meta[[pair_col]]) | meta[[pair_col]] == ""), ]
  meta$pair <- factor(meta[[pair_col]])
}
counts <- counts[, meta$pseudo_sample_id, drop = FALSE]
meta$condition <- factor(meta$condition, levels = c(condition_b, condition_a))  # condition_b = reference

cat(sprintf("Loaded %d genes x %d pseudo-samples across %d cell types (comparing '%s' vs '%s'%s)\n",
            nrow(counts), ncol(counts), length(unique(meta$cell_type)), condition_a, condition_b,
            if (is.null(pair_col)) "" else sprintf(", paired by '%s'", pair_col)))

all_results <- list()
vst_by_celltype <- list()

for (ct in sort(unique(meta$cell_type))) {
  ct_meta <- meta[meta$cell_type == ct, ]
  if (!is.null(pair_col)) {
    # Keep only pairs observed in both conditions for this cell type; a pair
    # seen in one condition only is absorbed entirely by its pair term and
    # contributes nothing to the within-pair condition estimate.
    has_a <- unique(as.character(ct_meta$pair[ct_meta$condition == condition_a]))
    has_b <- unique(as.character(ct_meta$pair[ct_meta$condition == condition_b]))
    complete <- intersect(has_a, has_b)
    dropped <- setdiff(unique(as.character(ct_meta$pair)), complete)
    if (length(dropped) > 0) {
      cat(sprintf("  '%s': dropping incomplete pair(s) %s\n", ct, paste(dropped, collapse = ", ")))
    }
    ct_meta <- ct_meta[as.character(ct_meta$pair) %in% complete, ]
    ct_meta$pair <- droplevels(ct_meta$pair)
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

  cat(sprintf("=== DESeq2: %s (%d vs %d pseudo-samples%s) ===\n", ct, n_a, n_b,
              if (is.null(pair_col)) "" else sprintf(", %d pairs", nlevels(ct_meta$pair))))
  ct_counts <- counts[, ct_meta$pseudo_sample_id, drop = FALSE]
  keep_genes <- rowSums(ct_counts) > 0
  ct_counts <- ct_counts[keep_genes, , drop = FALSE]

  design <- if (is.null(pair_col)) ~condition else ~pair + condition
  dds <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = design)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", condition_a, condition_b), alpha = padj_threshold)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$cell_type <- ct
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
heatmap_mat <- heatmap_mat[rowSums(!is.na(heatmap_mat)) > 0, , drop = FALSE]
# Pseudo-samples from skipped cell types (too few replicates/pairs) have no
# VST values at all; an all-NA column makes hclust fail, so drop them.
heatmap_mat <- heatmap_mat[, colSums(!is.na(heatmap_mat)) > 0, drop = FALSE]
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))  # z-score per gene across pseudo-samples
# Any NA left (a gene untested in one cell type) can still break distance
# computation for clustering -- fall back to unclustered order in that case.
can_cluster <- !anyNA(heatmap_mat_scaled)
if (!can_cluster) cat("  Heatmap has genes untested in some cell types -- drawing without clustering\n")

anno_idx <- match(colnames(heatmap_mat_scaled), meta$pseudo_sample_id)
col_anno <- if (is.null(pair_col)) {
  HeatmapAnnotation(condition = meta$condition[anno_idx], cell_type = meta$cell_type[anno_idx])
} else {
  HeatmapAnnotation(condition = meta$condition[anno_idx], cell_type = meta$cell_type[anno_idx],
                    pair = as.character(meta$pair[anno_idx]))
}

png(file.path(output_dir, "top_de_genes_heatmap.png"), width = 2400, height = max(1200, 40 * nrow(heatmap_mat_scaled)), res = 200)
draw(Heatmap(heatmap_mat_scaled, name = "z-score (VST)",
             top_annotation = col_anno,
             show_column_names = FALSE,
             row_names_gp = grid::gpar(fontsize = 8),
             cluster_rows = can_cluster, cluster_columns = can_cluster,
             na_col = "grey90"))
dev.off()

write.csv(combined_df, file.path(output_dir, "de_all_celltypes_combined.csv"), row.names = FALSE)

cat("=== Done. Outputs in", output_dir, "===\n")
cat("  de_<cell_type>.csv          -- per-cell-type DESeq2 results\n")
cat("  de_all_celltypes_combined.csv -- all cell types concatenated\n")
cat("  volcano_plots/volcano_<cell_type>.png\n")
cat("  top_de_genes_heatmap.png\n")

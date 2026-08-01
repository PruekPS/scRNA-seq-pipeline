#!/usr/bin/env Rscript
# Step 0b: Convert a Trailmaker-exported Seurat object (.rds) to Parse's own
# DGE_filtered flat-file format, so it can be picked up by
# scripts/00_load_parse_output.py (which already handles that format).
#
# Trailmaker's Insights module does not support uploading .h5ad, and its
# "Seurat object download" is the primary export path for code-based
# analysis -- so a Trailmaker export lands as an .rds, not a raw
# DGE_filtered/ folder.
#
# Deliberately NOT bridging live through reticulate/zellkonverter here:
# zellkonverter's writeH5AD() goes through Bioconductor's basilisk, which
# provisions its OWN isolated Python (via pyenv) rather than reusing the
# scrna conda env -- on this machine that meant trying to compile Python
# 3.14 from source, which failed on a missing Tcl/Tk dependency. Writing
# flat files and reusing the already-tested Python loader avoids needing
# any live R<->Python bridge at all, which is also more robust for a real
# HPC pipeline where the R and Python steps may run as separate job stages.
#
# Usage:
#   Rscript 00b_seurat_to_h5ad.R /path/to/trailmaker_export.rds /path/to/output_dir/
#   python 00_load_parse_output.py --sample-dir /path/to/output_dir/ <sample_name> --output raw.h5ad

suppressMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 00b_seurat_to_h5ad.R <input.rds> <output_dir>")
}
input_rds <- args[1]
output_dir <- args[2]

cat(sprintf("Loading Seurat object from %s...\n", input_rds))
seurat_obj <- readRDS(input_rds)

if (!inherits(seurat_obj, "Seurat")) {
  stop(sprintf("%s does not contain a Seurat object (got class: %s)",
               input_rds, paste(class(seurat_obj), collapse = ", ")))
}

cat(sprintf("Seurat object: %d cells x %d features, assays: %s\n",
            ncol(seurat_obj), nrow(seurat_obj),
            paste(names(seurat_obj@assays), collapse = ", ")))

DefaultAssay(seurat_obj) <- "RNA"
counts <- GetAssayData(seurat_obj, layer = "counts")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Writing count_matrix.mtx (cells x genes, matching split-pipe's own orientation)...\n")
# Seurat stores features x cells; split-pipe's DGE convention is cells x
# genes (rows = cells) -- transpose so the Python loader's assumptions hold.
Matrix::writeMM(Matrix::t(counts), file.path(output_dir, "count_matrix.mtx"))

cat("Writing all_genes.csv...\n")
genes_df <- data.frame(
  gene_id = rownames(seurat_obj),
  gene_name = rownames(seurat_obj),
  genome = "unknown"
)
write.csv(genes_df, file.path(output_dir, "all_genes.csv"), row.names = FALSE)

cat("Writing cell_metadata.csv...\n")
cell_meta <- seurat_obj@meta.data
cell_meta$bc_wells <- colnames(seurat_obj)
cell_meta <- cell_meta[, c("bc_wells", setdiff(colnames(cell_meta), "bc_wells"))]
write.csv(cell_meta, file.path(output_dir, "cell_metadata.csv"), row.names = FALSE)

cat(sprintf("Done. Load with:\n  python scripts/00_load_parse_output.py --sample-dir %s <sample_name> --output raw.h5ad\n", output_dir))

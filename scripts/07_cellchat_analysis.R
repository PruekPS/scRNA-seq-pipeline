#!/usr/bin/env Rscript
# Step 7: Cell-cell communication analysis with CellChat.
#
# Reads the flat MTX/CSV files written by scripts/06_export_for_cellchat.py
# (log-normalized matrix, genes, cell metadata with a `cell_type` grouping
# column) and runs CellChat's standard workflow: overexpressed-gene/
# interaction detection -> communication-probability inference -> pathway
# aggregation -> network summary. No live R<->Python bridge is used anywhere
# in this step (see 06_export_for_cellchat.py docstring for why).
#
# CellChat is NOT installable via conda/bioconda as a single package, and
# neither are two of its dependencies (NMF needs >=0.28, conda caps at
# 0.21.0; presto is GitHub-only and required, not optional -- without it
# identifyOverExpressedGenes() errors outright). Run
# `bash envs/install_r_extras.sh` once (as a SLURM job, not on the login
# node) to install all of it -- see that file for the full story, including
# why jinworks/CellChat (not the archived sqjin/CellChat) is the fork to use.
#
# Validated on the real 907,696-cell full-scale run (2026-08-02): 694,803
# cells across 10 cell types, 291 significant L-R interactions across 19
# pathways, completed in ~1.5h. --workers=1 is REQUIRED at this scale, not
# just a default -- future::multisession with >1 worker duplicates
# CellChat's large intermediate data per worker process and OOM-killed a
# 128gb job where --workers=1 completed fine in ~300gb. This is the
# slowest step in the pipeline at real cell counts (permutation-based
# significance testing) -- budget hours, not minutes.
#
# Usage (local):
#   Rscript 07_cellchat_analysis.R --input-dir cellchat_input/ \
#     --output-dir cellchat_results/ --species human
#
# Usage (HPC): see pipeline/run_downstream.slurm.

suppressMessages({
  library(CellChat)
  library(Matrix)
})

parse_args <- function(args) {
  out <- list(species = "human", min.cells = 10, workers = 1)
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
if (is.null(args$input.dir) || is.null(args$output.dir)) {
  stop("Usage: Rscript 07_cellchat_analysis.R --input-dir=DIR --output-dir=DIR [--species=human|mouse] [--min-cells=10] [--workers=1]")
}
input_dir <- args$input.dir
output_dir <- args$output.dir
species <- args$species
min_cells <- as.integer(args$min.cells)
workers <- as.integer(args$workers)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Loading exported data from", input_dir, "===\n")
mat <- Matrix::readMM(file.path(input_dir, "normalized_matrix.mtx"))
genes <- read.csv(file.path(input_dir, "genes.csv"))$gene
meta <- read.csv(file.path(input_dir, "cell_metadata.csv"), row.names = 1)

mat <- as(mat, "CsparseMatrix")
rownames(mat) <- genes
colnames(mat) <- rownames(meta)

if (!"cell_type" %in% colnames(meta)) {
  stop("cell_metadata.csv has no 'cell_type' column -- re-export with 06_export_for_cellchat.py")
}
meta$cell_type <- as.factor(meta$cell_type)
cat(sprintf("Loaded %d genes x %d cells across %d cell types\n",
            nrow(mat), ncol(mat), nlevels(meta$cell_type)))

cat("=== Creating CellChat object ===\n")
cellchat <- createCellChat(object = mat, meta = meta, group.by = "cell_type")

db <- if (species == "mouse") CellChatDB.mouse else CellChatDB.human
cellchat@DB <- db

cellchat <- subsetData(cellchat)

if (workers > 1) {
  suppressMessages(future::plan("multisession", workers = workers))
}

cat("=== Identifying overexpressed genes/interactions ===\n")
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cat("=== Computing communication probabilities (slowest step) ===\n")
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = min_cells)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cat("=== Writing results ===\n")
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file.path(output_dir, "significant_interactions.csv"), row.names = FALSE)
cat(sprintf("%d significant ligand-receptor interactions across %d pathways\n",
            nrow(df.net), length(unique(df.net$pathway_name))))

pdf(file.path(output_dir, "network_plots.pdf"), width = 8, height = 8)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                  label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                  label.edge = FALSE, title.name = "Interaction strength")
tryCatch({
  netVisual_bubble(cellchat, remove.isolate = FALSE)
}, error = function(e) cat("netVisual_bubble skipped:", conditionMessage(e), "\n"))
dev.off()

saveRDS(cellchat, file.path(output_dir, "cellchat_object.rds"))
cat("Wrote", output_dir, "/{significant_interactions.csv, network_plots.pdf, cellchat_object.rds}\n")
cat("=== Done ===\n")

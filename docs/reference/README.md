# Reference: design decisions and validation reports

The [guide](../01-getting-started.md) says what to run. These pages say *why it
is built this way* and *how well it works* — each validation report checks a
step against real or public data and against an alternative method, with the
numbers computed, not estimated.

| Page | What it covers | Data |
|---|---|---|
| [Design decisions](design-decisions.md) | scope, the pipeline as built, the decisions behind it, QC metrics and their biology, roadmap, all references | — |
| [Quality control](qc.md) | adaptive (MAD) thresholds vs. the standard fixed cutoff | PBMC3k (public) |
| [Normalization, HVG, PCA](normalization-pca.md) | two gene-selection methods compared | PBMC3k |
| [Clustering](clustering.md) | Leiden vs. Louvain at five resolutions | PBMC3k |
| [Markers and annotation](annotation.md) | Wilcoxon vs. t-test; how marker *specificity* makes or breaks annotation | PBMC3k |
| [Real-data validation](real-data-validation.md) | the full chain on Parse's real 1M-cell dataset; what only shows up at scale | Parse 1M PBMC (public) |
| [Cell-cell communication](cellchat.md) | CellChat results, and what its permutation test adds | Parse 1M PBMC, 694,803 cells |
| [Trajectory design](trajectory-design.md) | PAGA + CytoTRACE + CellRank: method, and what the experiment must contain | design only, not yet run |

**What has been validated where:**

| Steps | Status |
|---|---|
| 00–04b, 02b Harmony, 06–07 CellChat | validated on real Parse data at full scale (~1M cells) |
| 10–11 pseudobulk DE (incl. per-experiment design), 12–13 TCR | tested on synthetic data in split-pipe's output format with planted truth; not yet on real data |
| 08–09 trajectories | built and checked against CellRank's documentation; not yet run |

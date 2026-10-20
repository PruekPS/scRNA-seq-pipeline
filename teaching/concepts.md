# Concepts — why each step exists

One section per stage, in pipeline order. Each explains what problem the step
solves, the idea behind the method, and the judgement call it leaves to you.
Commands are in the [guide](../docs/01-getting-started.md).

---

## Testing small

A full run on ~1 million cells takes many hours and hundreds of gigabytes. A
bug that appears at the end costs all of that. While this pipeline was being
built, QC was first run at full scale with no small test: it ran for 24 hours
and then ran out of memory — something a 5-minute test on 50,000 cells would
have shown.

So: every new step, every new dataset, every changed setting — run it small
first, *look at the output*, then scale up. A run that finishes but produces
nonsense is worse than a crash, because it looks like a result. The pipeline's
synthetic test data goes one step further: it has a planted answer, so you can
check the result is *right*, not only that it ran.

---

## Quality control

**The problem.** Not every barcode is a healthy single cell. Some are empty or
near-empty, some are dying cells, some are two cells captured together.

**Three metrics per cell:**

- **Total counts** — how much RNA was captured. Too low: degraded or empty.
- **Genes detected** — complexity. Very low: empty barcode. Unusually high for
  the depth: possibly two cells.
- **Mitochondrial %** — when a cell's membrane is damaged, cytoplasmic RNA
  leaks out but mitochondria (with their own membranes) stay longer, so the
  mitochondrial share rises. High mitochondrial % flags stressed or dying
  cells.

**Adaptive vs. fixed thresholds.** A fixed rule ("remove cells with fewer than
200 genes") assumes every dataset has the same depth and composition. They do
not: tissue, cell type and chemistry all shift what "normal" looks like. The
pipeline instead flags cells far from *this dataset's* median — more than *k*
median absolute deviations (MAD) away. MAD is a spread measure that a few
extreme cells cannot distort. The trade-off: thresholds are dataset-specific,
so two datasets are not filtered by the same numbers. On PBMC3k this removed
every cell the standard fixed rule removed, plus low-complexity cells the fixed
rule cannot see ([report](../docs/reference/qc.md)).

**Doublets.** Two cells that end up with the same barcode look like one cell
with a mixed profile — for example, a "cell" expressing both T-cell and
monocyte genes, which could be mistaken for a new cell type. Scrublet simulates
artificial doublets by adding random pairs of real cells, then flags real cells
that look like those simulations. In split-pool chemistry, doublets are two
cells that by chance receive the same barcode combination, so their rate rises
with the number of cells per sublibrary.

**Why per sample.** Scrublet's model assumes one library. Pooling many samples
breaks that assumption — and on 987,585 cells it also ran out of memory at
400 GB. Running per sample fixed both. A fix that is right scientifically is
not always also the performance fix, but it is worth checking.

**Your judgement call:** *k*. Look at the QC violin plot: the flagged cells
should be a distinct tail, not a slice through the main population.

---

## Normalization, highly variable genes and PCA

**Normalization.** A cell sequenced twice as deeply shows roughly twice the
counts for every gene. Rescaling each cell to the same total (then taking logs)
removes that, so differences reflect biology, not depth.

**Highly variable genes (HVGs).** Of ~30,000–60,000 genes, most carry little
information between cell types — housekeeping genes, or genes barely detected.
Keeping the ~2,000 that vary most keeps the signal and drops most of the noise.
How "variable" is measured matters: two standard methods disagree on about one
gene in five ([report](../docs/reference/normalization-pca.md)).

**Scaling** gives every gene mean 0 and variance 1, so a highly expressed gene
does not dominate just because its numbers are larger.

**PCA** finds the few directions along which cells differ most and summarises
each cell by its position along them (50 numbers instead of 2,000). Everything
downstream — neighbours, clusters, UMAP — works in this compressed space.

**Your judgement call:** how many components carry signal. The elbow plot
should drop and then flatten; components after the flattening are mostly
noise.

---

## Clustering

The pipeline connects each cell to its nearest neighbours in PCA space, then the
Leiden algorithm finds groups of cells more connected to each other than to the
rest. A cluster is a stand-in for "probably the same cell type or state".

**Resolution** sets how finely to split: low merges related populations, high
can split one real population into noise-driven pieces. There is no correct
value — that is why the pipeline scans several. Leiden replaced the older
Louvain algorithm because Louvain can produce badly connected clusters; on
PBMC3k they agree closely ([report](../docs/reference/clustering.md)).

**UMAP** is only a picture: a 2-D layout that keeps neighbours close. Distances
between far-apart islands, and island sizes, do not mean much. Never base a
conclusion on UMAP alone.

**Your judgement call:** the resolution. Look at the top marker genes of each
cluster at two resolutions: does a split look like real biology (distinct
markers) or noise?

---

## Annotation

**How a label is computed.** First, a Wilcoxon test ranks each cluster's most
distinctive genes against all other cells. Then each cluster is scored against
marker gene sets for known cell types (CD3E for T cells, MS4A1 for B cells,
LYZ for monocytes, …).

**Why "Ambiguous" is a feature.** The first version always took the highest
score. On real data with 20 clusters it confidently called several clusters
"Platelet" because every score was negative and platelet's was merely the
least negative ([report](../docs/reference/real-data-validation.md#annotation-at-scale)).
The current annotator gives a label only when the best score is clearly ahead
of the runner-up; otherwise it says `Ambiguous (A / B)`. An honest "I am not
sure between these two" tells you where to look; a confident wrong label sends
you in the wrong direction. Remember this whenever you build a classifier.

**Specificity beats relevance.** Adding genes that genuinely mark dendritic
cells but are also expressed by monocytes made annotation *worse*: it pulled
two monocyte clusters into the dendritic-cell label
([report](../docs/reference/annotation.md)). Choose markers that are specific,
not merely associated.

**Your judgement call:** what the Ambiguous clusters are. Their top marker
genes usually tell you.

---

## Batch effects

A **batch effect** is a technical difference — processing run, lane, reagent
lot — big enough that cells cluster by batch instead of by cell type. Harmony
moves cells in PCA space so batches mix, while keeping cell types apart.

The hard part is deciding **what the batch is**:

- In split-pipe data every sublibrary contains cells from every sample, so the
  technical batch is the **sublibrary**.
- If each sample was processed separately, the sample is the batch.
- **Your biological comparison is never the batch.** Correcting by condition,
  time point or tissue removes exactly what you want to study. Correcting by
  sample in split-pipe data removes real differences between samples too.

**Your judgement call:** whether correction helped or over-corrected. Before:
do batches form their own islands within one cell type? After: they should
mix — but distinct cell types must stay distinct.

---

## Cell-cell communication

CellChat asks: for each pair of cell types, does one express a ligand and the
other its receptor, more than expected by chance? It turns average expression
into a **communication probability**, then tests it by shuffling cell-type
labels many times to see how often a probability that high arises at random.

It is an inference from expression, not a measurement: a ligand and receptor
both being expressed does not prove the cells signal to each other (proximity,
protein levels and timing all matter). Strong hits are hypotheses.

The permutation test matters most for **small groups**: a few high-expressing
cells in a rare cell type can produce an impressive-looking probability that is
no different from chance ([report](../docs/reference/cellchat.md)).

---

## Trajectories

Clustering gives categories; differentiation is a process. Trajectory methods
try to order cells along that process from a single snapshot:

- **PAGA** — which clusters are plausibly connected.
- **Pseudotime** — each cell's distance from a starting population you choose.
- **CytoTRACE** — cells that express more distinct genes tend to be less
  differentiated, giving a direction without time-course data.
- **CellRank** — combines these into each cell's probability of ending in each
  final state.

**The catches:** it is inference from one time point, not lineage tracing; the
result depends heavily on the starting population you pick; and a trajectory
between two tissues is only meaningful if both tissues were actually sampled
([design](../docs/reference/trajectory-design.md)).

---

## Differential expression

**Why not test cell by cell?** Cells from one person share that person's
genetics, environment and sample handling, so they are not independent. A
per-cell test treats 3 donors × 3,000 cells as 9,000 independent observations,
and p-values become wildly over-confident. **Pseudobulk** sums the counts of
each cell type within each sample, so each sample is one observation, then
uses DESeq2, a method built for exactly this kind of count data.

**Matching the design to the experiment.** The statistical model must reflect
how the samples relate:

- **Unpaired** (different people per group): compare the groups.
- **Paired** (the same patients before and after): include the patient in the
  model, so each patient is compared with themselves. People differ a lot from
  each other; pairing removes that variation. On synthetic data with a known
  answer, pairing found 62 planted changes instead of 38, with no false hits.
- **Batched**: include the batch, if it is spread across both groups.

A covariate that lines up exactly with the comparison (every patient in only
one group) cannot be separated from it; the pipeline refuses that design rather
than producing a meaningless answer.

**Volcano plots** show each gene's fold change against its significance: genes
that change a lot *and* confidently sit in the upper corners.

---

## T-cell receptors

Each T cell has a receptor of two chains, alpha and beta. The **CDR3** region
binds antigen and is almost unique per T-cell family. When a T cell meets its
antigen it divides, and its descendants keep the same receptor. So cells with
identical alpha and beta CDR3s form a **clonotype**, and a large clonotype has
**expanded** — usually a sign it is responding to something.

Three questions, and their traps:

- **Diversity.** A sample with more cells always shows more distinct
  clonotypes. Compare samples only after cutting each to the same number of
  cells (**rarefaction**).
- **Sharing.** The same clonotype in blood and tissue, or before and after
  treatment, suggests trafficking or persistence. The same clonotype in two
  different people ("public") can happen for common antigens.
- **Phenotype.** Which cell states (from annotation) the expanded clones
  occupy.

As for differential expression, cells of one sample are not independent: to
compare conditions, use per-sample fractions as the data points.

---

## Computing on a cluster

- **Login node vs. compute node.** The login node is shared by everyone, for
  light tasks. Real work — analysis *and* installs — goes to compute nodes via
  `sbatch` (or an interactive `srun` session).
- **Memory.** `sacct` shows the peak memory (`MaxRSS`) a job actually used.
  Ask for a bit more than that, not as much as possible: bigger requests wait
  longer.
- **More workers is not always faster.** Some parallel methods copy all the
  data into every worker, so 8 workers can need ~8× the memory.
- **Chaining.** `--dependency=afterok:<jobid>` starts a job only if the
  previous one succeeded.

---

Now try the [exercises](exercises.md).

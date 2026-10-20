# 7. Cell-cell communication (CellChat)

CellChat looks for pairs of cell types where one expresses a **ligand** and the
other its **receptor**, at levels unlikely to arise by chance, and reports a
communication probability for each. It is an inference from expression, not a
measurement of signalling: treat strong hits as hypotheses to check.

---

## Run it

Two steps — export from Python, analyse in R:

```bash
python scripts/06_export_for_cellchat.py --input 04b_annotated.h5ad \
  --output-dir cellchat_input/ --group-col cell_type_v2 --exclude-ambiguous
Rscript scripts/07_cellchat_analysis.R --input-dir=cellchat_input/ \
  --output-dir=cellchat_results/ --species=human --workers=1
```

| Option | Why |
|---|---|
| `--exclude-ambiguous` | cells 04b could not label confidently are left out; an uncertain label fed into a communication test only turns uncertainty into a false-looking signal |
| `--species` | `human` or `mouse` — picks CellChat's ligand-receptor database |
| `--workers=1` | **keep at 1.** More workers start separate R processes that each copy CellChat's large intermediate data: 8 workers ran out of memory where 1 finished |

The step exchanges flat files (counts + metadata) between Python and R rather
than using a live Python–R bridge; the bridge tool tried first failed to
install on the cluster.

---

## What you get

| Output | Shows |
|---|---|
| network plots | number and strength of significant interactions between each pair of cell types |
| pathway bar chart | which signalling pathways carry the most interactions |
| bubble plot | every significant ligand-receptor pair, by source → target |

On the real 694,803-cell PBMC run: 291 significant interactions across 19
pathways, about 1.5 hours, almost all of it in the permutation test. Full
results and a check of what the permutation test adds over simply ranking by
probability: [reference/cellchat.md](reference/cellchat.md).

---

## Limits to keep in mind

- **No donor structure.** All cells are treated as one population; this
  cannot yet show donor- or condition-specific signalling.
- **Excluded cells are invisible.** Signalling involving `Ambiguous` cells is
  not represented.
- **Finite database.** A pathway missing from the results may be missing from
  CellChat's database, not from the biology.
- **One snapshot.** Nothing about change over time.

---

Next: [8. Trajectories and fate mapping](08-trajectory.md)

# Test data

Scripts in `scripts/` are validated against the public 10x Genomics PBMC3k
dataset (the same one used by the canonical Seurat/scanpy tutorials) — not
project data. Files land in `tests/data/` and are gitignored (`*.h5ad`);
nothing here is committed.

Regenerate at any time:

```python
import scanpy as sc
adata = sc.datasets.pbmc3k()
adata.write_h5ad("tests/data/pbmc3k_raw.h5ad")
```

This is a smoke test, not a substitute for validating against real Parse
Evercode output once that's available.

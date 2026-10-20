#!/usr/bin/env python3
"""Build a small synthetic WT .h5ad + split-pipe-style TCR AIRR tables, for smoke-testing
scripts/12_tcr_attach.py and scripts/13_tcr_analysis.py (and the covariate design in 10/11).

The planted truth is known so the scripts can be checked against it:
  * 2 patients (P1, P2) x 2 time points (pre, 6h) = 4 samples, 3 sublibraries
  * TCR is detected mostly in T-lineage cells (planted ~85%) and rarely elsewhere (~3% noise)
  * P1 has 3 expanded CD8 clones present at BOTH time points (planted sharing)
  * 1 "public" clone: exactly 2 T cells in P1 and 2 in P2
  * every other T cell carries a unique clone; ~8% have only a beta chain (orphan)
  * cell ids follow split-pipe's combined form `<bc>__s<N>`

AIRR tables use standard AIRR rearrangement fields (sequence_id, locus, v_call, junction_aa, ...) plus
the per-cell columns split-pipe adds to its TCR output (cell_id, read_count, transcript_count,
isMultiplet, isPaired, sample_name; productive as "T"/"F").
This proves the plumbing, NOT the real file layout -- confirm against real output.

Expression also carries planted DE truth for scripts 10/11 (see make_wt): G0-G9 up 4x at
"6h", G10-G39 up 3x in patient P2. layers["counts"] holds the raw counts.

Usage: python tests/synthetic_tcr.py OUT_DIR
"""
import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd

AA = list("ACDEFGHIKLMNPQRSTVWY")
CELL_TYPES = {"CD4 T cell": 0.25, "CD8 T cell": 0.25, "Treg": 0.05, "NK cell": 0.08,
              "B cell": 0.12, "CD14+ Monocyte": 0.2, "cDC": 0.05}
T_LINEAGE = {"CD4 T cell", "CD8 T cell", "Treg"}
SAMPLES = {"P1_pre": ("P1", "pre"), "P1_6h": ("P1", "6h"), "P2_pre": ("P2", "pre"), "P2_6h": ("P2", "6h")}
CELLS_PER_SAMPLE = 300


def random_cdr3(rng, prefix, n=12):
    return prefix + "".join(rng.choice(AA, n)) + "F"


def random_nt(rng, aa_len):
    return "".join(rng.choice(list("ACGT"), 3 * aa_len))


def make_wt(rng):
    rows = []
    for s, (patient, tp) in SAMPLES.items():
        for i in range(CELLS_PER_SAMPLE):
            ct = rng.choice(list(CELL_TYPES), p=list(CELL_TYPES.values()))
            sub = f"s{rng.integers(1, 4)}"
            bc = f"{rng.integers(1, 49):02d}_{rng.integers(1, 97):02d}_{rng.integers(1, 97):02d}"
            rows.append({"cell": f"{bc}__{sub}", "sample": s, "patient": patient, "timepoint": tp,
                         "cell_type_v2": ct, "sublibrary": sub})
    obs = pd.DataFrame(rows).drop_duplicates("cell").set_index("cell")
    centers = {ct: rng.normal(0, 6, 2) for ct in CELL_TYPES}
    umap = np.vstack([centers[c] + rng.normal(0, 0.8, 2) for c in obs["cell_type_v2"]])
    # expression with planted DE truth for scripts 10/11:
    #   G0-G9   4x higher at "6h" than "pre" in every cell type (the condition effect)
    #   G10-G39 3x higher in patient P2 (a patient effect a paired design should absorb)
    n_genes = 300
    base = rng.lognormal(0.0, 1.0, n_genes)
    rate = np.tile(base, (len(obs), 1)) * rng.lognormal(0.0, 0.3, (len(obs), 1))
    rate[(obs["timepoint"] == "6h").to_numpy()[:, None] & (np.arange(n_genes) < 10)[None, :]] *= 4
    p2 = (obs["patient"] == "P2").to_numpy()[:, None] & ((np.arange(n_genes) >= 10) & (np.arange(n_genes) < 40))[None, :]
    rate[p2] *= 3
    counts = rng.poisson(rate).astype("float32")
    adata = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))
    adata.layers["counts"] = counts.copy()
    adata.obsm["X_umap"] = umap
    return adata


def make_airr(adata, rng):
    expanded = {f"EXP{i}": (random_cdr3(rng, "CAS"), random_cdr3(rng, "CAV")) for i in range(3)}
    public = (random_cdr3(rng, "CAS"), random_cdr3(rng, "CAV"))
    chains, pub_count = [], {"P1": 0, "P2": 0}
    for cell, row in adata.obs.iterrows():
        is_t = row["cell_type_v2"] in T_LINEAGE
        if rng.random() > (0.85 if is_t else 0.03):
            continue
        if row["patient"] == "P1" and row["cell_type_v2"] == "CD8 T cell" and rng.random() < 0.35:
            beta, alpha = expanded[f"EXP{rng.integers(0, 3)}"]
        elif is_t and pub_count[row["patient"]] < 2:
            beta, alpha = public  # exactly 2 cells per patient -> a guaranteed cross-patient (public) clone
            pub_count[row["patient"]] += 1
        else:
            beta, alpha = random_cdr3(rng, "CAS"), random_cdr3(rng, "CAV")
        orphan = rng.random() < 0.08
        for locus, cdr3, v, j in (("TRB", beta, "TRBV5-1", "TRBJ2-7"), ("TRA", alpha, "TRAV12-2", "TRAJ33")):
            if orphan and locus == "TRA":
                continue
            chains.append({
                "sequence_id": f"{cell}__{locus}", "sequence": "ACGT" * 20, "locus": locus, "rev_comp": "F",
                "productive": "T", "v_call": v, "d_call": "", "j_call": j, "c_call": locus[:2] + "C",
                "junction": random_nt(rng, len(cdr3)), "junction_aa": cdr3, "cdr3": random_nt(rng, len(cdr3)),
                "cdr3_aa": cdr3, "vj_in_frame": "T", "stop_codon": "F",
                "cell_id": cell, "read_count": int(rng.integers(5, 200)), "transcript_count": int(rng.integers(1, 30)),
                "isMultiplet": "0", "isPaired": "0" if orphan else "1", "sample_name": row["sample"]})
    return pd.DataFrame(chains)


def main(out: pathlib.Path):
    rng = np.random.default_rng(7)
    out.mkdir(parents=True, exist_ok=True)
    adata = make_wt(rng)
    adata.write_h5ad(out / "wt_annotated.h5ad")
    airr = make_airr(adata, rng)
    for s, df in airr.groupby("sample_name"):
        d = out / "airr_root" / s / "TCR_filtered"
        d.mkdir(parents=True, exist_ok=True)
        df.to_csv(d / "tcr_annotation_airr.tsv", sep="\t", index=False)
    print(f"wrote {adata.n_obs} WT cells, {airr['cell_id'].nunique()} TCR cells, {len(airr)} chains to {out}")


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))

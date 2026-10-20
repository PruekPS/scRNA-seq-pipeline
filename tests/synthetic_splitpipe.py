#!/usr/bin/env python3
"""Synthetic data in split-pipe v1.9.0's OUTPUT LAYOUT (WT + TCR), for an end-to-end
integration test of the whole pipeline: 00 load -> 01 QC -> 02 PCA -> 02b Harmony ->
03 cluster -> 04b annotate -> 12/13 TCR -> 10/11 pseudobulk DE.

Formats follow split-pipe's documented outputs and the real Parse output already
validated in this repo (1M PBMC dataset):
  WT  <wt_comb>/<sample>/DGE_filtered/{count_matrix.mtx, all_genes.csv, cell_metadata.csv}
      plus <wt_comb>/all-sample/DGE_filtered/ (the same cells again, all samples together)
      count_matrix.mtx = cells x genes; all_genes.csv = gene_id,gene_name,genome;
      cell_metadata.csv = bc_wells,sample,species,gene_count,tscp_count,read_count,
                          bc1_well..bc3_well,bc1_wind..bc3_wind
      bc_wells = "NN_NN_NN__s<N>" (well indexes + sublibrary position in the comb list)
  TCR <tcr_comb>/<sample>/TCR_filtered/tcr_annotation_airr.tsv (+ all-sample/), cell_id =
      the same bc_wells string; standard AIRR fields + split-pipe's per-cell columns.
  Also: tcr_persub/s<N>/all-sample/... with bare NN_NN_NN ids (per-sublibrary runs) and
        tcr_comb_misordered/ (TCR combined with sublibraries s2/s3 swapped -- must be refused).

Planted truth (truth_cells.csv, truth.json) so results can be checked, not just run:
  6 cell types with their real marker genes (from 04b's human dictionary);
  2 patients x 2 time points (pre, 6h); 3 sublibraries, every sample in every sublibrary;
  DEG1-DEG10 4x higher at 6h in every cell type; PATG1-PATG30 3x higher in patient P2;
  3 expanded CD8 clones in P1 present at BOTH time points; 1 public clone (2 T cells per patient).
This is synthetic: it proves the pipeline handles split-pipe's formats and recovers known
signal. It is not a substitute for checking the real outputs once they exist.

Usage: python tests/synthetic_splitpipe.py OUT_DIR
"""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

SAMPLES = {"P1_pre": ("P1", "pre"), "P1_6h": ("P1", "6h"), "P2_pre": ("P2", "pre"), "P2_6h": ("P2", "6h")}
CELLS_PER_SAMPLE = 500
N_SUBLIBS = 3
TYPES = {  # proportion, marker genes (04b human dictionary)
    "CD4 T cell": (0.28, ["IL7R", "CD3D", "CD3E", "CD4"]),
    "CD8 T cell": (0.22, ["CD8A", "CD8B", "CD3D", "CD3E"]),
    "NK cell": (0.10, ["NKG7", "GNLY", "KLRD1", "NCAM1"]),
    "B cell": (0.15, ["MS4A1", "CD79A", "CD79B", "CD19"]),
    "CD14+ Monocyte": (0.20, ["CD14", "LYZ", "S100A8", "S100A9"]),
    "cDC": (0.05, ["FCER1A", "CLEC10A", "CD1C"]),
}
T_TYPES = {"CD4 T cell", "CD8 T cell"}
OTHER_MARKERS = ["FOXP3", "IL2RA", "CTLA4", "MZB1", "JCHAIN", "XBP1", "FCGR3A", "MS4A7", "CDKN1C",
                 "LILRA4", "IL3RA", "GZMB", "PPBP", "PF4", "ITGA2B", "CD34", "SOX4", "SPINK2"]
MT = ["MT-CO1", "MT-CO2", "MT-ND1", "MT-ND4", "MT-ATP6", "MT-CYB"]
AA = list("ACDEFGHIKLMNPQRSTVWY")


def plate_well(wind: int, cols: int = 12) -> str:
    return "ABCDEFGHIJKLMNOP"[(wind - 1) // cols] + str((wind - 1) % cols + 1)


def make_cells(rng):
    rows, used = [], set()
    types, probs = list(TYPES), np.array([v[0] for v in TYPES.values()])
    for si, (s, (patient, tp)) in enumerate(SAMPLES.items()):
        for _ in range(CELLS_PER_SAMPLE):
            sub = int(rng.integers(1, N_SUBLIBS + 1))
            while True:  # round-1 well encodes the sample (WT kit: 48 wells, 12 per sample here)
                w1 = 12 * si + int(rng.integers(1, 13))
                w2, w3 = int(rng.integers(1, 97)), int(rng.integers(1, 97))
                if (sub, w1, w2, w3) not in used:
                    used.add((sub, w1, w2, w3))
                    break
            rows.append({"bc_wells": f"{w1:02d}_{w2:02d}_{w3:02d}__s{sub}", "sample": s, "patient": patient,
                         "timepoint": tp, "true_type": rng.choice(types, p=probs / probs.sum()), "sublib": sub,
                         "bc1_wind": w1, "bc2_wind": w2, "bc3_wind": w3})
    return pd.DataFrame(rows)


def make_counts(cells, rng):
    filler = [f"FILL{i}" for i in range(2500)]
    degs = [f"DEG{i}" for i in range(1, 11)]
    patg = [f"PATG{i}" for i in range(1, 31)]
    markers = sorted({g for _, m in TYPES.values() for g in m} | set(OTHER_MARKERS))
    genes = markers + MT + degs + patg + filler
    gi = {g: i for i, g in enumerate(genes)}
    n, g = len(cells), len(genes)
    base = np.full(g, 0.03)
    base[[gi[x] for x in filler]] = rng.lognormal(-1.2, 1.0, len(filler))
    base[[gi[x] for x in degs + patg]] = 1.5
    base[[gi[x] for x in MT]] = 4.0
    rate = np.tile(base, (n, 1))
    for ct, (_, mk) in TYPES.items():
        m = (cells["true_type"] == ct).to_numpy()
        rate[np.ix_(m, [gi[x] for x in mk])] = 6.0
    rate[np.ix_((cells["timepoint"] == "6h").to_numpy(), [gi[x] for x in degs])] *= 4
    rate[np.ix_((cells["patient"] == "P2").to_numpy(), [gi[x] for x in patg])] *= 3
    rate *= rng.lognormal(0, 0.35, (n, 1))
    X = sp.csr_matrix(rng.poisson(rate).astype(np.int32))
    genes_df = pd.DataFrame({"gene_id": [f"ENSGSYN{i:07d}" for i in range(g)], "gene_name": genes, "genome": "hg38"})
    return X, genes_df


def cell_metadata(cells, X):
    md = pd.DataFrame({
        "bc_wells": cells["bc_wells"], "sample": cells["sample"], "species": "hg38",
        "gene_count": np.asarray((X > 0).sum(axis=1)).ravel(), "tscp_count": np.asarray(X.sum(axis=1)).ravel()})
    md["read_count"] = (md["tscp_count"] * 1.6).astype(int)
    md["bc1_well"] = cells["bc1_wind"].map(lambda w: plate_well(w, 12))
    md["bc2_well"] = cells["bc2_wind"].map(lambda w: plate_well(w, 12))
    md["bc3_well"] = cells["bc3_wind"].map(lambda w: plate_well(w, 12))
    for k in (1, 2, 3):
        md[f"bc{k}_wind"] = cells[f"bc{k}_wind"]
    return md


def write_dge(d: pathlib.Path, X, genes_df, md):
    d.mkdir(parents=True, exist_ok=True)
    scipy.io.mmwrite(str(d / "count_matrix.mtx"), X)
    genes_df.to_csv(d / "all_genes.csv", index=False)
    md.to_csv(d / "cell_metadata.csv", index=False)


def make_tcr(cells, rng):
    def cdr3(prefix):
        return prefix + "".join(rng.choice(AA, 12)) + "F"
    expanded = [(cdr3("CAS"), cdr3("CAV")) for _ in range(3)]
    public = (cdr3("CAS"), cdr3("CAV"))
    pub_left = {"P1": 2, "P2": 2}
    rows, truth_clone = [], {}
    for _, c in cells.iterrows():
        is_t = c["true_type"] in T_TYPES
        if rng.random() > (0.85 if is_t else 0.02):
            continue
        if c["patient"] == "P1" and c["true_type"] == "CD8 T cell" and rng.random() < 0.3:
            k = int(rng.integers(0, 3)); beta, alpha = expanded[k]; truth_clone[c["bc_wells"]] = f"EXP{k}"
        elif is_t and pub_left[c["patient"]] > 0:
            beta, alpha = public; pub_left[c["patient"]] -= 1; truth_clone[c["bc_wells"]] = "PUBLIC"
        else:
            beta, alpha = cdr3("CAS"), cdr3("CAV")
        for locus, seq, v, j in (("TRB", beta, "TRBV5-1", "TRBJ2-7"), ("TRA", alpha, "TRAV12-2", "TRAJ33")):
            rows.append({
                "sequence_id": f"{c['bc_wells']}__CTG{locus}", "sequence": "".join(rng.choice(list("ACGT"), 300)),
                "locus": locus, "rev_comp": "F", "productive": "T", "v_call": v, "d_call": "", "j_call": j,
                "c_call": locus[:2] + "C", "junction": "".join(rng.choice(list("ACGT"), 3 * len(seq))),
                "junction_aa": seq, "cdr3": "".join(rng.choice(list("ACGT"), 3 * len(seq) - 6)), "cdr3_aa": seq[1:-1],
                "cell_id": c["bc_wells"], "read_count": int(rng.integers(5, 300)),
                "transcript_count": int(rng.integers(1, 40)), "isMultiplet": "0", "isPaired": "1",
                "sample_name": c["sample"]})
    return pd.DataFrame(rows), truth_clone


def write_tcr_tree(root: pathlib.Path, tcr: pd.DataFrame):
    for s, df in list(tcr.groupby("sample_name")) + [("all-sample", tcr)]:
        d = root / s / "TCR_filtered"
        d.mkdir(parents=True, exist_ok=True)
        df.to_csv(d / "tcr_annotation_airr.tsv", sep="\t", index=False)


def main(out: pathlib.Path):
    rng = np.random.default_rng(11)
    cells = make_cells(rng)
    X, genes_df = make_counts(cells, rng)
    md = cell_metadata(cells, X)
    for s in SAMPLES:
        m = (cells["sample"] == s).to_numpy()
        write_dge(out / "wt_comb" / s / "DGE_filtered", X[m], genes_df, md[m])
    write_dge(out / "wt_comb" / "all-sample" / "DGE_filtered", X, genes_df, md)

    tcr, truth_clone = make_tcr(cells, rng)
    write_tcr_tree(out / "tcr_comb", tcr)
    bare = tcr.assign(cell_id=tcr["cell_id"].str.replace(r"__s\d+$", "", regex=True),
                      sub=tcr["cell_id"].str.extract(r"__s(\d+)$")[0])
    for sub, df in bare.groupby("sub"):
        write_tcr_tree(out / "tcr_persub" / f"s{sub}", df.drop(columns="sub"))
    swap = {"__s2": "__s3", "__s3": "__s2"}
    mis = tcr.assign(cell_id=tcr["cell_id"].str.replace(r"__s[23]$", lambda m: swap[m.group(0)], regex=True))
    write_tcr_tree(out / "tcr_comb_misordered", mis)

    pd.DataFrame([{"sample": s, "patient": p, "timepoint": t} for s, (p, t) in SAMPLES.items()]).to_csv(
        out / "sample_sheet.csv", index=False)
    cells[["bc_wells", "sample", "true_type"]].assign(true_clone=cells["bc_wells"].map(truth_clone)).to_csv(
        out / "truth_cells.csv", index=False)
    json.dump({"planted_de_genes": [f"DEG{i}" for i in range(1, 11)],
               "patient_effect_genes": [f"PATG{i}" for i in range(1, 31)],
               "n_cells": int(len(cells)), "n_tcr_cells": int(tcr["cell_id"].nunique())},
              open(out / "truth.json", "w"), indent=2)
    print(f"wrote {len(cells)} cells x {X.shape[1]} genes, {tcr['cell_id'].nunique()} TCR cells -> {out}")


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))

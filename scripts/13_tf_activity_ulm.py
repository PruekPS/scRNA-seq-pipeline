"""Step 13 (optional): transcription-factor activity with a univariate linear
model (ULM) over a signed TF -> target network (e.g. CollecTRI).

Method = decoupler's `run_ulm` (Badia-i-Mompel et al., Bioinformatics Advances
2022): for each TF and each sample/contrast, regress the gene-level values y
(a DE statistic, or centered log-expression) on the TF's signed target weights
x (+1 activation / -1 repression for targets, 0 for every other gene), with an
intercept, across all genes. The TF activity score is the t-value of the
slope -- positive = its targets move in the direction the TF would push them.
Re-implemented here (vectorised numpy) because installing decoupler 1.x into
this environment would downgrade numpy 2.x -> 1.26 via its numba pin.

Network CSV: columns source, target, weight (one row per TF-target edge).
Fetch CollecTRI once (Muller-Dott et al., Nucleic Acids Res 2023) e.g. from
OmniPath's REST service -- not inside this analysis step.

Two input modes (either or both):
  --de-dir DIR          every de_*.csv from scripts/11 (uses its 'stat' column)
                        -> one activity per TF per contrast
  --pseudobulk-dir DIR  pseudobulk_counts.mtx / genes.csv / pseudobulk_metadata.csv
                        from scripts/10 -> log-CPM, centred per gene across
                        pseudo-samples -> one activity per TF per pseudo-sample

Outputs (in --output-dir): tf_activity_de.csv / tf_activity_pseudobulk.csv
(long format: tf, sample_or_contrast, score, pvalue, n_targets).

Usage:
    python 13_tf_activity_ulm.py --net collectri.csv --de-dir de/ --output-dir tf/
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy import stats


def ulm(Y, genes, net, min_targets=5):
    """Y: genes x samples array (finite), genes: list. Returns (scores, pvals, ntargets) DataFrames TF x samples."""
    gidx = pd.Index(genes)
    net = net[net["target"].isin(gidx)]
    sizes = net.groupby("source").size()
    tfs = sizes.index[sizes >= min_targets]
    W = np.zeros((len(gidx), len(tfs)))
    tpos = {t: i for i, t in enumerate(tfs)}
    for s, t, w in net[net["source"].isin(tfs)][["source", "target", "weight"]].itertuples(index=False):
        W[gidx.get_loc(t), tpos[s]] = w
    n = Y.shape[0]
    Xc = W - W.mean(axis=0)                      # centred predictors (genes x TFs)
    Yc = Y - Y.mean(axis=0)                      # centred responses (genes x samples)
    sxx = (Xc ** 2).sum(axis=0)                  # per TF
    beta = (Xc.T @ Yc) / sxx[:, None]            # TFs x samples
    syy = (Yc ** 2).sum(axis=0)                  # per sample
    rss = syy[None, :] - beta ** 2 * sxx[:, None]
    se = np.sqrt(np.maximum(rss, 1e-300) / (n - 2) / sxx[:, None])
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), df=n - 2)
    return (pd.DataFrame(t, index=tfs), pd.DataFrame(p, index=tfs), sizes.reindex(tfs))


def to_long(T, P, sizes, cols):
    T.columns = cols; P.columns = cols
    out = T.stack().rename("score").to_frame().join(P.stack().rename("pvalue"))
    out.index.names = ["tf", "sample_or_contrast"]
    out = out.reset_index()
    out["n_targets"] = out["tf"].map(sizes)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--net", required=True)
    ap.add_argument("--de-dir", dest="de_dir")
    ap.add_argument("--pseudobulk-dir", dest="pb_dir")
    ap.add_argument("--output-dir", required=True, dest="output_dir")
    ap.add_argument("--min-targets", type=int, default=5, dest="min_targets")
    args = ap.parse_args()
    if not (args.de_dir or args.pb_dir):
        raise SystemExit("give --de-dir and/or --pseudobulk-dir")
    out = pathlib.Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    net = pd.read_csv(args.net)
    print(f"network: {net['source'].nunique()} TFs, {len(net):,} edges")

    if args.de_dir:
        files = sorted(f for f in pathlib.Path(args.de_dir).glob("de_*.csv") if f.name != "de_all_celltypes_combined.csv")
        series = {}
        for f in files:
            d = pd.read_csv(f).dropna(subset=["stat"]).drop_duplicates("gene").set_index("gene")["stat"]
            name = str(pd.read_csv(f, nrows=1)["cell_type"].iloc[0]) if "cell_type" in pd.read_csv(f, nrows=1) else f.stem[3:]
            series[name] = d
        M = pd.DataFrame(series)
        M = M.dropna()  # genes tested in every contrast
        T, P, sizes = ulm(M.values, M.index.tolist(), net, args.min_targets)
        res = to_long(T, P, sizes, M.columns.tolist())
        res.to_csv(out / "tf_activity_de.csv", index=False)
        print(f"DE mode: {M.shape[1]} contrasts x {M.shape[0]:,} shared genes -> {T.shape[0]} TFs scored")

    if args.pb_dir:
        d = pathlib.Path(args.pb_dir)
        X = scipy.io.mmread(d / "pseudobulk_counts.mtx")
        X = (X.toarray() if sp.issparse(X) else np.asarray(X)).astype(float)   # genes x pseudo-samples
        genes = pd.read_csv(d / "genes.csv")["gene"].tolist()
        meta = pd.read_csv(d / "pseudobulk_metadata.csv")
        keep = X.sum(axis=1) >= 10
        X, genes = X[keep], [g for g, k in zip(genes, keep) if k]
        cpm = X / X.sum(axis=0, keepdims=True) * 1e6
        L = np.log2(cpm + 1)
        L = L - L.mean(axis=1, keepdims=True)
        T, P, sizes = ulm(L, genes, net, args.min_targets)
        res = to_long(T, P, sizes, meta["pseudo_sample_id"].tolist())
        res = res.merge(meta.rename(columns={"pseudo_sample_id": "sample_or_contrast"}), on="sample_or_contrast", how="left")
        res.to_csv(out / "tf_activity_pseudobulk.csv", index=False)
        print(f"pseudobulk mode: {L.shape[1]} pseudo-samples x {L.shape[0]:,} genes -> {T.shape[0]} TFs scored")


if __name__ == "__main__":
    main()

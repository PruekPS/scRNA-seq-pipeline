"""Step 16 (optional): a publication overview figure of a dataset, Nature style.

Double-column (183 mm) figure:
  a  UMAP coloured by cell type, labelled on the map
  b  cell-type composition per sample, with the number of cells
  c  the same UMAP split by group (only with --group-col)
  d  marker genes on the UMAP (only with --genes)
Each panel is also saved alone in panels/, plus the numbers behind the figure in
source_data/ (journals ask for this), and Fig_atlas_legend.txt, a draft legend.

Run this on a computer with Arial for the final version, or accept the look-alike font
(see figstyle.py). --labels-csv lets you rename, reorder and colour the cell types:
columns `label` (as in the data), `display` (optional new name), `color` (optional hex).
Rows are the display order. Keep cell types of one lineage in one hue family.

Check before using the figure:
  - a sample that is really a copy of another (e.g. duplicated libraries) must not be
    counted twice: pool them with --pool-samples 'A+B=A and B'.
  - remove low-quality cells and doublets first, and recompute the UMAP on the cells you
    show. An embedding made before cleanup keeps their shape in the map.

Usage:
    python 16_atlas_figure.py --input 04b_annotated.h5ad --output-dir figures/atlas \
        --label-col cell_type_v2 --sample-col sample --group-col condition \
        --genes CD3E CD8A MS4A1 LYZ NKG7 PPBP
"""
import argparse
import pathlib

import anndata as ad
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from matplotlib.colors import LinearSegmentedColormap

import figstyle

CMAP = LinearSegmentedColormap.from_list("expr", ["#e9e9e9", "#c6dbef", "#4292c6", "#08306b"])


def load_table(args):
    a = ad.read_h5ad(args.input)
    if args.sample_col not in a.obs:
        print(f"  note: no '{args.sample_col}' column: treating all cells as one sample")
        a.obs[args.sample_col] = "All cells"
    for col in [args.label_col] + ([args.group_col] if args.group_col else []):
        if col not in a.obs:
            raise SystemExit(f"'{col}' not in adata.obs (available: {list(a.obs.columns)})")
    if args.embedding not in a.obsm:
        raise SystemExit(f"--embedding '{args.embedding}' not in adata.obsm ({list(a.obsm.keys())})")
    keep = np.ones(a.n_obs, dtype=bool)
    if args.exclude_labels:
        keep = ~a.obs[args.label_col].astype(str).isin(args.exclude_labels).to_numpy()
        print(f"Excluded {int((~keep).sum()):,} cells with labels {args.exclude_labels}")
    d = pd.DataFrame({"cell": a.obs_names, "x": a.obsm[args.embedding][:, 0], "y": a.obsm[args.embedding][:, 1],
                      "label": a.obs[args.label_col].astype(str).values,
                      "sample": a.obs[args.sample_col].astype(str).values})
    if args.group_col:
        d["group"] = a.obs[args.group_col].astype(str).values
    genes = []
    for g in args.genes or []:
        if g not in a.var_names and (a.raw is None or g not in a.raw.var_names):
            print(f"  note: gene {g} not in the data, skipped")
            continue
        genes.append(g)
        d[f"gene:{g}"] = expression(a, g)
    return d[keep].reset_index(drop=True), genes


def expression(a, gene):
    """log1p(counts per 10,000), from raw counts when available."""
    if "counts" in a.layers:
        counts = a.layers["counts"]
        library = np.maximum(np.asarray(counts.sum(axis=1)).ravel(), 1)
        col = counts[:, a.var_names.get_loc(gene)]
        col = col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()
        return np.log1p(col / library * 1e4)
    source = a.raw if a.raw is not None and gene in a.raw.var_names else a
    col = source[:, gene].X
    return col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()


def order_and_colors(d, labels_csv):
    present = list(d["label"].value_counts().index)
    if labels_csv:
        table = pd.read_csv(labels_csv, dtype=str, keep_default_na=False)
        if "label" not in table.columns:
            raise SystemExit(f"{labels_csv} needs a 'label' column (has {list(table.columns)})")
        unknown = sorted(set(present) - set(table["label"]))
        if unknown:
            raise SystemExit(f"Labels missing from {labels_csv}: {unknown}")
        names = dict(zip(table["label"], table["display"])) if "display" in table.columns else {}
        names = {k: v.replace("\\n", " ") for k, v in names.items() if v}
        d["label"] = d["label"].map(lambda v: names.get(v, v))
        order = [names.get(v, v) for v in table["label"] if names.get(v, v) in set(d["label"])]
        colors = figstyle.categorical_colors(order)
        if "color" in table.columns:
            colors.update({names.get(k, k): v for k, v in zip(table["label"], table["color"]) if v})
        return order, colors
    return present, figstyle.categorical_colors(present)


def place_labels(ax, d, order, small_fraction=0.015):
    """Large populations are labelled on their cluster; small ones are pushed outward with a thin leader
    line. Labels are then nudged apart. Check the result by eye: a few may need --labels-csv short names."""
    span = np.array([d.x.max() - d.x.min(), d.y.max() - d.y.min()])
    centre = np.array([d.x.median(), d.y.median()])
    anchors, positions = {}, {}
    for lab in order:
        s = d[d.label == lab]
        if len(s) < 10:
            continue
        cx, cy = s.x.median(), s.y.median()
        k = np.argmin((s.x - cx) ** 2 + (s.y - cy) ** 2)
        anchors[lab] = np.array([s.x.iloc[k], s.y.iloc[k]])
        positions[lab] = anchors[lab].copy()
        if len(s) < small_fraction * len(d):
            direction = (anchors[lab] - centre) / span
            positions[lab] = anchors[lab] + direction / (np.linalg.norm(direction) + 1e-9) * 0.09 * span
    labs = list(positions)
    P = np.array([positions[lab] for lab in labs], dtype=float)
    for _ in range(200):
        moved = False
        for i in range(len(P)):
            for j in range(i + 1, len(P)):
                gap = (P[i] - P[j]) / span
                if abs(gap[0]) < 0.16 and abs(gap[1]) < 0.035:
                    push = np.array([0, 0.012 if gap[1] >= 0 else -0.012]) * span
                    P[i] += push
                    P[j] -= push
                    moved = True
        if not moved:
            break
    for lab, (x, y) in zip(labs, P):
        if np.linalg.norm((np.array([x, y]) - anchors[lab]) / span) > 0.02:
            ax.plot([anchors[lab][0], x], [anchors[lab][1], y], color="#4d4d4d", lw=0.35, zorder=3)
        text = ax.text(x, y, lab, fontsize=5.5, ha="center", va="center", zorder=4)
        text.set_path_effects([path_effects.withStroke(linewidth=1.5, foreground="white")])


def panel_umap(ax, d, order, colors, size, labels=True):
    ax.scatter(d.x, d.y, c=d.label.map(colors), s=size, linewidths=0, rasterized=True)
    figstyle.clean_embedding(ax, arrows=True, x=d.x, y=d.y)
    if labels:
        place_labels(ax, d, order)


def panel_composition(ax, d, order, colors, samples, groups):
    table = pd.crosstab(d["sample"], d["label"], normalize="index").mul(100).reindex(index=samples, columns=order).fillna(0)
    counts = d["sample"].value_counts().reindex(samples)
    y = np.arange(len(samples))[::-1]
    left = np.zeros(len(samples))
    for lab in order:
        ax.barh(y, table[lab], left=left, color=colors[lab], height=0.78, linewidth=0)
        left += table[lab].values
    ax.set_yticks(y, samples)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("Cells (%)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    for yi, s in zip(y, samples):
        ax.text(101.5, yi, f"{counts[s]:,}", va="center", ha="left", fontsize=5)
    ax.text(101.5, 1.0, "n", transform=ax.get_xaxis_transform(), va="bottom", ha="left", fontsize=5, style="italic")
    return table, counts


def panel_groups(axes, d, colors, size, groups):
    for ax, g in zip(axes, groups):
        s = d[d.group == g]
        ax.scatter(d.x, d.y, c="#e6e6e6", s=size * 0.6, linewidths=0, rasterized=True)
        ax.scatter(s.x, s.y, c=s.label.map(colors), s=size * 1.2, linewidths=0, rasterized=True)
        figstyle.clean_embedding(ax)
        ax.set_title(f"{g} (n = {len(s):,})", fontsize=5.5, pad=2)


def panel_genes(axes, d, genes, size):
    for ax, g in zip(axes, genes):
        o = d.sort_values(f"gene:{g}")
        v = o[f"gene:{g}"].to_numpy()
        vmax = np.quantile(v[v > 0], 0.99) if (v > 0).any() else 1
        ax.scatter(o.x, o.y, c=v, cmap=CMAP, vmin=0, vmax=vmax, s=size * 0.6, linewidths=0, rasterized=True)
        figstyle.clean_embedding(ax)
        ax.set_title(g, fontsize=5.5, style="italic", pad=1.5)
    for ax in axes[len(genes):]:
        ax.axis("off")


def gene_colorbar(fig, cax):
    bar = fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP), cax=cax, orientation="horizontal")
    bar.set_ticks([0, 1])
    bar.set_ticklabels(["0", "max"])
    bar.outline.set_linewidth(0.4)
    bar.ax.tick_params(labelsize=5, length=1.5, width=0.4)
    cax.set_title("Expression (log CP10k, scaled per gene)", fontsize=5, pad=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Annotated .h5ad with a UMAP")
    parser.add_argument("--output-dir", required=True, help="Where to write the figure (created if needed)")
    parser.add_argument("--label-col", default="cell_type_v2", help="obs column with cell types (default cell_type_v2)")
    parser.add_argument("--sample-col", default="sample", help="obs column with samples (default sample)")
    parser.add_argument("--group-col", default=None, help="obs column with groups for panel c (e.g. condition)")
    parser.add_argument("--genes", nargs="*", default=None, help="Marker genes for panel d (up to 12 fit)")
    parser.add_argument("--embedding", default="X_umap", help="obsm key to plot (default X_umap)")
    parser.add_argument("--labels-csv", default=None, help="CSV to rename, order and colour labels: label,display,color")
    parser.add_argument("--exclude-labels", nargs="*", default=None, help="Labels to leave out, e.g. 'Doublet'")
    parser.add_argument("--sample-order", nargs="*", default=None, help="Sample order in panel b (default: by group, name)")
    parser.add_argument("--group-order", nargs="*", default=None, help="Group order in panel c (default: sorted)")
    parser.add_argument("--pool-samples", nargs="*", default=None, metavar="A+B=NAME",
                        help="Count these samples as one, e.g. 'S1+S2=S1 and S2' when they are not independent")
    args = parser.parse_args()

    out = pathlib.Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    figstyle.apply_style()
    d, genes = load_table(args)
    genes = genes[:12]
    order, colors = order_and_colors(d, args.labels_csv)

    for spec in args.pool_samples or []:
        if "=" not in spec:
            raise SystemExit(f"--pool-samples needs the form A+B=NAME, got {spec!r}")
        members, name = spec.split("=", 1)
        members = members.split("+")
        unknown = sorted(set(members) - set(d["sample"]))
        if unknown:
            raise SystemExit(f"--pool-samples: {unknown} not in the data")
        d["sample"] = d["sample"].replace({m: name for m in members})

    groups = args.group_order or (sorted(d["group"].unique()) if args.group_col else [])
    if args.sample_order:
        samples = [s for s in args.sample_order if s in set(d["sample"])]
    elif args.group_col:
        sample_group = d.groupby("sample")["group"].first()
        samples = sorted(sample_group.index, key=lambda s: (groups.index(sample_group[s]) if sample_group[s] in groups else 99, s))
    else:
        samples = sorted(d["sample"].unique())

    d = d.iloc[np.random.default_rng(0).permutation(len(d))].reset_index(drop=True)  # no group painted on top
    size = 0.25 if len(d) > 50000 else figstyle.point_size(len(d))

    # ---- full figure: rows are laid out in mm so an absent panel leaves no gap
    n_gene_rows = int(np.ceil(len(genes) / 6))
    heights = {"top": 88, "groups": 42 if groups else 0, "genes": 28 * n_gene_rows + (10 if genes else 0)}
    total = sum(heights.values()) + 4
    fig = figstyle.new_figure(figstyle.DOUBLE, total)
    y_top = total - 2 - heights["top"]
    ax_a = figstyle.add_axes_mm(fig, 2, y_top, 112, heights["top"] - 2)
    ax_b = figstyle.add_axes_mm(fig, 136, y_top + 14, 36, heights["top"] - 22)
    panel_umap(ax_a, d, order, colors, size)
    table, counts = panel_composition(ax_b, d, order, colors, samples, groups)
    fig.text(2 / figstyle.DOUBLE, (y_top + heights["top"] - 2) / total, "a", fontsize=8, fontweight="bold", va="top")
    fig.text(118 / figstyle.DOUBLE, (y_top + heights["top"] - 2) / total, "b", fontsize=8, fontweight="bold", va="top")
    y_cursor = y_top
    letters = iter("cd")
    if groups:
        y_cursor -= heights["groups"]
        n = max(len(groups), 1)
        w = min(45, 176 / n)
        axes_c = [figstyle.add_axes_mm(fig, 2 + i * (w + 1), y_cursor + 2, w, heights["groups"] - 8) for i in range(n)]
        panel_groups(axes_c, d, colors, size, groups)
        fig.text(2 / figstyle.DOUBLE, (y_cursor + heights["groups"] - 1) / total, next(letters), fontsize=8,
                 fontweight="bold", va="top")
    if genes:
        y_cursor -= heights["genes"]
        axes_d = [figstyle.add_axes_mm(fig, 2 + (i % 6) * 30, y_cursor + 10 + (n_gene_rows - 1 - i // 6) * 28, 28, 23)
                  for i in range(len(genes))]
        panel_genes(axes_d, d, genes, size)
        gene_colorbar(fig, figstyle.add_axes_mm(fig, 70, y_cursor + 3, 40, 1.5))
        fig.text(2 / figstyle.DOUBLE, (y_cursor + heights["genes"] - 1) / total, next(letters), fontsize=8,
                 fontweight="bold", va="top")
    legend = (f"Overview of {len(d):,} cells from {len(samples)} samples. a, UMAP coloured by cell type. b, Cell-type "
              "composition per sample; n, cells per sample."
              + (" c, The same UMAP split by group." if groups else "")
              + (f" {'d' if groups else 'c'}, Marker genes: log1p(counts per 10,000), colour capped at each gene's "
                 "99th percentile among expressing cells." if genes else ""))
    figstyle.save_figure(fig, out / "Fig_atlas", legend=legend, check=False)

    # ---- single panels at one-column width
    panels = out / "panels"
    fig, ax = plt.subplots(figsize=(figstyle.SINGLE * figstyle.MM, 80 * figstyle.MM))
    panel_umap(ax, d, order, colors, size)
    figstyle.save_figure(fig, panels / "panel_umap_labelled", check=False)
    fig, ax = plt.subplots(figsize=(figstyle.SINGLE * figstyle.MM, 80 * figstyle.MM))
    panel_umap(ax, d, order, colors, size, labels=False)
    handles = [plt.Line2D([], [], marker="o", ls="", color=colors[lab], markersize=3) for lab in order]
    ax.legend(handles, order, loc="center left", bbox_to_anchor=(1.0, 0.5), handletextpad=0.2, labelspacing=0.35)
    figstyle.save_figure(fig, panels / "panel_umap_with_legend", check=False)
    fig, ax = plt.subplots(figsize=(60 * figstyle.MM, (10 + 5 * len(samples)) * figstyle.MM))
    panel_composition(ax, d, order, colors, samples, groups)
    figstyle.save_figure(fig, panels / "panel_composition", check=False)

    # ---- source data
    source = out / "source_data"
    source.mkdir(exist_ok=True)
    gene_cols = [f"gene:{g}" for g in genes]
    cells = d[["cell", "x", "y", "label", "sample"] + (["group"] if groups else []) + gene_cols]
    cells.rename(columns={"x": "UMAP1", "y": "UMAP2"}).to_csv(source / "cells.csv.gz", index=False)
    table.round(3).assign(n_cells=counts).to_csv(source / "composition_percent.csv")
    d["label"].value_counts().reindex(order).rename("cells").to_csv(source / "cell_type_counts.csv")
    print(f"Wrote {out / 'Fig_atlas'}.* ({len(d):,} cells, {len(order)} cell types, {len(samples)} samples)")


if __name__ == "__main__":
    main()

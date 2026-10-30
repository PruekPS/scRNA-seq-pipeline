"""Step 17 (optional): compare a per-sample value between groups, Nature style.

One panel per value: bar = group mean, dot = one sample. The unit of comparison is the
sample (mouse, patient), never the cell. Cells from one sample are not independent.
A test is drawn only if you ask (--test), and only where it is possible: exact
Mann-Whitney for 2 groups, Kruskal-Wallis with Dunn's comparisons for 3 or more. With
few samples a test may be unable to reach significance: an exact test on 2 vs 2 samples
can never give P below 0.33. The legend file records n and the test.

Two ways to give the values:
  --values-csv  a table with one row per sample: sample column, group column, value columns
  --input       an .h5ad: for each sample, the percent of cells with each label in --labels
                (within the cells matching --within-col/--within-values, if given)

Usage:
    python 17_group_figure.py --values-csv scores.csv --sample-col sample --group-col condition \
        --values "Effector Treg" "Tissue Treg" --ylabel "% of Treg" --output figures/treg
    python 17_group_figure.py --input 04b_annotated.h5ad --sample-col sample --group-col condition \
        --label-col cell_type_v2 --labels "B cell" "NK cell" --ylabel "% of cells" --output figures/composition
"""
import argparse
import pathlib

import numpy as np
import pandas as pd

import figstyle


def table_from_h5ad(args):
    import anndata as ad
    obs = ad.read_h5ad(args.input, backed="r").obs
    for col in (args.sample_col, args.group_col, args.label_col):
        if col is None or col not in obs.columns:
            raise SystemExit(f"Column {col!r} not in adata.obs (available: {list(obs.columns)})")
    if not args.labels:
        raise SystemExit("--input needs --labels")
    if args.within_col:
        if args.within_col not in obs.columns or not args.within_values:
            raise SystemExit("--within-col needs an obs column and --within-values")
        obs = obs[obs[args.within_col].astype(str).isin(args.within_values)]
    percent = pd.crosstab(obs[args.sample_col].astype(str), obs[args.label_col].astype(str), normalize="index") * 100
    missing = sorted(set(args.labels) - set(percent.columns))
    if missing:
        print(f"  note: labels with no cells, shown as 0: {missing}")
    percent = percent.reindex(columns=args.labels, fill_value=0.0)
    groups = obs.groupby(args.sample_col, observed=True)[args.group_col].first().astype(str)
    groups.index = groups.index.astype(str)
    return percent.assign(**{args.group_col: groups.reindex(percent.index)}).reset_index(names=args.sample_col), args.labels


def draw_panel(ax, values, order, colors, ylabel, test):
    positions = {g: i for i, g in enumerate(order)}
    y_max = max(float(np.nanmax(v)) for v in values.values() if len(v))
    y_min = min(0.0, min(float(np.nanmin(v)) for v in values.values() if len(v)))
    for g, i in positions.items():
        v = np.asarray(values[g], float)
        if len(v):
            ax.bar(i, v.mean(), 0.62, color=colors[g], alpha=0.25, linewidth=0, zorder=1)
            jitter = np.linspace(-0.14, 0.14, len(v)) if len(v) > 1 else np.zeros(1)
            ax.scatter(i + jitter, v, s=7, color=colors[g], edgecolor="white", linewidth=0.3, zorder=3, clip_on=False)
    ax.set_xticks(range(len(order)), order, rotation=45 if max(len(g) for g in order) > 6 else 0,
                  ha="right" if max(len(g) for g in order) > 6 else "center", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_xlim(-0.55, len(order) - 0.45)
    ax.set_ylabel(ylabel)
    result = figstyle.group_test(values) if test else None
    top = y_max
    if result:
        drawn = {k: p for k, p in result["pairs"].items() if len(order) == 2 or p < 0.05}
        top = figstyle.add_stat_brackets(ax, drawn, positions, y_max, (y_max - y_min) or 1.0)
    ax.set_ylim(y_min, max(top, y_max * 1.05) if y_max > 0 else 1)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--values-csv", default=None, help="Table with one row per sample")
    parser.add_argument("--input", default=None, help="Annotated .h5ad, as an alternative to --values-csv")
    parser.add_argument("--sample-col", default="sample", help="Column with the sample name (default sample)")
    parser.add_argument("--group-col", required=True, help="Column with the group")
    parser.add_argument("--values", nargs="*", default=None, help="--values-csv: columns to plot, one panel each")
    parser.add_argument("--label-col", default=None, help="--input: obs column with cell types")
    parser.add_argument("--labels", nargs="*", default=None, help="--input: cell types to plot, one panel each")
    parser.add_argument("--within-col", default=None, help="--input: compute percent only within these cells")
    parser.add_argument("--within-values", nargs="*", default=None, help="--input: values of --within-col")
    parser.add_argument("--group-order", nargs="*", default=None, help="Group order and colours (default: sorted)")
    parser.add_argument("--colors", nargs="*", default=None, help="Hex colours, one per group in --group-order")
    parser.add_argument("--ylabel", default="Value", help="y-axis label; a list of one per panel is not supported")
    parser.add_argument("--test", action="store_true", help="Add a test across groups (see the top of this file)")
    parser.add_argument("--output", required=True, help="Output path without extension; its folder is created")
    args = parser.parse_args()
    if bool(args.values_csv) == bool(args.input):
        raise SystemExit("Give exactly one of --values-csv or --input")

    if args.values_csv:
        table = pd.read_csv(args.values_csv)
        value_cols = args.values or []
        if not value_cols:
            raise SystemExit("--values-csv needs --values")
        missing = [c for c in [args.sample_col, args.group_col] + value_cols if c not in table.columns]
        if missing:
            raise SystemExit(f"Columns not in {args.values_csv}: {missing} (has {list(table.columns)})")
    else:
        table, value_cols = table_from_h5ad(args)
    table[args.group_col] = table[args.group_col].astype(str)
    order = args.group_order or sorted(table[args.group_col].unique())
    unknown = sorted(set(table[args.group_col]) - set(order))
    if unknown:
        raise SystemExit(f"Groups missing from --group-order: {unknown}")
    if args.colors and len(args.colors) != len(order):
        raise SystemExit(f"{len(args.colors)} --colors for {len(order)} groups")
    colors = dict(zip(order, args.colors)) if args.colors else figstyle.categorical_colors(order)

    figstyle.apply_style()
    n = len(value_cols)
    ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))
    panel_w = max(32, 8 * len(order) + 14)
    fig, axes = figstyle.subplots(nrows, ncols, w_mm=min(183, panel_w * ncols), h_mm=46 * nrows)
    notes = []
    for i, (ax, col) in enumerate(zip(axes.flat, value_cols)):
        values = {g: table.loc[table[args.group_col] == g, col].dropna().to_numpy(float) for g in order}
        if n > 1:
            ax.set_title(col)
            figstyle.panel_letter(ax, "abcdefghijklmnopqrstuvwxyz"[i], dx=-0.28, dy=1.03)
        result = draw_panel(ax, values, order, colors, args.ylabel, args.test)
        if args.test:
            notes.append(f"{col}: {result['test']}, P = {result['p']:.3g}." if result else f"{col}: no test possible (too few samples).")
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    fig.tight_layout()

    sizes = table.groupby(args.group_col)[args.sample_col].nunique().reindex(order).fillna(0).astype(int)
    legend = ("Bars, group mean; dots, individual samples. n per group: "
              + ", ".join(f"{g} {sizes[g]}" for g in order) + ". " + " ".join(notes))
    output = pathlib.Path(args.output)
    figstyle.save_figure(fig, output, legend=legend, check=False,
                         source_data=table[[args.sample_col, args.group_col] + value_cols].set_index(args.sample_col))
    print(f"Wrote {output}.*")
    print(legend)


if __name__ == "__main__":
    main()

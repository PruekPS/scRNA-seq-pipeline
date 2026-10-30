"""Shared Nature-style figure settings and helpers, used by every plotting step.

Style (Nature's figure guidelines): Arial or Helvetica at 5-7 pt, bold lowercase
panel letters, 0.5 pt lines, widths of 89 mm (single column), 120 mm (1.5 column)
or 183 mm (double column), editable text in vector files, cell scatter layers
rasterized inside them, no titles or explanatory sentences inside the figure --
those belong in the legend. Colours are colour-blind safe (Okabe-Ito first).

Every figure saved with save_figure() becomes a small package:
    name.pdf, name.svg      vector, editable text
    name.png, name.tiff     600 dpi
    name_legend.txt         what is plotted, the n, the test (if legend= is given)
    name_source_data.csv    the numbers behind the figure (if source_data= is given)

Arial is not installed on most clusters. The fallback is Liberation Sans, which
has the same letter widths, then DejaVu Sans. To redraw with real Arial, run the
plotting step on a computer that has it, or install the font on the cluster.

Set the environment variable FIG_FORMATS to write fewer files during a test run,
e.g. FIG_FORMATS=png.

Not a pipeline step: nothing to run. Other scripts do `import figstyle`.
"""
import itertools
import math
import os
import pathlib
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

MM = 1 / 25.4  # inches per mm
SINGLE, ONE_AND_HALF, DOUBLE = 89, 120, 183  # Nature column widths in mm
INK = "#1a1a1a"
FONT_PREFERENCE = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
# Okabe-Ito (Nature Methods 2011) first, then lighter partners for many categories.
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#4D4D4D",
           "#882255", "#44AA99", "#AA4499", "#117733", "#999933", "#88CCEE", "#CC6677", "#DDCC77",
           "#332288", "#B3B3B3", "#661100", "#6699CC"]
SEQUENTIAL = "viridis"


def apply_style():
    """Set Nature-style defaults for all matplotlib figures in this process."""
    installed = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in FONT_PREFERENCE if f in installed), "DejaVu Sans")
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [font] + FONT_PREFERENCE,
        "mathtext.fontset": "custom", "mathtext.rm": font, "mathtext.it": f"{font}:italic",
        "mathtext.bf": f"{font}:bold", "mathtext.default": "regular",
        "font.size": 6, "axes.titlesize": 6, "axes.titleweight": "normal", "axes.labelsize": 6,
        "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5, "legend.title_fontsize": 5.5,
        "figure.titlesize": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.major.size": 2, "ytick.major.size": 2,
        "xtick.major.pad": 1.5, "ytick.major.pad": 1.5, "axes.labelpad": 2, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": False, "legend.frameon": False, "lines.linewidth": 0.75,
        "patch.linewidth": 0.5, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.dpi": 600, "figure.dpi": 150, "text.color": INK, "axes.labelcolor": INK,
        "axes.edgecolor": INK, "xtick.color": INK, "ytick.color": INK,
    })
    return font


def new_figure(w_mm=SINGLE, h_mm=60):
    return plt.figure(figsize=(w_mm * MM, h_mm * MM))


def subplots(nrows=1, ncols=1, w_mm=SINGLE, h_mm=60, **kwargs):
    kwargs.setdefault("squeeze", False)
    return plt.subplots(nrows, ncols, figsize=(w_mm * MM, h_mm * MM), **kwargs)


def add_axes_mm(fig, x, y, w, h):
    """Axes placed in mm from the bottom-left corner of the figure."""
    W, H = fig.get_size_inches() / MM
    return fig.add_axes([x / W, y / H, w / W, h / H])


def panel_letter(ax, letter, dx=-0.15, dy=1.02):
    """Bold lowercase panel letter, placed in axes-fraction coordinates."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", ha="left", va="bottom")


def categorical_colors(names):
    """name -> colour, in the order given; colours repeat only past 20 categories."""
    return {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(names)}


def rasterize_points(ax_or_fig):
    """Rasterize scatter layers so a vector file with many cells stays small and editable."""
    axes = ax_or_fig.axes if hasattr(ax_or_fig, "axes") and not hasattr(ax_or_fig, "transData") else [ax_or_fig]
    for ax in axes:
        for coll in ax.collections:
            coll.set_rasterized(True)


def clean_embedding(ax, arrows=False, x=None, y=None):
    """Remove the frame around a UMAP; arrows=True draws short UMAP1/UMAP2 axis arrows
    (x, y = the plotted coordinates)."""
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(""); ax.set_ylabel("")
    for s in ax.spines.values():
        s.set_visible(False)
    if arrows and x is not None:
        x0, x1, y0, y1 = np.min(x), np.max(x), np.min(y), np.max(y)
        L = 0.12 * max(x1 - x0, y1 - y0)
        kw = dict(arrowstyle="-|>,head_length=0.25,head_width=0.12", lw=0.5, color=INK, shrinkA=0, shrinkB=0)
        ax.annotate("", xy=(x0 + L, y0), xytext=(x0, y0), arrowprops=kw)
        ax.annotate("", xy=(x0, y0 + L), xytext=(x0, y0), arrowprops=kw)
        ax.text(x0 + L / 2, y0 - 0.02 * (y1 - y0), "UMAP1", ha="center", va="top", fontsize=5)
        ax.text(x0 - 0.02 * (x1 - x0), y0 + L / 2, "UMAP2", ha="right", va="center", rotation=90, fontsize=5)


def point_size(n_cells):
    """Scatter dot area (pt^2) that keeps a UMAP readable on a small Nature panel."""
    return float(min(6.0, max(0.3, 3000 / max(n_cells, 1))))


def umap(adata, color, ax, title="", legend_loc="right margin", **kwargs):
    """scanpy UMAP panel with the frame off, a suitable dot size and rasterized dots."""
    import scanpy as sc  # imported here so figures that need no AnnData load faster
    kwargs.setdefault("size", point_size(adata.n_obs))
    kwargs.setdefault("legend_fontsize", 5)
    sc.pl.umap(adata, color=color, ax=ax, show=False, frameon=False, title=title, legend_loc=legend_loc, **kwargs)
    rasterize_points(ax)


def marker_dotplot(adata, markers, groupby, path, legend):
    """scanpy dotplot of {group: genes}, at Nature width with thin lines and italic gene names."""
    import scanpy as sc
    present = {k: [g for g in v if g in adata.raw.var_names] for k, v in markers.items()}
    present = {k: v for k, v in present.items() if v}
    n_groups = adata.obs[groupby].nunique()
    dp = sc.pl.dotplot(adata, present, groupby=groupby, use_raw=True, show=False, return_fig=True, dot_max=1.0,
                       figsize=(DOUBLE * MM, (30 + 4.5 * n_groups) * MM))
    dp.style(dot_edge_lw=0.2).legend(width=0.8)
    dp.make_figure()
    fig = plt.gcf()
    for name, ax in dp.ax_dict.items():
        for patch in ax.patches:
            patch.set_linewidth(0.5)
    for t in dp.ax_dict["mainplot_ax"].get_xticklabels():
        t.set_fontstyle("italic")
    save_figure(fig, path, legend=legend, check=False)


def _visible_texts(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    hidden = set()  # tick labels whose tick lies outside the view limits are never drawn
    for ax in fig.axes:
        for k, (labels, lim) in enumerate(((ax.get_xticklabels(), ax.get_xlim()), (ax.get_yticklabels(), ax.get_ylim()))):
            lo, hi = min(lim), max(lim)
            for t in labels:
                if not ax.axison or not lo - 1e-9 <= t.get_position()[k] <= hi + 1e-9:
                    hidden.add(id(t))
    out = []
    for t in fig.findobj(matplotlib.text.Text):
        if id(t) in hidden or not t.get_visible() or not t.get_text().strip():
            continue
        try:
            box = t.get_window_extent(renderer)
        except Exception:
            continue
        if box.width >= 1 and box.height >= 1:
            out.append((t, box))
    return out, renderer


def _text_polygon(t, box, renderer, shrink=0.5):
    """The real (rotated) rectangle a text occupies, in display pixels."""
    angle = math.radians(t.get_rotation() % 180)
    if abs(math.sin(angle)) < 1e-6 or abs(math.cos(angle)) < 1e-6:
        return np.array([[box.x0 + shrink, box.y0 + shrink], [box.x1 - shrink, box.y0 + shrink],
                         [box.x1 - shrink, box.y1 - shrink], [box.x0 + shrink, box.y1 - shrink]])
    rotation = t.get_rotation(); t.set_rotation(0)
    flat = t.get_window_extent(renderer); t.set_rotation(rotation)
    w, h = flat.width / 2 - shrink, flat.height / 2 - shrink
    cx, cy = (box.x0 + box.x1) / 2, (box.y0 + box.y1) / 2
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[cx + dx * c - dy * s, cy + dx * s + dy * c] for dx, dy in ((-w, -h), (w, -h), (w, h), (-w, h))])


def _convex_overlap(p, q):
    """Separating-axis test for two convex polygons."""
    for poly in (p, q):
        for i in range(len(poly)):
            edge = poly[(i + 1) % len(poly)] - poly[i]
            normal = np.array([-edge[1], edge[0]])
            a, b = p @ normal, q @ normal
            if a.max() <= b.min() or b.max() <= a.min():
                return False
    return True


def layout_problems(fig, min_pt=4.9):
    """Plain-language list of overlapping text and text too small for Nature's 5 pt minimum.
    Text running off the figure is not listed: saving uses a tight box that keeps it."""
    texts, renderer = _visible_texts(fig)
    shapes = [(t, box, _text_polygon(t, box, renderer)) for t, box in texts]
    problems = []
    for (t1, b1, p1), (t2, b2, p2) in itertools.combinations(shapes, 2):
        if b1.overlaps(b2) and _convex_overlap(p1, p2):
            problems.append(f"overlap: {t1.get_text()!r} and {t2.get_text()!r}")
    small = sorted({t.get_text() for t, _ in texts if t.get_fontsize() < min_pt})
    if small:
        problems.append(f"text under {min_pt:g} pt: {small[:5]}")
    return problems


def save_figure(fig, path, legend=None, source_data=None, check=True):
    """Write the figure package described at the top of this file. `path` has no extension;
    its folder is created if needed. Layout problems are printed, never fatal."""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if check:
        for p in layout_problems(fig):
            print(f"  layout note ({path.name}): {p}")
    formats = [f.strip() for f in os.environ.get("FIG_FORMATS", "pdf,svg,png,tiff").split(",") if f.strip()]
    for ext in formats:
        kwargs = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        fig.savefig(f"{path}.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.02, **kwargs)
    if legend:
        pathlib.Path(f"{path}_legend.txt").write_text(textwrap.dedent(legend).strip() + "\n")
    if source_data is not None:
        source_data.to_csv(f"{path}_source_data.csv")
    plt.close(fig)


# ---------------------------------------------------------------- statistics used by figures
def fmt_p(p):
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "n/a"
    if p < 0.0001:
        return "P < 0.0001"
    if p < 0.001:
        return f"P = {p:.4f}"
    if p >= 0.9995:
        return "P > 0.99"
    return f"P = {p:.3g}" if p < 0.1 else f"P = {p:.2f}"


def _n_assignments(ns):
    out = math.factorial(sum(ns))
    for k in ns:
        out //= math.factorial(k)
    return out


def kruskal_exact(groups, max_enum=400000):
    """Kruskal-Wallis H with an exact permutation P when the design is small enough to
    enumerate, otherwise the asymptotic chi-square P. Returns (H, P, method)."""
    groups = [np.asarray(g, float) for g in groups]
    ns = [len(g) for g in groups]
    allv = np.concatenate(groups)
    if len(set(allv)) < 2:
        return 0.0, 1.0, "exact"
    H = stats.kruskal(*groups).statistic
    ranks = stats.rankdata(allv)
    N = len(allv)
    tie = 1 - (np.sum([(c ** 3 - c) for c in np.unique(allv, return_counts=True)[1]]) / (N ** 3 - N))

    def h_of(labels):
        s = sum(ranks[labels == k].sum() ** 2 / n for k, n in enumerate(ns))
        return (12 / (N * (N + 1)) * s - 3 * (N + 1)) / tie

    if _n_assignments(ns) > max_enum:
        return H, float(stats.chi2.sf(H, len(ns) - 1)), "asymptotic"
    h0 = h_of(np.repeat(np.arange(len(ns)), ns))
    count = total = 0

    def assign(remaining, k, labels):
        nonlocal count, total
        if k == len(ns) - 1:
            labels[list(remaining)] = k
            total += 1
            count += h_of(labels) >= h0 - 1e-9
            return
        for comb in itertools.combinations(sorted(remaining), ns[k]):
            next_labels = labels.copy()
            next_labels[list(comb)] = k
            assign(remaining - set(comb), k + 1, next_labels)

    assign(set(range(N)), 0, np.full(N, -1))
    return H, count / total, "exact"


def dunn(groups, names):
    """Dunn's multiple comparisons after Kruskal-Wallis, Bonferroni-adjusted over all pairs.
    Returns {(a, b): (z, p_raw, p_adjusted)}."""
    groups = [np.asarray(g, float) for g in groups]
    allv = np.concatenate(groups)
    N = len(allv)
    ranks = stats.rankdata(allv)
    _, counts = np.unique(allv, return_counts=True)
    var0 = N * (N + 1) / 12 - np.sum(counts ** 3 - counts) / (12 * (N - 1))
    start = np.cumsum([0] + [len(g) for g in groups])
    mean_rank = [ranks[start[i]:start[i + 1]].mean() for i in range(len(groups))]
    pairs = list(itertools.combinations(range(len(groups)), 2))
    out = {}
    for i, j in pairs:
        se = math.sqrt(var0 * (1 / len(groups[i]) + 1 / len(groups[j])))
        z = (mean_rank[i] - mean_rank[j]) / se if se > 0 else 0.0
        p = 2 * stats.norm.sf(abs(z))
        out[(names[i], names[j])] = (z, p, min(1.0, p * len(pairs)))
    return out


def mannwhitney_exact(a, b, max_enum=2_000_000):
    """Two-sided Mann-Whitney U; exact permutation P (ties handled) when feasible. Returns (U, P, method)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    U = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    allv = np.concatenate([a, b])
    n1, N = len(a), len(allv)
    if math.comb(N, n1) > max_enum:
        return U, stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic").pvalue, "asymptotic"
    ranks = stats.rankdata(allv)
    mu = n1 * (N + 1) / 2
    observed = abs(ranks[:n1].sum() - mu)
    count = total = 0
    for c in itertools.combinations(range(N), n1):
        total += 1
        count += abs(ranks[list(c)].sum() - mu) >= observed - 1e-9
    return U, count / total, "exact"


def group_test(values_by_group):
    """values_by_group: {group: values}. Kruskal-Wallis + Dunn for 3+ groups, exact Mann-Whitney
    for 2. Groups with one value cannot be tested. Returns a dict, or None if no test is possible."""
    testable = {g: np.asarray(v, float) for g, v in values_by_group.items() if len(v) >= 2}
    if len(testable) < 2 or sum(len(v) for v in testable.values()) < 4:
        return None
    names, vals = list(testable), list(testable.values())
    if len(names) == 2:
        U, p, method = mannwhitney_exact(*vals)
        return {"test": f"two-sided Mann-Whitney U test ({method} P)", "U": U, "p": p, "groups": names,
                "pairs": {(names[0], names[1]): p}}
    H, p, method = kruskal_exact(vals)
    pairs = {k: v[2] for k, v in dunn(vals, names).items()}
    return {"test": f"Kruskal-Wallis ({method} P) with Dunn's comparisons, Bonferroni-adjusted", "H": H, "p": p,
            "groups": names, "pairs": pairs}


def add_stat_brackets(ax, pairs, positions, y_top, y_span):
    """Draw a bracket and P value for each {(group_a, group_b): p}."""
    ordered = sorted(pairs.items(), key=lambda kv: abs(positions[kv[0][1]] - positions[kv[0][0]]))
    for k, ((a, b), p) in enumerate(ordered):
        y = y_top + (0.08 + 0.13 * k) * y_span
        h = 0.025 * y_span
        ax.plot([positions[a], positions[a], positions[b], positions[b]], [y, y + h, y + h, y], color=INK, lw=0.5,
                clip_on=False)
        ax.text((positions[a] + positions[b]) / 2, y + h + 0.01 * y_span, fmt_p(p).replace("P", r"$\mathit{P}$"),
                ha="center", va="bottom", fontsize=5)
    return y_top + (0.08 + 0.13 * len(ordered)) * y_span


__all__ = ["apply_style", "new_figure", "subplots", "add_axes_mm", "panel_letter", "categorical_colors",
           "rasterize_points", "clean_embedding", "point_size", "umap", "marker_dotplot", "layout_problems", "save_figure", "fmt_p", "kruskal_exact",
           "dunn", "mannwhitney_exact", "group_test", "add_stat_brackets", "MM", "SINGLE", "ONE_AND_HALF", "DOUBLE",
           "INK", "PALETTE", "SEQUENTIAL"]

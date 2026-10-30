"""Step 9 (optional): CellRank fate mapping without RNA velocity.

Uses CellRank 2 (Weiler, Lange, Klein, Pe'er & Theis, Nature Methods 2024,
doi:10.1038/s41592-024-02303-9; code: github.com/theislab/cellrank) built
on two velocity-free kernels:
  - CytoTRACEKernel: a pure-Python reimplementation, inside CellRank, of
    CytoTRACE's transcriptional-diversity score (Gulati et al., Science
    2020, doi:10.1126/science.aax0249) -- cells expressing more distinct
    genes score as less differentiated. No R dependency, no separate
    CytoTRACE install.
  - ConnectivityKernel: the existing neighbor graph, included to stabilize
    the directed CytoTRACE signal against noise (CellRank's own documented
    pattern: combine kernels as a weighted sum, e.g. `0.8 * ctk + 0.2 * ck`).
  - PseudotimeKernel (optional, --use-dpt-pseudotime): blends in diffusion
    pseudotime from scripts/08_trajectory_paga.py as a third signal, if
    that step was run first.

The estimator is GPCCA (Generalized Perron Cluster Cluster Analysis) --
CellRank's recommended estimator, which computes coarse-grained
"macrostates", automatically classifies some as terminal/initial states,
and computes each cell's fate probability toward each terminal state.

Deliberately NO RNA velocity in this script -- see
docs/reference/trajectory-design.md for why: standard Parse
split-pipe/Trailmaker output doesn't separate spliced/unspliced reads,
so true velocity needs raw FASTQ reprocessing, which is out of scope
until that's confirmed feasible and worth doing.

Tested on synthetic branching data (tests/synthetic_trajectory.py), which
proves the code runs, not that a real trajectory is meaningful: see
docs/reference/trajectory-design.md for what the experiment must contain.

Input: an annotated, clustered .h5ad -- either scripts/04b_better_annotation.py's
direct output, or scripts/08_trajectory_paga.py's output if using
--use-dpt-pseudotime.

Usage:
    python 09_cellrank_fate_mapping.py --input 04b_annotated.h5ad --output 09_cellrank.h5ad \
        --plot-dir cellrank_plots/ --annotation-col cell_type_v2 --n-states 8
    python 09_cellrank_fate_mapping.py --input 08_paga.h5ad --output 09_cellrank.h5ad \
        --use-dpt-pseudotime --annotation-col cell_type_v2 --n-states 8
"""
import argparse
import functools
import pathlib

import matplotlib.pyplot as plt
import cellrank as cr
import scanpy as sc

import figstyle


def build_combined_kernel(adata, args):
    print("Computing CytoTRACE kernel (pure Python, no R dependency)...")
    ctk = cr.kernels.CytoTRACEKernel(adata)
    ctk.compute_cytotrace(layer=args.cyto_layer)
    ctk.compute_transition_matrix(threshold_scheme=args.threshold_scheme)

    print("Computing connectivity kernel (stabilizes the directed signal)...")
    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()

    weights = {"cyto": args.cyto_weight, "conn": args.conn_weight}
    kernels = {"cyto": ctk, "conn": ck}

    if args.use_dpt_pseudotime:
        if "dpt_pseudotime" not in adata.obs:
            raise SystemExit(
                "--use-dpt-pseudotime given but adata.obs['dpt_pseudotime'] not found -- "
                "run scripts/08_trajectory_paga.py first and pass its output as --input."
            )
        print("Computing pseudotime kernel from dpt_pseudotime (scripts/08_trajectory_paga.py)...")
        ptk = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
        ptk.compute_transition_matrix(threshold_scheme=args.threshold_scheme)
        weights["dpt"] = args.dpt_weight
        kernels["dpt"] = ptk

    total = sum(weights.values())
    if total <= 0:
        raise SystemExit(f"Kernel weights must sum to a positive number, got {weights}")
    weights = {k: v / total for k, v in weights.items()}
    print(f"Combining kernels with normalized weights: {weights}")

    # functools.reduce, not sum(): sum() starts from the number 0, which CellRank kernels cannot add to.
    return functools.reduce(lambda a, b: a + b, (weights[name] * kernels[name] for name in kernels))


def save_plots(adata, g, annotation_col, plot_dir):
    plot_dir = pathlib.Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # CellRank draws its own figures, so a plot it cannot make (e.g. no initial state found) is
    # reported and skipped rather than losing hours of computed results.
    plots = {
        "eigenvalue_spectrum": lambda: g.plot_spectrum(real_only=True),
        "macrostates_all": lambda: g.plot_macrostates(which="all", legend_loc="right", s=6),
        "macrostates_terminal": lambda: g.plot_macrostates(which="terminal", legend_loc="right", s=6),
        "macrostates_initial": lambda: g.plot_macrostates(which="initial", s=6),
        "coarse_transition_matrix": lambda: g.plot_coarse_T(annotate=False),
        "fate_probabilities": lambda: g.plot_fate_probabilities(same_plot=True),
        "circular_projection": lambda: cr.pl.circular_projection(adata, keys=[annotation_col], legend_loc="right"),
    }
    for name, draw in plots.items():
        try:
            draw()
            figstyle.rasterize_points(plt.gcf())
            figstyle.save_figure(plt.gcf(), plot_dir / name, check=False)
        except Exception as e:
            print(f"{name} plot skipped: {e}")
        plt.close("all")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to annotated, clustered .h5ad")
    parser.add_argument("--output", required=True, help="Path to write the .h5ad with fate probabilities")
    parser.add_argument("--plot-dir", default="cellrank_plots", help="Directory for CellRank plots")
    parser.add_argument("--drivers-csv", default=None, dest="drivers_csv",
                         help="Optional path to write per-lineage driver genes as CSV")
    parser.add_argument("--annotation-col", default="cell_type_v2", dest="annotation_col",
                         help="adata.obs column macrostates are compared against (default: cell_type_v2)")
    parser.add_argument("--n-states", type=int, default=8, dest="n_states",
                         help="Number of macrostates for GPCCA to compute -- should comfortably "
                              "exceed the expected number of distinct initial+terminal populations")
    parser.add_argument("--cyto-layer", default="counts", dest="cyto_layer",
                         help="adata.layers key CytoTRACE reads (default: 'counts', raw counts -- "
                              "CytoTRACE's gene-diversity signal is designed for unnormalized data)")
    parser.add_argument("--threshold-scheme", default="soft", dest="threshold_scheme",
                         choices=["soft", "hard"],
                         help="How kernels convert pseudotime/CytoTRACE ordering into directed "
                              "transitions (default: soft -- down-weights backward edges rather "
                              "than removing them outright, gentler than 'hard')")
    parser.add_argument("--cyto-weight", type=float, default=0.7, dest="cyto_weight")
    parser.add_argument("--conn-weight", type=float, default=0.3, dest="conn_weight")
    parser.add_argument("--dpt-weight", type=float, default=0.3, dest="dpt_weight",
                         help="Only used with --use-dpt-pseudotime; cyto/conn weights are "
                              "renormalized alongside it, not added on top")
    parser.add_argument("--use-dpt-pseudotime", action="store_true", dest="use_dpt_pseudotime",
                         help="Blend in diffusion pseudotime from scripts/08_trajectory_paga.py "
                              "(requires --input to be that script's output)")
    args = parser.parse_args()

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.input)

    if args.annotation_col not in adata.obs:
        raise SystemExit(f"--annotation-col '{args.annotation_col}' not in adata.obs "
                          f"(available: {list(adata.obs.columns)})")
    if args.cyto_layer not in adata.layers:
        raise SystemExit(f"--cyto-layer '{args.cyto_layer}' not in adata.layers "
                          f"(available: {list(adata.layers.keys())})")

    combined_kernel = build_combined_kernel(adata, args)

    print(f"Fitting GPCCA estimator (n_states={args.n_states})...")
    g = cr.estimators.GPCCA(combined_kernel)
    g.compute_schur()
    g.compute_macrostates(n_states=args.n_states, cluster_key=args.annotation_col)
    g.predict_terminal_states()
    # A real dataset can have a state that looks both initial and terminal. That is worth knowing
    # but is no reason to lose the run, so it is reported below rather than raised as an error.
    g.predict_initial_states(allow_overlap=True)

    print("Terminal states:", list(g.terminal_states.cat.categories))
    print("Initial states:", list(g.initial_states.cat.categories))

    g.compute_fate_probabilities()

    print("Computing lineage driver genes...")
    driver_df = g.compute_lineage_drivers()

    figstyle.apply_style()
    save_plots(adata, g, args.annotation_col, args.plot_dir)

    if args.drivers_csv:
        driver_df.to_csv(args.drivers_csv)
        print(f"Wrote driver genes to {args.drivers_csv}")

    estimator_path = pathlib.Path(args.output).with_suffix(".gpcca.pickle")
    g.write(str(estimator_path))
    print(f"Wrote fitted GPCCA estimator to {estimator_path} (avoids recomputing macrostates/fate "
          f"probabilities if only downstream plotting/reporting changes)")

    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()

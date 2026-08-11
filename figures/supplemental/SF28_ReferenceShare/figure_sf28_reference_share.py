"""SF28 — Reference-isoform share of gene expression across the floored pop_BC
(n = 1,548 genes; retained after the 25% reference-share floor, 4-CT re-scope 2026-07-11).

Standalone rebuild of the reference-share panel that was previously bundled
as Panel B of PairSetDescriptives. Split per Yul-era paper numbering so
the paper's SF28 reference resolves to one figure.

Value plotted per gene: reference-isoform DMSO mean expression (all_samples
basis) divided by the parent gene's total expression across all isoforms
(both from the Isopair pipeline; see Methods, pop_BC / Reference isoform).
By construction every retained gene has share >= 25%; the distribution starts
at the floor and is centered near the true dominant share (median ~64%).

Style: matplotlib rendered with ggplot-mimic theme (grey panel + white
gridlines) so the panel visually matches SF1-SF23. See
figures/lib/ggplot_style.py.
"""

from __future__ import annotations
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
LIB  = HERE.parents[1] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    render_and_validate,
    docx_body_fs,
    BAR_COLOR,
    TITLE_C,
    HEADER_FS,
)

apply_ggplot_rcparams()

NATIVE_W = 6.5
# Single wide panel with ample whitespace — target 13 pt (below the ~14 pt
# ceiling) so axis titles read well; the default 10 pt floor looked
# undersized against the rcParams-default tick labels (Pete 2026-07-10).
BODY_FS = docx_body_fs(NATIVE_W, target_pt=13)

DATA = HERE / "data"


def main():
    fra = pd.read_csv(DATA / "ref_expression_fraction.tsv", sep="\t")
    summary = pd.read_csv(DATA / "descriptives_summary.tsv", sep="\t").set_index("metric")["value"]
    n_genes = int(summary.loc["n_genes_in_pop_BC"])

    vals = (fra["ref_fraction_of_gene"].dropna() * 100).values  # → %
    med = float(np.median(vals))
    frac_ge_50 = float((vals >= 50).mean() * 100)

    fig, ax = plt.subplots(figsize=(NATIVE_W, 4.0))
    fig.subplots_adjust(left=0.11, right=0.95, top=0.88, bottom=0.14)

    style_axes_ggplot(ax)

    bins = np.arange(0, 105, 5)
    ax.hist(vals, bins=bins, color=BAR_COLOR, edgecolor="white", linewidth=0.6, zorder=3)

    floor_pct = float(summary.loc["floor_pct"])
    ax.axvline(med, linestyle="--", color="#c0392b", linewidth=1.5, zorder=4)
    ax.axvline(floor_pct, linestyle="-", color="#2c3e50", linewidth=1.2, zorder=4)
    # The bars are ~uniform in height across the whole range, so there is no
    # clear space inside the axes for the line labels. Add headroom above the
    # bars and place both labels in it (preserving the original y-ticks so the
    # headroom doesn't introduce a stray empty tick).
    yticks = [t for t in ax.get_yticks() if 0 <= t <= ax.get_ylim()[1]]
    data_ymax = ax.get_ylim()[1]
    ax.set_ylim(0, data_ymax * 1.20)
    ax.set_yticks(yticks)
    ytext = data_ymax * 1.10
    ax.text(
        med + 1.5, ytext,
        f"median = {med:.1f}%",
        ha="left", va="center",
        fontsize=BODY_FS,
        color="#c0392b",
    )
    ax.text(
        floor_pct + 1.5, ytext,
        f"{floor_pct:.0f}% floor",
        ha="left", va="center",
        fontsize=BODY_FS,
        color="#2c3e50",
    )

    ax.set_xlabel("Reference share of gene expression (%)", fontsize=BODY_FS)
    ax.set_ylabel("Number of genes", fontsize=BODY_FS)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    # Match tick-label size to the local BODY_FS (else ticks fall back to the
    # module rcParams default, dwarfing the axis titles).
    ax.tick_params(axis="both", labelsize=BODY_FS)

    # No overall figure title — caption carries the title role (Yul-style).
    render_and_validate(fig, HERE / "figure_sf28_reference_share",
                        native_width_in=NATIVE_W)
    print(f"  n_genes = {n_genes:,}  median = {med:.1f}%  ≥50% = {frac_ge_50:.1f}%")


if __name__ == "__main__":
    main()

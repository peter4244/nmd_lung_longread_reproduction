"""SF27 — Isoform count per gene in the Isopair pair-set cohort (n = 1,585).

Standalone rebuild of the isoform-count panel that was previously bundled
as Panel A of PairSetDescriptives. Split per Yul-era paper numbering so
the paper's SF27 reference resolves to one figure.

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
BODY_FS = docx_body_fs(NATIVE_W)

DATA = HERE / "data"


def main():
    iso = pd.read_csv(DATA / "isoforms_per_gene.tsv", sep="\t")
    summary = pd.read_csv(DATA / "descriptives_summary.tsv", sep="\t").set_index("metric")["value"]
    n_genes = int(summary.loc["n_genes_in_pop_BC"])
    med     = int(summary.loc["median_isoforms_per_gene"])

    fig, ax = plt.subplots(figsize=(NATIVE_W, 4.0))
    fig.subplots_adjust(left=0.11, right=0.90, top=0.88, bottom=0.14)

    style_axes_ggplot(ax)

    # Histogram with 1-count bins
    counts = iso["n_isoforms"].astype(int)
    max_bin = min(int(counts.quantile(0.995)) + 1, 30)
    bins = np.arange(1, max_bin + 2) - 0.5
    ax.hist(counts, bins=bins, color=BAR_COLOR, edgecolor="white", linewidth=0.6, zorder=3)

    # Explicit y-headroom + capped ticks so the top y-tick label never clips the
    # figure top edge (4-CT re-scope: median dropped to 6, tallest bin grew and
    # the auto-limit placed a 250 tick that clipped).
    bin_heights, _ = np.histogram(counts, bins=bins)
    top_h = int(bin_heights.max())
    ax.set_ylim(0, top_h * 1.12)
    ax.set_yticks(np.arange(0, top_h, 50))

    # Median reference line
    ax.axvline(med, linestyle="--", color="#c0392b", linewidth=1.5, zorder=4)
    ymax = ax.get_ylim()[1]
    ax.text(
        med + 0.4, ymax * 0.94,
        f"median = {med}",
        ha="left", va="top",
        fontsize=BODY_FS,
        color="#c0392b",
    )

    ax.set_xlabel("Number of expressed isoforms per gene", fontsize=BODY_FS)
    ax.set_ylabel("Number of genes", fontsize=BODY_FS)
    ax.set_xlim(0.5, max_bin + 0.5)
    ax.set_xticks([1, 5, 10, 15, 20, 25, 30])

    # No overall figure title — caption carries the title role (Yul-style).
    render_and_validate(fig, HERE / "figure_sf27_isoforms_per_gene",
                        native_width_in=NATIVE_W)
    print(f"  n_genes = {n_genes:,}  median = {med}  max bin shown = {max_bin}")


if __name__ == "__main__":
    main()

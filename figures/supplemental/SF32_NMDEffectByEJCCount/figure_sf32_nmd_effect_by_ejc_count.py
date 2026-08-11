"""SF32 — NMD effect size by downstream EJC count (PTC+ comparators).

Violin + box distributions of mashr posterior-mean logFC (SMG1i vs DMSO,
averaged across AT2/LAE/FB/MV) as a function of the comparator's number of
downstream EJCs. Restricted to comparators with at least one downstream
EJC from SF31's population. EJC counts ≥ 7 collapsed into a "7+" bin.

Downstream-EJC count is measured against the reference AUG-traced stop
codon (unbiased) — same anchor used by SF31. This avoids the TD2-CDS
attenuation that would arise from anchoring on the TD2-called stop.

Data source: `sf32_ejc_count_logfc.tsv` (produced by data_export.R).
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
    CONTROL_COLOR,
    AXIS_C,
    TITLE_C,
    HEADER_FS,
)

apply_ggplot_rcparams()

NATIVE_W = 9.0
BODY_FS = docx_body_fs(NATIVE_W)

DATA = HERE / "data"


def main():
    df = pd.read_csv(DATA / "sf32_ejc_count_logfc.tsv", sep="\t")
    df = df.dropna(subset=["mean_logFC", "n_downstream_ejc"])

    # 7+ bin
    df["ejc_bin"] = np.where(df["n_downstream_ejc"] >= 7, "7+",
                              df["n_downstream_ejc"].astype(int).astype(str))
    order = [str(i) for i in range(1, 7)] + ["7+"]
    df = df[df["ejc_bin"].isin(order)].copy()

    fig, ax = plt.subplots(figsize=(NATIVE_W, 5.0))
    fig.subplots_adjust(left=0.09, right=0.78, top=0.86, bottom=0.18)
    style_axes_ggplot(ax, xgrid=False, ygrid=True)

    data = [df.loc[df["ejc_bin"] == b, "mean_logFC"].to_numpy() for b in order]
    positions = np.arange(len(order))

    parts = ax.violinplot(data, positions=positions, widths=0.75,
                           showmeans=False, showmedians=False, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor(CONTROL_COLOR)
        body.set_edgecolor(AXIS_C)
        body.set_alpha(0.55)
        body.set_zorder(3)

    # Overlay narrow boxplots
    ax.boxplot(data, positions=positions, widths=0.14,
                showfliers=False, patch_artist=True,
                boxprops=dict(facecolor="white", edgecolor=TITLE_C, linewidth=1.0),
                medianprops=dict(color=TITLE_C, linewidth=1.4),
                whiskerprops=dict(color=TITLE_C, linewidth=1.0),
                capprops=dict(color=TITLE_C, linewidth=1.0),
                zorder=5)

    # Median line
    medians = [np.median(d) for d in data]
    ax.plot(positions, medians, color="#c0392b", linewidth=1.2,
            marker="o", markersize=4, zorder=6, label="Median trend")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.00), frameon=True, facecolor="white",
              edgecolor="none", fontsize=BODY_FS - 1)

    # Sample sizes folded into the x-tick labels (below the axis) — avoids
    # the collision with violin outlines that above-violin placement causes.
    n_by_bin = [len(df.loc[df["ejc_bin"] == b]) for b in order]
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{b}\n(n = {n:,})" for b, n in zip(order, n_by_bin)],
                        fontsize=BODY_FS, color=TITLE_C)
    ax.set_xlabel("Number of downstream EJCs (reference AUG-traced stop)",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Mean logFC (mashr posterior mean)",
                   fontsize=BODY_FS, color=TITLE_C)
    # Set the y-limits from the ACTUAL violin outline extents (the KDE, which
    # extends past the data by the default `cut`), so no violin tail is clipped
    # at the top/bottom of the axes -- the quantile-based limit used before sat
    # below the tallest violin (bin 5) and cut it off.
    vy = np.concatenate([b.get_paths()[0].vertices[:, 1] for b in parts["bodies"]])
    ax.set_ylim(vy.min() - 0.1, vy.max() + 0.1)
    ax.set_xlim(-0.6, len(order) - 0.4)

    # No overall figure title — caption carries the title role (Yul-style).
    render_and_validate(fig, HERE / "figure_sf32_nmd_effect_by_ejc_count",
                        native_width_in=NATIVE_W)
    print(f"  medians: {[round(m,2) for m in medians]}")
    print(f"  n_by_bin: {[len(d) for d in data]}")


if __name__ == "__main__":
    main()

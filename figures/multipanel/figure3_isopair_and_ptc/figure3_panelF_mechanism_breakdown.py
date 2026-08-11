"""
Figure 3 Panel F — PTC-causing splice events stacked by mechanism.

Horizontal stacked bar chart: for each splice event type attributed as the
cause of a PTC in an NMD comparator, the number of pairs is split by
mechanism (Frameshift / In-frame stop / 3'UTR splice). Event types are
ordered by total count descending — Skipped exon (SE) at top.

Mirrors the original Figure 4 Panel D ("PTC-Causing Splice Events (PTC+
Pairs)" / fig_sankey_ptcpos in the Rmd).

Scope: ENST-only effectively_ptc with direct event attribution (n=1,434
attributed of 1,489 effectively_ptc NMD pairs). Attribution is ref-AUG-derived
uniformly; mechanism populated by 05r_ref_atg_analysis.R for ALL effectively_ptc
pairs (no TD2-PTC-stop path).

Data:
  ./data/panelF_mechanism_breakdown.tsv
    Long-form: event_type, mechanism, n
    Source: all_attr from panel_e_compute.R (ENST-only pop_ptc_plus,
            1,434 attributed PTC+ pairs).

Style: HEADER_FS=18 bold, BODY_FS=14 regular. Mechanism palette from the
Rmd: 3'UTR splice teal, In-frame stop blue, Frameshift coral.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

HEADER_FS = 18
BODY_FS = 14

# Mechanism palette — matches Rmd's mech_colors dict.
# Stack order in each row from left to right: 3'UTR splice, In-frame stop, Frameshift.
MECH_ORDER = ["3'UTR splice", "In-frame stop", "Frameshift"]
MECH_COLORS = {
    "3'UTR splice":  "#66c2a5",
    "In-frame stop": "#4393c3",
    "Frameshift":    "#d6604d",
}
C_TITLE = "#222222"
C_AXIS = "#555555"
# At 6×4, the smaller bars don't have room for BODY_FS labels — restrict
# in-segment value labels to substantial segments and shrink their font.
LABEL_MIN = 100
INNER_FS = 9
LEGEND_FS = 11


def load_data():
    df = pd.read_csv(HERE.parent / "data" / "panelF_mechanism_breakdown.tsv", sep="\t")
    pivoted = df.pivot_table(index="event_type", columns="mechanism",
                              values="n", fill_value=0)
    # Restore the descending-by-total ordering and the mechanism column order
    pivoted = pivoted.reindex(columns=[m for m in MECH_ORDER if m in pivoted.columns],
                                fill_value=0)
    pivoted = pivoted.assign(_total=pivoted.sum(axis=1)).sort_values("_total",
                                                                     ascending=True)
    pivoted = pivoted.drop(columns="_total")
    return pivoted


def build_figure(piv):
    # Composite cell target: 6.0 × 4.0 in (1.5:1). y-axis category labels + bottom legend.
    fig, ax = plt.subplots(figsize=(6, 4))

    y = np.arange(len(piv.index))
    left = np.zeros(len(piv.index))
    for mech in [m for m in MECH_ORDER if m in piv.columns]:
        vals = piv[mech].to_numpy()
        ax.barh(y, vals, left=left, height=0.65,
                color=MECH_COLORS[mech], edgecolor="none", label=mech)
        # In-segment value labels (white bold) for stacks >= LABEL_MIN
        for i, v in enumerate(vals):
            if v >= LABEL_MIN:
                ax.text(left[i] + v / 2, y[i], f"{int(v):,}",
                        ha="center", va="center",
                        fontsize=INNER_FS, fontweight="bold", color="white")
        left += vals

    ax.set_yticks(y)
    ax.set_yticklabels(piv.index, fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="x", labelsize=BODY_FS, colors=C_TITLE)

    ax.set_xlabel("Number of pairs", fontsize=BODY_FS, color=C_TITLE)
    # No in-panel title / subtitle — caption / composite letter label carries it.
    # The "1,812 attributed pairs | stacked by mechanism" provenance line moves
    # to the figure caption per the grant titleless convention.

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    legend = ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.42),
                       ncol=3, frameon=False, fontsize=LEGEND_FS,
                       handlelength=1.2, handletextpad=0.5, columnspacing=1.0)
    for text in legend.get_texts():
        text.set_color(C_TITLE)

    fig.subplots_adjust(left=0.22, right=0.97, bottom=0.30, top=0.96)
    return fig


def main():
    piv = load_data()
    fig = build_figure(piv)
    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)
    out_dir = HERE.parent
    fig.savefig(out_dir / "figure3_panelF_mechanism_breakdown.pdf",
                facecolor="white")
    fig.savefig(out_dir / "figure3_panelF_mechanism_breakdown.png",
                dpi=300, facecolor="white")
    print("Saved: figure3_panelF_mechanism_breakdown.pdf and .png")


if __name__ == "__main__":
    main()

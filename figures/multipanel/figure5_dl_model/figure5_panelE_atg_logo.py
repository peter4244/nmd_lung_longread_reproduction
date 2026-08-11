"""
Figure 5 Panel E — Per-nucleotide signed SHAP × input around the start
codon, NMD class. Built natively from the joint DeepSHAP motif logo TSV.

Data:
  data/motif_logo_atg_joint_atg500_stop500.tsv
    Columns: relative_position, channel, nmd_mean_contrib,
             ctrl_mean_contrib, nmd_freq, ctrl_freq, n_nmd, n_ctrl
    Source: scripts/export_joint_motif_logos.py in the model repo
            (5 DeepSHAP runs pooled; n_test = 2,268 NMD + 7,863 Control).

Coordinate system in the TSV: relative_position 0 is the middle nucleotide
of the start codon. Mapping to canonical Kozak coordinates (A of ATG = +1):
   TSV pos      Kozak pos      Identity
   −1           +1             A (always)
    0           +2             T (always)  ← relative_position = 0
   +1           +3             G (always)
   −2           −1             variable
   −3           −2             variable
   −4           −3             variable (canonical purine)
   +2           +4             variable (canonical G)

Letter height = nmd_mean_contrib (signed SHAP × input averaged across
NMD test isoforms). Positive (above zero) = pushes prediction toward
NMD; negative (below zero) = pushes away from NMD. The three invariant
ATG positions (−1, 0, +1 in TSV coords) have zero attribution by
construction (one-hot inputs identical across all isoforms) and are
shaded with a light box marking the start codon.

Style overlay: publication composite cell, no in-panel title, validator-
clean, 6 × 4 cell at 1.5:1 aspect.
"""

import sys
from pathlib import Path

import logomaker
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    assert_text_within_canvas,
    TITLE_C,
    AXIS_C,
)
apply_ggplot_rcparams()

# Panel-local font sizes bumped from the shim defaults for composite readability
# (the composite is 22 in wide; body ticks below ~13 pt read too small).
AXIS_FS = 14
TICK_FS = 13

DATA = HERE.parent / "data" / "panelE_motif_logo_atg.tsv"

WINDOW_LEFT  = -10        # TSV-coordinate window for the panel
WINDOW_RIGHT = +13
NUC_COLORS = {
    "A": "#109648",   # green
    "C": "#255C99",   # blue
    "G": "#F7B32B",   # gold
    "U": "#D62839",   # red — RNA notation; data column is "T", display as "U"
}


def load_letter_heights(class_col: str = "nmd_mean_contrib") -> pd.DataFrame:
    """Build a logomaker-compatible matrix (rows = positions, cols = A/C/G/U)
    of signed SHAP × input contributions for the chosen class column.

    The TSV stores the channel as DNA letter "T"; we relabel to RNA letter
    "U" here so all downstream figure annotations use RNA notation.
    """
    df = pd.read_csv(DATA, sep="\t")
    df = df[df["channel"].isin(["A", "C", "G", "T"])]
    df = df[(df["relative_position"] >= WINDOW_LEFT) &
            (df["relative_position"] <= WINDOW_RIGHT)]
    matrix = df.pivot(index="relative_position", columns="channel",
                      values=class_col).fillna(0.0)
    matrix = matrix.rename(columns={"T": "U"})
    matrix = matrix[["A", "C", "G", "U"]]
    return matrix


def build_figure():
    mat = load_letter_heights("nmd_mean_contrib")

    fig, ax = plt.subplots(figsize=(6, 4))
    # Plain white axes with a bottom + left spine — no panel bg, no gridlines.
    ax.set_facecolor("white")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(AXIS_C)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(axis="both", length=3, colors=TITLE_C)

    logo = logomaker.Logo(
        mat, ax=ax,
        color_scheme={"A": NUC_COLORS["A"], "C": NUC_COLORS["C"],
                      "G": NUC_COLORS["G"], "U": NUC_COLORS["U"]},
        flip_below=True,
        font_name="Arial",
        stack_order="big_on_top",
    )

    # Shade the three start-codon positions (-1, 0, +1 in TSV coords)
    ax.axvspan(-1.5, 1.5, color="#EAEAEA", alpha=0.55, zorder=1)

    # Zero baseline
    ax.axhline(0, color=TITLE_C, linewidth=0.6, zorder=2)

    # X-axis: convert TSV positions to codon-relative coordinates. Same
    # convention as Panel F: no position 0, first nt of codon = +1.
    #   TSV ≥ −1 (start codon and downstream): pos = TSV + 2
    #   TSV ≤ −2 (upstream):                   pos = TSV + 1
    def tsv_to_pos(p):
        return p + 2 if p >= -1 else p + 1

    visible_ticks_tsv = list(range(WINDOW_LEFT, WINDOW_RIGHT + 1, 2))
    ax.set_xticks(visible_ticks_tsv)
    ax.set_xticklabels([tsv_to_pos(p) for p in visible_ticks_tsv],
                       fontsize=TICK_FS, color=TITLE_C)
    ax.set_xlabel("Position relative to start codon (A of AUG = +1)",
                  fontsize=AXIS_FS, color=TITLE_C)

    ax.set_ylabel("Signed SHAP × input", fontsize=AXIS_FS, color=TITLE_C)
    ax.tick_params(axis="y", labelsize=TICK_FS, colors=TITLE_C)

    # Tight margins for composite use
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.19, top=0.95)
    return fig


def main():
    fig = build_figure()

    try:
        from validate_figure_layout import validate_figure_layout
        validate_figure_layout(fig, fig.axes[0], verbose=True)
    except Exception as e:
        print(f"[validator] non-fatal: {e}")

    assert_text_within_canvas(fig)

    out = HERE.parent / "figure5_panelE_atg_logo"
    fig.savefig(f"{out}.pdf", facecolor="white")
    fig.savefig(f"{out}.png", dpi=300, facecolor="white")
    print(f"Saved: {out}.pdf and .png")


if __name__ == "__main__":
    main()

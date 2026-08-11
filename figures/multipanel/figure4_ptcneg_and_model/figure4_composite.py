"""
Figure 4 composite — 2×2 uniform layout (4 panels A–D).

Layout (2 cols × 2 rows, identical 1.5:1 cells, 12"×8" landscape):
  Row 1 — Section A. Cleanest small-n GENCODE-only evidence
          (all-3-ENST + coding-CDS; own-GENCODE-stop PTC determination):
    A 5'UTR length         | B Longest 5'UTR ORF length
  Row 2 — Section C. ENST-reference + ref-AUG-traceable scope:
    C 5'UTR length         | D Longest 5'UTR ORF length

The §4 TD2-bias justification (TD2 vs ref-AUG ORF length, paired Kozak,
TD2-ATG position relative to ref AUG) has been moved to the consolidated
TD2BiasEvidence supplemental figure.

All panels render at 6.0 × 4.0 in (1.5:1) without bbox_inches="tight",
so the composite is a uniform `GridSpec(2, 2)` with no per-row scaling.

Pre-registered in ./RATIONALE.md.
"""

import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
LIB = HERE.parents[1] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

LBL_FS = 22
LBL_COLOR = "#111111"

CELL_W = 6.0
CELL_H = 4.0
N_COLS = 2
N_ROWS = 2

PANELS = {
    "A": "figure4_panelA_5utr_length_all3enst.png",
    "B": "figure4_panelB_longest_5utr_orf_all3enst.png",
    "C": "figure4_panelC_5utr_length_refaug.png",
    "D": "figure4_panelD_longest_5utr_orf_refaug.png",
}
ROWS = [["A", "B"], ["C", "D"]]


def build_figure():
    paths = {k: HERE / v for k, v in PANELS.items()}

    fig_w = CELL_W * N_COLS
    fig_h = CELL_H * N_ROWS
    fig = plt.figure(figsize=(fig_w, fig_h))

    gs = GridSpec(
        N_ROWS, N_COLS, figure=fig,
        hspace=0.15, wspace=0.02,
        left=0.025, right=0.995, top=0.93, bottom=0.03,
    )

    for ri, row in enumerate(ROWS):
        for ci, letter in enumerate(row):
            ax = fig.add_subplot(gs[ri, ci])
            img = mpimg.imread(str(paths[letter]))
            ax.imshow(img, aspect="auto")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            ax.text(
                0.0, 1.03, letter,
                transform=ax.transAxes,
                fontsize=LBL_FS, fontweight="bold", color=LBL_COLOR,
                ha="left", va="bottom",
            )

    return fig


def main():
    fig = build_figure()
    pdf = HERE / "figure4_composite.pdf"
    png = HERE / "figure4_composite.png"
    fig.savefig(pdf, facecolor="white")
    fig.savefig(png, dpi=300, facecolor="white")
    print(f"Saved: {pdf.name} and {png.name}")


if __name__ == "__main__":
    main()

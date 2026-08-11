"""
Figure 3 composite — 2x3 uniform multipanel layout of the Isopair + PTC analysis.

Layout (2 cols × 3 rows, identical cells):
  Row 1: A pair concept schematic | B sequence similarity
  Row 2: C event prevalence       | D stop-to-EJC distance
  Row 3: E PTC-causing events     | F mechanism breakdown

Design principle: every panel script renders at the same target aspect
(CELL_W × CELL_H = 6.0 × 4.0 in, 1.5:1) and saves WITHOUT bbox_inches="tight"
— so the saved PNG dimensions are exactly the target. The composite is then
a plain GridSpec(3, 2) with equal cells and no per-row scaling, which
guarantees zero aspect-ratio distortion of any panel.

Source PNGs come from the per-panel render scripts in this folder. Re-run
each panel script (or `data_export.R` + all panel scripts) if any source
has changed.
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

# Per-grant convention: panel labels at 22pt bold
LBL_FS = 22
LBL_COLOR = "#111111"

# Cell geometry — must match each panel's figsize.
CELL_W = 6.0
CELL_H = 4.0
N_COLS = 2
N_ROWS = 3

# Panel PNG paths
PANELS = {
    "A": "figure3_panelA_pair_concept.png",
    "B": "figure3_panelB_sequence_similarity.png",
    "C": "figure3_panelC_event_prevalence.png",
    "D": "figure3_panelD_stop_codon_distance.png",
    "E": "figure3_panelE_ptc_event_attribution.png",
    "F": "figure3_panelF_mechanism_breakdown.png",
}
ROWS = [["A", "B"], ["C", "D"], ["E", "F"]]


def build_figure():
    paths = {k: HERE / v for k, v in PANELS.items()}

    fig_w = CELL_W * N_COLS
    fig_h = CELL_H * N_ROWS
    fig = plt.figure(figsize=(fig_w, fig_h))

    # top < 1 and hspace > 0 leave headroom for the panel letter labels
    # (positioned at axes-coords y=1.02; needs ~3% of cell height of clear
    # space above each cell or the letter ascenders get clipped at the
    # figure edge or by the row above).
    gs = GridSpec(
        N_ROWS, N_COLS, figure=fig,
        hspace=0.10, wspace=0.02,
        left=0.005, right=0.995, top=0.96, bottom=0.005,
    )

    for ri, row in enumerate(ROWS):
        for ci, letter in enumerate(row):
            ax = fig.add_subplot(gs[ri, ci])
            img = mpimg.imread(str(paths[letter]))
            # aspect="auto" stretches the image to fill the cell — since
            # every panel was rendered at CELL_W:CELL_H and the cell matches
            # that aspect, the stretch ratio is 1:1 (no distortion).
            ax.imshow(img, aspect="auto")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            # Panel label in top-left corner of the cell.
            ax.text(
                -0.005, 1.02, letter,
                transform=ax.transAxes,
                fontsize=LBL_FS, fontweight="bold", color=LBL_COLOR,
                ha="left", va="bottom",
            )

    return fig


def main():
    fig = build_figure()
    pdf = HERE / "figure3_composite.pdf"
    png = HERE / "figure3_composite.png"
    fig.savefig(pdf, facecolor="white")
    fig.savefig(png, dpi=300, facecolor="white")
    print(f"Saved: {pdf.name} and {png.name}")


if __name__ == "__main__":
    main()

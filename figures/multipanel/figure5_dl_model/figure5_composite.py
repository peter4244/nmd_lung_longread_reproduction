"""
Figure 5 composite — DL model panels A–G with aspect-preserving layout.

Layout (~22.32 in × 12.00 in):

  ┌────────┬─────────────────────────┐
  │        │ B  ROC + PR     │  C    │
  │        │ inset           │       │
  │   A    ├─────────────────────────┤
  │ Model  │ D  Structural   │  E    │
  │ arch.  │ features        │       │
  │        ├─────────────────────────┤
  │        │ F  STOP SHAP    │  G    │
  └────────┴─────────────────────────┘

Right side: 3 rows × 2 cols, alphabetical order B C / D E / F G. Each
cell is 6×4 in (1.5:1, native for B/C/D/G). E and F (native 2:1) fit
width-wise with small top/bottom whitespace, preserving their native
aspect.

Panel A (native aspect 0.86) occupies the left column at full composite
height (12 in), preserving its native shape (width = 0.86 × 12 ≈ 10.32 in).

Cell width W = 6.0 in (1.5:1 → cell height 4.0 in).
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

# Cell geometry (matches B/C/D/G native 1.5:1)
CELL_W = 6.0
CELL_H = 4.0
N_RIGHT_COLS = 2
N_RIGHT_ROWS = 3

RIGHT_W = CELL_W * N_RIGHT_COLS              # = 12.0
RIGHT_H = CELL_H * N_RIGHT_ROWS              # = 12.0

# Panel A column preserves native aspect 0.86 at full composite height
A_ASPECT = 14.0 / 16.33                      # ≈ 0.857
A_H = RIGHT_H
A_W = A_ASPECT * A_H                         # ≈ 5.14

FIG_W = A_W + RIGHT_W
FIG_H = RIGHT_H

PANELS = {
    "A": "figure5_panelA_architecture.png",
    "B": "figure5_panelB_performance.png",
    "C": "figure5_panelC_branch_importance.png",
    "D": "figure5_panelD_structural_features.png",
    "E": "figure5_panelE_atg_logo.png",
    "F": "figure5_panelF_stop_logo.png",
    "G": "figure5_panelG_uorf_attention.png",
}


def add_panel(fig, gs_slot, letter, paths, aspect="auto"):
    ax = fig.add_subplot(gs_slot)
    img = mpimg.imread(str(paths[letter]))
    ax.imshow(img, aspect=aspect)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(
        0.0, 1.01, letter,
        transform=ax.transAxes,
        fontsize=LBL_FS, fontweight="bold", color=LBL_COLOR,
        ha="left", va="bottom",
    )
    return ax


# A COMPOSITE THAT READS PNGs OFF DISK CANNOT TELL A FRESH PANEL FROM A COMMITTED ONE, and on
# 2026-08-10 that stopped being hypothetical. In clean-room job 9049007 panels A, E and F failed to
# render -- A needs the sibling model repo, E and F need logomaker, which is not in the image -- and
# this script assembled the COMMITTED copies of all three and reported ok. The run said Figure 5
# reproduced while three of its seven panels came from a laptop.
#
# The guard is opt-in and takes a reference FILE rather than a duration. mtime-against-a-clock does
# not work here: a git checkout stamps every file with the checkout time, so "modified recently" is
# true of a stale panel and a fresh one alike. "Newer than this marker, which the stage touched
# before it started rendering" is true only of a panel this run actually wrote.
def check_panel_freshness(paths, reference):
    ref_mtime = reference.stat().st_mtime
    stale = [(k, paths[k]) for k in sorted(paths)
             if paths[k].stat().st_mtime <= ref_mtime]
    if stale:
        lines = "\n".join(
            f"    {k}  {p.name}  (mtime {p.stat().st_mtime:.0f} <= marker {ref_mtime:.0f})"
            for k, p in stale)
        raise SystemExit(
            f"STALE PANELS -- {len(stale)} of {len(paths)} were not written by this run:\n"
            f"{lines}\n"
            f"  Marker: {reference}\n"
            f"  Their producers did not run, or ran and failed. Assembling them would report a\n"
            f"  composite this run did not build. Fix the producers, or drop --require-fresh to\n"
            f"  assemble deliberately from what is on disk.")


def build_figure():
    paths = {k: HERE / v for k, v in PANELS.items()}

    fig = plt.figure(figsize=(FIG_W, FIG_H))

    # Outer split: A column on left, right grid on right.
    # top=0.955 leaves just enough headroom for the panel letters
    # (y=1.01 in axes coords, 22-pt bold).
    outer = GridSpec(
        1, 2, figure=fig,
        width_ratios=[A_W, RIGHT_W],
        left=0.005, right=0.998, top=0.955, bottom=0.015,
        wspace=0.015,
    )

    # Panel A: aspect="auto" since the cell is sized to A's native aspect.
    add_panel(fig, outer[0, 0], "A", paths, aspect="auto")

    # Right grid: 3 rows × 2 cols of 1.5:1 cells (matches B/C/D/G native).
    right_gs = outer[0, 1].subgridspec(
        N_RIGHT_ROWS, N_RIGHT_COLS,
        hspace=0.07, wspace=0.04,
    )

    # Row 1: B, C — both native 1.5:1.
    add_panel(fig, right_gs[0, 0], "B", paths, aspect="auto")
    add_panel(fig, right_gs[0, 1], "C", paths, aspect="auto")

    # Row 2: D, E. D is 1.5:1 (auto); E is 2:1 (equal → whitespace).
    add_panel(fig, right_gs[1, 0], "D", paths, aspect="auto")
    add_panel(fig, right_gs[1, 1], "E", paths, aspect="equal")

    # Row 3: F, G. F is 2:1 (equal → whitespace); G is 1.5:1 (auto).
    add_panel(fig, right_gs[2, 0], "F", paths, aspect="equal")
    add_panel(fig, right_gs[2, 1], "G", paths, aspect="auto")

    return fig


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--require-fresh", metavar="MARKER_FILE", default=None,
                    help="refuse to assemble any panel PNG not newer than this file")
    args = ap.parse_args()

    paths = {k: HERE / v for k, v in PANELS.items()}
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        raise SystemExit("missing panel PNG(s):\n  " + "\n  ".join(missing))

    # Always say WHAT was assembled, whether or not the guard is on. A composite is the one figure
    # whose inputs are invisible in its output.
    print("panels assembled:")
    for k in sorted(paths):
        st = paths[k].stat()
        print(f"    {k}  {paths[k].name:44s} {st.st_size:>9,d} B  mtime {st.st_mtime:.0f}")

    if args.require_fresh:
        marker = Path(args.require_fresh)
        if not marker.exists():
            raise SystemExit(f"--require-fresh marker does not exist: {marker}")
        check_panel_freshness(paths, marker)
        print(f"    all {len(paths)} panels newer than {marker}")

    fig = build_figure()

    from validate_figure_layout import validate_multipanel_layout
    validate_multipanel_layout(fig, verbose=True)

    out_dir = HERE
    fig.savefig(out_dir / "figure5_composite.pdf", facecolor="white")
    fig.savefig(out_dir / "figure5_composite.png", dpi=200, facecolor="white")
    print(f"Saved: figure5_composite.pdf and .png  ({FIG_W:.2f} × {FIG_H:.2f} in)")


if __name__ == "__main__":
    main()

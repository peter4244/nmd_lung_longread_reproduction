"""
Figure 3 Panel A — Isopair Isoform Pair Framework schematic.

Conceptual diagram showing the three-isoform pair construction:
  Reference (dominant, non-NMD)   — 6 exons, blue
  NMD comparator (SMG1i-responsive) — 5 exons (exon 4 skipped), coral
  Control comparator (non-NMD)     — 6 exons, A3SS lengthens exon 3, light blue

Each gene contributes paired comparisons sharing the same reference isoform:
  NMD pair    = Reference vs NMD comparator
  Control pair = Reference vs Control comparator

Splice events distinguishing the comparator from the reference are detected
by Isopair and validated through reference reconstruction.

Ported from:
  ../make_pair_concept_figure.R (Mar 2026 R/ggplot2 version)
Pattern reference:
  grants/nmd_2026/figure_aim1.py (matplotlib + figures/lib/ style)
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

# Typography (matches matplotlib.md convention)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

HEADER_FS = 15
BODY_FS = 11

# Colors — matched to R source for visual consistency with prior renderings
C_REF_FILL = "#4393c3"
C_REF_LBL = "#2166ac"
C_NMD_FILL = "#ef8a62"
C_NMD_LBL = "#d6604d"
C_CTRL_FILL = "#67a9cf"
C_CTRL_LBL = "#2b8cbe"
C_INTRON = "#777777"
C_EXON_EDGE = "#4a4a4a"
C_BRACKET = "#333333"
C_TITLE = "#222222"

# Gene model coordinates (genomic, + strand). Match R source exactly.
REF_STARTS = [200, 600, 1000, 1400, 1800, 2200]
REF_ENDS = [350, 800, 1200, 1550, 2000, 2500]

NMD_STARTS = [200, 600, 1000, 1800, 2200]            # exon 4 skipped
NMD_ENDS = [350, 800, 1200, 2000, 2500]

CTRL_STARTS = [200, 600, 940, 1400, 1800, 2200]      # A3SS lengthens exon 3
CTRL_ENDS = [350, 800, 1200, 1550, 2000, 2500]

REF_Y = 5.0
NMD_Y = 3.0
CTRL_Y = 1.0
EXON_H = 0.42
EVENT_BAR_OFFSET = 0.55
EVENT_TEXT_OFFSET = 0.95

BRACKET_NMD_X = 2600
BRACKET_CTRL_X = 3050  # wider gap so the "NMD pair" label fits between
TICK_LEN = 40


def draw_track(ax, starts, ends, y, fill_color):
    """Draw one isoform track: introns + chevrons + exon rectangles."""
    # Intron lines
    for i in range(len(starts) - 1):
        ax.plot(
            [ends[i], starts[i + 1]],
            [y, y],
            color=C_INTRON,
            linewidth=1.0,
            solid_capstyle="butt",
            zorder=1,
        )

    # Chevron arrows on introns (1–3 per intron based on length)
    arm = 0.11
    for i in range(len(starts) - 1):
        intron_s = ends[i]
        intron_e = starts[i + 1]
        intron_len = intron_e - intron_s
        n_chev = max(1, min(3, intron_len // 120))
        step = intron_len / (n_chev + 1)
        for j in range(1, int(n_chev) + 1):
            cx = intron_s + step * j
            w = min(25, intron_len * 0.12)
            ax.plot(
                [cx - w, cx], [y + arm, y],
                color=C_INTRON, linewidth=1.0, solid_capstyle="butt", zorder=1,
            )
            ax.plot(
                [cx - w, cx], [y - arm, y],
                color=C_INTRON, linewidth=1.0, solid_capstyle="butt", zorder=1,
            )

    # Exon rectangles
    for s, e in zip(starts, ends):
        rect = Rectangle(
            (s, y - EXON_H / 2),
            e - s,
            EXON_H,
            facecolor=fill_color,
            edgecolor=C_EXON_EDGE,
            linewidth=0.8,
            zorder=2,
        )
        ax.add_patch(rect)


def draw_bracket(ax, x, y_top, y_bot, tick_len, color):
    """Vertical bracket with left-pointing ticks at top and bottom."""
    ax.plot([x, x], [y_bot, y_top], color=color, linewidth=1.4, solid_capstyle="butt")
    ax.plot([x, x - tick_len], [y_top, y_top], color=color, linewidth=1.4, solid_capstyle="butt")
    ax.plot([x, x - tick_len], [y_bot, y_bot], color=color, linewidth=1.4, solid_capstyle="butt")


def build_figure():
    # Composite cell target: 6.0 × 4.0 in (1.5:1). All Figure-3 panels share
    # this aspect so the 2×3 grid composes without per-row scaling.
    # Font sizes scaled down (HEADER_FS=15, BODY_FS=11) so left-of-track
    # labels still fit; xlim extends leftward to give them horizontal room.
    fig, ax = plt.subplots(figsize=(6, 4))

    # Three isoform tracks
    draw_track(ax, REF_STARTS, REF_ENDS, REF_Y, C_REF_FILL)
    draw_track(ax, NMD_STARTS, NMD_ENDS, NMD_Y, C_NMD_FILL)
    draw_track(ax, CTRL_STARTS, CTRL_ENDS, CTRL_Y, C_CTRL_FILL)

    # Pair brackets on right
    draw_bracket(ax, BRACKET_NMD_X, REF_Y, NMD_Y, TICK_LEN, C_BRACKET)
    draw_bracket(ax, BRACKET_CTRL_X, REF_Y, CTRL_Y, TICK_LEN, C_BRACKET)

    # Track labels (left of gene models) — body, never bold per principles.md
    ax.text(140, REF_Y, "Reference\n(dominant, non-NMD)",
            ha="right", va="center", fontsize=BODY_FS, color=C_REF_LBL)
    ax.text(140, NMD_Y, "NMD comparator\n(SMG1i-responsive)",
            ha="right", va="center", fontsize=BODY_FS, color=C_NMD_LBL)
    ax.text(140, CTRL_Y, "Control comparator\n(non-NMD, next-best)",
            ha="right", va="center", fontsize=BODY_FS, color=C_CTRL_LBL)

    # Bracket labels (right of brackets) — body, never bold
    ax.text(BRACKET_NMD_X + 40, (REF_Y + NMD_Y) / 2, "NMD\npair",
            ha="left", va="center", fontsize=BODY_FS, color=C_NMD_LBL)
    ax.text(BRACKET_CTRL_X + 40, (REF_Y + CTRL_Y) / 2, "Control\npair",
            ha="left", va="center", fontsize=BODY_FS, color=C_CTRL_LBL)

    # Event annotations: solid colored bar + italic body label, BELOW the track
    # Skipped exon (NMD track) — lands in the gap between NMD and Control rows
    ax.plot([1370, 1580], [NMD_Y - EVENT_BAR_OFFSET, NMD_Y - EVENT_BAR_OFFSET],
            color=C_NMD_LBL, linewidth=2.0, solid_capstyle="butt")
    ax.text(1475, NMD_Y - EVENT_TEXT_OFFSET, "Skipped exon",
            ha="center", va="top", fontsize=BODY_FS,
            fontstyle="italic", color=C_NMD_LBL)

    # A3SS (Control track) — lands below the Control row
    ax.plot([910, 1030], [CTRL_Y - EVENT_BAR_OFFSET, CTRL_Y - EVENT_BAR_OFFSET],
            color=C_REF_LBL, linewidth=2.0, solid_capstyle="butt")
    ax.text(970, CTRL_Y - EVENT_TEXT_OFFSET, "A3SS",
            ha="center", va="top", fontsize=BODY_FS,
            fontstyle="italic", color=C_REF_LBL)

    # No in-panel title — caption / composite letter label carries it.
    # Per grant convention (figure_nmd_prelim_panel*.py): panels are titleless.

    # Axes setup — xlim extended leftward to fit left-of-track labels at
    # the smaller 6×4 canvas. Right side keeps room for two brackets + labels.
    ax.set_xlim(-1200, 3600)
    ax.set_ylim(-0.5, 6.0)
    ax.set_aspect("auto")
    ax.axis("off")

    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)
    return fig


def main():
    fig = build_figure()

    from validate_figure_layout import validate_figure_layout
    result = validate_figure_layout(fig, fig.axes[0], verbose=True)
    if result['summary']['n_errors'] > 0:
        raise SystemExit(
            f"Panel A validator: {result['summary']['n_errors']} errors — fix before saving"
        )

    out_dir = HERE.parent
    pdf_path = out_dir / "figure3_panelA_pair_concept.pdf"
    png_path = out_dir / "figure3_panelA_pair_concept.png"

    fig.savefig(pdf_path, facecolor="white")
    fig.savefig(png_path, dpi=300, facecolor="white")

    print(f"Saved: {pdf_path.name} and {png_path.name}")


if __name__ == "__main__":
    main()

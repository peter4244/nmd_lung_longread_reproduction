"""
Figure 3 Panel B — Sequence similarity NMD vs Control.

Density plot of exonic sequence shared between the comparator isoform and
its reference. NMD pairs shift right (more similar to reference) than
Control pairs, consistent with targeted splicing changes that perturb a
small fraction of the transcript rather than reshape it.

Data:
  ./data/panelB_sequence_similarity.tsv
    Columns: comparison ("NMD" | "Control"), pct_shared (numeric 0-100)
    Source: Isopair::quantifyPairDivergence() → div_c2/c4$pct_shared.
    Exported by ./data_export.R.

Style: matches Panel A — HEADER_FS=18 bold, BODY_FS=14 regular per
figures/lib/principles.md. NMD = coral, Control = light blue (consistent
with Panel A's NMD comparator / Control comparator track colors).
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

HEADER_FS = 18
BODY_FS = 14

# Colors — match Panel A
C_NMD_FILL = "#ef8a62"
C_NMD_LINE = "#d6604d"
C_CTRL_FILL = "#67a9cf"
C_CTRL_LINE = "#2b8cbe"
C_TITLE = "#222222"


def load_data():
    tsv = HERE.parent / "data" / "panelB_sequence_similarity.tsv"
    df = pd.read_csv(tsv, sep="\t")
    return df


def build_figure(df):
    # Sized to the composite cell: 6.0 × 4.0 in = 1.5:1 aspect.
    # All Figure-3 panels target this aspect so the 2×3 grid composes
    # without per-row scaling — see figure3_composite_methodology.md.
    fig, ax = plt.subplots(figsize=(6, 4))

    x_grid = np.linspace(0, 100, 500)

    for label, fill, line in [
        ("NMD", C_NMD_FILL, C_NMD_LINE),
        ("Control", C_CTRL_FILL, C_CTRL_LINE),
    ]:
        vals = df.loc[df["comparison"] == label, "pct_shared"].to_numpy()
        kde = gaussian_kde(vals, bw_method="scott")
        y = kde(x_grid)
        ax.fill_between(x_grid, y, alpha=0.45, color=fill, linewidth=0,
                        label=f"{label} (n={len(vals):,})")
        ax.plot(x_grid, y, color=line, linewidth=1.6)

    # No in-panel title — caption / composite letter label carries it.

    # Axis labels — body
    ax.set_xlabel("Sequence shared (%)", fontsize=BODY_FS, color=C_TITLE)
    # No y-axis label or tick labels — absolute KDE density values aren't
    # meaningful for shape comparison; remove for visual cleanliness.
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=BODY_FS, colors=C_TITLE)

    # Hide top/right/left spines (clean look — no y-axis exposed)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#555555")
    ax.spines["bottom"].set_linewidth(0.8)

    ax.set_xlim(0, 100)
    ax.set_ylim(bottom=0)

    legend = ax.legend(fontsize=BODY_FS, loc="upper left", frameon=False)
    for text in legend.get_texts():
        text.set_color(C_TITLE)

    fig.subplots_adjust(left=0.04, right=0.97, bottom=0.16, top=0.96)
    return fig


def main():
    df = load_data()
    fig = build_figure(df)

    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)

    out_dir = HERE.parent
    pdf_path = out_dir / "figure3_panelB_sequence_similarity.pdf"
    png_path = out_dir / "figure3_panelB_sequence_similarity.png"

    fig.savefig(pdf_path, facecolor="white")
    fig.savefig(png_path, dpi=300, facecolor="white")

    print(f"Saved: {pdf_path.name} and {png_path.name}")


if __name__ == "__main__":
    main()

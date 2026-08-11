"""
Figure 4 Panel C — 5'UTR length, ENST-reference + ref-AUG-traceable scope.

3-group boxplot (log y): NMD+/PTC+ (n=756) vs NMD+/PTC- (n=63; ORF match
+ ORF diff merged) vs Control (n=819). 5'UTR boundary anchored on the
ref-AUG-projected start codon (utr5_features_refaug.rds).

Section C of Figure 4. Pre-registered in ./RATIONALE.md.

Data:
  ./data/panelC_5utr_length_long.tsv
  ./data/panelC_5utr_length_descriptives.tsv
  ./data/panelC_5utr_length_pairwise.tsv
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

BODY_FS = 14
LABEL_FS = 12
C_TITLE = "#222222"
C_AXIS = "#555555"

GROUP_ORDER = ["NMD+/PTC+", "NMD+/PTC-", "Control"]
COLORS = {"NMD+/PTC+": "#ef8a62", "NMD+/PTC-": "#fee08b", "Control": "#92c5de"}


from p_to_stars import p_to_stars as sig_marker  # noqa: E402 — unified 4-tier scheme


def build_figure():
    long = pd.read_csv(HERE.parent / "data" / "panelC_5utr_length_long.tsv", sep="\t")
    stats = pd.read_csv(HERE.parent / "data" / "panelC_5utr_length_pairwise.tsv", sep="\t")
    desc = pd.read_csv(HERE.parent / "data" / "panelC_5utr_length_descriptives.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(6, 4))
    plot_df = long.assign(
        y=np.log10(np.clip(long["utr5_length_nt"], 1, None))
    )
    data = [plot_df.loc[plot_df["group"] == g, "y"].to_numpy() for g in GROUP_ORDER]

    sns.violinplot(data=plot_df, x="group", y="y", order=GROUP_ORDER,
                   hue="group", palette=COLORS, legend=False,
                   cut=0, inner="quart", density_norm="width",
                   linewidth=0.9, ax=ax)
    for patch in ax.collections:
        patch.set_edgecolor(C_TITLE); patch.set_alpha(0.78)
    for line in ax.lines:
        line.set_color(C_TITLE); line.set_linewidth(1.0)

    ymax = max(d.max() for d in data)
    y_top = ymax + 0.18
    spacing = 0.50
    pairs = [("NMD+/PTC-", "NMD+/PTC+", 0, 1, y_top),
             ("NMD+/PTC-", "Control",  1, 2, y_top + spacing),
             ("NMD+/PTC+", "Control",  0, 2, y_top + 2 * spacing)]
    for gx, gy, ix, iy, yy in pairs:
        row = stats[((stats["group_x"] == gx) & (stats["group_y"] == gy)) |
                    ((stats["group_x"] == gy) & (stats["group_y"] == gx))]
        if row.empty: continue
        ax.plot([ix, ix, iy, iy], [yy - 0.04, yy, yy, yy - 0.04],
                color=C_TITLE, linewidth=0.6)
        ax.text((ix + iy) / 2, yy + 0.02, sig_marker(row.iloc[0]["wilcox_p"]),
                ha="center", va="bottom", fontsize=BODY_FS, color=C_TITLE)

    n_per_group = {g: int(desc.loc[desc["group"] == g, "n"].iloc[0]) for g in GROUP_ORDER}
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels([f"{g}\n(n={n_per_group[g]:,})" for g in GROUP_ORDER],
                       fontsize=LABEL_FS, color=C_TITLE)
    ax.set_xlabel("")
    ax.set_ylabel("5'UTR length (nt, log10)", fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=C_TITLE)
    ax.set_ylim(bottom=0, top=y_top + 3 * spacing)

    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_color(C_AXIS); ax.spines[s].set_linewidth(0.8)

    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.20, top=0.96)
    return fig


def main():
    fig = build_figure()
    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)
    out_dir = HERE.parent
    fig.savefig(out_dir / "figure4_panelC_5utr_length_refaug.pdf", facecolor="white")
    fig.savefig(out_dir / "figure4_panelC_5utr_length_refaug.png", dpi=300, facecolor="white")
    print("Saved: figure4_panelC_5utr_length_refaug.{pdf,png}")


if __name__ == "__main__":
    main()

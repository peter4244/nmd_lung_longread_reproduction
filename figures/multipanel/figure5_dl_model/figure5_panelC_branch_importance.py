"""
Figure 5 Panel C — Branch-level Shapley importance (Structural / Stop / ATG).

Mirrors orf_model_report fig2a_branch_importance: KernelSHAP with 500
background, exact 2^3 coalitions, on NMD test transcripts only. Bars show
mean |SHAP| per branch; labels show the share-of-total percentage.

Data:
  ./data/panelC_branch_importance.tsv
    Columns: branch, mean_abs_shap, pct, n_nmd
    Source: kernel_shap_branch_atg500_stop500.tsv (4ct model, results_4ct/).
    Exported by ./data_export.py.

Style: publications overlay — 6×4 cell at 1.5:1, no in-panel title,
no bbox_inches="tight" on save, validator wired into main().
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
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

# Branch palette — matches orf_model_report fig2a_branch_importance colours.
BRANCH_COLORS = {
    "Structural": "#d95f02",
    "Stop":       "#FF6B6B",
    "ATG":        "#4ECDC4",
}
C_TITLE = "#222222"
C_AXIS = "#555555"


def load_data():
    return pd.read_csv(HERE.parent / "data" / "panelC_branch_importance.tsv", sep="\t")


def build_figure(df):
    fig, ax = plt.subplots(figsize=(6, 4))

    branches = df["branch"].tolist()  # Structural / Stop / ATG (already ranked)
    vals = df["mean_abs_shap"].to_numpy()
    pcts = df["pct"].to_numpy()
    colors = [BRANCH_COLORS[b] for b in branches]

    x = range(len(branches))
    bars = ax.bar(x, vals, color=colors, edgecolor="none", width=0.62)

    # In-bar percentage labels (white, bold)
    for i, (v, p) in enumerate(zip(vals, pcts)):
        # Show pct above the bar; numeric value just below the top inside.
        ax.text(i, v + max(vals) * 0.03, f"{p:.1f}%",
                ha="center", va="bottom",
                fontsize=BODY_FS, fontweight="bold", color=C_TITLE)
        ax.text(i, v / 2, f"{v:.2f}",
                ha="center", va="center",
                fontsize=BODY_FS - 2, fontweight="bold", color="white")

    ax.set_xticks(list(x))
    ax.set_xticklabels(branches, fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=C_TITLE)

    ax.set_ylabel("Mean |SHAP|", fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylim(0, max(vals) * 1.20)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.16, top=0.96)
    return fig


def main():
    df = load_data()
    fig = build_figure(df)

    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)

    out_dir = HERE.parent
    fig.savefig(out_dir / "figure5_panelC_branch_importance.pdf", facecolor="white")
    fig.savefig(out_dir / "figure5_panelC_branch_importance.png", dpi=300, facecolor="white")
    print("Saved: figure5_panelC_branch_importance.pdf and .png")


if __name__ == "__main__":
    main()

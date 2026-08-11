"""
Figure 5 Panel D — Per-feature DeepSHAP importance for the structural branch.

Mirrors orf_model_report fig2b_structural_importance: mean |SHAP| per
structural feature on NMD test transcripts, averaged across 5 DeepSHAP
runs. Downstream EJC count dominates by 15×, consistent with the
manuscript's "EJC count was by far the most important" claim.

Data:
  ./data/panelD_structural_features.tsv
    Columns: feature (friendly label), channel (raw), mean_abs_shap
    Source: deepshap_summary_atg500_stop500_run{1..5}.tsv (4ct).
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

# Structural-branch colour matches Panel E ("Structural") bar.
C_BAR = "#d95f02"
C_TITLE = "#222222"
C_AXIS = "#555555"


def load_data():
    return pd.read_csv(HERE.parent / "data" / "panelD_structural_features.tsv", sep="\t")


def build_figure(df):
    fig, ax = plt.subplots(figsize=(6, 4))

    # Reverse for matplotlib's bottom-up barh ordering (largest at top)
    df_plot = df.iloc[::-1].reset_index(drop=True)
    y = range(len(df_plot))

    ax.barh(y, df_plot["mean_abs_shap"], color=C_BAR, edgecolor="none", height=0.62)

    # End-of-bar value labels
    xmax = df_plot["mean_abs_shap"].max()
    for i, v in enumerate(df_plot["mean_abs_shap"]):
        ax.text(v + xmax * 0.015, i, f"{v:.3f}",
                ha="left", va="center",
                fontsize=BODY_FS, color=C_TITLE)

    ax.set_yticks(list(y))
    ax.set_yticklabels(df_plot["feature"], fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="x", labelsize=BODY_FS, colors=C_TITLE)

    ax.set_xlabel("Mean |SHAP|", fontsize=BODY_FS, color=C_TITLE)
    ax.set_xlim(0, xmax * 1.18)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    fig.subplots_adjust(left=0.42, right=0.97, bottom=0.18, top=0.96)
    return fig


def main():
    df = load_data()
    fig = build_figure(df)

    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)

    out_dir = HERE.parent
    fig.savefig(out_dir / "figure5_panelD_structural_features.pdf", facecolor="white")
    fig.savefig(out_dir / "figure5_panelD_structural_features.png", dpi=300, facecolor="white")
    print("Saved: figure5_panelD_structural_features.pdf and .png")


if __name__ == "__main__":
    main()

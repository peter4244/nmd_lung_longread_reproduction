"""
Figure 3 Panel C — Splice event prevalence NMD vs Control.

Paired barplot: percent of pairs containing each splice event type in NMD
(c2) vs Control (c4) gene-matched pair sets, ordered by NMD prevalence
descending. Significance stars above each pair come from a per-event-type
Fisher's exact test (precomputed in 05_final_report_mashr.Rmd table 1b).

Headline result: skipped exons (SE) are ~2× as frequent in NMD pairs
(44.2%) as in Control (21.2%, p ~ 10⁻⁸¹).

Data:
  ./data/panelC_event_prevalence.tsv
    Columns: event_type, n_pairs_with_event_NMD, pct_of_pairs_NMD,
             n_pairs_with_event_Ctrl, pct_of_pairs_Ctrl, fisher_OR,
             fisher_p, direction
    Source: isopair_wrapper/tables/table1b_event_enrichment.csv, copied by
    ./data_export.R.

Style: HEADER_FS=18 bold, BODY_FS=14 regular. NMD coral, Control light blue.
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

C_NMD = "#ef8a62"
C_CTRL = "#67a9cf"
C_TITLE = "#222222"
C_AXIS = "#555555"

EVENT_LABELS = {
    "SE": "Skipped exon",
    "Alt_TSS": "Alt TSS",
    "Alt_TES": "Alt TES",
    "IR": "Intron retention",
    "A5SS": "A5SS",
    "A3SS": "A3SS",
    "Partial_IR_5": "Partial IR 5′",
    "Partial_IR_3": "Partial IR 3′",
    "Missing_Internal": "Missing internal",
}


from p_to_stars import p_to_stars as sig_stars  # noqa: E402 — unified 4-tier scheme


def load_data():
    df = pd.read_csv(HERE.parent / "data" / "panelC_event_prevalence.tsv", sep="\t")
    df = df.sort_values("pct_of_pairs_NMD", ascending=False).reset_index(drop=True)
    df["label"] = df["event_type"].map(EVENT_LABELS).fillna(df["event_type"])
    df["stars"] = df["fisher_p"].apply(sig_stars)
    return df


def build_figure(df):
    # Composite cell target: 6.0 × 4.0 in (1.5:1). 10 rotated categories on x.
    fig, ax = plt.subplots(figsize=(6, 4))

    n = len(df)
    x = np.arange(n)
    bw = 0.38

    ax.bar(x - bw / 2, df["pct_of_pairs_NMD"], bw,
           color=C_NMD, edgecolor="none", label="NMD")
    ax.bar(x + bw / 2, df["pct_of_pairs_Ctrl"], bw,
           color=C_CTRL, edgecolor="none", label="Control")

    # Significance labels above the tallest bar in each pair
    y_max = max(df["pct_of_pairs_NMD"].max(), df["pct_of_pairs_Ctrl"].max())
    star_pad = y_max * 0.03
    for i, row in df.iterrows():
        top = max(row["pct_of_pairs_NMD"], row["pct_of_pairs_Ctrl"]) + star_pad
        ax.text(i, top, row["stars"], ha="center", va="bottom",
                fontsize=BODY_FS, color=C_TITLE)

    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], rotation=35, ha="right",
                       fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=C_TITLE)

    ax.set_ylabel("Pairs with event (%)", fontsize=BODY_FS, color=C_TITLE)
    # No in-panel title — caption / composite letter label carries it.

    ax.set_ylim(0, y_max * 1.15)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    legend = ax.legend(fontsize=BODY_FS, loc="upper right", frameon=False)
    for text in legend.get_texts():
        text.set_color(C_TITLE)

    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.38, top=0.96)
    return fig


def main():
    df = load_data()
    fig = build_figure(df)
    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)
    out_dir = HERE.parent
    fig.savefig(out_dir / "figure3_panelC_event_prevalence.pdf",
                facecolor="white")
    fig.savefig(out_dir / "figure3_panelC_event_prevalence.png",
                dpi=300, facecolor="white")
    print("Saved: figure3_panelC_event_prevalence.pdf and .png")


if __name__ == "__main__":
    main()

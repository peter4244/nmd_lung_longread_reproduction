"""
Figure 3 Panel E — PTC-introducing event attribution.

Paired barplot: among events Isopair attributes as the cause of a PTC in
NMD comparators (NMD-causing), and among ALL splice events in Control
comparators (the background), what fraction is each event type?

Headline: skipped exon (SE) accounts for 57.0% of PTC-causing events vs
13.2% of Control events (~4× enrichment, p ≈ 2×10⁻²⁰³). A3SS (2.2×) and
A5SS (2.5×) are also significantly enriched among PTC-causing events;
terminal events (Alt TES) and 5'UTR-side partial-IR are depleted.

Scope: ENST-only effectively_ptc (n=1,489 NMD comparators; 1,434 with direct
event attribution) vs all detailed_events in ENST-only pop_traceable Control
(3,218 events across 1,286 c4 pairs). Attribution is ref-AUG-derived
uniformly — populated by 05r_ref_atg_analysis.R for ALL effectively_ptc
pairs (no TD2-PTC-stop path).

Data:
  ./data/panelE_ptc_event_attribution.tsv
    Columns: event_type, n_ptc_events, pct_of_ptc, n_ctrl_events, pct_ctrl,
             enrichment, fisher_p, direction
    Source: panel_e_compute.R; output of ref-AUG-derived attribution.

Style: HEADER_FS=18 bold, BODY_FS=14 regular. PTC-causing (NMD-attributed)
coral, Control all-events light blue — matches Panels A/B/C/D.
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
    df = pd.read_csv(HERE.parent / "data" / "panelE_ptc_event_attribution.tsv", sep="\t")
    df = df.sort_values("pct_of_ptc", ascending=False).reset_index(drop=True)
    df["label"] = df["event_type"].map(EVENT_LABELS).fillna(df["event_type"])
    df["stars"] = df["fisher_p"].apply(sig_stars)
    return df


def build_figure(df):
    # Composite cell target: 6.0 × 4.0 in (1.5:1). 11 rotated categories + title+subtitle.
    fig, ax = plt.subplots(figsize=(6, 4))

    n = len(df)
    x = np.arange(n)
    bw = 0.38

    ax.bar(x - bw / 2, df["pct_of_ptc"], bw,
           color=C_NMD, edgecolor="none", label="PTC-causing (NMD)")
    ax.bar(x + bw / 2, df["pct_ctrl"], bw,
           color=C_CTRL, edgecolor="none", label="All events (Control)")

    y_max = max(df["pct_of_ptc"].max(), df["pct_ctrl"].max())
    star_pad = y_max * 0.03
    for i, row in df.iterrows():
        top = max(row["pct_of_ptc"], row["pct_ctrl"]) + star_pad
        ax.text(i, top, row["stars"], ha="center", va="bottom",
                fontsize=BODY_FS, color=C_TITLE)

    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], rotation=35, ha="right",
                       fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=C_TITLE)

    ax.set_ylabel("% of events", fontsize=BODY_FS, color=C_TITLE)

    # No in-panel title / subtitle — caption / composite letter label carries it.
    # The "PTCs identified: 1,489 NMD | 207 Control" (ENST-only) provenance
    # line moves to the figure caption per the grant titleless convention.

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
    fig.savefig(out_dir / "figure3_panelE_ptc_event_attribution.pdf",
                facecolor="white")
    fig.savefig(out_dir / "figure3_panelE_ptc_event_attribution.png",
                dpi=300, facecolor="white")
    print("Saved: figure3_panelE_ptc_event_attribution.pdf and .png")


if __name__ == "__main__":
    main()

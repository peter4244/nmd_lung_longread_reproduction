"""SF30 — GAIN direction by splice event type: NMD vs Control.

Grouped bar chart, one bar-pair per event type. The bar height is the
percentage of that event type's occurrences (across all detected events in
the pop_BC pair set) where the comparator is the GAIN direction relative to
its reference isoform (i.e., comparator has the additional sequence).
Sorted by NMD-side %GAIN, descending.

Significance stars: Fisher's exact test on the 2×2 (GAIN/LOSS × NMD/Control)
count table for each event type.

Data source: `table1c_gain_loss_events.csv` (rebuilt by this dir's data_export.R
on floored pop_BC-matched profiles; the legacy-Rmd producer is deprecated)
(from `05_final_report_mashr.Rmd` § "Gain/Loss Direction by Event Type").
"""

from __future__ import annotations
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
LIB  = HERE.parents[1] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    render_and_validate,
    docx_body_fs,
    NMD_COLOR,
    CONTROL_COLOR,
    TITLE_C,
    HEADER_FS,
)

apply_ggplot_rcparams()

NATIVE_W = 11.0
BODY_FS = docx_body_fs(NATIVE_W)

DATA = HERE / "data"


def sig_marker(p):
    if pd.isna(p): return ""
    if p < 1e-10:  return "****"
    if p < 1e-4:   return "***"
    if p < 1e-3:   return "**"
    if p < 0.05:   return "*"
    return "n.s."


# Cosmetic labels — replace underscores with spaces for readability
DISPLAY_LABEL = {
    "Alt_TSS":          "Alt TSS",
    "Alt_TES":          "Alt TES",
    "A5SS":             "A5SS",
    "A3SS":             "A3SS",
    "SE":               "SE",
    "Missing_Internal": "Missing Internal",
    "Partial_IR_5":     "Partial IR 5′",
    "Partial_IR_3":     "Partial IR 3′",
    "IR":               "Full IR",
    "IR_diff_5":        "IR diff 5′",
    "IR_diff_3":        "IR diff 3′",
    "IR_diff_5_3":      "IR diff 5′+3′",
}


def main():
    df = pd.read_csv(DATA / "table1c_gain_loss_events.csv")
    # Drop event types with no events in one arm (uninformative empty bar-pair;
    # e.g. IR diff 5'+3' has zero Control events).
    df = df[(df["n_NMD_GAIN"] + df["n_NMD_LOSS"] > 0) &
            (df["n_Control_GAIN"] + df["n_Control_LOSS"] > 0)].copy()
    df = df.sort_values("pct_NMD_GAIN", ascending=False).reset_index(drop=True)

    events = df["event_type"].tolist()
    x = np.arange(len(events))
    bar_width = 0.38

    fig, ax = plt.subplots(figsize=(NATIVE_W, 5.4))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.85, bottom=0.33)
    style_axes_ggplot(ax, xgrid=False, ygrid=True)

    bars_nmd = ax.bar(x - bar_width / 2, df["pct_NMD_GAIN"], width=bar_width,
                       color=NMD_COLOR, edgecolor="white", linewidth=0.6,
                       label="NMD", zorder=3)
    bars_ctrl = ax.bar(x + bar_width / 2, df["pct_Control_GAIN"], width=bar_width,
                       color=CONTROL_COLOR, edgecolor="white", linewidth=0.6,
                       label="Control", zorder=3)

    # Significance annotations above each bar-pair
    for i, row in df.iterrows():
        peak = max(row["pct_NMD_GAIN"], row["pct_Control_GAIN"])
        marker = sig_marker(row["fisher_p"])
        if not marker:
            continue
        ax.text(x[i], peak + 3, marker,
                ha="center", va="bottom",
                fontsize=BODY_FS, color=TITLE_C, fontweight="bold",
                zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels([DISPLAY_LABEL.get(e, e) for e in events],
                       fontsize=BODY_FS - 1, color=TITLE_C,
                       rotation=90, ha="right", va="center",
                       rotation_mode="anchor")
    ax.set_xlabel("")
    ax.set_ylabel("Events with Gain direction (%)", fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(-0.7, len(events) - 0.3)

    ax.legend(loc="upper right", frameon=True, facecolor="white",
              edgecolor="none", fontsize=BODY_FS)

    # No overall figure title or subtitle — caption carries both roles (Yul-style).
    render_and_validate(fig, HERE / "figure_sf30_gain_direction_by_event",
                        native_width_in=NATIVE_W)
    print(f"  events: {len(events)}")


if __name__ == "__main__":
    main()

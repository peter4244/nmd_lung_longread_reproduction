"""SF37 — Absolute Shapley values across the AUG and stop-codon windows.

Two-panel supplemental figure showing per-position |SHAP| aggregated across
all input channels, comparing NMD-susceptible vs Control transcript classes,
across the two 1000-nt sequence windows the DL model consumes.

Data source: the DEPOSITED shap_profile_{atg,stop}_joint_atg1000_stop1000_seed42.tsv
(from 09b_export_subgroup_profiles.py, 5-run pooled DeepSHAP).

Style: ggplot-mimic theme (grey panel + white gridlines) to match SF1-SF23.
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
    facet_header,
    render_and_validate,
    docx_body_fs,
    NMD_COLOR,
    CONTROL_COLOR,
    TITLE_C,
    HEADER_FS,
    panel_letter,
)

apply_ggplot_rcparams()

NATIVE_W = 11.5
BODY_FS = docx_body_fs(NATIVE_W)

# REPOINTED 2026-08-06 onto the deposited 1000nt model. This read HERE/"data", a directory that
# did not exist, naming the superseded 500nt files -- so SF37 has not been renderable. It now reads
# the deposit directly (same four-column schema, twice the window) and STOPS if the file is absent
# rather than falling back to results_4ct, which would redraw the published figure.
DEPOSIT = HERE.parents[2] / "data_deposit" / "source_data" / "model"
TAG = "atg1000_stop1000_seed42"


def load_window(branch):
    path = DEPOSIT / f"shap_profile_{branch}_joint_{TAG}.tsv"
    if not path.exists():
        raise SystemExit(f"missing deposited SHAP profile: {path}\n"
                         f"  No fallback to results_4ct: that is the superseded 500nt vintage.")
    df = pd.read_csv(path, sep="\t")
    # Aggregate across channels: total |SHAP| per position, per class.
    agg = df.groupby("relative_position", as_index=False)[["nmd_abs_shap", "ctrl_abs_shap"]].sum()
    agg = agg.sort_values("relative_position").reset_index(drop=True)
    return agg


def smooth(series, window=11):
    """Rolling-mean smoothing — the per-position raw signal has trinucleotide-
    frame periodicity that reads as noise at figure scale; a small rolling
    average clarifies the class-level pattern without hiding real features."""
    return series.rolling(window=window, center=True, min_periods=1).mean()


def draw_panel(ax, agg, *, codon_label):
    style_axes_ggplot(ax)
    ax.plot(agg["relative_position"], smooth(agg["nmd_abs_shap"]),
            color=NMD_COLOR, linewidth=1.6, label="NMD susceptible", zorder=3)
    ax.plot(agg["relative_position"], smooth(agg["ctrl_abs_shap"]),
            color=CONTROL_COLOR, linewidth=1.6, label="Control", zorder=3)
    ax.axvline(0, color="#555555", linewidth=0.8, linestyle=":", zorder=2)
    ax.set_xlabel(f"Position relative to {codon_label} (nt)", fontsize=BODY_FS)
    ax.set_ylabel("Total |SHAP|", fontsize=BODY_FS)
    ax.set_xlim(agg["relative_position"].min(), agg["relative_position"].max())
    # Ticks span the ACTUAL window. These were [-200..200], sized for the 500nt model, and would
    # have labelled only the middle two-fifths of a 1000nt window without erroring.
    ax.set_xticks([-400, -200, 0, 200, 400])
    ax.set_ylim(bottom=0)
    return ax


def main():
    atg = load_window("atg")
    stop = load_window("stop")

    fig, axes = plt.subplots(1, 2, figsize=(NATIVE_W, 4.8))
    fig.subplots_adjust(left=0.09, right=0.96, top=0.82, bottom=0.22, wspace=0.36)

    draw_panel(axes[0], atg, codon_label="start codon")
    draw_panel(axes[1], stop, codon_label="stop codon")

    facet_header(axes[0], "Start codon window", height=0.08)
    facet_header(axes[1], "Stop codon window",   height=0.08)

    # x=-0.12, not -0.07: the 1000nt window's |SHAP| totals need a wider tick label, which
    # collided with panel B's letter at the old offset.
    panel_letter(axes[0], "A", x=-0.18, y=1.14)
    panel_letter(axes[1], "B", x=-0.18, y=1.14)

    # Shared legend in the bottom margin — the |SHAP| curves fill both
    # panels, so an in-panel legend can't clear the data; both panels share
    # the same two series, so one figure-level legend below is cleaner.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False,
               fontsize=BODY_FS, bbox_to_anchor=(0.5, 0.01))

    # Facet headers extend past the axis edges by design — the per-axis
    # text-vs-spine check fires a false positive on them.
    render_and_validate(fig, HERE / "figure_sf37_shap_across_windows",
                        native_width_in=NATIVE_W, strict_per_axis=False)
    print(f"  n positions ATG: {len(atg)}   STOP: {len(stop)}")


if __name__ == "__main__":
    main()

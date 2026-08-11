"""
Supplemental Figure — Head-to-head comparison of our NMD predictor against
two variant-level NMD prediction models (NMDetective-B, Lindeboom 2019;
NMDEP-rule baseline, Saadat 2025) on the same cohort of long-read-observed
isoforms, using the mashr posterior-mean log fold-change under SMG1 inhibition
as the continuous gold standard.

Two panels:

  (A) Predicted score vs gold-standard log fold-change, one subplot per
      model, points coloured by PTC subclass (NMD+/PTC+, NMD+/PTC−, Control).
      Spearman correlation annotated per panel.

  (B) Within-subgroup Spearman correlation per model (grouped bar chart;
      bars = subclass × model). Shows where each model's predictions hold
      up vs break down.

Data: analysis/predictor_comparison/ -> $OUT/predictor_comparison/
      per_isoform_scores_2026.8.6.tsv and metrics_summary_2026.8.6.tsv

Style: canonical NMD palette (peach/yellow/blue for PTC+/PTC−/Control);
       Arial; validator-clean.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    panel_letter,
    docx_body_fs,
    render_and_validate,
)
apply_ggplot_rcparams()

NATIVE_W = 12.0
BODY_FS = docx_body_fs(NATIVE_W)

# 2026-07-26: was HERE.parents[3]/"code"/"nmd_predictor_comparison" -- the LEGACY layout, which
# does not exist in this repo (it is analysis/predictor_comparison), so this figure could not have
# found its inputs at all. The chain writes to .P$OUT/predictor_comparison.
import sys as _s
_r = HERE
while not (_r / "config" / "paths.yml").exists():
    if _r.parent == _r:
        raise RuntimeError("repo root not found")
    _r = _r.parent
_s.path.insert(0, str(_r / "python"))
from nmd_paths import P as _P
ANALYSIS = _P["OUT"] / "predictor_comparison"
# 2026.8.6: the predictor-comparison chain was rerun against the DEPOSITED 1000nt model. The
# 2026.7.11 files it used to name were built from the superseded 500nt run and do not exist here.
SCORES = ANALYSIS / "per_isoform_scores_2026.8.6.tsv"
METRICS = ANALYSIS / "metrics_summary_2026.8.6.tsv"
for _f in (SCORES, METRICS):
    if not _f.exists():
        raise SystemExit(f"missing {_f}\n  Run analysis/predictor_comparison/01..04 first.")

# Canonical palette (project_nmd_figure_palette)
GROUP_LABEL = {
    "NMD+/PTC+": "NMD+/PTC+",
    "NMD+/PTC-": "NMD+/PTC−",
    "Control":   "Control",
}
# Local color deviation from the project palette: project canonical PTC-
# is #d95f02 (dark orange), but in scatter plots with ~1k PTC+ orange
# points overlapping ~100 PTC- points the two oranges blur. Swap PTC-
# to a contrasting purple just in this figure.
GROUP_COLOR = {
    "NMD+/PTC+": "#ef8a62",
    "NMD+/PTC-": "#762a83",
    "Control":   "#67a9cf",
}
# Plot order: heaviest-overplotted group first (Control), most-distinctive
# group last so it lands on top of the overplot stack (PTC-).
SUBCLASS_ORDER = ["NMD+/PTC+", "NMD+/PTC-", "Control"]
SUBCLASS_PLOT_ORDER = ["Control", "NMD+/PTC+", "NMD+/PTC-"]

MODEL_LABEL = {
    "nmdetective_b": "NMDetective-B\n(Lindeboom 2019)",
    "nmdep_rule":    "NMDEP rule baseline\n(Saadat 2025)",
    "our_model":     "Our model\n(this work)",
}
MODEL_ORDER = ["nmdetective_b", "nmdep_rule", "our_model"]
MODEL_SCORE_COL = {
    "nmdetective_b": "nmdetective_b_score",
    "nmdep_rule":    "nmdep_rule_score",
    "our_model":     "our_model_prob",
}


def load():
    s = pd.read_csv(SCORES, sep="\t")
    m = pd.read_csv(METRICS, sep="\t")
    # Head-to-head intersection restricted to H5 split == "test" (chr1/3/5/7
    # holdout the deep-learning model never saw at training time). This is
    # the methodologically correct frame for performance comparison: the
    # rule-based models also see these isoforms for the first time.
    mask = (
        s["nmdetective_b_score"].notna()
        & s["nmdep_rule_score"].notna()
        & s["our_model_prob"].notna()
        & (s["h5_split"] == "test")
    )
    return s[mask].copy(), m


def panel_a_scatter(ax, df, model_key):
    """One subplot in Panel A: predicted score vs gold-standard log-FC."""
    style_axes_ggplot(ax)
    col = MODEL_SCORE_COL[model_key]
    for sub in SUBCLASS_PLOT_ORDER:
        d = df[df["subclass"] == sub]
        # PTC- is the smallest, most-overplotted group — draw with a slight
        # edge + higher alpha so the 95 dots stay visible against PTC+.
        if sub == "NMD+/PTC-":
            ax.scatter(
                d[col], d["mashr_posterior_mean_logfc"],
                color=GROUP_COLOR[sub], alpha=0.85, s=22,
                edgecolor="white", linewidth=0.3, label=GROUP_LABEL[sub],
            )
        else:
            ax.scatter(
                d[col], d["mashr_posterior_mean_logfc"],
                color=GROUP_COLOR[sub], alpha=0.5, s=16,
                edgecolor="none", label=GROUP_LABEL[sub],
            )
    # Pooled Spearman
    sp = df[[col, "mashr_posterior_mean_logfc"]].corr(method="spearman").iloc[0, 1]
    ax.text(
        0.04, 0.95, f"Spearman = {sp:.2f}",
        transform=ax.transAxes, fontsize=BODY_FS, color="#222222",
        va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#999", lw=0.5),
    )
    # No subpanel title — fold the model identity into the x-axis label so the
    # reader can still tell the three sub-panels apart (Yul-style).
    ax.set_xlabel(
        f"Predicted NMD score\n{MODEL_LABEL[model_key]}",
        fontsize=BODY_FS, color="#222222",
    )
    if model_key == "our_model":
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    else:
        ax.set_xticks([0.0, 0.2, 0.4, 0.6])
    ax.tick_params(axis="both", labelsize=BODY_FS - 2, colors="#222222")


def panel_b_bars(ax, metrics, head_to_head_n):
    """Grouped bar chart of per-subclass Spearman within head-to-head set."""
    style_axes_ggplot(ax, xgrid=False, ygrid=True)
    rows = []
    for model in MODEL_ORDER:
        for sub in SUBCLASS_ORDER:
            key = f"head-to-head:test:{sub}"
            r = metrics[(metrics["model"] == model) & (metrics["stratum"] == key)]
            if r.empty:
                continue
            rows.append({
                "model":   MODEL_LABEL[model].replace("\n", " "),
                "subclass": sub,
                "spearman": r.iloc[0]["spearman"],
                "n":         int(r.iloc[0]["n"]),
            })
    d = pd.DataFrame(rows)

    x = np.arange(len(MODEL_ORDER))
    w = 0.27
    for i, sub in enumerate(SUBCLASS_ORDER):
        vals_by_model = []
        for m in MODEL_ORDER:
            row = d[(d["model"] == MODEL_LABEL[m].replace("\n", " "))
                    & (d["subclass"] == sub)]
            vals_by_model.append(row.iloc[0]["spearman"] if not row.empty else np.nan)
        bars = ax.bar(
            x + (i - 1) * w, vals_by_model, width=w,
            color=GROUP_COLOR[sub], edgecolor="black", linewidth=0.5,
            label=f"{GROUP_LABEL[sub]} (n = {int(d.loc[d['subclass']==sub, 'n'].iloc[0])})",
        )
        # Bar labels
        for rect, v in zip(bars, vals_by_model):
            if np.isnan(v):
                continue
            ax.text(
                rect.get_x() + rect.get_width() / 2,
                v + (0.03 if v >= 0 else -0.06),
                f"{v:.2f}",
                ha="center", va="bottom" if v >= 0 else "top",
                fontsize=BODY_FS - 3, color="#222222",
            )

    ax.axhline(0, color="#555555", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(
        ["NMDetective-B\n(Lindeboom 2019)",
         "NMDEP rule baseline\n(Saadat 2025)",
         "Our model\n(this work)"],
        fontsize=BODY_FS - 2, color="#222222",
    )
    ax.set_ylabel("Spearman correlation with\nmashr log fold-change (SMG1i)",
                  fontsize=BODY_FS, color="#222222")
    ax.set_ylim(-0.35, 0.75)
    ax.set_yticks([-0.2, 0.0, 0.2, 0.4, 0.6])
    ax.tick_params(axis="y", labelsize=BODY_FS - 2, colors="#222222")
    ax.legend(
        loc="upper left", fontsize=BODY_FS - 2, frameon=False,
        handlelength=1.5, handletextpad=0.5, borderaxespad=0.4,
        title=f"Within-subclass (n = {head_to_head_n} total)",
        title_fontsize=BODY_FS - 2,
    )


def build_figure():
    df, metrics = load()
    n_hh = len(df)

    # 2-row layout matches the docx-approved SF43: Panel A (3 scatter
    # subpanels) on TOP row spanning full width, Panel B (grouped bar
    # chart) on BOTTOM row also full width. Restored 2026-07-10 after
    # an earlier session inadvertently flattened to 1-row wide layout.
    fig = plt.figure(figsize=(12, 12))
    gs_outer = fig.add_gridspec(
        nrows=2, ncols=1,
        height_ratios=[1.0, 1.0],
        hspace=0.35, left=0.12, right=0.96, bottom=0.06, top=0.94,
    )
    gs_a = gs_outer[0, 0].subgridspec(nrows=1, ncols=3, wspace=0.45)

    axes_a = [fig.add_subplot(gs_a[0, i]) for i in range(3)]
    ax_b = fig.add_subplot(gs_outer[1, 0])

    # Panel A: scatter per model
    for ax, m in zip(axes_a, MODEL_ORDER):
        panel_a_scatter(ax, df, m)
    axes_a[0].set_ylabel(
        "Gold standard: mashr posterior mean\nlog fold-change under SMG1i",
        fontsize=BODY_FS, color="#222222",
    )
    # Shared horizontal legend ABOVE Panel A (between the panel label
    # row and the subplot title row). Keeps the legend out of the
    # data area entirely.
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color="w",
                   markerfacecolor=GROUP_COLOR[s], markersize=7,
                   markeredgecolor="white" if s == "NMD+/PTC-" else "none",
                   label=GROUP_LABEL[s])
        for s in SUBCLASS_ORDER
    ]
    fig.legend(
        handles=handles, loc="upper center",
        bbox_to_anchor=(0.5, 0.98),
        ncol=3, fontsize=BODY_FS - 2, frameon=False,
        handlelength=1.2, handletextpad=0.4, columnspacing=1.6,
        title=None,
    )

    # Panel B: stratified bars
    panel_b_bars(ax_b, metrics, n_hh)

    # Panel labels — Panel A[0] is a narrow subplot in a 3-col subgrid,
    # so -0.15 axes-fraction is a small absolute distance; Panel B spans
    # the full figure width, so its letter needs a smaller axes-fraction
    # offset to stay inside the canvas.
    panel_letter(axes_a[0], "A", x=-0.15, y=1.08)
    panel_letter(ax_b, "B", x=-0.05, y=1.08)

    return fig, axes_a + [ax_b]


def main():
    fig, axes = build_figure()
    render_and_validate(fig, HERE.parent / "figure_s_model_comparison",
                        native_width_in=NATIVE_W)


if __name__ == "__main__":
    main()

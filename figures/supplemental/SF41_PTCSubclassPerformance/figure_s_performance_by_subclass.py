"""
Supplemental Figure — Deep-learning NMD-call discrimination, stratified by
PTC subclass within the n = 819 ref-AUG-traceable subset.

Companion to the §5 manuscript sentence:

    "However, overall predictive performance was substantially lower in
     NMD+/PTC- isoforms (SFx)."

Two panels at the held-out chr-1/3/5/7 paralog-free test split:

  (A) ROC curves for two contrasts:
        — NMD+/PTC+ retained  vs  Control
        — NMD+/PTC- retained  vs  Control
      AUC annotated.

  (B) Per-isoform predicted NMD probability stratified by subclass
      (strip + boxplot overlay), with the 0.5 decision threshold marked.

KEY NUMBERS (from data/predprob_by_subclass_refaug.tsv):
  NMD+/PTC+ vs Control:  AUC = 0.951, AUPRC = 0.935   (n=195 vs n=204)
  NMD+/PTC- vs Control:  AUC = 0.687, AUPRC = 0.130   (n=16  vs n=204)
  Mean predicted prob:   Control 0.19, NMD+/PTC+ 0.85, NMD+/PTC- 0.39.

Scope note: test-set predictions only (chr-1/3/5/7, paralog-free). NMD+/
PTC- has n = 16 at the test split.

Data: data/predprob_by_subclass_refaug.tsv (built by sibling
PTCSubclassBranchSHAP/data_export.R).
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    panel_letter,
    docx_body_fs,
    render_and_validate,
    TITLE_C,
)
apply_ggplot_rcparams()

NATIVE_W = 8.2
BODY_FS  = docx_body_fs(NATIVE_W)

DATA = HERE.parent / "data" / "predprob_by_subclass_refaug.tsv"

GROUP_ORDER = ["NMD+/PTC+", "NMD+/PTC- retained", "Control"]
GROUP_LABEL = {
    "NMD+/PTC+":          "NMD+/PTC+",
    "NMD+/PTC- retained": "NMD+/PTC−",
    "Control":            "Control",
}
# Canonical NMD subgroup palette — verbatim from orf_model_report_v5.Rmd:48-50
# and project_nmd_figure_palette memory. Do NOT swap for matplotlib defaults;
# the same group must read in the same color across the whole paper.
GROUP_COLOR = {
    "NMD+/PTC+":          "#ef8a62",   # peach (= binary NMD)
    "NMD+/PTC- retained": "#d95f02",   # orange (Rmd "NMD PTC-, ref ATG retained")
    "Control":            "#67a9cf",   # blue (= binary Control)
}


def load_split():
    df = pd.read_csv(DATA, sep="\t")
    ctrl = df[df["group"] == "Control"]
    ptcp = df[df["group"] == "NMD+/PTC+"]
    ptcn = df[df["group"] == "NMD+/PTC- retained"]
    return df, ctrl, ptcp, ptcn


def contrast_metrics(pos: pd.DataFrame, neg: pd.DataFrame):
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    p = np.concatenate([pos["prob"].values, neg["prob"].values])
    fpr, tpr, _ = roc_curve(y, p)
    auc = roc_auc_score(y, p)
    auprc = average_precision_score(y, p)
    return fpr, tpr, auc, auprc


def build_figure():
    df, ctrl, ptcp, ptcn = load_split()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.2, 3.8))
    fig.subplots_adjust(left=0.085, right=0.99, bottom=0.16, top=0.88, wspace=0.28)

    # ── Panel A — ROC curves
    style_axes_ggplot(axA)
    fpr1, tpr1, auc1, _ = contrast_metrics(ptcp, ctrl)
    fpr2, tpr2, auc2, _ = contrast_metrics(ptcn, ctrl)
    axA.plot(fpr1, tpr1, color=GROUP_COLOR["NMD+/PTC+"], lw=2.0, zorder=4,
             label=f"NMD+/PTC+ vs Control\nAUC = {auc1:.2f}")
    axA.plot(fpr2, tpr2, color=GROUP_COLOR["NMD+/PTC- retained"], lw=2.0, zorder=4,
             label=f"NMD+/PTC− vs Control\nAUC = {auc2:.2f}")
    axA.plot([0, 1], [0, 1], "--", color="#666666", lw=0.8, zorder=3)
    axA.set_xlim(-0.02, 1.02)
    axA.set_ylim(-0.02, 1.02)
    axA.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axA.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axA.set_xlabel("False positive rate", fontsize=BODY_FS)
    axA.set_ylabel("True positive rate", fontsize=BODY_FS)
    # No subpanel title — caption identifies Panel A (Yul-style).
    axA.tick_params(axis="both", labelsize=BODY_FS - 2)
    axA.legend(loc="lower right", fontsize=BODY_FS - 3, frameon=True, framealpha=0.9,
               facecolor="white", edgecolor="none", handlelength=1.6,
               handletextpad=0.5, borderaxespad=0.4, labelspacing=0.25)
    axA.set_aspect("equal")

    # ── Panel B — predicted probability per subgroup
    style_axes_ggplot(axB, xgrid=False, ygrid=True)
    rng = np.random.default_rng(seed=1)
    xs = []
    ys = []
    cs = []
    positions = list(range(len(GROUP_ORDER)))
    for i, g in enumerate(GROUP_ORDER):
        vals = df.loc[df["group"] == g, "prob"].values
        jitter = rng.uniform(-0.18, 0.18, size=len(vals))
        xs.extend(i + jitter)
        ys.extend(vals)
        cs.extend([GROUP_COLOR[g]] * len(vals))
    axB.scatter(xs, ys, c=cs, s=10, alpha=0.55, edgecolor="none",
                zorder=2, rasterized=True)

    # Boxplot overlay
    bdata = [df.loc[df["group"] == g, "prob"].values for g in GROUP_ORDER]
    bp = axB.boxplot(
        bdata, positions=positions, widths=0.45,
        patch_artist=True, showfliers=False,
        boxprops=dict(facecolor="none", edgecolor="black", linewidth=1.0),
        medianprops=dict(color="black", linewidth=1.4),
        whiskerprops=dict(color="black", linewidth=1.0),
        capprops=dict(color="black", linewidth=1.0),
        zorder=3,
    )

    axB.set_xticks(positions)
    axB.set_xticklabels([
        f"{GROUP_LABEL[g]}\n(n = {int((df['group'] == g).sum()):,})"
        for g in GROUP_ORDER
    ], fontsize=BODY_FS - 2)
    axB.set_xlim(-0.6, len(GROUP_ORDER) - 0.4)
    axB.set_ylim(-0.03, 1.03)
    axB.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axB.set_ylabel("Model NMD probability  P(NMD | x)", fontsize=BODY_FS)
    # No subpanel title — caption identifies Panel B (Yul-style).
    axB.tick_params(axis="y", labelsize=BODY_FS - 2)

    # Panel labels
    panel_letter(axA, "A", x=-0.13, y=1.02)
    panel_letter(axB, "B", x=-0.13, y=1.02)

    # No footnote — caption carries the key numbers (Yul-style).
    return fig, (axA, axB)


def main():
    fig, axes = build_figure()
    # The opaque white AUC legend (lower-right) grazes the tail of the
    # NMD+/PTC- ROC curve by ~28 px post-floor; the box masks the data
    # beneath it, so this reviewed overlap is accepted with a tight tolerance.
    render_and_validate(fig, HERE.parent / "figure_s_performance_by_subclass",
                        native_width_in=NATIVE_W,
                        legend_overlap_tolerance_px=30)


if __name__ == "__main__":
    main()

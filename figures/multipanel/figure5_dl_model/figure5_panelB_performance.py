"""
Figure 5 Panel B — Held-out model performance: ROC curve with PR-curve inset.

Best NMD-prediction model (ATG=500, STOP=500), evaluated on the
chr-1/3/5/7 paralog-free holdout test set. AUC and AUPRC annotated; PR
curve shown as an inset in the lower-right quadrant of the ROC.

Data:
  ./data/panelB_predictions.tsv   — per-test-transcript (label, prob)
  ./data/panelB_metrics.tsv       — n_test, n_nmd, auprc (and auc for cross-check)

Style: publications overlay — drops bbox_inches="tight"; validator wired
into main().
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def roc_curve(y_true, y_score):
    """ROC. Returns fpr, tpr (both monotonically increasing)."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    order = np.argsort(-y_score)
    yt = y_true[order]
    P = yt.sum()
    N = len(yt) - P
    tps = np.cumsum(yt)
    fps = np.cumsum(1 - yt)
    fpr = np.concatenate(([0.0], fps / N))
    tpr = np.concatenate(([0.0], tps / P))
    return fpr, tpr, None  # third item matches sklearn signature


def precision_recall_curve(y_true, y_score):
    """PR. Returns precision, recall arrays (precision monotonically updated)."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    order = np.argsort(-y_score)
    yt = y_true[order]
    P = yt.sum()
    tps = np.cumsum(yt)
    fps = np.cumsum(1 - yt)
    precision = tps / np.maximum(tps + fps, 1)
    recall = tps / P
    # Prepend (recall=0, precision=1) for plotting
    return np.concatenate(([1.0], precision)), np.concatenate(([0.0], recall)), None


def roc_auc_score(fpr, tpr):
    """Trapezoidal AUC over arrays sorted by fpr."""
    return float(np.trapezoid(tpr, fpr) if hasattr(np, "trapezoid") else np.trapz(tpr, fpr))

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

BODY_FS = 14
LABEL_FS = 12
INSET_FS = 10

# Curves: NMD-coral for both ROC and PR; subtle gray for the diagonal /
# baseline reference lines.
C_ROC = "#d6604d"
C_PR = "#d6604d"
C_DIAG = "#888888"
C_BASE = "#888888"
C_TITLE = "#222222"
C_AXIS = "#555555"


def load_data():
    preds = pd.read_csv(HERE.parent / "data" / "panelB_predictions.tsv", sep="\t")
    metrics = pd.read_csv(HERE.parent / "data" / "panelB_metrics.tsv", sep="\t").iloc[0]
    return preds, metrics


def build_figure(preds, metrics):
    fig, ax = plt.subplots(figsize=(6, 4))

    # ── Main: ROC ─────────────────────────────────────────────
    fpr, tpr, _ = roc_curve(preds["label"], preds["prob"])
    auc_val = roc_auc_score(fpr, tpr)

    ax.plot(fpr, tpr, color=C_ROC, linewidth=2.0, solid_capstyle="round")
    ax.plot([0, 1], [0, 1], color=C_DIAG, linestyle="--", linewidth=1.0, alpha=0.7)
    ax.fill_between(fpr, tpr, color=C_ROC, alpha=0.08)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("False positive rate", fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylabel("True positive rate", fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="both", labelsize=LABEL_FS, colors=C_TITLE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    # ── AUC / AUPRC annotation (two lines in the safe pocket above the
    # diagonal and below the curve; sample counts → caption to avoid
    # crowding) ──
    auprc = float(metrics["auprc"])
    txt = f"AUC = {auc_val:.2f}\nAUPRC = {auprc:.2f}"
    ax.text(
        0.10, 0.65, txt,
        transform=ax.transAxes, ha="left", va="top",
        fontsize=BODY_FS, color=C_TITLE,
    )

    # ── PR-curve inset (lower-right of ROC) ──
    # axes coords [x, y, w, h]: pushed to the lower-right so it sits in the
    # natural empty quadrant of a ROC plot. Slightly inset from the edges
    # to keep tick labels readable.
    pr_ax = ax.inset_axes([0.58, 0.10, 0.39, 0.39])

    precision, recall, _ = precision_recall_curve(preds["label"], preds["prob"])
    pr_ax.plot(recall, precision, color=C_PR, linewidth=1.6, solid_capstyle="round")
    # Random-classifier baseline = positive-class prevalence
    prev = float(np.mean(preds["label"]))
    pr_ax.axhline(prev, color=C_BASE, linestyle="--", linewidth=0.8, alpha=0.7)

    pr_ax.set_xlim(0, 1)
    pr_ax.set_ylim(0, 1.02)
    pr_ax.set_xlabel("Recall", fontsize=INSET_FS, color=C_TITLE, labelpad=2)
    pr_ax.set_ylabel("Precision", fontsize=INSET_FS, color=C_TITLE, labelpad=2)
    pr_ax.tick_params(axis="both", labelsize=INSET_FS - 1, colors=C_TITLE,
                      length=2, pad=2)
    pr_ax.set_xticks([0, 0.5, 1.0])
    pr_ax.set_yticks([0, 0.5, 1.0])
    for side in ("top", "right"):
        pr_ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        pr_ax.spines[side].set_color(C_AXIS)
        pr_ax.spines[side].set_linewidth(0.6)

    # NO CORNER LABEL. There was a bold "PR" at (0.05, 0.95), placed on the assumption that a
    # corner is empty space. That holds for a ROC and is exactly wrong for a precision-recall curve:
    # precision starts near 1.0 at low recall, so the top-left IS the data, and the better the model
    # the worse the collision. At AUPRC 0.82 the curve enters at ~0.95 and the label sat on it.
    #
    # It is deleted rather than moved, because it carried no information: this inset's axes already
    # read "Precision" and "Recall". Its only effect was to compete with the curve and to invite the
    # question of why it did not say AUPRC -- which is a different quantity, annotated in the main
    # panel. A label a careful reader has to decode has failed even when it is defensible.

    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.16, top=0.96)
    return fig


def main():
    preds, metrics = load_data()
    fig = build_figure(preds, metrics)

    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)

    # THE INSET IS VALIDATED TOO, for OVERLAPS only. It used to be skipped outright -- "the inset is
    # intentionally crowded" -- and that skip is why a label printed on top of the curve shipped:
    # the one axis nobody checked was the one with the defect. The premise was sound and the scope
    # was too wide. crowding_distance=0 keeps the density WARNINGS off, which is the part that would
    # be noise in a small inset, while hard text-on-data collisions still raise.
    # Resolve the inset from the PARENT's child_axes, not from fig.axes -- an inset created with
    # Axes.inset_axes does not appear in fig.axes at all (measured: len(fig.axes) == 1 here). The
    # first attempt indexed fig.axes[1] and the assertion caught it immediately, which is the only
    # reason this validates the inset rather than something else while reporting success.
    insets = fig.axes[0].child_axes
    assert len(insets) == 1, f"expected exactly one inset, found {len(insets)}"
    validate_figure_layout(fig, insets[0], crowding_distance=0, verbose=True)

    out_dir = HERE.parent
    fig.savefig(out_dir / "figure5_panelB_performance.pdf", facecolor="white")
    fig.savefig(out_dir / "figure5_panelB_performance.png", dpi=300, facecolor="white")
    print("Saved: figure5_panelB_performance.pdf and .png")


if __name__ == "__main__":
    main()

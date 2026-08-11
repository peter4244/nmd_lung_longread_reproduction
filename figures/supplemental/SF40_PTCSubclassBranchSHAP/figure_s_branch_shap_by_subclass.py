"""
Supplemental Figure — Branch-level KernelSHAP stratified by PTC subclass
across the ref-AUG-traceable subset (not test-set-restricted; post 25%
reference-share floor 2026-07-10). Shown: the 1,570 of 1,638 ref-AUG-traceable
isoforms that have a computed model branch attribution (68 fall outside the
model's scored cohort — see data_export.R join coverage).

Style: matches Figure 5 Panel C exactly — branch order Structural / Stop /
ATG, same colour palette, percent label above each bar, mean |SHAP|
inside each bar, no in-panel title. Three side-by-side mini-panels (one
per subgroup) share a common y-axis so heights compare directly.

Companion to the §5 manuscript sentence:

    "In NMD+/PTC- isoforms ... the model ... assigned twice as much
     importance to the start window compared to the NMD+/PTC+ isoforms."
     (manuscript 2026.7.17, §5)

KEY NUMBER (full cohort, no test-set restriction; post-floor):
  Mean within-isoform ATG-branch share:
    NMD+/PTC+          (n=735) →  9.1%
    NMD+/PTC- retained (n= 54) → 17.0%
    Control            (n=781) → 10.3%
  Ratio NMD+/PTC- vs NMD+/PTC+ ≈ 1.99×.

The Panel-C-style bar labels show share-of-subgroup-total (sum of
subgroup mean |SHAP|), giving:
    NMD+/PTC+   Structural 61.6% / Stop 29.4% / ATG  9.1%
    NMD+/PTC-   Structural 49.6% / Stop 33.5% / ATG 17.0%
    Control     Structural 63.0% / Stop 26.7% / ATG 10.3%

Manuscript ↔ figure now agree: the current §5 text says "twice as much
importance", matching the post-floor share-of-total ratio of 1.87×
(NMD+/PTC- 17.0% vs NMD+/PTC+ 9.1%). An earlier draft said "roughly three
times"; that was corrected to "twice" and this figure is consistent with the
corrected text. (Pre-floor the ratio was 2.21×.)

Data: data/branch_shap_by_subclass_refaug.tsv (built by data_export.R).
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    facet_header,
    panel_letter,
    render_and_validate,
    docx_body_fs,
    TITLE_C as C_TITLE,
    AXIS_C as C_AXIS,
)
apply_ggplot_rcparams()

NATIVE_W = 15.0
BODY_FS = docx_body_fs(NATIVE_W)
HEADER_FS = BODY_FS + 4

BRANCH_COLORS = {
    "Structural": "#d95f02",
    "Stop":       "#FF6B6B",
    "ATG":        "#4ECDC4",
}

BRANCH_ORDER = ["Structural", "Stop", "ATG"]

# Display labels: the start-codon branch is stored as "ATG" in the data but
# shown as "AUG" (RNA, not DNA).
BRANCH_DISPLAY = {"Structural": "Structural", "Stop": "Stop", "ATG": "Start"}

DATA = HERE.parent / "data" / "branch_shap_by_subclass_refaug.tsv"

GROUP_ORDER = ["NMD+/PTC+", "NMD+/PTC- retained", "Control"]
GROUP_LABEL = {
    "NMD+/PTC+":          "NMD+/PTC+",
    "NMD+/PTC- retained": "NMD+/PTC−",
    "Control":            "Control",
}


def per_group_summary(df: pd.DataFrame, group: str):
    """Return a DataFrame matching Figure 5 Panel C's per-branch table:
       columns branch, mean_abs_shap, pct (within-subgroup share-of-total),
       with branch order = Structural / Stop / ATG."""
    sub = df[df["group"] == group]
    means = (
        sub.groupby("branch")["abs_shap"].mean().reindex(BRANCH_ORDER).reset_index()
    )
    total = means["abs_shap"].sum()
    means["pct"] = 100 * means["abs_shap"] / total
    means.columns = ["branch", "mean_abs_shap", "pct"]
    means["n"] = int((df["group"] == group).shape[0]) if False else \
        int(df.loc[df["group"] == group, "comparator_isoform_id"].nunique())
    return means


def draw_panel(ax, summary: pd.DataFrame, group_label: str, n: int, y_top: float):
    """Render one Panel-C-style mini-panel on `ax` for one subgroup.

    Y-axis is within-subgroup share of total mean |SHAP|, in %. The mean
    |SHAP| absolute value is retained inside each bar (white bold) so the
    Panel-C dual-label pattern is preserved.
    """
    style_axes_ggplot(ax, xgrid=False, ygrid=True)

    branches = summary["branch"].tolist()
    pcts = summary["pct"].to_numpy()                # within-subgroup share, %
    vals = summary["mean_abs_shap"].to_numpy()      # absolute mean |SHAP|
    colors = [BRANCH_COLORS[b] for b in branches]

    x = range(len(branches))
    ax.bar(x, pcts, color=colors, edgecolor="white", width=0.62,
            linewidth=0.6, zorder=3)

    # Percent ABOVE each bar (bold black, 1 decimal).
    for i, (v, p) in enumerate(zip(vals, pcts)):
        ax.text(i, p + y_top * 0.025, f"{p:.1f}%",
                ha="center", va="bottom",
                fontsize=BODY_FS, fontweight="bold", color=C_TITLE)
        # Mean |SHAP| INSIDE each bar (white bold), or above if bar is too
        # short to host a centered white label legibly.
        if p > y_top * 0.10:
            ax.text(i, p / 2, f"{v:.2f}",
                    ha="center", va="center",
                    fontsize=BODY_FS - 2, fontweight="bold", color="white")
        else:
            ax.text(i, p + y_top * 0.09, f"{v:.2f}",
                    ha="center", va="bottom",
                    fontsize=BODY_FS - 2, fontweight="bold", color=C_TITLE)

    ax.set_xticks(list(x))
    ax.set_xticklabels([BRANCH_DISPLAY[b] for b in branches],
                       fontsize=BODY_FS, color=C_TITLE)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=C_TITLE)
    ax.set_ylim(0, y_top)

    # Facet-style grey strip header naming the subgroup + its n.
    facet_header(ax, f"{group_label}   (n = {n:,})", height=0.08)


def build_figure():
    df = pd.read_csv(DATA, sep="\t")
    summaries = {g: per_group_summary(df, g) for g in GROUP_ORDER}
    n_by_group = {
        g: int(df.loc[df["group"] == g, "comparator_isoform_id"].nunique())
        for g in GROUP_ORDER
    }
    # Y-axis is now percentage (within-subgroup share). The largest bar
    # across all subgroups is the Control Structural at ~62%; cap at 75%
    # to leave headroom for the percent-above-bar labels without colliding
    # with the subgroup header.
    y_top = 75.0

    # Three Panel-C-style mini-panels side by side. Width per panel
    # matches Panel C's 6×4 sizing; total width 18 in.
    fig, axes = plt.subplots(1, 3, figsize=(NATIVE_W, 5.4), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.16, top=0.80,
                        wspace=0.16)

    for ax, group in zip(axes, GROUP_ORDER):
        draw_panel(ax, summaries[group], GROUP_LABEL[group],
                   n_by_group[group], y_top)

    # Y-axis label only on leftmost panel (matches Panel C).
    axes[0].set_ylabel("Within-subgroup share of |SHAP| (%)",
                       fontsize=BODY_FS - 2, color=C_TITLE)
    # Tick marks every 20% — clean for a 0–75% range.
    for ax in axes:
        ax.set_yticks([0, 20, 40, 60])

    # Panel labels (top-left of each axes)
    for ax, letter in zip(axes, ["A", "B", "C"]):
        panel_letter(ax, letter, x=-0.04, y=1.15)

    # No footnote — caption carries the key numbers (Yul-style).
    return fig, axes


def main():
    fig, axes = build_figure()
    # Facet headers extend past axis edges by design.
    render_and_validate(fig, HERE.parent / "figure_s_branch_shap_by_subclass",
                        native_width_in=NATIVE_W, strict_per_axis=False)


if __name__ == "__main__":
    main()

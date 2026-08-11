"""
Figure 5 Panel G — attention on uORFs, by NMD/PTC subgroup, full deposited cohort.

WHAT CHANGED AND WHY (2026-08-06). This panel used to plot `strict_attn` over the ~1,600-isoform
Subset 2 scope. W348 measured that statistic and it does not survive: `strict_attn` came from a
LEFT JOIN against uorf_features_in_priority_slots.tsv followed by `is.na -> 0`, and NOT ONE of the
832 Control or 767 NMD+/PTC+ isoforms appears in that slot file. Every zero in both groups was the
fill. The two "significant" p-values (~1e-90, ~1e-97) were therefore testing whether the slot file
has coverage, not whether attention differs between groups, and the third comparison's NaN was the
only honest result of the three. The old docstring read the zero mass as biology — "they shouldn't
need uORFs to explain NMD trajectory" — which is exactly the misreading a fill-as-data-point
invites.

The panel now plots `uorf_attention_frac` from the deposited uorf_attention_metrics.tsv, which has
neither defect: 41,776 isoforms, zero blanks, and every zero carries n_uorf_slots == 0. That is
CHECKED in data_export_deposit.py:panel_g() on every build rather than assumed, because the entire
lesson of W348 is that a structural zero and a filled zero are indistinguishable once both are 0.
The comparison groups now carry real coverage too — 2,197 of 6,344 PTC+ and 15,698 of 32,455
Control isoforms are non-zero — so the tests compare two measured distributions.

NO STRIP OVERLAY. The superseded panel drew one over ~1,600 points. At 41,776 — 32,455 of them in
one group — a strip plot is a solid block that reports nothing. The zero inflation is real and
still has to be legible, so it is stated as the non-zero percentage under each group rather than
implied by ink density.

Data:
  ./data/panelG_uorf_attention_v1.tsv           per-isoform: subgroup, uorf_attention_frac
  ./data/panelG_uorf_attention_v1_pairwise.tsv  Mann-Whitney, computed by the producer
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

BODY_FS = 14
LABEL_FS = 12
C_TITLE = "#222222"
C_AXIS = "#555555"

GROUP_ORDER = ["NMD_PTC_neg_lost", "NMD_PTC_neg_retained", "NMD_PTC_pos", "non_NMD"]
GROUP_DISPLAY = {
    "NMD_PTC_neg_lost": "NMD+/PTC−\nstart lost",
    "NMD_PTC_neg_retained": "NMD+/PTC−\nstart retained",
    "NMD_PTC_pos": "NMD+/PTC+",
    "non_NMD": "Control",
}
COLORS = {
    "NMD_PTC_neg_lost": "#fdb863",
    "NMD_PTC_neg_retained": "#fee08b",
    "NMD_PTC_pos": "#ef8a62",
    "non_NMD": "#92c5de",
}

from p_to_stars import p_to_stars as sig_marker  # noqa: E402 — unified 4-tier scheme


def lookup_p(stats, ga, gb):
    row = stats[((stats["group_x"] == ga) & (stats["group_y"] == gb)) |
                ((stats["group_x"] == gb) & (stats["group_y"] == ga))]
    if row.empty:
        raise SystemExit(f"no pairwise row for {ga} vs {gb} -- the producer and the panel disagree "
                         f"about the group set, which is how a bracket gets drawn over the wrong "
                         f"comparison")
    return float(row.iloc[0]["p"])


def build_figure():
    df = pd.read_csv(HERE.parent / "data" / "panelG_uorf_attention_v1.tsv", sep="\t")
    stats = pd.read_csv(HERE.parent / "data" / "panelG_uorf_attention_v1_pairwise.tsv", sep="\t")
    df = df[df["subgroup"].isin(GROUP_ORDER)].copy()
    for g in GROUP_ORDER:
        if not (df["subgroup"] == g).any():
            raise SystemExit(f"group {g} is empty -- an empty factor level renders without error")

    fig, ax = plt.subplots(figsize=(6, 4))

    sns.violinplot(
        data=df, x="subgroup", y="uorf_attention_frac", order=GROUP_ORDER,
        hue="subgroup", palette=COLORS, legend=False,
        cut=0, inner=None, density_norm="width", linewidth=0.9, ax=ax,
    )
    for patch in ax.collections:
        patch.set_edgecolor(C_TITLE)
        patch.set_alpha(0.55)

    n_per = df.groupby("subgroup").size().to_dict()
    medians = df.groupby("subgroup")["uorf_attention_frac"].median().to_dict()
    nonzero = (df[df["uorf_attention_frac"] > 0].groupby("subgroup").size()
               .reindex(GROUP_ORDER).fillna(0).to_dict())

    for i, g in enumerate(GROUP_ORDER):
        ax.hlines(medians[g], i - 0.30, i + 0.30, colors=C_TITLE, linewidth=1.8, zorder=5)

    # Only the two comparisons the claim rests on: each PTC-negative arm against PTC-positive.
    # Control-vs-PTC+ is deliberately not bracketed. Both are median-zero, and at n=32,455 vs
    # 6,344 a Mann-Whitney returns p ~ 1e-120 on a difference in zero fraction that no reader
    # would call an effect -- a star there would be a statement about sample size.
    # spacing must clear the RENDERED HEIGHT of the star glyph on the lower bracket, not just the
    # offset to its baseline: at 0.13 the first '****' intersected the second bracket's segment.
    y_top, spacing, tick, text_off = 1.06, 0.24, 0.018, 0.028
    for (ga, gb, xa, xb, yy) in [
        ("NMD_PTC_neg_lost", "NMD_PTC_pos", 0, 2, y_top),
        ("NMD_PTC_neg_retained", "NMD_PTC_pos", 1, 2, y_top + spacing),
    ]:
        marker = sig_marker(lookup_p(stats, ga, gb))
        ax.plot([xa, xa, xb, xb], [yy - tick, yy, yy, yy - tick], color=C_TITLE, linewidth=0.7)
        ax.text((xa + xb) / 2, yy + text_off, marker, ha="center", va="bottom",
                fontsize=BODY_FS, color=C_TITLE)

    ax.set_xticks(range(len(GROUP_ORDER)))
    # n and the non-zero share go on SEPARATE lines. On one line the last two labels
    # ("(n=6,344; 35% >0)(n=32,455; 48% >0)") ran together with no gap between them -- a text
    # overflow none of the six validators catches, so it has to be designed out rather than
    # checked for.
    ax.set_xticklabels(
        [f"{GROUP_DISPLAY[g]}\nn={n_per[g]:,}\n{100*nonzero[g]/n_per[g]:.0f}% > 0"
         for g in GROUP_ORDER],
        fontsize=LABEL_FS - 2, color=C_TITLE,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Attention fraction on uORFs", fontsize=BODY_FS, color=C_TITLE)
    ax.set_ylim(-0.02, y_top + spacing + 5 * text_off)
    ax.set_yticks([0.0, 0.25, 0.50, 0.75, 1.00])
    ax.tick_params(axis="y", labelsize=LABEL_FS, colors=C_TITLE)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(C_AXIS)
        ax.spines[side].set_linewidth(0.8)

    fig.subplots_adjust(left=0.15, right=0.97, bottom=0.30, top=0.96)
    return fig


def main():
    fig = build_figure()

    # THE VALIDATOR MUST RAISE. validate_figure_layout RETURNS a dict; it does not throw. The bare
    # call this replaces printed "Result: 1 error(s) must be fixed" and then saved the figure
    # anyway, so the gate was decorative -- feedback_validator_surgical_skip's exact failure mode,
    # found here by the collision it let through.
    from validate_figure_layout import validate_figure_layout
    result = validate_figure_layout(fig, fig.axes[0], verbose=True)
    if result["summary"]["n_errors"]:
        raise SystemExit(f"layout validation failed with {result['summary']['n_errors']} error(s); "
                         f"figure NOT saved")

    out_dir = HERE.parent
    fig.savefig(out_dir / "figure5_panelG_uorf_attention.pdf", facecolor="white")
    fig.savefig(out_dir / "figure5_panelG_uorf_attention.png", dpi=300, facecolor="white")
    print("Saved: figure5_panelG_uorf_attention.pdf and .png")


if __name__ == "__main__":
    main()

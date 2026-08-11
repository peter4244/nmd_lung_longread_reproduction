"""
SF36 — CDS length and two 3'UTR-length measures in the reference AUG-traceable
cohort (n = 756 / 63 / 819). Reference-AUG-projected stop anchor.
1×3 panels A/B/C.

Data reused from the combined SF33+SF36 export in the sibling dir; no new
data_export needed.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    panel_letter,
    docx_body_fs,
    place_violin_median_label,
    render_and_validate,
    SUBGROUP_PAL,
    TITLE_C,
    AXIS_C,
)
apply_ggplot_rcparams()

# Native figure width — passed through to docx_body_fs so BODY_FS lands at
# the 10 pt docx-scale target regardless of native render size.
NATIVE_W = 12.0
BODY_FS  = docx_body_fs(NATIVE_W)

# Shared export dir data_export.R writes to (was SF33's local ./data, now stale).
# REPOINTED 2026-08-08. The path named a PRE-RENAME directory -- these figures were
# CDSand3UTR_GENCODEonly/ and TD2BiasEvidence/ before the SF-numbering rename, and the scripts
# kept the old names while their data_export.R wrote to the new ones. So the export succeeded,
# the render failed on a missing file, and the two halves had been disagreeing since the rename.
# SF36 has no exporter of its own; it shares SF33's feature tables.
DATA = HERE.parents[1] / "SF33_CdsAnd3UTR_GENCODE" / "data"

# Canonical SUBGROUP_PAL uses Unicode minus (U+2212) in the "NMD+/PTC−"
# key. Data files store the ASCII hyphen form; we remap to the Unicode
# form after loading so SUBGROUP_PAL keys match without a local override.
GROUP_ORDER = ["NMD+/PTC+", "NMD+/PTC−", "Control"]
GROUP_LABEL_REMAP = {"NMD+/PTC-": "NMD+/PTC−"}

PANELS = [
    {"letter": "A", "slug": "panelD_cds_length_refaug",
     "value_col": "cds_nt",
     "ylabel": "CDS length (ref AUG ORF; nt, log10)"},
    {"letter": "B", "slug": "panelE_utr3_translation_refaug",
     "value_col": "utr3_to_tx_end_nt",
     "ylabel": "3'UTR length from annotated stop (nt, log10)"},
    {"letter": "C", "slug": "panelF_utr3_non_ptc_stop_refaug",
     "value_col": "utr3_via_non_ptc_stop_nt",
     "ylabel": "3'UTR length, non-PTC stop (nt, log10)"},
]


def sig_marker(p):
    if pd.isna(p): return "n.s."
    if p < 1e-10:  return "****"
    if p < 1e-4:   return "***"
    if p < 1e-3:   return "**"
    if p < 0.05:   return "*"
    return "n.s."


def render_panel(ax, long, desc, stats, value_col, ylabel):
    style_axes_ggplot(ax, xgrid=False, ygrid=True)
    plot_df = long.copy()
    plot_df["y"] = np.log10(np.clip(plot_df[value_col], 1, None))
    plot_df = plot_df.dropna(subset=["y"])

    sns.violinplot(data=plot_df, x="group", y="y", order=GROUP_ORDER,
                   hue="group", palette=SUBGROUP_PAL, legend=False,
                   cut=0, inner=None, density_norm="width",
                   linewidth=0.9, ax=ax)
    for patch in ax.collections:
        patch.set_edgecolor(TITLE_C); patch.set_alpha(0.85); patch.set_zorder(3)
    # Median label placed inside the violin via helper — automatic
    # font-fit + fallback to "median = X" below the xtick. See
    # project_nmd_sf33_sf36_companion memory.
    data = [plot_df.loc[plot_df["group"] == g, "y"].to_numpy() for g in GROUP_ORDER]
    median_fallbacks = {}
    for i, g in enumerate(GROUP_ORDER):
        v_log = plot_df.loc[plot_df["group"] == g, "y"]
        if len(v_log) == 0:
            continue
        m_log = float(np.median(v_log))
        m_raw = float(desc.loc[desc["group"] == g, "median"].iloc[0])
        result, _fs = place_violin_median_label(
            ax, i, v_log.to_numpy(), m_raw, m_log, BODY_FS, NATIVE_W)
        if result == "fallback":
            ax.plot([i - 0.20, i + 0.20], [m_log, m_log],
                    color=TITLE_C, linewidth=1.6, zorder=5)
            median_fallbacks[g] = m_raw

    ymax = max(d.max() for d in data)
    y_top = ymax + 0.18
    spacing = 0.50
    # Use Unicode minus (U+2212) to match GROUP_LABEL_REMAP-applied data.
    # ASCII hyphen here silently drops 2 of 3 brackets.
    pairs = [("NMD+/PTC−", "NMD+/PTC+", 0, 1, y_top),
             ("NMD+/PTC−", "Control",  1, 2, y_top + spacing),
             ("NMD+/PTC+", "Control",  0, 2, y_top + 2 * spacing)]
    for gx, gy, ix, iy, yy in pairs:
        row = stats[((stats["group_x"] == gx) & (stats["group_y"] == gy)) |
                    ((stats["group_x"] == gy) & (stats["group_y"] == gx))]
        if row.empty: continue
        p = row.iloc[0]["wilcox_p"]
        ax.plot([ix, ix, iy, iy], [yy - 0.04, yy, yy, yy - 0.04],
                color=TITLE_C, linewidth=0.6)
        ax.text((ix + iy) / 2, yy + 0.02, sig_marker(p),
                ha="center", va="bottom", fontsize=BODY_FS, color=TITLE_C)

    n_per = {g: int(desc.loc[desc["group"] == g, "n"].iloc[0]) for g in GROUP_ORDER}

    ax.set_xticks(range(len(GROUP_ORDER)))
    # Tilted 30° 2-row labels at BODY_FS-4 — matches SF33 companion.
    ax.set_xticklabels(
        [f"{g}\nn = {n_per[g]:,}" for g in GROUP_ORDER],
        fontsize=BODY_FS - 4, color=TITLE_C,
        rotation=30, ha="right", va="top")
    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=BODY_FS - 2, color=TITLE_C)
    ax.tick_params(axis="y", labelsize=BODY_FS - 2, colors=TITLE_C)
    ax.set_ylim(bottom=0, top=y_top + 3 * spacing)

    # Fallback median annotation below xtick for any group whose median
    # label didn't fit inside the violin at readable font.
    for i, g in enumerate(GROUP_ORDER):
        if g in median_fallbacks:
            m_raw = median_fallbacks[g]
            ax.annotate(
                f"median = {int(round(m_raw)):,}",
                xy=(i, 0), xycoords=("data", "axes fraction"),
                xytext=(0, -60), textcoords="offset points",
                fontsize=BODY_FS - 1, color=TITLE_C,
                ha="center", va="top",
            )


def build_figure():
    # 2×2 layout with Panel C bottom-center — matches SF33 companion.
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(NATIVE_W, 11.0))
    gs = GridSpec(2, 4, figure=fig,
                  left=0.08, right=0.96, top=0.92, bottom=0.10,
                  wspace=0.9, hspace=0.55)
    axes = [
        fig.add_subplot(gs[0, 0:2]),   # A — top-left
        fig.add_subplot(gs[0, 2:4]),   # B — top-right
        fig.add_subplot(gs[1, 1:3]),   # C — bottom-center
    ]

    long = pd.read_csv(DATA / "features_1166_subset_long.tsv", sep="\t")
    long["group"] = long["group"].replace(GROUP_LABEL_REMAP)
    for spec, ax in zip(PANELS, axes):
        desc  = pd.read_csv(DATA / f"{spec['slug']}_descriptives.tsv", sep="\t")
        stats = pd.read_csv(DATA / f"{spec['slug']}_pairwise.tsv", sep="\t")
        desc["group"] = desc["group"].replace(GROUP_LABEL_REMAP)
        stats["group_x"] = stats["group_x"].replace(GROUP_LABEL_REMAP)
        stats["group_y"] = stats["group_y"].replace(GROUP_LABEL_REMAP)
        render_panel(ax, long, desc, stats, spec["value_col"], spec["ylabel"])
        panel_letter(ax, spec["letter"], x=-0.14, y=1.08)
    return fig


def main():
    fig = build_figure()
    render_and_validate(fig, HERE.parent / "figure_sf36",
                        native_width_in=NATIVE_W)


if __name__ == "__main__":
    main()

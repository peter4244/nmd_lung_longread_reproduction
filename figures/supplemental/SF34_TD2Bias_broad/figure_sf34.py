"""
SF34 — TransDecoder2 vs reference AUG-traced open reading frames on the broad
reference AUG-traceable cohort (n = 819 gene-matched NMD comparator pairs).
1×3: length KDE / Kozak violin / position histogram.

Data reused from the combined SF34+SF35 export in the sibling dir; no new
data_export needed.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams, style_axes_ggplot, panel_letter, docx_body_fs,
    render_and_validate, TITLE_C, AXIS_C,
)
apply_ggplot_rcparams()
# Native figure width — passed through to docx_body_fs so BODY_FS lands at
# the 10 pt docx-scale target regardless of native render size, and to
# assert_docx_readable so the 9 pt floor is enforced against the right scale.
NATIVE_W = 11.0
BODY_FS  = docx_body_fs(NATIVE_W)
# Read from the shared TD2BiasEvidence/data dir that data_export.R writes to
# (the local ./data was a pre-reorg leftover; see REFERENCE_FLOOR_PLAN.md).
# REPOINTED 2026-08-08. The path named a PRE-RENAME directory -- these figures were
# CDSand3UTR_GENCODEonly/ and TD2BiasEvidence/ before the SF-numbering rename, and the scripts
# kept the old names while their data_export.R wrote to the new ones. So the export succeeded,
# the render failed on a missing file, and the two halves had been disagreeing since the rename.
DATA = HERE.parent / "data"

# Panel A uses the Panel B colour regime: Reference = blue, TD2 = salmon.
C_REF_FILL  = "#4393c3"; C_REF_LINE  = "#2b8cbe"
C_TD2_FILL  = "#ef8a62"; C_TD2_LINE  = "#d6604d"
C_REF_ATG = "#4393c3"; C_TD2_ATG = "#ef8a62"
C_HIST_FILL = "#67a9cf"; C_HIST_LINE = "#2b8cbe"
C_REF_VLINE = "#d6604d"


def render_length_kde(ax, tsv_name):
    style_axes_ggplot(ax, xgrid=True, ygrid=False)
    df = pd.read_csv(DATA / tsv_name, sep="\t")
    x_grid = np.linspace(np.log10(50), np.log10(15000), 400)
    peaks = {}
    # Data-side column value → display label (data uses "ATG" wording; RNA-side
    # display uses "AUG" per project convention).
    for data_key, display, fill, line in [
        ("Reference ATG ORF", "Reference AUG ORF", C_REF_FILL, C_REF_LINE),
        ("TD2-selected CDS",  "TD2-selected CDS",  C_TD2_FILL, C_TD2_LINE),
    ]:
        vals = np.log10(df.loc[df["ORF_source"] == data_key, "length_nt"].to_numpy())
        kde = gaussian_kde(vals, bw_method="scott")
        y = kde(x_grid)
        ax.fill_between(10 ** x_grid, y, alpha=0.55, color=fill, linewidth=0, zorder=3)
        ax.plot(10 ** x_grid, y, color=line, linewidth=1.6, zorder=4)
        peak_i = int(np.argmax(y))
        peaks[display] = (10 ** x_grid[peak_i], y[peak_i], line, len(vals))
    ax.set_xscale("log")
    ax.set_xlim(50, 15000)
    ax.set_xticks([100, 1000, 10000])
    y_max = max(p[1] for p in peaks.values())
    ax.set_ylim(0, y_max * 1.65)
    ax.set_xlabel("ORF length (nt, log scale)", fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Density", fontsize=BODY_FS, color=TITLE_C)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=BODY_FS, colors=TITLE_C)
    # Curve labels — one line each, no (n=...) tag. Sample size lives in
    # Panel B xticks; repeating it here just costs space and creates
    # overlaps at close panel widths.
    label_list = list(peaks.items())
    y_top_pos = [0.94, 0.85]
    for (label, (_, _, line_c, n)), y_pos in zip(label_list, y_top_pos):
        ax.text(0.03, y_pos, label,
                transform=ax.transAxes,
                ha="left", va="top", fontsize=BODY_FS,
                color=line_c, fontweight="bold", zorder=5)


def render_kozak_violin(ax, tsv_name, pvalue_text):
    style_axes_ggplot(ax, xgrid=False, ygrid=True)
    df = pd.read_csv(DATA / tsv_name, sep="\t")
    # Data-side column uses "ATG" wording; display uses "AUG" for the reference.
    sources = ["Reference ATG", "TD2-predicted ATG"]
    display = {"Reference ATG": "Reference\nstart codon", "TD2-predicted ATG": "TD2-predicted\nstart codon"}
    palette = {"Reference ATG": C_REF_ATG, "TD2-predicted ATG": C_TD2_ATG}
    n_per = {s: int((df["ATG_source"] == s).sum()) for s in sources}

    ref = df[df["ATG_source"] == "Reference ATG"][["comparator_isoform_id", "score"]]\
            .rename(columns={"score": "ref_score"})
    td2 = df[df["ATG_source"] == "TD2-predicted ATG"][["comparator_isoform_id", "score"]]\
            .rename(columns={"score": "td2_score"})
    paired = ref.merge(td2, on="comparator_isoform_id")
    n_pairs = len(paired)
    n_ref_higher = int((paired["ref_score"] > paired["td2_score"]).sum())
    pct_ref_higher = 100 * n_ref_higher / n_pairs

    sns.violinplot(data=df, x="ATG_source", y="score", order=sources,
                   hue="ATG_source", palette=palette, legend=False,
                   cut=0, inner=None, density_norm="width",
                   linewidth=0.9, ax=ax)
    for patch in ax.collections:
        patch.set_edgecolor(TITLE_C); patch.set_alpha(0.85); patch.set_zorder(3)
    # Median tick + numeric label per violin (SF29 style, no IQR marks).
    for i, s in enumerate(sources):
        v = df.loc[df["ATG_source"] == s, "score"].to_numpy()
        if v.size == 0:
            continue
        m = float(np.median(v))
        ax.plot([i - 0.22, i + 0.22], [m, m],
                color=TITLE_C, linewidth=1.6, zorder=5)
        ax.text(i, m + 0.15, f"{m:.2f}",
                ha="center", va="bottom", fontsize=BODY_FS,
                fontweight="bold", color=TITLE_C, zorder=6)

    data = [df.loc[df["ATG_source"] == s, "score"].to_numpy() for s in sources]
    ymax = max(np.max(d) for d in data); ymin = min(np.min(d) for d in data)
    y_br = ymax + 0.6
    ax.plot([0, 0, 1, 1], [y_br - 0.15, y_br, y_br, y_br - 0.15],
            color=TITLE_C, linewidth=0.8)
    ax.text(0.5, y_br + 0.1, f"paired Wilcoxon\n{pvalue_text}",
            ha="center", va="bottom", fontsize=BODY_FS, color=TITLE_C)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"{display[s]}\n(n={n_per[s]:,})" for s in sources],
                       fontsize=BODY_FS, color=TITLE_C)
    ax.set_xlabel("")
    ax.set_ylabel("Kozak PWM score (log-odds)", fontsize=BODY_FS, color=TITLE_C)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=TITLE_C)
    ax.set_ylim(bottom=ymin - 0.3, top=y_br + 1.4)

    # Directional-comparison count ("Ref AUG > TD2 Kozak in N/N pairs")
    # moved to the legend — showing it in-panel forces a bbox that fights
    # the violins for space at close panel widths.


def render_position_hist(ax, tsv_name):
    style_axes_ggplot(ax, xgrid=True, ygrid=True)
    df = pd.read_csv(DATA / tsv_name, sep="\t")
    n = len(df)
    n_eq   = int((df["distance_nt"] == 0).sum())
    n_down = int((df["distance_nt"] > 0).sum())
    pct_eq   = 100 * n_eq / n
    pct_down = 100 * n_down / n
    pos = df.loc[df["distance_nt"] > 0, "distance_nt"]
    median_pos = float(np.median(pos)) if len(pos) > 0 else 0.0

    x_min, x_max = -500, 2500
    clipped = int(((df["distance_nt"] < x_min) | (df["distance_nt"] > x_max)).sum())
    vals = df["distance_nt"].clip(x_min, x_max)
    bins = np.linspace(x_min, x_max, 51)
    counts, edges, _ = ax.hist(vals, bins=bins, color=C_HIST_FILL,
                                edgecolor=C_HIST_LINE, linewidth=0.6, alpha=0.90,
                                zorder=3)

    ymax = counts.max() * 1.42
    ax.set_xlim(x_min, x_max); ax.set_ylim(0, ymax)
    # Drop the -500 clip-floor tick: it crowds the origin against the 0 tick
    # (500 nt apart vs 1000 nt for the rest). The red dashed line marks 0, so
    # even [0, 1000, 2000] spacing reads cleaner without losing the anchor.
    ax.set_xticks([0, 1000, 2000])
    # Explicit yticks — auto-ticks otherwise include values above ymax
    # (e.g. "1000" when ymax=820), whose label bbox overlaps the panel
    # letter placed above the axis.
    _ymax_int = int(ymax)
    _step = 200 if _ymax_int > 300 else 20
    ax.set_yticks([y for y in range(0, _ymax_int + 1, _step) if y < ymax])
    ax.axvline(0, color=C_REF_VLINE, linewidth=1.4, linestyle="--", zorder=3)

    # Category counts (TD2 = Ref / TD2 downstream / median offset) and the
    # clipped-outlier note moved to the legend. The dashed red line at x=0
    # visually anchors the reference AUG, and the median dotted line
    # marks the median offset — no in-panel captions needed.

    ax.set_xlabel("TD2 start codon position\nrelative to reference AUG (nt)",
                  fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Number of isoforms",
                  fontsize=BODY_FS, color=TITLE_C)
    ax.tick_params(axis="both", labelsize=BODY_FS, colors=TITLE_C)


def build_figure():
    from matplotlib.gridspec import GridSpec
    # 2×2 with Panel C centered on the bottom row. Implemented as a 2×4
    # grid where each top panel takes 2 columns and Panel C takes the
    # middle 2 columns (cols 1:3), which centers it below A/B while giving
    # it the same width as either.
    fig = plt.figure(figsize=(NATIVE_W, 9.0))
    gs = GridSpec(2, 4, figure=fig,
                  left=0.075, right=0.96, top=0.92, bottom=0.13,
                  wspace=0.9, hspace=0.60)
    axes = [
        fig.add_subplot(gs[0, 0:2]),   # A — top-left
        fig.add_subplot(gs[0, 2:4]),   # B — top-right
        fig.add_subplot(gs[1, 1:3]),   # C — bottom-center
    ]

    render_length_kde(axes[0], "panelA_td2_vs_refaug_length.tsv")
    render_kozak_violin(axes[1], "panelB_kozak.tsv", "p = 1.5×10$^{-36}$")
    render_position_hist(axes[2], "panelC_td2_position.tsv")

    for ax, letter in zip(axes, ["A", "B", "C"]):
        panel_letter(ax, letter, x=-0.15, y=1.10)
    return fig


def main():
    fig = build_figure()
    render_and_validate(fig, HERE.parent / "figure_sf34",
                        native_width_in=NATIVE_W)
    print("Saved: figure_sf34.{pdf,png}")


if __name__ == "__main__":
    main()

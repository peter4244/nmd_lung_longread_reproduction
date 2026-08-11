"""SF31 — NMD-response magnitude vs stop-codon distance to the last EJC.

Scatter of mashr posterior-mean logFC (SMG1i vs DMSO, averaged across AT2,
LAE, FB, MV) against distance from the comparator's stop codon to the last
exon-junction complex. Loess overlay + 95% CI band. Vertical dashed line at
50 nt marks the operational PTC threshold.

Stop anchor: reference AUG-traced (unbiased). The comparator's stop codon
is defined as the first in-frame stop reached by projecting the reference
gene's AUG onto the comparator, and the distance to the last EJC is
computed from the comparator's exon structure. This avoids the TD2-CDS
attenuation that would arise from anchoring on the TD2-called stop for
occult-PTC isoforms.

No rank-correlation statistic is annotated: the phenomenon is a
threshold-like ~50 nt EJC-rule response, not a monotone-scaled quantity.
Summary badge reports the Wilcoxon comparison of logFC above vs at/below
the 50 nt threshold.

Data source: `sf31_ptc_distance_logfc.tsv` (produced by data_export.R from
`ref_atg_analysis.rds$c2` + `structures.rds` + mashr posterior means).
"""

from __future__ import annotations
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
LIB  = HERE.parents[1] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    render_and_validate,
    docx_body_fs,
    NMD_COLOR,
    TITLE_C,
    HEADER_FS,
)

NATIVE_W = 8.0
BODY_FS = docx_body_fs(NATIVE_W)

apply_ggplot_rcparams()

DATA = HERE / "data"


def _tricube(u):
    u = np.abs(u)
    w = np.where(u < 1, (1 - u ** 3) ** 3, 0.0)
    return w


def loess_with_ci(x, y, frac=0.30, n_boot=200, seed=42, n_grid=200):
    """Simple local-linear (LOESS-like) fit with bootstrap 95% CI.

    Iterates over a fixed x-grid; at each grid point, fits a weighted linear
    regression to the nearest `k = ceil(frac * n)` observations, weighted by
    tricube of scaled distance. This reproduces the shape of ggplot's
    `geom_smooth(method="loess")` closely enough for a §4 supplemental
    figure — the point is the dose-response curve, not the exact smoother.
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    k = max(int(np.ceil(frac * n)), 8)

    x_grid = np.linspace(np.quantile(x, 0.005), np.quantile(x, 0.995), n_grid)

    def _fit(xx, yy):
        out = np.empty_like(x_grid)
        for i, xg in enumerate(x_grid):
            d = np.abs(xx - xg)
            # local-window bandwidth = distance to k-th nearest
            if k >= len(xx):
                h = d.max() if d.max() > 0 else 1.0
            else:
                h = np.partition(d, k - 1)[k - 1]
                if h == 0:
                    h = np.max(d[d > 0]) if np.any(d > 0) else 1.0
            w = _tricube(d / h)
            if w.sum() == 0:
                out[i] = np.nan
                continue
            # Weighted linear regression via normal equations
            xw = xx - xg  # center on grid point
            sw   = w.sum()
            swx  = (w * xw).sum()
            swy  = (w * yy).sum()
            swxx = (w * xw * xw).sum()
            swxy = (w * xw * yy).sum()
            det = sw * swxx - swx * swx
            if det == 0:
                out[i] = swy / sw
            else:
                # intercept at xg
                out[i] = (swxx * swy - swx * swxy) / det
        return out

    fit_main = _fit(x, y)
    boots = np.empty((n_boot, n_grid))
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = _fit(x[idx], y[idx])
    lo = np.nanquantile(boots, 0.025, axis=0)
    hi = np.nanquantile(boots, 0.975, axis=0)
    return x_grid, fit_main, lo, hi


def main():
    df = pd.read_csv(DATA / "sf31_ptc_distance_logfc.tsv", sep="\t")
    df = df.dropna(subset=["stop_to_last_ejc_via_ref_nt", "mean_logFC"])

    # x-axis clip 1st/99th centile
    x_lo, x_hi = df["stop_to_last_ejc_via_ref_nt"].quantile([0.01, 0.99]).values
    plot_df = df[(df["stop_to_last_ejc_via_ref_nt"] >= x_lo)
                 & (df["stop_to_last_ejc_via_ref_nt"] <= x_hi)].copy()
    n_plot = len(plot_df)

    x = plot_df["stop_to_last_ejc_via_ref_nt"].to_numpy()
    y = plot_df["mean_logFC"].to_numpy()

    # Threshold comparison: > 50 nt vs ≤ 50 nt (the operational EJC-rule cutoff)
    above = df.loc[df["stop_to_last_ejc_via_ref_nt"] >  50, "mean_logFC"].to_numpy()
    at_or_below = df.loc[df["stop_to_last_ejc_via_ref_nt"] <= 50, "mean_logFC"].to_numpy()
    med_above = float(np.median(above)) if above.size else np.nan
    med_below = float(np.median(at_or_below)) if at_or_below.size else np.nan
    if above.size and at_or_below.size:
        _, p_thr = mannwhitneyu(above, at_or_below, alternative="greater")
    else:
        p_thr = np.nan

    fig, ax = plt.subplots(figsize=(NATIVE_W, 5.0))
    fig.subplots_adjust(left=0.12, right=0.92, top=0.86, bottom=0.14)
    style_axes_ggplot(ax)

    ax.scatter(x, y, s=6, alpha=0.18, color="#3a3a3a",
                edgecolors="none", zorder=3, rasterized=True)

    xg, yfit, lo, hi = loess_with_ci(x, y, frac=0.30)
    ax.fill_between(xg, lo, hi, color=NMD_COLOR, alpha=0.30, linewidth=0, zorder=4)
    ax.plot(xg, yfit, color=NMD_COLOR, linewidth=2.0, zorder=5)

    ax.axvline(50, linestyle="--", color="#c0392b", linewidth=1.0, zorder=2)
    # Place label a visible gap right of the threshold line (~3% of x-axis span
    # in data coords, y in axes coords so it stays pinned near the top).
    ax.text(50 + 0.03 * (x_hi - x_lo), 0.97,
            "PTC threshold (50 nt)",
            transform=ax.get_xaxis_transform(),
            ha="left", va="top", fontsize=BODY_FS - 1, color="#c0392b")

    ax.set_xlabel("Distance from stop codon (reference AUG-traced) to last EJC (nt)",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Mean logFC (mashr posterior mean)",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.set_xlim(x_lo, x_hi)
    ax.set_xticks([-500, 0, 500, 1000, 1500, 2000, 2500])
    # Add top headroom so the top y-tick doesn't touch the canvas edge —
    # post-floor the logFC range reaches ~6 and the auto '6' tick clipped.
    _ymn, _ymx = ax.get_ylim()
    ax.set_ylim(_ymn, _ymx + 0.10 * (_ymx - _ymn))

    # Threshold-honest summary: median logFC above vs at/below 50 nt + Wilcoxon p.
    if np.isnan(p_thr):
        p_str = "n/a"
    elif p_thr < 1e-160:
        p_str = "< 10⁻¹⁶⁰"
    else:
        p_str = f"= {p_thr:.2g}"
    ax.text(0.98, 0.04,
             f"median logFC > 50 nt = {med_above:.2f}  (n = {above.size:,})\n"
             f"median logFC ≤ 50 nt = {med_below:.2f}  (n = {at_or_below.size:,})\n"
             f"Wilcoxon p (one-sided) {p_str}\n"
             f"total n = {len(df):,}  (shown: {n_plot:,})",
             transform=ax.transAxes, ha="right", va="bottom",
             fontsize=BODY_FS - 1, color=TITLE_C,
             bbox=dict(facecolor="white", edgecolor="none", alpha=0.85,
                       boxstyle="square,pad=0.4"))

    # No overall figure title — caption carries the title role (Yul-style).
    render_and_validate(fig, HERE / "figure_sf31_ptc_distance_dose_response",
                        native_width_in=NATIVE_W)
    print(f"  n comparators total: {len(df):,}")
    print(f"  n shown (1st-99th centile): {n_plot:,}")
    print(f"  median logFC >50 / ≤50 nt = {med_above:.3f} / {med_below:.3f}   p = {p_thr:.2e}")


if __name__ == "__main__":
    main()

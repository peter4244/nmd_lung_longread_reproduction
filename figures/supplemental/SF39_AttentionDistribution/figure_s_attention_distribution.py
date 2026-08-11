"""
SF39 — Attention distribution across the model's five candidate ORFs,
NMD susceptible vs non-NMD.

Renders natively from the current 4CT model's per-isoform attention weights
(`uorf_attention_predictions.tsv`, columns attn_0…attn_4). Replaces the
earlier legacy-PNG placeholder that carried stale 10-ORF ranks, "SQANTI CDS"
subtitle, and a pre-4CT-switch n = 15,584 cohort.

Population: held-out test split (chr 1/3/5/7), excluding the "test_paralog"
partition per the training-time paralog dedup convention.

Panels:
  (A) Attention weight per candidate ORF rank, by class (boxplot).
  (B) Shannon entropy of the per-isoform attention vector across the five
      candidates, KDE per class.
"""

# `from __future__` MUST be the first statement after the docstring. It was below the path
# bootstrap, which is a hard SyntaxError -- so this file has never been importable and SF39 could
# not be rendered at all, independently of any data problem. Pre-existing in HEAD; found 2026-08-06
# while repointing the figure at the deposited model.
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while not (_r / "config" / "paths.yml").exists():
    if _r.parent == _r:
        raise RuntimeError("repo root not found")
    _r = _r.parent
_sys.path.insert(0, str(_r / "python"))
from nmd_paths import P

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    render_and_validate,
    docx_body_fs,
    panel_letter,
    TITLE_C,
    AXIS_C,
)
apply_ggplot_rcparams()

NATIVE_W = 11.0
BODY_FS = docx_body_fs(NATIVE_W)
LABEL_FS = BODY_FS - 2

# REPOINTED 2026-08-06. BOTH inputs now come from the deposited 1000nt model, and that is the point
# of the change rather than a tidy-up: these two paths were reading DIFFERENT MODEL VINTAGES.
# uorf_attention_predictions.tsv resolved to the deposit-era file while predictions_atg500_stop500
# .tsv is the superseded 500nt run, so the figure joined attention from one model to labels from
# another and would have rendered without complaint. Same failure class as W348 -- a join that
# silently mixes two populations.
DEPOSIT     = P["ROOT"] / "data_deposit" / "source_data" / "model" if "ROOT" in P else None
_DEP        = (Path(__file__).resolve().parents[3] / "data_deposit" / "source_data" / "model")
DATA_ATTN   = _DEP / "uorf_attention_predictions.tsv"
DATA_LABELS = _DEP / "predictions_atg1000_stop1000_seed42_all.tsv"
for _p in (DATA_ATTN, DATA_LABELS):
    if not _p.exists():
        raise SystemExit(f"missing deposited input: {_p}\n"
                         f"  No fallback to results_4ct -- that is the superseded 500nt vintage, "
                         f"and mixing it with deposit attention is the defect this line fixes.")

# Main-paper Fig 4 / Fig 5G palette
COLOR_NMD  = "#ef8a62"   # coral
COLOR_CTRL = "#92c5de"   # light blue

N_ORFS = 5

# POPULATION: the full cohort, per D77. This is an INTERPRETABILITY quantity -- how attention is
# distributed across candidate ORFs -- and D77 puts every section 5 interpretability quantity over
# all scored isoforms, train and held-out pooled. Pooling is legitimate for interpretation and
# never for a performance number; nothing here is a performance number.
#
# TWO DEFECTS REMOVED, both of which rendered cleanly:
#  1. The old code selected the test split by filtering chr to 1/3/5/7. A `split` column now exists
#     (evaluate.py writes one per row), so the chromosome filter was a workaround for a limitation
#     that no longer applies -- the same one 10_export_stop_codon_freq_sf37.py's header retired on
#     2026-08-03. Worse, chr 1/3/5/7 catches `test_paralog` too, so the figure silently included
#     122 paralog isoforms that test_clean deliberately excludes.
#  2. Its comment promised "the same 2,268 / 7,863 test-set n's" as SF37/38/41/42. Those are the
#     500nt model's counts and no longer describe anything in the deposit.
#
# Labels still come from the predictions file rather than the attention TSV's own label column:
# a handful of isoforms disagreed between the two, and one source for one quantity is the rule.


def load():
    attn = pd.read_csv(DATA_ATTN, sep="\t")
    labels = pd.read_csv(DATA_LABELS, sep="\t")[["isoform_id", "label"]]
    attn = attn.drop(columns=[c for c in ["label"] if c in attn.columns])
    df = attn.merge(labels, on="isoform_id", how="inner")
    if len(df) != len(attn):
        raise SystemExit(f"attention rows {len(attn):,} but joined {len(df):,} -- the two deposited "
                         f"files disagree about the isoform set, which must not be silently dropped")
    return df


def shannon_entropy(row):
    p = np.array([row[f"attn_{k}"] for k in range(N_ORFS)], dtype=float)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    return float(-(p * np.log2(p)).sum())


def render_panel_A(ax, df):
    style_axes_ggplot(ax, xgrid=False, ygrid=True)

    positions_ctrl = np.arange(N_ORFS) * 3 - 0.55
    positions_nmd  = np.arange(N_ORFS) * 3 + 0.55

    data_ctrl = [df.loc[df["label"] == 0, f"attn_{k}"].to_numpy() for k in range(N_ORFS)]
    data_nmd  = [df.loc[df["label"] == 1, f"attn_{k}"].to_numpy() for k in range(N_ORFS)]

    def draw_boxes(data, positions, color):
        ax.boxplot(
            data, positions=positions, widths=0.85,
            patch_artist=True, showfliers=False,
            boxprops=dict(facecolor=color, edgecolor=TITLE_C, linewidth=0.9),
            medianprops=dict(color=TITLE_C, linewidth=1.4),
            whiskerprops=dict(color=TITLE_C, linewidth=0.9),
            capprops=dict(color=TITLE_C, linewidth=0.9),
            zorder=3,
        )

    draw_boxes(data_ctrl, positions_ctrl, COLOR_CTRL)
    draw_boxes(data_nmd,  positions_nmd,  COLOR_NMD)

    ax.set_xticks(np.arange(N_ORFS) * 3)
    ax.set_xticklabels([str(k) for k in range(N_ORFS)],
                        fontsize=BODY_FS, color=TITLE_C)
    ax.set_xlim(-2, (N_ORFS - 1) * 3 + 2)
    ax.set_ylim(-0.03, 1.03)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Candidate ORF rank (0 = main ORF)",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Attention weight",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.tick_params(axis="y", labelsize=BODY_FS, colors=TITLE_C)


def render_panel_B(ax, df):
    style_axes_ggplot(ax, xgrid=False, ygrid=True)

    df["entropy"] = df.apply(shannon_entropy, axis=1)
    e_nmd  = df.loc[df["label"] == 1, "entropy"].to_numpy()
    e_ctrl = df.loc[df["label"] == 0, "entropy"].to_numpy()

    max_ent = float(np.log2(N_ORFS))
    x_grid = np.linspace(0, max_ent, 300)

    for e, color in [(e_nmd, COLOR_NMD), (e_ctrl, COLOR_CTRL)]:
        kde = gaussian_kde(e, bw_method="scott")
        y = kde(x_grid)
        ax.fill_between(x_grid, y, alpha=0.55, color=color,
                        edgecolor=TITLE_C, linewidth=0.9, zorder=3)

    ax.set_xlim(0, max_ent * 1.02)
    ax.set_xlabel("Per-isoform attention entropy (bits)",
                   fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Density", fontsize=BODY_FS, color=TITLE_C)
    ax.tick_params(axis="both", labelsize=BODY_FS, colors=TITLE_C)


def build_figure():
    df = load()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(NATIVE_W, 4.8))
    fig.subplots_adjust(left=0.09, right=0.95, top=0.80, bottom=0.22,
                        wspace=0.32)

    render_panel_A(axA, df)
    render_panel_B(axB, df)

    panel_letter(axA, "A", x=-0.14, y=1.22)
    panel_letter(axB, "B", x=-0.14, y=1.22)

    # Shared legend in the bottom margin — both panels share the same two series; per-panel
    # legends overlapped the boxes/density.
    #
    # THE n's ARE IN THE LEGEND, not only the caption. They used to live in the caption alone, on
    # the reasoning that the legend should carry colors only -- but captions get separated from
    # figures, and this figure's population just changed (test split -> full cohort, D77), which is
    # exactly the kind of change a reader cannot detect from a caption they are not holding.
    n_nmd = int((df["label"] == 1).sum())
    n_ctrl = int((df["label"] == 0).sum())
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_NMD,  edgecolor=TITLE_C,
                      linewidth=0.9, label=f"NMD susceptible (n = {n_nmd:,})"),
        plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_CTRL, edgecolor=TITLE_C,
                      linewidth=0.9, label=f"Control (n = {n_ctrl:,})"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=BODY_FS, bbox_to_anchor=(0.5, 0.01))

    return fig, df


def main():
    fig, df = build_figure()
    render_and_validate(fig, HERE.parent / "figure_s_attention_distribution",
                        native_width_in=NATIVE_W)
    print(f"Saved figure_s_attention_distribution.{{pdf,png}}")
    print(f"  NMD:     n = {int((df['label'] == 1).sum()):,}")
    print(f"  Control: n = {int((df['label'] == 0).sum()):,}")


if __name__ == "__main__":
    main()

"""SF29 — Transcript length distribution by role in the Isopair pair set.

Standalone rebuild of the transcript-length-by-role panel that was previously
bundled as Panel C of PairSetDescriptives. Split per Yul-era paper numbering
so the paper's SF29 reference resolves to one figure.

Three roles within each pair-set triplet:
  - NMD comparator: NMD-classed isoform
  - Reference: highest-expressed non-NMD isoform in DMSO
  - Control comparator: non-NMD, non-reference isoform

Style: matplotlib rendered with ggplot-mimic theme (grey panel + white
gridlines) so the panel visually matches SF1-SF23. See
figures/lib/ggplot_style.py.
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
    render_and_validate,
    docx_body_fs,
    NMD_COLOR,
    CONTROL_COLOR,
    REF_COLOR,
    AXIS_C,
    TITLE_C,
    HEADER_FS,
)

apply_ggplot_rcparams()

NATIVE_W = 6.5
BODY_FS = docx_body_fs(NATIVE_W)

DATA = HERE / "data"

ROLES = ["NMD comparator", "Reference", "Control comparator"]
COLORS = {
    "NMD comparator":     NMD_COLOR,
    "Reference":          REF_COLOR,
    "Control comparator": CONTROL_COLOR,
}
XTICK_LABELS = ["NMD\ncomparator", "Reference", "Control\ncomparator"]


def main():
    txl = pd.read_csv(DATA / "tx_length_by_role_long.tsv", sep="\t")
    summary = pd.read_csv(DATA / "descriptives_summary.tsv", sep="\t").set_index("metric")["value"]
    n_triplets = int(summary.loc["n_pop_bc_triplets"])

    fig, ax = plt.subplots(figsize=(NATIVE_W, 4.6))
    fig.subplots_adjust(left=0.19, right=0.97, top=0.92, bottom=0.19)

    style_axes_ggplot(ax)

    data = [np.log10(txl.loc[txl["role"] == r, "length_nt"].values + 1) for r in ROLES]
    parts = ax.violinplot(
        data, positions=range(len(ROLES)), widths=0.7,
        showmeans=False, showmedians=False, showextrema=False,
    )
    for body, role in zip(parts["bodies"], ROLES):
        body.set_facecolor(COLORS[role])
        body.set_edgecolor(AXIS_C)
        body.set_alpha(0.85)
        body.set_zorder(3)

    for i, r in enumerate(ROLES):
        vals = txl.loc[txl["role"] == r, "length_nt"].values
        med = float(np.median(vals))
        ax.plot([i - 0.20, i + 0.20], [np.log10(med + 1)] * 2,
                color=TITLE_C, lw=1.6, zorder=5)
        ax.text(i, np.log10(med + 1) + 0.045, f"{int(med):,}",
                ha="center", va="bottom", fontsize=BODY_FS,
                fontweight="bold", color=TITLE_C, zorder=6)

    ax.set_xticks(range(len(ROLES)))
    ax.set_xticklabels(XTICK_LABELS, fontsize=BODY_FS, color=TITLE_C)
    ax.set_ylabel("Transcript length (nt, log scale)", fontsize=BODY_FS)
    ax.set_yticks([np.log10(10), np.log10(100), np.log10(1_000),
                   np.log10(10_000), np.log10(100_000)])
    ax.set_yticklabels(["10", "100", "1,000", "10,000", "100,000"])
    ax.set_xlim(-0.6, len(ROLES) - 0.4)

    # No overall figure title — caption carries the title role (Yul-style).
    render_and_validate(fig, HERE / "figure_sf29_transcript_length_by_role",
                        native_width_in=NATIVE_W)
    for r in ROLES:
        med = int(np.median(txl.loc[txl["role"] == r, "length_nt"].values))
        print(f"  median {r}: {med:,}")


if __name__ == "__main__":
    main()

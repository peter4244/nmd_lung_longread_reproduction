"""SF42 — Raw GC content across the AUG and stop-codon windows.

Two-panel supplemental figure showing mean rolling-GC content of transcript
sequences per position, comparing NMD-susceptible vs Control classes, across
the AUG and stop-codon windows the deep-learning model consumes.

The paper's claim (§5, at SF42): GC content differentiates NMD from Control
transcripts in the stop-codon window but not in the AUG window, consistent
with PTCs turning exon sequence (GC-rich) into 3'UTR sequence (GC-poor)
downstream of a premature stop. Post-stop, Control transcripts show the
expected sharp drop in GC content (real 3'UTR); NMD transcripts stay
GC-rich because their "post-stop" sequence is still coding.

Data source: model:results_deposit/gc_content_across_{atg,stop}_window_refaug_only_atg1000_stop1000_seed42.tsv
(from 09_export_gc_content.py --branch {atg,stop}). 50-bp sliding window
with 10-bp step. Shading is ±1 standard error of the mean.

Style: ggplot-mimic theme (grey panel + white gridlines) to match SF1-SF23.
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

# paths.yml is the single resolver (NMD_<KEY> > paths.yml > repo-relative); walk up to it rather
# than counting parents, which is what the old sibling-repo path did and why it was unportable.
_r = HERE
while not (_r / "config" / "paths.yml").exists():
    if _r.parent == _r:
        raise RuntimeError("repo root not found")
    _r = _r.parent
sys.path.insert(0, str(_r / "python"))
from nmd_paths import P  # noqa: E402

from ggplot_style import (
    apply_ggplot_rcparams,
    style_axes_ggplot,
    facet_header,
    render_and_validate,
    docx_body_fs,
    NMD_COLOR,
    CONTROL_COLOR,
    TITLE_C,
    HEADER_FS,
    panel_letter,
)

apply_ggplot_rcparams()

NATIVE_W = 11.5
BODY_FS = docx_body_fs(NATIVE_W)

# REPOINTED 2026-08-06 onto the DEPOSITED 1000nt model. HERE/"data" did not exist and the files
# it named were the superseded 500nt run, so SF42 has not been renderable. The tables were
# produced on Explorer from the DeepSHAP npz (job 8968572 for the full-cohort pair; the
# refaug_only pair this figure uses was already there from 2026-08-05) and fetched into the model
# repo's results_deposit/.
Y_RANGE = (0.10, 0.62)   # shared by both panels; asserted in draw_panel
# REPOINTED 2026-08-10 off the sibling repo: a path built from parents[3].parent resolves only
# where the model repo happens to sit beside this one, and the two GC tables are now deposited.
DATA = P["MODEL_RESULTS"]
TAG = "atg1000_stop1000_seed42"


def _require(path):
    if not path.exists():
        raise SystemExit(f"missing GC table: {path}\n"
                         f"  Produced by model:09_export_gc_content.py from the DeepSHAP npz on\n"
                         f"  Explorer (results_deposit_h5_2026-08-04). No 500nt fallback.")
    return path


def load(branch):
    # ref-AUG-only cohort (see 09d_export_gc_content_refaug_only.py in the
    # model repo): NMD isoforms restricted to is_ref_cds=1 on ORF0; Controls
    # unchanged. Drops the ~23% of NMD test transcripts that fall through to
    # the TD2 fallback, whose "stop" would be a downstream non-PTC codon.
    df = pd.read_csv(
        _require(DATA / f"gc_content_across_{branch}_window_refaug_only_{TAG}.tsv"),
        sep="\t",
    )
    # A GC FRACTION CANNOT EXCEED 1, AND A BLANK PANEL DOES NOT RAISE. Measured 2026-08-06 on
    # gc_content_across_stop_window_refaug_only_atg1000_stop1000_seed42.tsv: 11 of 92 rows carry
    # mean_gc between 1.000 and 1.066, and the whole stop series sits at 0.719-1.066 -- entirely
    # above this figure's 0.35-0.60 axis. The panel therefore rendered COMPLETELY EMPTY and all six
    # validators passed it, because they check layout, not whether anything was drawn.
    #
    # The same producer's ATG output is in range (0.306-0.548) and the standard full-cohort export
    # is in range for both branches (stop 0.131-0.502, job 8968572), so the defect is specific to
    # 09d_export_gc_content_refaug_only.py on the stop branch, not to GC content generally.
    #
    # Refusing rather than clamping, rescaling, or quietly switching to the full-cohort file: that
    # file is a DIFFERENT POPULATION (all isoforms, not the ref-AUG-only cohort this figure is
    # about), and substituting the nearest convenient artifact for the thing itself is the standing
    # rule this project keeps paying to relearn.
    gc = df["mean_gc"].astype(float)
    if not gc.between(0.0, 1.0).all():
        n_bad = int((~gc.between(0.0, 1.0)).sum())
        raise SystemExit(
            f"{branch}: {n_bad} of {len(gc)} mean_gc values are outside [0, 1] "
            f"(range {gc.min():.3f}-{gc.max():.3f}). A GC fraction cannot exceed 1, so this table "
            f"is wrong and the figure must not be drawn from it.\n"
            f"  Fix model:09d_export_gc_content_refaug_only.py for the stop branch and re-export.\n"
            f"  Do NOT substitute gc_content_across_{branch}_window_{TAG}.tsv -- that is the full\n"
            f"  cohort, a different population from the ref-AUG-only one this figure reports.")
    return df.sort_values(["class", "rel_mid"]).reset_index(drop=True)


def draw_panel(ax, df, *, codon_label):
    style_axes_ggplot(ax)
    ax.axvline(0, color="#555555", linewidth=0.8, linestyle=":", zorder=2)
    for cls, color, label in [
        ("NMD",     NMD_COLOR,     "NMD susceptible"),
        ("Control", CONTROL_COLOR, "Control"),
    ]:
        sub = df[df["class"] == cls]
        x   = sub["rel_mid"].values
        y   = sub["mean_gc"].values
        se  = sub["se_gc"].values
        ax.fill_between(x, y - se, y + se, color=color, alpha=0.15, zorder=3)
        # n IN THE LEGEND, not only the caption. The NMD arm is the ref-AUG-only subset --
        # 5,905 of 9,321 NMD isoforms -- while Control is unfiltered, so the two curves have
        # very different denominators and a reader holding the figure alone cannot tell.
        n_cls = int(sub["n_transcripts"].iloc[0]) if "n_transcripts" in sub else len(sub)
        ax.plot(x, y, color=color, linewidth=1.8, zorder=4,
                label=f"{label} (n = {n_cls:,})")
    ax.set_xlabel(f"Position relative to {codon_label} (nt)", fontsize=BODY_FS)
    # TWO LINES, measured rather than styled. On one line this rotated label spans figure y
    # 0.1378-0.9022 -- 0.76 of the canvas against an axes height of 0.60 -- so it overflowed
    # its own axis at both ends and left 0.0018 to the panel letter above it. That is a
    # WARNING on this laptop and a hard validator ERROR in the container, whose fonts are
    # wider (job 9049007). Wrapping halves the height; nothing about the wording changes.
    ax.set_ylabel("Mean GC content\n(50-nt window)", fontsize=BODY_FS)
    ax.set_xlim(df["rel_mid"].min(), df["rel_mid"].max())
    # BOTH OF THESE WERE SIZED FOR THE 500nt MODEL and neither would have errored.
    # Ticks stopped at +/-200 on a window that now runs to +/-475, labelling under half of it;
    # and the 0.35 floor clipped the tails, which fall to 0.112. A clipped curve and a
    # mislabelled axis both render perfectly cleanly.
    ax.set_xticks([-400, -200, 0, 200, 400])
    #
    # ONE SHARED Y RANGE ACROSS BOTH PANELS, not per-panel autoscale. The two panels are read
    # against each other -- the point is that post-stop GC stays high in NMD while the start
    # window peaks and falls -- and two different y scales would make that comparison a
    # visual artifact. Asserted rather than assumed so a future export cannot quietly fall
    # outside it and be clipped, which is the defect this figure just had.
    lo, hi = float(df["mean_gc"].min()), float(df["mean_gc"].max())
    if lo < Y_RANGE[0] or hi > Y_RANGE[1]:
        raise SystemExit(f"data spans {lo:.3f}-{hi:.3f}, outside the shared axis {Y_RANGE}. "
                         f"Widen Y_RANGE deliberately rather than letting the curve clip.")
    ax.set_ylim(*Y_RANGE)
    return ax


def main():
    atg  = load("atg")
    stop = load("stop")

    fig, axes = plt.subplots(1, 2, figsize=(NATIVE_W, 4.8), sharey=True)
    # left=0.09 gave the y-axis label 0.0020 of clearance from panel letter 'A' on this laptop --
    # a WARNING here and a hard validator ERROR inside the container (job 9049007), because the
    # image resolves a different font and the wider glyphs turn near-touching into overlapping.
    # A layout whose validity depends on which fonts are installed is not a layout; widened so the
    # gap is real rather than marginal.
    fig.subplots_adjust(left=0.115, right=0.96, top=0.82, bottom=0.22, wspace=0.28)

    draw_panel(axes[0], atg,  codon_label="start codon")
    draw_panel(axes[1], stop, codon_label="stop codon")
    axes[1].set_ylabel("")

    facet_header(axes[0], "Start codon window", height=0.08)
    facet_header(axes[1], "Stop codon window",  height=0.08)

    panel_letter(axes[0], "A", x=-0.14, y=1.14)
    panel_letter(axes[1], "B", x=-0.14, y=1.14)

    # Shared legend in the bottom margin — the GC curves fill both panels,
    # so an in-panel legend can't clear the data.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False,
               fontsize=BODY_FS, bbox_to_anchor=(0.5, 0.01))

    # Facet headers extend past axis edges by design.
    render_and_validate(fig, HERE / "figure_sf42_gc_content_stop_window",
                        native_width_in=NATIVE_W, strict_per_axis=False)
    print(f"  ATG rows: {len(atg)}   STOP rows: {len(stop)}")


if __name__ == "__main__":
    main()

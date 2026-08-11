"""
Figure 3 Panel D — Stop-codon distance from last EJC (NMD vs Control).

Density of distance from the comparator's stop codon to its last exon-exon
junction. NMD comparators concentrate to the RIGHT of the 50-nt PTC threshold
(stops downstream of the last EJC); Control comparators concentrate LEFT
(normal stops in the last exon).

Stop source: ref-AUG-traced `comp_stop_tx_pos` from ref_atg_analysis.rds.
Scope is pop_traceable ∩ ENST-only references (categories effectively_ptc,
no_downstream_ejc, truncated_no_ejc), matching the §4 framework
(figures/multipanel/figure4_ptcneg_and_model/RATIONALE.md §2). No TD2 fallback
is used — pairs without a ref-AUG-traceable stop fall outside this panel's
scope by construction.

Data:
  ./data/panelD_stop_codon_distance.tsv
    Columns: comparator_isoform_id, comparison ("NMD" | "Control"),
             distance (int, nt), category (ref-AUG class)
    Source: pop_traceable_ENST (NMD c2 = 1,659 / Control c4 = 1,286)
            + structures.rds (last EJC tx position) + ref_atg_analysis.rds
            (ref-AUG stop). Exported by ./data_export.R.

Sample size: 1,659 NMD / 1,286 Control (ENST-only pop_traceable).

X-axis clipped to [-1000, 1500] nt; clipped counts printed to stdout for
methodology file disclosure.

Style: HEADER_FS=18 bold, BODY_FS=14 regular. NMD coral, Control light blue.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

HERE = Path(__file__).resolve()
LIB = HERE.parents[2] / "lib"
sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

HEADER_FS = 18
BODY_FS = 14

C_NMD_FILL = "#ef8a62"
C_NMD_LINE = "#d6604d"
C_CTRL_FILL = "#67a9cf"
C_CTRL_LINE = "#2b8cbe"
C_TITLE = "#222222"
C_AXIS = "#555555"
C_THRESHOLD = "#cc0000"

X_MIN, X_MAX = -1000, 1500


def load_data():
    return pd.read_csv(HERE.parent / "data" / "panelD_stop_codon_distance.tsv", sep="\t")


def build_figure(df):
    # Composite cell target: 6.0 × 4.0 in (1.5:1).
    fig, ax = plt.subplots(figsize=(6, 4))

    x_grid = np.linspace(X_MIN, X_MAX, 800)
    n_clipped = {}
    for label, fill, line in [
        ("NMD", C_NMD_FILL, C_NMD_LINE),
        ("Control", C_CTRL_FILL, C_CTRL_LINE),
    ]:
        vals = df.loc[df["comparison"] == label, "distance"].dropna().to_numpy()
        kde = gaussian_kde(vals, bw_method="scott")
        y = kde(x_grid)
        ax.fill_between(x_grid, y, alpha=0.45, color=fill, linewidth=0,
                        label=f"{label} (n={len(vals):,})")
        ax.plot(x_grid, y, color=line, linewidth=1.6)
        n_clipped[label] = int(((vals < X_MIN) | (vals > X_MAX)).sum())

    # PTC threshold — short label sits just RIGHT of the line, top of plot.
    # Multi-line "50-nt PTC threshold" was overlapping the upper-left legend
    # at 6×4; compact "PTC ≥50 nt" + legend moved to upper-right resolves it.
    ax.axvline(50, color=C_THRESHOLD, linestyle="--", linewidth=1.2, alpha=0.7)
    ymax = ax.get_ylim()[1]
    ax.text(50 + 30, ymax * 0.96, "PTC ≥50 nt",
            ha="left", va="top", fontsize=BODY_FS, color=C_THRESHOLD)

    # No in-panel title — caption / composite letter label carries it.
    ax.set_xlabel("Distance from last EJC (nt)", fontsize=BODY_FS, color=C_TITLE)
    # No y-axis label or tick labels — absolute KDE density values aren't
    # meaningful for shape comparison; remove for visual cleanliness.
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=BODY_FS, colors=C_TITLE)

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(bottom=0)

    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(C_AXIS)
    ax.spines["bottom"].set_linewidth(0.8)

    # Legend at upper-RIGHT — the NMD-coral tail at x≳400 sits well below
    # y=0.001, so the upper-right quadrant is mostly empty. Anchoring
    # there prevents the legend text from running into the Control-blue
    # peak (which sits at x≈-100, dominating the upper-LEFT region).
    legend = ax.legend(fontsize=BODY_FS, loc="upper right",
                       bbox_to_anchor=(1.0, 0.96), frameon=False)
    for text in legend.get_texts():
        text.set_color(C_TITLE)

    print(f"  Clipped at x∈[{X_MIN}, {X_MAX}]: NMD={n_clipped['NMD']}, Control={n_clipped['Control']}")

    fig.subplots_adjust(left=0.04, right=0.97, bottom=0.16, top=0.96)
    return fig


def main():
    df = load_data()
    fig = build_figure(df)
    from validate_figure_layout import validate_figure_layout
    validate_figure_layout(fig, fig.axes[0], verbose=True)
    out_dir = HERE.parent
    fig.savefig(out_dir / "figure3_panelD_stop_codon_distance.pdf",
                facecolor="white")
    fig.savefig(out_dir / "figure3_panelD_stop_codon_distance.png",
                dpi=300, facecolor="white")
    print("Saved: figure3_panelD_stop_codon_distance.pdf and .png")


if __name__ == "__main__":
    main()

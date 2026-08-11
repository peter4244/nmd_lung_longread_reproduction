"""Matplotlib style helpers that mimic ggplot2's default theme.

Used by all §4 and §5 Supplemental Figures so they visually match the
R/ggplot-generated SF1–SF23 while keeping matplotlib as the tooling.

The three functions to call from a panel script:

  - `apply_ggplot_rcparams()` at import time to set font + fonttype.
  - `style_axes_ggplot(ax)` inside the panel render, right after building
    the axes. Applies the grey panel background and white gridlines.
  - `facet_header(ax, text)` to place a plain bold panel header
    above an axes (matches the "AT2 / LAE / FB / MV" strip you see in
    SF1, SF3, SF19, SF22, etc.).

Paper-canonical colours are re-exported here so panels can `from
ggplot_style import CT_COLORS, NMD_COLOR, CONTROL_COLOR, ...` rather than
each rebuilding the palette.
"""

from __future__ import annotations
import matplotlib.pyplot as plt

# ── palette (matches project_nmd_figure_palette memory + Panel C py) ──
NMD_COLOR     = "#ef8a62"   # NMD comparator; also used for NMD/Susceptible class
CONTROL_COLOR = "#67a9cf"   # Control comparator; also used for non-NMD class
REF_COLOR     = "#7f7f7f"   # Reference isoform (structurally distinct)
BAR_COLOR     = "#4292c6"   # blue for descriptive histograms

# Cell-type palette used by the R/ggplot SF20 + SF21 pages.
CT_COLORS = {
    "AT2": "#F19797",   # peachy pink
    "LAE": "#A0C954",   # olive/green
    "FB":  "#3ABAB4",   # teal
    "MV":  "#B08CD3",   # lavender/purple
}
CT_ORDER = ["AT2", "LAE", "FB", "MV"]

# PTC-subgroup palette (matches project_nmd_figure_palette memory).
# Keys use the Unicode minus (U+2212) in "NMD+/PTC−" to match Yul-style
# canonical labelling; figure scripts that load ASCII-hyphen data should
# remap "-" → "−" before joining against this palette.
SUBGROUP_PAL = {
    "NMD+/PTC+": "#ef8a62",   # peach — matches binary NMD_COLOR
    "NMD+/PTC−": "#d95f02",   # orange — ColorBrewer Dark2 (ref-ATG-retained)
    "Control":   "#67a9cf",   # blue  — matches binary CONTROL_COLOR
}

# Structural colours for the theme itself.
PANEL_BG   = "#EBEBEB"   # ggplot theme_grey panel background
GRID_COLOR = "#FFFFFF"   # ggplot theme_grey major grid
STRIP_BG   = "#D9D9D9"   # facet strip background
TITLE_C    = "#222222"   # near-black for chart text
AXIS_C     = "#555555"   # spine + tick lines

# Font sizes: matched to Fig 3/4/5 (HEADER_FS=18, BODY_FS=14) so SFs read at
# the same scale as the main paper figures. Two sizes only per
# figures_principles.md.
HEADER_FS = 18
BODY_FS   = 14
PANEL_LETTER_FS = 20  # bold A/B/C labels top-left of each panel

# ─── Docx-scale readability accounting ─────────────────────────────────
# The supplemental-figures Word document scales every embedded figure to a
# fixed 6.5 " content width, so effective font size in the reader's docx is
#   effective_pt = native_pt × (6.5 / native_width_in)
# A figure rendered at 16 " with native BODY_FS = 14 pt reads at only 5.7 pt
# in the docx — unreadable. Calibrated empirically against Pete's SF25 /
# SF29 / SF30 review (2026-07-09): target = 10 pt effective, floor = 9 pt,
# ceiling ~14 pt before text starts to dominate the panel.
PAPER_CONTENT_W_IN = 6.5
PAPER_TARGET_PT    = 10.0
# Floor lowered from 9 → 7 on 2026-07-10 to accommodate Yul-style
# xtick labels on 3-group violins (SF32/SF35), which need ~7 pt
# effective to fit tilted 2-row labels without overlap. Body text
# still targets 10 pt via docx_body_fs(NATIVE_W); the floor is only
# for tick labels and secondary annotations.
PAPER_FLOOR_PT     = 7.0


def docx_body_fs(native_width_in, *, target_pt=PAPER_TARGET_PT):
    """Native BODY_FS that renders at `target_pt` in the 6.5 " docx.

    Every panel that will be embedded in the supplemental Word doc SHOULD
    set its own BODY_FS via this helper at the top of the module, e.g.

        NATIVE_W = 11.0
        BODY_FS  = docx_body_fs(NATIVE_W)

    rather than accepting the shared 14-pt constant blindly — the shared
    constant is calibrated for a 6.5 "-wide render, not a 16 "-wide one.
    """
    return int(round(target_pt * native_width_in / PAPER_CONTENT_W_IN))


def docx_effective_pt(native_pt, native_width_in):
    """Effective docx-rendered point size of `native_pt` at scale."""
    return native_pt * PAPER_CONTENT_W_IN / native_width_in


def render_and_validate(fig, out_path, *, native_width_in, dpi=300,
                        run_data_free=True, run_tick_disjoint=True,
                        strict_per_axis=True,
                        run_legend_clear=True,
                        legend_clearance_px=4,
                        legend_overlap_tolerance_px=0):
    """Run the full seven-validator gate + save the figure in one call.

    This is the SSOT for the figure-shipping contract. Every figure's
    main() should reduce to::

        fig = build_figure()
        render_and_validate(fig, HERE.parent / "figure_xx",
                            native_width_in=NATIVE_W)

    instead of copy-pasting 30 lines of individual validator calls.
    Adding a new check in the future is a one-line edit here — every
    figure that uses this helper picks it up automatically.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    out_path : str | Path
        Base path WITHOUT extension. Two files are written:
        ``{out_path}.pdf`` and ``{out_path}.png``.
    native_width_in : float
        Passed to `assert_docx_readable` so the docx-scale floor is
        checked against the correct render width.
    dpi : int
        PNG DPI. PDF is always vector.
    run_data_free : bool
        Skip `assert_data_free` if the figure's non-data content is
        intentionally on top of data (rare — only opt out with a
        comment explaining why).
    run_tick_disjoint : bool
        Skip `assert_tick_labels_disjoint` if a figure has an
        acceptable within-axis tick overlap that has been reviewed
        (rare).
    run_legend_clear : bool
        Skip `assert_legend_clear` (the zero-tolerance legend-vs-data
        clearance check). Only set False for figures where a
        specific, reviewed legend overlap has been accepted.
    legend_clearance_px : int
        Keep-out zone around every data pixel for the legend check
        (default 4 ≈ 1pt at DPI=300). Increase for stricter breathing
        room, decrease if a tight-corner legend is acceptable.
    legend_overlap_tolerance_px : int
        Pixel count of intersection that still passes. Default 0 —
        strict. Bump above 0 only after visual review of a specific
        acceptable overlap, and leave a comment explaining why.

    Raises
    ------
    AssertionError
        If any validator fails. Message names the failing check and
        the offending element.
    """
    # Local imports to avoid a top-level ggplot_style ↔ validate cycle.
    from validate_figure_layout import (
        assert_style_symmetric,
        assert_docx_readable,
        assert_data_free,
        assert_tick_labels_disjoint,
        assert_legend_clear,
        validate_figure_layout,
        validate_multipanel_layout,
    )

    assert_text_within_canvas(fig)
    assert_style_symmetric(fig)
    assert_docx_readable(fig, native_width_in=native_width_in)
    if run_tick_disjoint:
        assert_tick_labels_disjoint(fig)

    _mp = validate_multipanel_layout(fig)
    if _mp["summary"]["n_errors"] > 0:
        msgs = "\n".join(f"  [{e.check}] {e.message}" for e in _mp["errors"])
        raise AssertionError(
            f"validate_multipanel_layout: "
            f"{_mp['summary']['n_errors']} errors:\n{msgs}"
        )

    for ax in fig.axes:
        if strict_per_axis:
            r = validate_figure_layout(fig, ax, verbose=False)
            if r["summary"]["n_errors"] > 0:
                msgs = "\n".join(
                    f"  [{e.check}] {e.message}" for e in r["errors"]
                )
                raise AssertionError(
                    f"validate_figure_layout: "
                    f"{r['summary']['n_errors']} errors on {ax}:\n{msgs}"
                )
        if run_data_free:
            assert_data_free(fig, ax)
        if run_legend_clear:
            assert_legend_clear(
                fig, ax,
                clearance_px=legend_clearance_px,
                overlap_tolerance_px=legend_overlap_tolerance_px,
            )

    from pathlib import Path
    p = Path(out_path)
    fig.savefig(f"{p}.pdf", facecolor="white")
    fig.savefig(f"{p}.png", dpi=dpi, facecolor="white")
    return str(p)


def place_violin_median_label(ax, i, v_log, m_raw, m_log, body_fs, native_w,
                              *, min_effective_pt=PAPER_FLOOR_PT,
                              peak_width=0.8, width_margin=0.85):
    """Place a numeric median label inside a violin at position `i`.

    The label serves as the median marker (no separate tick line needed).
    Font is shrunk progressively — from `body_fs` down to the effective-
    readability floor at `native_w` — until the label bbox fits within
    the violin's width.

    **Width check is at the label's TOP edge**, not at the median.
    Violins ALWAYS taper above the median, so the top of the label bbox
    is where the width constraint bites. Same at the bottom edge in case
    the KDE peak sits above the median. The narrower of the two applies.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Panel axes.
    i : int
        Group index (x-position of the violin).
    v_log : array-like
        Y-values (e.g., log10 of the raw values) for the group's data —
        the KDE is fit to these.
    m_raw : float
        Raw (un-logged) median to display.
    m_log : float
        Y-position of the median in the plotting scale (typically
        median of `v_log`).
    body_fs : int
        Native font size to try first.
    native_w : float
        Figure native width in inches (for docx-scale floor calc).
    min_effective_pt : float
        Readability floor at docx scale (default: PAPER_FLOOR_PT).
    peak_width : float
        Seaborn violinplot peak width (default 0.8 matches seaborn's
        default `width=0.8`).
    width_margin : float
        Safety factor — label width must be < this fraction of violin
        width at the narrower of top/bottom of the label bbox.

    Returns
    -------
    ('inside', native_pt) if the label fits inside the violin.
    ('fallback', None)    if no font from body_fs down to the readability
                           floor fits — caller places the label elsewhere
                           (below xtick as 'median = X' first choice,
                           above violin second choice).
    """
    import numpy as np
    from scipy.stats import gaussian_kde

    label = f"{int(round(m_raw)):,}"
    v_arr = np.asarray(v_log, dtype=float)
    if len(v_arr) < 3:
        return "fallback", None

    kde = gaussian_kde(v_arr)
    y_min, y_max = float(v_arr.min()), float(v_arr.max())
    y_grid = np.linspace(y_min, y_max, 200)
    max_d = float(kde.evaluate(y_grid).max())

    def _violin_w_at(y):
        d = float(kde.evaluate([y])[0])
        return peak_width * (d / max_d)

    floor_native = int(np.ceil(min_effective_pt * native_w / PAPER_CONTENT_W_IN))
    fig = ax.figure

    for fs in range(int(round(body_fs)), floor_native - 1, -1):
        tmp = ax.text(0, 0, label, fontsize=fs, fontweight="bold",
                      ha="center", va="center", alpha=0)
        fig.canvas.draw()
        bb = tmp.get_window_extent(renderer=fig.canvas.get_renderer())
        tmp.remove()
        inv = ax.transData.inverted()
        (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        label_w = abs(x1 - x0)
        label_h_half = abs(y1 - y0) / 2.0

        # Width at TOP of label (violin narrower above median — Pete 2026-07-10).
        y_top = min(m_log + label_h_half, y_max)
        w_top = _violin_w_at(y_top)
        # Width at BOTTOM of label (in case KDE peak sits above median).
        y_bot = max(m_log - label_h_half, y_min)
        w_bot = _violin_w_at(y_bot)
        min_w = min(w_top, w_bot) * width_margin

        if label_w < min_w:
            ax.text(i, m_log, label, ha="center", va="center",
                    fontsize=fs, fontweight="bold", color=TITLE_C, zorder=6)
            return "inside", fs

    return "fallback", None


# ─── NIH grant-figure readability accounting ───────────────────────────
# NIH R01 single-column page content width is 7.5 " after the 0.5 "
# margins. Figures placed edge-to-edge span 8.5 " (rare). Calibrated
# empirically 2026-07-09 against the submitted 2026 NMD grant figures:
#   figure_aim1.py     (8.5 ", BODY_FS=10)   → 10 pt effective — reads clean
#   figure_overview.py (15.5 ", BODY_FS=21)  → 10 pt effective — reads clean
#   figure_saec_lineage.py (15 ", BODY_FS=17) → 8.5 pt — clean but tight
# Target 11 pt body / 13 pt header at page scale. Floor 8 pt (looser than
# publications' 9 pt because grant figures are schematic — text anchors
# to labeled regions, not competing with dense data plotting).
GRANT_PAGE_W_IN     = 7.5      # NIH single-column content width
GRANT_TARGET_PT     = 11.0     # body text target
GRANT_HEADER_PT     = 13.0     # header text target
GRANT_FLOOR_PT      = 8.0      # readability floor


def grant_body_fs(native_width_in, *, target_pt=GRANT_TARGET_PT,
                   page_w_in=GRANT_PAGE_W_IN):
    """Native BODY_FS that renders at `target_pt` on the NIH grant page.

    Use inside grant figures:

        NATIVE_W = 15.5
        BODY_FS  = grant_body_fs(NATIVE_W)       # → 23 at 15.5"

    If the figure spans the full 8.5" page edge-to-edge (rare), pass
    page_w_in=8.5.
    """
    return int(round(target_pt * native_width_in / page_w_in))


def grant_effective_pt(native_pt, native_width_in,
                        *, page_w_in=GRANT_PAGE_W_IN):
    """Effective NIH-page point size of `native_pt` at grant page scale."""
    return native_pt * page_w_in / native_width_in


def apply_ggplot_rcparams():
    """Set matplotlib rcParams to fonts/embedding compatible with SF1-SF23.

    Call this once at import time in a panel script (before creating any
    Figure). Sets Arial-family fonts and PDF/PS fonttype=42 so text is
    embedded as vector glyphs (searchable, non-rasterised).
    """
    plt.rcParams.update({
        "font.family": "sans-serif",
        # THE LAST ENTRY WAS THE ONLY REACHABLE ONE IN THE CONTAINER, and it is the one face
        # here that is NOT metric-compatible with Arial. The image has no Arial, no Helvetica
        # Neue and no Helvetica, so every container render fell through to DejaVu Sans, which
        # is wider -- that is why SF42 passed the layout validator on the laptop and failed it
        # in clean-room job 9049007. Liberation Sans and Nimbus Sans are Arial/Helvetica
        # metric clones, so inserting them ahead of DejaVu makes the fallback preserve advance
        # widths instead of changing them. Nimbus Sans is already in the image; Liberation is
        # closer in glyph shape and is added to it. Nothing changes on a machine with Arial.
        "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica",
                            "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.titlesize": HEADER_FS,
        "axes.labelsize": BODY_FS,
        "xtick.labelsize": BODY_FS,
        "ytick.labelsize": BODY_FS,
        "legend.fontsize": BODY_FS,
    })


def style_axes_ggplot(ax, *, xgrid=True, ygrid=True, minor=False):
    """Apply the ggplot theme_grey look to `ax`.

    - Grey panel background (`PANEL_BG`).
    - White major gridlines behind data.
    - No visible spines (ggplot's default hides them; the panel bg
      provides the visual frame).
    - Tick marks removed, only labels kept.
    - Text drawn in `TITLE_C`.

    Call this immediately after creating the axes and before drawing.
    """
    ax.set_facecolor(PANEL_BG)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    # Grid behind data. axisbelow=True is critical so bars don't cover the grid.
    ax.set_axisbelow(True)
    ax.grid(
        axis="both" if (xgrid and ygrid) else ("x" if xgrid else "y"),
        which="major",
        color=GRID_COLOR,
        linewidth=0.8,
        zorder=0,
    )
    if minor:
        ax.minorticks_on()
        ax.grid(which="minor", color=GRID_COLOR, linewidth=0.4, zorder=0)
    # Suppress tick marks; keep labels in TITLE_C.
    ax.tick_params(
        axis="both",
        which="both",
        length=0,
        colors=TITLE_C,
        labelsize=BODY_FS,
    )
    # Axis labels — set colour only; leave the actual label text to caller.
    ax.xaxis.label.set_color(TITLE_C)
    ax.yaxis.label.set_color(TITLE_C)
    return ax


def facet_header(ax, text, *, height=0.06):
    """Place a panel header above an axes as plain bold text.

    Renders `text` in bold above the axes with no background box (the
    old grey strip read as clutter). `height` sets the header's vertical
    offset above the axes in axes-relative coordinates.
    """
    ax.text(
        0.5, 1.0 + height / 2,
        text,
        transform=ax.transAxes,
        ha="center", va="center",
        fontsize=HEADER_FS - 1,
        fontweight="bold",
        color=TITLE_C,
    )


def panel_letter(ax, letter, *, x=-0.02, y=1.02):
    """Place a bold A/B/C letter in the top-left of a panel axes.

    Matches the sequential letter-label convention used across this manuscript's composite
    figures: bold, top-left of each panel, above the axes.

    NB: this text sits ABOVE the axes (y > 1). If you also have a
    facet_header on the axes, increase y further (typically 1.10-1.14)
    AND make sure fig.subplots_adjust(top=...) leaves ≥ 0.15 headroom
    above the axes so the letter's ascender doesn't clip at the canvas
    top. Run assert_text_within_canvas(fig) before savefig to catch
    clipping automatically.
    """
    ax.text(
        x, y, letter,
        transform=ax.transAxes,
        fontsize=PANEL_LETTER_FS,
        fontweight="bold",
        color=TITLE_C,
        ha="left", va="bottom",
    )


def assert_text_within_canvas(fig, *, tolerance_px=1.0):
    """Raise AssertionError if any Text artist extends outside the figure canvas.

    This is the systematic guard for the "panel letter clipped at top"
    class of bug — panel letters and facet-header text sit at y > 1 in
    axes coords, and if fig.subplots_adjust(top=…) doesn't leave enough
    headroom their ascenders get chopped by the canvas boundary.

    Call this after all text is placed and BEFORE savefig. Raises with a
    detailed message pointing at the offending Text object so the user
    knows what to move.
    """
    # Force a draw so text bounding boxes are populated.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_bbox = fig.bbox  # in display pixels
    offenders = []
    for ax in fig.axes:
        for txt in list(ax.texts) + [ax.xaxis.label, ax.yaxis.label, ax.title]:
            s = txt.get_text() if hasattr(txt, "get_text") else ""
            if not s or not s.strip():
                continue
            try:
                bb = txt.get_window_extent(renderer=renderer)
            except Exception:
                continue
            if (bb.x0 < fig_bbox.x0 - tolerance_px
                    or bb.y0 < fig_bbox.y0 - tolerance_px
                    or bb.x1 > fig_bbox.x1 + tolerance_px
                    or bb.y1 > fig_bbox.y1 + tolerance_px):
                offenders.append(
                    f"  '{s}' extends outside canvas "
                    f"(text bbox x0={bb.x0:.0f} y0={bb.y0:.0f} x1={bb.x1:.0f} y1={bb.y1:.0f}; "
                    f"canvas x=[0..{fig_bbox.x1:.0f}] y=[0..{fig_bbox.y1:.0f}])"
                )
    if offenders:
        raise AssertionError(
            "Text elements would be clipped by the figure canvas — adjust "
            "fig.subplots_adjust(top=…) / panel_letter(y=…) / margins:\n"
            + "\n".join(offenders)
        )

#!/usr/bin/env python3
"""
validate_figure_layout.py

Validates spatial layout of matplotlib annotation-based figures (architecture
diagrams, flowcharts, schematics).  Checks for text-arrow collisions, text
overlap, crowding, centering, alignment, and padding.

Python port of validate_figure_layout.R — same 7 checks, adapted for
matplotlib's artist model.

Usage:
    from validate_figure_layout import validate_figure_layout
    fig, ax = plt.subplots(...)
    # ... build figure ...
    result = validate_figure_layout(fig, ax)

All thresholds are in data-coordinate units and can be tuned via arguments.
"""

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.text as mtext
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass
class TextElement:
    idx: int
    x: float
    y: float
    label: str
    fontsize: float
    ha: str
    va: str
    alpha: float
    bb_xmin: float = 0.0
    bb_xmax: float = 0.0
    bb_ymin: float = 0.0
    bb_ymax: float = 0.0
    # True when the text was placed with transform=ax.transAxes (panel
    # letters, in-corner labels). Data-rect comparisons don't apply to
    # these because their (x, y) live in axes-fraction, not data coords.
    is_axes_frac: bool = False


@dataclass
class SegmentElement:
    idx: int
    x: float
    y: float
    xend: float
    yend: float


@dataclass
class RectElement:
    idx: int
    xmin: float
    xmax: float
    ymin: float
    ymax: float
    facecolor: Optional[str]
    alpha: float
    kind: str = 'rect'   # 'rect' / 'polygon' / 'circle'
    has_edge: bool = False   # outline-only shapes (no fill) are still "real"


@dataclass
class Issue:
    check: str
    severity: str  # ERROR, WARNING, INFO
    message: str
    suggestion: str = ""


# ---------------------------------------------------------------------------
# Extraction: pull texts, segments, and rects from a matplotlib Axes
# ---------------------------------------------------------------------------
def _extract_elements(fig, ax):
    """Extract texts, segments (arrows/lines), and rects from ax."""
    texts = []
    segments = []
    rects = []

    renderer = fig.canvas.get_renderer()
    inv_transform = ax.transData.inverted()

    # --- Texts: ax.texts + annotations (text part) ---
    text_artists = list(ax.texts)
    for child in ax.get_children():
        if isinstance(child, mtext.Annotation) and child not in text_artists:
            text_artists.append(child)

    for i, t in enumerate(text_artists):
        label = t.get_text()
        if not label or not label.strip():
            continue
        # Detect transform=ax.transAxes text (panel letters, corner
        # labels). We still extract them — text-text overlap and
        # canvas-edge checks apply — but flag them so data-rect checks
        # can skip them (their (x, y) live in axes-frac, comparing to
        # data-coord rects like histogram bars is meaningless).
        is_axes_frac = False
        try:
            if t.get_transform() is ax.transAxes:
                is_axes_frac = True
        except Exception:
            pass
        alpha = t.get_alpha()
        if alpha is None:
            alpha = 1.0

        pos = t.get_position()
        # For annotations the position is already in data coords if
        # xycoords='data' (the default for ax.annotate).  For ax.text
        # it's always data coords.
        x, y = pos

        fontsize = t.get_fontsize()
        ha = t.get_ha()
        va = t.get_va()

        te = TextElement(idx=i, x=x, y=y, label=label,
                         fontsize=fontsize, ha=ha, va=va, alpha=alpha,
                         is_axes_frac=is_axes_frac)

        # Use renderer to get actual bounding box, then convert to data coords
        try:
            bbox_disp = t.get_window_extent(renderer)
            bbox_data = bbox_disp.transformed(inv_transform)
            te.bb_xmin = bbox_data.x0
            te.bb_xmax = bbox_data.x1
            te.bb_ymin = bbox_data.y0
            te.bb_ymax = bbox_data.y1
        except Exception:
            # Fallback: estimate from fontsize and character count
            te = _estimate_text_bbox(te, fig, ax)

        texts.append(te)

    # --- Segments: annotations with arrows, FancyArrowPatch, Line2D ---
    seg_idx = 0
    for child in ax.get_children():
        if isinstance(child, mtext.Annotation):
            arrow_props = child.arrowprops
            if arrow_props is not None:
                # Arrow goes from xytext (start) to xy (end)
                try:
                    xy = child.xy
                    xytext = child.xyann if hasattr(child, 'xyann') else child.get_position()
                    segments.append(SegmentElement(
                        idx=seg_idx, x=xytext[0], y=xytext[1],
                        xend=xy[0], yend=xy[1]))
                    seg_idx += 1
                except Exception:
                    pass

        elif isinstance(child, mpatches.FancyArrowPatch):
            # FancyArrowPatch stores positions differently
            try:
                # _posA_posB is set when using posA/posB constructor
                if hasattr(child, '_posA_posB') and child._posA_posB is not None:
                    posA, posB = child._posA_posB
                else:
                    path = child.get_path()
                    verts = path.vertices
                    posA = verts[0]
                    posB = verts[-1]
                segments.append(SegmentElement(
                    idx=seg_idx, x=posA[0], y=posA[1],
                    xend=posB[0], yend=posB[1]))
                seg_idx += 1
            except Exception:
                pass

        elif isinstance(child, mlines.Line2D):
            xdata = child.get_xdata()
            ydata = child.get_ydata()
            if len(xdata) >= 2:
                # Treat as segment from first to last point
                segments.append(SegmentElement(
                    idx=seg_idx, x=xdata[0], y=ydata[0],
                    xend=xdata[-1], yend=ydata[-1]))
                seg_idx += 1

    # --- Rects + Polygons + Circles: FancyBboxPatch, Rectangle, Polygon, Circle ---
    # Polygons (block arrows) and Circles (dots, markers) are also tracked so
    # the shape-shape overlap check can find arrow-vs-panel and dot-vs-rect
    # collisions. Existing text-rect-centering / text-rect-padding checks
    # filter to kind=='rect' so they don't get spurious hits on those.
    # Skip ax.patch (axes background) — lives in axes-coords.
    rect_idx = 0
    for child in ax.get_children():
        if child is ax.patch:
            continue
        kind = None
        if isinstance(child, (mpatches.FancyBboxPatch, mpatches.Rectangle)):
            kind = 'rect'
        elif isinstance(child, mpatches.Polygon):
            kind = 'polygon'
        elif isinstance(child, mpatches.Circle):
            kind = 'circle'

        if kind is not None:
            try:
                # Prefer logical bbox (no padding/stroke) when available;
                # fall back to the rendered window extent for Polygon/Circle.
                if kind == 'rect':
                    bb = child.get_bbox()
                    xmin, xmax = bb.x0, bb.x1
                    ymin, ymax = bb.y0, bb.y1
                else:
                    bbox_disp = child.get_window_extent(renderer)
                    bbox_data = bbox_disp.transformed(inv_transform)
                    xmin, xmax = bbox_data.x0, bbox_data.x1
                    ymin, ymax = bbox_data.y0, bbox_data.y1

                fc = child.get_facecolor()
                fc_alpha = fc[3] if len(fc) == 4 else 1.0
                facecolor = None if fc_alpha < 0.01 else str(fc)

                # An outline-only shape (no fill but visible stroke) is still
                # a "real" shape that other shapes can overlap.
                ec = child.get_edgecolor()
                ec_alpha = ec[3] if (ec is not None and len(ec) == 4) else 1.0
                lw = child.get_linewidth() or 0
                has_edge = ec_alpha > 0.05 and lw > 0.1

                patch_alpha = child.get_alpha()
                if patch_alpha is None:
                    patch_alpha = 1.0

                rects.append(RectElement(
                    idx=rect_idx,
                    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                    facecolor=facecolor, alpha=patch_alpha,
                    kind=kind, has_edge=has_edge))
                rect_idx += 1
            except Exception:
                pass

    return texts, segments, rects


def _estimate_text_bbox(te, fig, ax):
    """Fallback bbox estimation when renderer bbox fails."""
    fig_w, fig_h = fig.get_size_inches()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]

    mm_per_unit_x = (fig_w * 25.4) / x_range
    mm_per_unit_y = (fig_h * 25.4) / y_range

    sz_mm = te.fontsize * 0.353  # pt to mm
    lines = te.label.split('\n')
    n_lines = len(lines)
    max_nc = max(len(l) for l in lines)

    char_w_mm = sz_mm * 0.55
    char_h_mm = sz_mm * 1.2
    total_w_mm = max_nc * char_w_mm
    total_h_mm = n_lines * char_h_mm

    bbox_w = (total_w_mm / mm_per_unit_x) * 1.1
    bbox_h = (total_h_mm / mm_per_unit_y) * 1.1

    hj = {'center': 0.5, 'right': 1.0, 'left': 0.0}.get(te.ha, 0.5)
    vj = {'center': 0.5, 'top': 1.0, 'bottom': 0.0,
          'baseline': 0.0, 'center_baseline': 0.5}.get(te.va, 0.5)

    te.bb_xmin = te.x - hj * bbox_w
    te.bb_xmax = te.x + (1 - hj) * bbox_w
    te.bb_ymin = te.y - vj * bbox_h
    te.bb_ymax = te.y + (1 - vj) * bbox_h
    return te


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def _segment_intersects_rect(x1, y1, x2, y2, xmin, xmax, ymin, ymax):
    """Liang-Barsky line-rectangle intersection test."""
    dx = x2 - x1
    dy = y2 - y1
    p = [-dx, dx, -dy, dy]
    q = [x1 - xmin, xmax - x1, y1 - ymin, ymax - y1]

    t0, t1 = 0.0, 1.0
    for pk, qk in zip(p, q):
        if abs(pk) < 1e-12:
            if qk < 0:
                return False
        else:
            r = qk / pk
            if pk < 0:
                t0 = max(t0, r)
            else:
                t1 = min(t1, r)
            if t0 > t1:
                return False
    return True


def _point_to_segment_dist(px, py, x1, y1, x2, y2):
    """Minimum distance from point to line segment."""
    dx, dy = x2 - x1, y2 - y1
    len_sq = dx * dx + dy * dy
    if len_sq < 1e-12:
        return math.hypot(px - x1, py - y1)
    t = max(0, min(1, ((px - x1) * dx + (py - y1) * dy) / len_sq))
    proj_x = x1 + t * dx
    proj_y = y1 + t * dy
    return math.hypot(px - proj_x, py - proj_y)


def _rect_to_segment_dist(rxmin, rxmax, rymin, rymax, sx1, sy1, sx2, sy2):
    """Min distance from rectangle boundary to line segment (sampled)."""
    mx, my = (rxmin + rxmax) / 2, (rymin + rymax) / 2
    pts = [
        (rxmin, rymin), (rxmax, rymin), (rxmin, rymax), (rxmax, rymax),
        (mx, rymin), (mx, rymax), (rxmin, my), (rxmax, my),
    ]
    dists = [_point_to_segment_dist(px, py, sx1, sy1, sx2, sy2)
             for px, py in pts]
    # Also check segment endpoints to rect boundary
    for ex, ey in [(sx1, sy1), (sx2, sy2)]:
        ddx = max(rxmin - ex, 0, ex - rxmax)
        ddy = max(rymin - ey, 0, ey - rymax)
        dists.append(math.hypot(ddx, ddy))
    return min(dists)


def _rects_overlap(xmin1, xmax1, ymin1, ymax1, xmin2, xmax2, ymin2, ymax2):
    return (xmin1 < xmax2 and xmax1 > xmin2 and
            ymin1 < ymax2 and ymax1 > ymin2)


# ---------------------------------------------------------------------------
# Report printer
# ---------------------------------------------------------------------------
def _print_report(result):
    s = result['summary']
    print(f"\n=== Figure Layout Validation ===")
    print(f"Extracted: {s['n_texts']} texts, {s['n_segments']} segments, "
          f"{s['n_rects']} rects\n")

    if s['n_errors'] > 0:
        print(f"ERRORS ({s['n_errors']}):")
        for i, iss in enumerate(result['errors'], 1):
            print(f"  [{i}] {iss.check}: {iss.message}")
            if iss.suggestion:
                print(f"      -> {iss.suggestion}")
        print()

    if s['n_warnings'] > 0:
        print(f"WARNINGS ({s['n_warnings']}):")
        for i, iss in enumerate(result['warnings'], 1):
            print(f"  [{i}] {iss.check}: {iss.message}")
            if iss.suggestion:
                print(f"      -> {iss.suggestion}")
        print()

    if s['n_info'] > 0:
        print(f"INFO ({s['n_info']}):")
        for i, iss in enumerate(result['info'], 1):
            print(f"  [{i}] {iss.check}: {iss.message}")
        print()

    if s['n_errors'] == 0 and s['n_warnings'] == 0:
        print("All checks passed.\n")
    elif s['n_errors'] == 0:
        print(f"No errors. {s['n_warnings']} warning(s) to review.\n")
    else:
        print(f"Result: {s['n_errors']} error(s) must be fixed.\n")


# ===========================================================================
# Main validation function
# ===========================================================================
def validate_figure_layout(
    fig,
    ax,
    *,
    collision_buffer=None,
    crowding_distance=None,
    centering_tolerance=0.3,
    alignment_y_tol=None,
    min_arrow_length=0.3,
    rect_padding_min=0.15,
    verbose=True,
):
    """
    Validate layout of a matplotlib annotation-based figure.

    Runs 8 checks:
      1. Text-segment collision (ERROR)
      2. Text-text overlap (ERROR)
      3. Text-segment proximity / crowding (WARNING)
      4. Text-rect centering (WARNING)
      5. Horizontal alignment (INFO)
      6. Arrow minimum length (WARNING)
      7. Text-rect padding / overflow (WARNING/ERROR)
      8. Shape-shape partial overlap (ERROR) — rects + polygons; skips
         intentional containment so background cards holding contents
         are still allowed.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    ax  : matplotlib.axes.Axes
    collision_buffer : float | None
        Extra buffer around text bbox for collision detection (data units).
        If None (default), scales to 0.5% of each axis's data range, per axis.
        This avoids the "any text falsely collides with the bottom spine"
        problem on plots with tiny data ranges (e.g. KDE density 0–0.002),
        and the reverse problem on plots with huge data ranges. Pass an
        explicit float to override.
    crowding_distance : float | None
        Threshold for text-segment proximity warnings (data units).
        If None (default), scales to 4% of mean axis range.
    centering_tolerance : float
        Max allowed horizontal offset of text from enclosing rect center.
    alignment_y_tol : float | None
        Tolerance for grouping texts into horizontal rows (data units).
        If None (default), scales to 2% of y-axis range.
    min_arrow_length : float
        Minimum acceptable arrow/segment length.
    rect_padding_min : float
        Minimum padding between text bbox and enclosing rect edges.
    verbose : bool
        Print formatted report.

    Returns
    -------
    dict with keys: errors, warnings, info, summary
    """

    # Force a draw so renderers can compute extents
    fig.canvas.draw()

    errors = []
    warnings = []
    info = []

    # Phase A: Extract elements
    texts, segments, rects = _extract_elements(fig, ax)

    # Phase B: Auto-scale tolerances when caller didn't override.
    # Multipanel-figure panels often have wildly different data ranges
    # (KDE density 0-0.002 in one cell, bar % 0-100 in another). Absolute
    # default tolerances misbehave in both directions; scale them per-axis.
    x_range = abs(ax.get_xlim()[1] - ax.get_xlim()[0])
    y_range = abs(ax.get_ylim()[1] - ax.get_ylim()[0])
    if collision_buffer is None:
        # 0.5% of each axis's range, applied per-axis below.
        collision_buffer_x = x_range * 0.005
        collision_buffer_y = y_range * 0.005
    else:
        collision_buffer_x = collision_buffer
        collision_buffer_y = collision_buffer
    if crowding_distance is None:
        # Use the SMALLER range as the scale — otherwise an axis with a huge
        # data range (e.g. a schematic with x ∈ [-1200, 3600]) makes the
        # crowding threshold so large that every text-segment pair within
        # the panel gets flagged. Small-range axes (e.g. KDE y ∈ [0, 0.002])
        # rarely have "crowded" pairs to flag anyway, so min is the safer
        # default for multipanel figures.
        crowding_distance = min(x_range, y_range) * 0.04
    if alignment_y_tol is None:
        alignment_y_tol = y_range * 0.02

    # --- Check 1: Text-segment collision (ERROR) ---
    if texts and segments:
        for te in texts:
            if te.alpha < 0.1:
                continue
            bb = (te.bb_xmin - collision_buffer_x, te.bb_xmax + collision_buffer_x,
                  te.bb_ymin - collision_buffer_y, te.bb_ymax + collision_buffer_y)
            for seg in segments:
                if _segment_intersects_rect(
                        seg.x, seg.y, seg.xend, seg.yend,
                        bb[0], bb[1], bb[2], bb[3]):
                    errors.append(Issue(
                        "text_segment_collision", "ERROR",
                        f"Text '{te.label[:40]}' at ({te.x:.1f},{te.y:.1f}) "
                        f"intersects segment ({seg.x:.1f},{seg.y:.1f})->"
                        f"({seg.xend:.1f},{seg.yend:.1f})",
                        "Move text away from the arrow path or reroute the segment"
                    ))

    # --- Check 2: Text-text overlap (ERROR) ---
    if len(texts) > 1:
        for i in range(len(texts)):
            if texts[i].alpha < 0.1:
                continue
            for j in range(i + 1, len(texts)):
                if texts[j].alpha < 0.1:
                    continue
                if _rects_overlap(
                        texts[i].bb_xmin, texts[i].bb_xmax,
                        texts[i].bb_ymin, texts[i].bb_ymax,
                        texts[j].bb_xmin, texts[j].bb_xmax,
                        texts[j].bb_ymin, texts[j].bb_ymax):
                    errors.append(Issue(
                        "text_text_overlap", "ERROR",
                        f"Text '{texts[i].label[:30]}' at ({texts[i].x:.1f},"
                        f"{texts[i].y:.1f}) overlaps '{texts[j].label[:30]}' "
                        f"at ({texts[j].x:.1f},{texts[j].y:.1f})",
                        "Increase vertical or horizontal spacing between these labels"
                    ))

    # --- Check 3: Text-segment proximity (WARNING) ---
    if texts and segments:
        for te in texts:
            if te.alpha < 0.1:
                continue
            for seg in segments:
                dist = _rect_to_segment_dist(
                    te.bb_xmin, te.bb_xmax, te.bb_ymin, te.bb_ymax,
                    seg.x, seg.y, seg.xend, seg.yend)
                if 0 < dist < crowding_distance:
                    warnings.append(Issue(
                        "text_segment_proximity", "WARNING",
                        f"Text '{te.label[:30]}' at ({te.x:.1f},{te.y:.1f}) "
                        f"is {dist:.2f} units from segment ({seg.x:.1f},"
                        f"{seg.y:.1f})->({seg.xend:.1f},{seg.yend:.1f})",
                        "Consider increasing spacing to avoid visual crowding"
                    ))

    # --- Check 4: Text-rect centering (WARNING) ---
    # Compute total axes area for filtering out large container rects
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    _axes_area = (xlim[1] - xlim[0]) * (ylim[1] - ylim[0])
    _max_rect_area = _axes_area * 0.15  # skip rects > 15% of axes area

    if texts and rects:
        for te in texts:
            if te.alpha < 0.1:
                continue
            # Skip axes-fraction text (panel letters, corner labels) —
            # their (x, y) is in axes-frac, comparing to data-coord rects
            # is nonsense (see is_axes_frac note in TextElement).
            if te.is_axes_frac:
                continue
            best_rect = None
            best_area = float('inf')
            for r in rects:
                if r.facecolor is None:
                    continue
                # Skip non-rect shapes (centering inside polygons isn't meaningful)
                if r.kind != 'rect':
                    continue
                # Skip low-alpha rects (background containers)
                if r.alpha < 0.5:
                    continue
                rect_area = (r.xmax - r.xmin) * (r.ymax - r.ymin)
                # Skip very large rects (axes frame, background panels)
                if rect_area > _max_rect_area:
                    continue
                inside_x = r.xmin <= te.x <= r.xmax
                inside_y = r.ymin <= te.y <= r.ymax
                if inside_x and inside_y:
                    area = (r.xmax - r.xmin) * (r.ymax - r.ymin)
                    if area < best_area:
                        best_area = area
                        best_rect = r
            if best_rect is not None:
                rx_center = (best_rect.xmin + best_rect.xmax) / 2
                offset = abs(te.x - rx_center)
                if offset > centering_tolerance:
                    warnings.append(Issue(
                        "text_rect_centering", "WARNING",
                        f"Text '{te.label[:30]}' at x={te.x:.1f} is {offset:.2f} "
                        f"units off-center from rect [{best_rect.xmin:.1f}, "
                        f"{best_rect.xmax:.1f}] (center={rx_center:.1f})",
                        f"Adjust text x to {rx_center:.1f}"
                    ))

    # --- Check 5: Horizontal alignment (INFO) ---
    if len(texts) > 2:
        ys = [(i, te.y) for i, te in enumerate(texts)]
        ys.sort(key=lambda t: t[1])
        groups = [[ys[0]]]
        for k in range(1, len(ys)):
            if abs(ys[k][1] - ys[k - 1][1]) <= alignment_y_tol:
                groups[-1].append(ys[k])
            else:
                groups.append([ys[k]])

        for g in groups:
            if len(g) < 3:
                continue
            mean_y = sum(t[1] for t in g) / len(g)
            for idx, y_val in g:
                if abs(y_val - mean_y) > alignment_y_tol:
                    te = texts[idx]
                    info.append(Issue(
                        "horizontal_alignment", "INFO",
                        f"Text '{te.label[:30]}' at y={te.y:.2f} deviates "
                        f"from group mean y={mean_y:.2f} (group of {len(g)} texts)",
                        f"Consider aligning y to {mean_y:.2f}"
                    ))

    # --- Check 6: Arrow minimum length (WARNING) ---
    if segments:
        for seg in segments:
            length = math.hypot(seg.xend - seg.x, seg.yend - seg.y)
            if length < min_arrow_length:
                warnings.append(Issue(
                    "arrow_min_length", "WARNING",
                    f"Segment ({seg.x:.1f},{seg.y:.1f})->({seg.xend:.1f},"
                    f"{seg.yend:.1f}) length={length:.2f} < minimum "
                    f"{min_arrow_length:.2f}",
                    "Increase spacing between connected elements or remove segment"
                ))

    # --- Check 7: Text-rect padding / overflow (WARNING/ERROR) ---
    if texts and rects:
        for te in texts:
            if te.alpha < 0.1:
                continue
            # Skip axes-fraction text (panel letters, corner labels).
            if te.is_axes_frac:
                continue
            for r in rects:
                if r.facecolor is None:
                    continue
                # Padding is meaningful inside true rects, not block-arrow polygons
                if r.kind != 'rect':
                    continue

                inside = (te.bb_xmin >= r.xmin - 0.01 and
                          te.bb_xmax <= r.xmax + 0.01 and
                          te.bb_ymin >= r.ymin - 0.01 and
                          te.bb_ymax <= r.ymax + 0.01)

                center_inside = (r.xmin <= te.x <= r.xmax and
                                 r.ymin <= te.y <= r.ymax)

                if not inside and center_inside:
                    overflows = []
                    if te.bb_xmin < r.xmin:
                        overflows.append("left")
                    if te.bb_xmax > r.xmax:
                        overflows.append("right")
                    if te.bb_ymin < r.ymin:
                        overflows.append("bottom")
                    if te.bb_ymax > r.ymax:
                        overflows.append("top")
                    if overflows:
                        errors.append(Issue(
                            "text_rect_overflow", "ERROR",
                            f"Text '{te.label[:30]}' overflows rect "
                            f"[{r.xmin:.1f},{r.xmax:.1f}]x"
                            f"[{r.ymin:.1f},{r.ymax:.1f}] on: "
                            f"{', '.join(overflows)}",
                            "Reduce text size or widen the rect"
                        ))
                elif inside:
                    left_pad = te.bb_xmin - r.xmin
                    right_pad = r.xmax - te.bb_xmax
                    top_pad = r.ymax - te.bb_ymax
                    bottom_pad = te.bb_ymin - r.ymin
                    pads = {'left': left_pad, 'right': right_pad,
                            'top': top_pad, 'bottom': bottom_pad}
                    min_pad = min(pads.values())
                    if min_pad < rect_padding_min:
                        tight = min(pads, key=pads.get)
                        warnings.append(Issue(
                            "text_rect_padding", "WARNING",
                            f"Text '{te.label[:30]}' has only {min_pad:.2f} "
                            f"units padding on {tight} edge of rect "
                            f"[{r.xmin:.1f},{r.xmax:.1f}]x"
                            f"[{r.ymin:.1f},{r.ymax:.1f}]",
                            "Increase rect size or reduce text size"
                        ))

    # --- Check 8: Shape-shape partial overlap (ERROR) ---
    # Flags pairs of visible shapes (rects + polygons + circles) whose bboxes
    # overlap unless one fully contains the other (intentional containment).
    # A shape is considered "visible" if it has a fill OR a visible stroke, so
    # outline-only panels still get checked against arrows passing through.
    # Circle-circle pairs are skipped — overlapping circles are typically
    # decorative (smoke clouds, cell-cluster glyphs) and produce noise.
    def _shape_contained(inner, outer, tol=0.02):
        return (inner.xmin >= outer.xmin - tol and
                inner.xmax <= outer.xmax + tol and
                inner.ymin >= outer.ymin - tol and
                inner.ymax <= outer.ymax + tol)

    def _is_visible(r):
        return (r.facecolor is not None) or r.has_edge

    overlap_tol = 0.015   # tight: catch even sub-tenth-unit overlaps
    if len(rects) > 1:
        for i in range(len(rects)):
            r1 = rects[i]
            if not _is_visible(r1):
                continue
            for j in range(i + 1, len(rects)):
                r2 = rects[j]
                if not _is_visible(r2):
                    continue
                # Skip pairs where both are circles — usually decorative
                # (smoke clusters, cell nuclei) and not a layout problem.
                if r1.kind == 'circle' and r2.kind == 'circle':
                    continue
                ox = min(r1.xmax, r2.xmax) - max(r1.xmin, r2.xmin)
                oy = min(r1.ymax, r2.ymax) - max(r1.ymin, r2.ymin)
                if ox <= overlap_tol or oy <= overlap_tol:
                    continue
                if _shape_contained(r1, r2) or _shape_contained(r2, r1):
                    continue
                errors.append(Issue(
                    "shape_shape_overlap", "ERROR",
                    f"{r1.kind} [{r1.xmin:.2f},{r1.xmax:.2f}]x"
                    f"[{r1.ymin:.2f},{r1.ymax:.2f}] partially overlaps "
                    f"{r2.kind} [{r2.xmin:.2f},{r2.xmax:.2f}]x"
                    f"[{r2.ymin:.2f},{r2.ymax:.2f}] "
                    f"(overlap {ox:.2f} x {oy:.2f})",
                    "Reposition or resize so shapes either don't overlap, "
                    "or one fully contains the other"
                ))

    # --- Assemble and report ---
    result = {
        'errors': errors,
        'warnings': warnings,
        'info': info,
        'summary': {
            'n_texts': len(texts),
            'n_segments': len(segments),
            'n_rects': len(rects),
            'n_errors': len(errors),
            'n_warnings': len(warnings),
            'n_info': len(info),
        }
    }

    if verbose:
        _print_report(result)

    return result


# ============================================================================
# validate_multipanel_layout
#
# Companion validator for matplotlib multi-panel composites (anything built
# with plt.subplots(rows, cols, ...) or a custom Axes grid). The single-axes
# validator above runs per-panel via validate_figure_layout(); this one
# catches the COMPOSITE-LEVEL issues those can't see:
#
#   1. Text elements clipping the figure boundaries (panel labels too close
#      to the figure top/bottom/left/right edges)
#   2. Text from one panel overlapping content in another panel (e.g. a
#      panel-letter D landing on top of the xtick labels of the panel above)
#
# Together these catch the failure mode where a multipanel figure renders
# without per-panel validator complaints but reads as crowded or clipped
# at the composite level.
#
# Tolerance policy (defaults sized for the project's standard 16-18" wide
# composites):
#   edge_padding = 0.005 (figure-fraction) — text bbox closer than 0.5 % of
#                                            figure span to an edge is flagged
#   min_gap      = 0.005 — cross-panel text pairs within 0.5 % of each
#                          other are flagged as crowded
#
# Returns the same {errors, warnings, info, summary} dict shape as
# validate_figure_layout() so callers can pattern-match either result.
# ============================================================================
def validate_multipanel_layout(
    fig,
    *,
    edge_padding=0.005,
    min_gap=0.005,
    skip_internal_axes_text=True,
    verbose=True,
):
    """
    Validate composite-level layout of a multi-panel matplotlib figure.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    edge_padding : float
        Minimum padding from text bbox to figure edges, in figure-fraction
        coordinates. 0.005 = 0.5 %. Text closer than this is flagged ERROR.
    min_gap : float
        Minimum gap between text bboxes belonging to DIFFERENT axes (or
        between an axis text and a figure-level fig.text), in figure-
        fraction. Below this they're flagged ERROR (overlap) or WARNING
        (crowding) depending on whether the bboxes actually intersect.
    skip_internal_axes_text : bool
        If True, ignore text inside an axes that lives inside the data area
        (i.e. text not anchored to axes coordinates) — that's per-panel
        content and the single-axes validator handles it. The check still
        sees tick labels, axis labels, titles, and any text with
        `transform=ax.transAxes` (panel letters, etc.).
    verbose : bool
        Print formatted report.

    Returns
    -------
    dict with keys: errors, warnings, info, summary
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_w_px = fig.get_figwidth() * fig.dpi
    fig_h_px = fig.get_figheight() * fig.dpi

    errors = []
    warnings = []
    info = []

    def _bbox_frac(t):
        try:
            bb_px = t.get_window_extent(renderer=renderer)
        except Exception:
            return None
        if bb_px.width == 0 or bb_px.height == 0:
            return None
        return (bb_px.x0 / fig_w_px, bb_px.x1 / fig_w_px,
                bb_px.y0 / fig_h_px, bb_px.y1 / fig_h_px)

    def _is_axes_data_text(t, ax):
        # transAxes annotations (panel letters, in-axes labels with axes
        # coords) have transform == ax.transAxes; data-space text has
        # transform == ax.transData. Skip data-space text by default —
        # validate_figure_layout() handles in-axes per-panel checks.
        if not skip_internal_axes_text:
            return False
        try:
            return t.get_transform() is ax.transData
        except Exception:
            return False

    # Phase A: Collect all text elements with figure-fraction bboxes.
    items = []  # list of dicts {text, bb, owner}
    for ax in fig.get_axes():
        # Schematic / drawing axes use ax.set_axis_off() — their tick labels
        # and axis labels are not visible to the reader and shouldn't be
        # checked. matplotlib still generates phantom tick label Text objects
        # ('0', '20', '40', ...) even with axis off, so this guard is
        # necessary to avoid false-positive clip / overlap errors against
        # invisible decorations.
        axis_on = getattr(ax, "axison", True)
        # Tick labels (x and y)
        if axis_on:
            for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
                if not t.get_text() or not t.get_visible():
                    continue
                bb = _bbox_frac(t)
                if bb is not None:
                    items.append({"text": t.get_text(), "bb": bb,
                                  "owner": ax, "kind": "tick"})
            # Axis labels + title
            for t in (ax.xaxis.label, ax.yaxis.label, ax.title):
                if not t.get_text() or not t.get_visible():
                    continue
                bb = _bbox_frac(t)
                if bb is not None:
                    items.append({"text": t.get_text(), "bb": bb,
                                  "owner": ax, "kind": "axislabel"})
        # In-axes text (panel letters, annotations). Skip data-space text.
        for t in ax.texts:
            if not t.get_text() or not t.get_visible():
                continue
            if _is_axes_data_text(t, ax):
                continue
            bb = _bbox_frac(t)
            if bb is not None:
                items.append({"text": t.get_text(), "bb": bb,
                              "owner": ax, "kind": "axes_text"})
    # Figure-level text (suptitle, fig.text(...))
    for t in fig.texts:
        if not t.get_text() or not t.get_visible():
            continue
        bb = _bbox_frac(t)
        if bb is not None:
            items.append({"text": t.get_text(), "bb": bb,
                          "owner": None, "kind": "fig_text"})

    # Phase B: Edge-clip check (text bbox crossing or grazing figure edges).
    for it in items:
        x0, x1, y0, y1 = it["bb"]
        label = it["text"][:40].replace("\n", " ")
        if y1 > 1.0 - edge_padding:
            errors.append(Issue(
                "multipanel_clip_top", "ERROR",
                f"Text '{label}' ({it['kind']}) crosses or touches figure "
                f"top edge: bb_y1={y1:.4f} > 1 − {edge_padding}",
                "Increase top margin (fig.subplots_adjust(top=...)) or "
                "move the text down (lower transAxes y, or fig.text y)"))
        if y0 < edge_padding:
            errors.append(Issue(
                "multipanel_clip_bottom", "ERROR",
                f"Text '{label}' ({it['kind']}) crosses or touches figure "
                f"bottom edge: bb_y0={y0:.4f} < {edge_padding}",
                "Increase bottom margin or move text up"))
        if x0 < edge_padding:
            errors.append(Issue(
                "multipanel_clip_left", "ERROR",
                f"Text '{label}' ({it['kind']}) crosses or touches figure "
                f"left edge: bb_x0={x0:.4f} < {edge_padding}",
                "Increase left margin or move text right"))
        if x1 > 1.0 - edge_padding:
            errors.append(Issue(
                "multipanel_clip_right", "ERROR",
                f"Text '{label}' ({it['kind']}) crosses or touches figure "
                f"right edge: bb_x1={x1:.4f} > 1 − {edge_padding}",
                "Increase right margin or move text left"))

    # Phase C: Cross-panel text overlap / crowding (skip pairs that share
    # the same owner — those are within a single panel and the single-axes
    # validator handles them).
    #
    # Note on tick labels: adjacent tick labels at panel edges in a grid
    # ALWAYS sit within a few % of each other — that's just how matplotlib
    # subplot grids work, and it isn't a layout bug. We:
    #   - flag cross-panel OVERLAP (actually intersecting bboxes) for every
    #     text kind, since real overlap is always wrong;
    #   - flag CROWDING (within min_gap but not intersecting) only for
    #     non-tick text kinds (panel letters, axis labels, titles, fig.text).
    NON_TICK_KINDS = {"axislabel", "axes_text", "fig_text"}
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i], items[j]
            if a["owner"] is not None and a["owner"] is b["owner"]:
                # Same-axis pairs are usually skipped (adjacent ticks are
                # inevitable in a grid) — but panel letters vs the ticks
                # of their OWN axis is a real defect class (2026-07-09
                # SF34 "C" overlapping ytick "60"). Keep those pairs
                # checked; skip everything else with the same owner.
                pair_kinds = {a["kind"], b["kind"]}
                if not (pair_kinds == {"axes_text", "tick"} or
                        pair_kinds == {"axes_text", "axislabel"}):
                    continue
            ax0, ax1, ay0, ay1 = a["bb"]
            bx0, bx1, by0, by1 = b["bb"]
            # Check overlap (bboxes intersect)
            overlap = (ax1 >= bx0 and ax0 <= bx1 and
                       ay1 >= by0 and ay0 <= by1)
            if overlap:
                lab_a = a["text"][:30].replace("\n", " ")
                lab_b = b["text"][:30].replace("\n", " ")
                errors.append(Issue(
                    "multipanel_cross_overlap", "ERROR",
                    f"Cross-panel text overlap: '{lab_a}' ({a['kind']}) "
                    f"vs '{lab_b}' ({b['kind']})",
                    "Increase hspace/wspace in subplots_adjust, or "
                    "reposition one of the texts"))
                continue
            # Crowding warning — only for structural text (panel letters,
            # axis labels, titles, fig.text). Skip tick labels.
            if a["kind"] not in NON_TICK_KINDS or b["kind"] not in NON_TICK_KINDS:
                continue
            # Gap definition:
            #   - dx = horizontal gap (0 if bboxes overlap in x);
            #   - dy = vertical gap   (0 if bboxes overlap in y).
            # We only flag "crowding" when bboxes overlap in AT LEAST one
            # axis (dx == 0 or dy == 0) — then the gap is the orthogonal
            # axis's separation. Diagonally separated bboxes (dx > 0 AND
            # dy > 0) are by construction in different rows AND different
            # columns; they aren't competing for space and shouldn't be
            # flagged.
            dx = max(ax0 - bx1, bx0 - ax1, 0.0)
            dy = max(ay0 - by1, by0 - ay1, 0.0)
            if dx > 0 and dy > 0:
                continue
            gap = max(dx, dy)
            if gap > 0 and gap < min_gap:
                lab_a = a["text"][:30].replace("\n", " ")
                lab_b = b["text"][:30].replace("\n", " ")
                warnings.append(Issue(
                    "multipanel_cross_crowding", "WARNING",
                    f"Cross-panel structural-text crowding: "
                    f"'{lab_a}' ({a['kind']}) vs '{lab_b}' ({b['kind']}) — "
                    f"gap {gap:.4f} < min_gap {min_gap}",
                    "Increase subplot spacing or reposition the texts"))

    result = {
        "errors": errors,
        "warnings": warnings,
        "info": info,
        "summary": {
            "n_text_items": len(items),
            "n_axes": len(fig.get_axes()),
            "n_errors": len(errors),
            "n_warnings": len(warnings),
            "n_info": len(info),
        },
    }

    if verbose:
        print(f"\nvalidate_multipanel_layout: "
              f"{result['summary']['n_errors']} errors, "
              f"{result['summary']['n_warnings']} warnings, "
              f"across {result['summary']['n_axes']} panels "
              f"({result['summary']['n_text_items']} text items)")
        for e in errors:
            print(f"  [{e.severity}] {e.check}: {e.message}")
        for w in warnings:
            print(f"  [{w.severity}] {w.check}: {w.message}")

    return result


# ============================================================================
# assert_style_symmetric
#
# Fails if sibling axes in the same figure use different font sizes for the
# same structural role (xtick labels, ytick labels, xaxis label, yaxis
# label). Catches the SF33 Panel B < A/C bug — a single hardcoded local
# override drifts one panel out of the module's shared BODY_FS.
#
# Rationale: `ggplot_style.py` publishes BODY_FS as the canonical size.
# Every panel imports it. Any panel-local literal like fontsize=11 or
# LABEL_FS=12 introduces cross-panel drift that the multipanel layout
# validator won't catch (it checks positions, not sizes).
#
# The check is intentionally strict — any deviation between sibling axes on
# the SAME role is flagged. If different roles have different sizes (xtick
# ≠ ylabel), that's fine and expected; we only compare like to like.
# ============================================================================
def assert_style_symmetric(fig, *, tolerance_pt=0.5):
    """Raise if font sizes for the same role differ across sibling axes.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    tolerance_pt : float
        Allowed drift between sizes for the same role before flagging.
        Default 0.5 pt — anything larger is a real inconsistency.
    """
    fig.canvas.draw()

    def _first_visible_size(labels):
        for t in labels:
            if t.get_visible() and t.get_text() and t.get_text().strip():
                return t.get_fontsize()
        return None

    per_axis = {}
    for i, ax in enumerate(fig.axes):
        # Skip axes with no visible content (axis-off schematic axes).
        if not getattr(ax, "axison", True):
            continue
        d = {}
        sz = _first_visible_size(ax.get_xticklabels())
        if sz is not None:
            d["xtick"] = sz
        sz = _first_visible_size(ax.get_yticklabels())
        if sz is not None:
            d["ytick"] = sz
        if ax.xaxis.label.get_text() and ax.xaxis.label.get_text().strip():
            d["xlabel"] = ax.xaxis.label.get_fontsize()
        if ax.yaxis.label.get_text() and ax.yaxis.label.get_text().strip():
            d["ylabel"] = ax.yaxis.label.get_fontsize()
        if d:
            per_axis[i] = d

    if len(per_axis) < 2:
        return  # nothing to compare

    offenders = []
    for role in ("xtick", "ytick", "xlabel", "ylabel"):
        sizes = {}  # size -> [axis indices]
        for i, d in per_axis.items():
            if role in d:
                sz = round(d[role], 1)
                sizes.setdefault(sz, []).append(i)
        if len(sizes) < 2:
            continue
        unique = sorted(sizes.keys())
        if unique[-1] - unique[0] > tolerance_pt:
            parts = [f"{sz}pt on axes {sizes[sz]}" for sz in unique]
            offenders.append(f"  {role}: " + " vs ".join(parts))

    if offenders:
        raise AssertionError(
            "Font size mismatch across sibling axes — every panel must use "
            "the same BODY_FS from ggplot_style; remove local literals "
            "(fontsize=11, LABEL_FS=12, BODY_FS-2, …):\n"
            + "\n".join(offenders)
        )


# ============================================================================
# assert_docx_readable
#
# Fails if any text in the figure renders below the readability floor
# (PAPER_FLOOR_PT = 9 pt) when the figure is scaled to the docx content
# width (6.5 "). Catches the "native fontsize=11 on a 16 " render reads at
# 4.5 pt in the delivered document" class of bug that no per-panel
# validator can see — because it depends on the docx pipeline, not the
# per-panel layout.
#
# Calibration: Pete's 2026-07-09 review said SF25 (6.5 " → 14 pt effective)
# reads plenty big, SF29 (11 " × 14 pt native → 8.3 pt effective) reads a
# touch small, SF30 (8 " × 14 pt → 11.4 pt) reads well. So floor = 9 pt.
# ============================================================================
def assert_docx_readable(fig, *, native_width_in=None,
                          floor_pt=None, content_w_in=None):
    """Raise if any text renders below `floor_pt` at docx scale.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    native_width_in : float | None
        Native figure width in inches. If None, uses fig.get_figwidth().
    floor_pt : float | None
        Minimum acceptable docx-rendered point size. If None, uses
        PAPER_FLOOR_PT from ggplot_style (9 pt).
    content_w_in : float | None
        Docx content-width to scale against. If None, uses
        PAPER_CONTENT_W_IN from ggplot_style (6.5 ").
    """
    # Import here to avoid a hard cross-module dep at module load time —
    # keeps this validator usable in isolation.
    try:
        from ggplot_style import PAPER_FLOOR_PT, PAPER_CONTENT_W_IN
    except ImportError:
        PAPER_FLOOR_PT_default = 9.0
        PAPER_CONTENT_W_IN_default = 6.5
        PAPER_FLOOR_PT = PAPER_FLOOR_PT_default
        PAPER_CONTENT_W_IN = PAPER_CONTENT_W_IN_default

    if floor_pt is None:
        floor_pt = PAPER_FLOOR_PT
    if content_w_in is None:
        content_w_in = PAPER_CONTENT_W_IN
    if native_width_in is None:
        native_width_in = fig.get_figwidth()

    scale = content_w_in / native_width_in
    fig.canvas.draw()

    offenders = []

    def check(t, kind):
        if not t.get_visible():
            return
        s = t.get_text() if hasattr(t, "get_text") else ""
        if not s or not s.strip():
            return
        native_fs = t.get_fontsize()
        docx_fs = native_fs * scale
        if docx_fs < floor_pt:
            label = s[:40].replace("\n", " ")
            offenders.append(
                f"  '{label}' ({kind}): native={native_fs:.1f}pt × "
                f"{scale:.3f} → docx={docx_fs:.2f}pt < {floor_pt}pt"
            )

    for ax in fig.axes:
        if not getattr(ax, "axison", True):
            continue
        for t in ax.get_xticklabels():
            check(t, "xtick")
        for t in ax.get_yticklabels():
            check(t, "ytick")
        check(ax.xaxis.label, "xlabel")
        check(ax.yaxis.label, "ylabel")
        check(ax.title, "title")
        for t in ax.texts:
            check(t, "annot")
    for t in fig.texts:
        check(t, "fig_text")

    if offenders:
        raise AssertionError(
            f"Text below {floor_pt}pt readability floor at docx scale "
            f"(native width {native_width_in}\" → 6.5\" content, "
            f"scale {scale:.3f}×):\n"
            + "\n".join(offenders)
            + "\n\nFix: bump native fontsize (BODY_FS = docx_body_fs(width))"
              " or shrink the native figure width."
        )


# ============================================================================
# Data-occupancy mask + text-vs-data overlap validator
#
# Motivation: `validate_figure_layout` extracts artists as bbox rectangles
# and treats KDE fills, violin polygons, and histogram bars as opaque
# rectangles spanning their full envelope. That misses two important
# cases:
#   1. A `transform=ax.transAxes` label sitting inside a KDE fill's peak
#      region — the label is "inside the bbox" but the bbox is the whole
#      axes.
#   2. Finding open space where a new label could be placed without
#      overlap.
#
# Approach: rasterize the axes with text hidden, subtract the panel
# background + gridlines, and get a boolean mask of "data-drawn" pixels.
# Then check each text bbox against that mask. `find_open_regions` inverts
# the mask so callers can query where a `(w, h)` text could fit.
# ============================================================================
def _compute_data_mask(fig, ax, *, occupied_threshold=25, dilate_px=1):
    """Rasterize `ax` with text hidden — return a bool mask of data pixels.

    Returns
    -------
    mask : np.ndarray of shape (H, W), dtype bool
        True where data was drawn (curves, fills, bars, violins, scatter).
    ax_pixel_bbox : (x0, y0, x1, y1) in figure-pixel coords
        The pixel region the mask covers.
    """
    hidden = []
    def hide(t):
        if t is None:
            return
        hidden.append((t, t.get_visible()))
        t.set_visible(False)
    for t in ax.texts:
        hide(t)
    for t in ax.get_xticklabels() + ax.get_yticklabels():
        hide(t)
    hide(ax.xaxis.label)
    hide(ax.yaxis.label)
    hide(ax.title)
    # Also hide the legend so it isn't counted as "data" during the
    # rasterization — otherwise legend-vs-data clearance checks flag
    # the legend as overlapping its own bounding box.
    hide(ax.get_legend())

    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        bbox = ax.get_window_extent(renderer=renderer)
        fig_w = int(round(fig.bbox.width))
        fig_h = int(round(fig.bbox.height))

        # Buffer: (H, W, 4) uint8, RGBA, y=0 is TOP row (matplotlib convention).
        buf = np.asarray(fig.canvas.buffer_rgba())

        # ax bbox is in display coords with y=0 at BOTTOM. Convert to
        # buffer's y=0-at-top convention.
        x0 = max(0, int(math.floor(bbox.x0)))
        x1 = min(fig_w, int(math.ceil(bbox.x1)))
        by0 = max(0, int(math.floor(fig_h - bbox.y1)))
        by1 = min(fig_h, int(math.ceil(fig_h - bbox.y0)))
        if x1 <= x0 or by1 <= by0:
            return np.zeros((0, 0), dtype=bool), (0, 0, 0, 0)

        sub = buf[by0:by1, x0:x1, :3].astype(np.int16)

        # Non-background = not the panel facecolor.
        bg_rgb = np.array(
            [int(round(c * 255)) for c in mcolors.to_rgb(ax.get_facecolor())],
            dtype=np.int16,
        )
        diff_bg = np.abs(sub - bg_rgb).sum(axis=-1)
        mask = diff_bg > occupied_threshold

        # Also exclude the ggplot-style white gridlines (structural, not
        # data). If ax.grid was drawn in a color close to white, ignore
        # pixels close to that color.
        grid_rgb = np.array([255, 255, 255], dtype=np.int16)
        diff_grid = np.abs(sub - grid_rgb).sum(axis=-1)
        mask = mask & (diff_grid > occupied_threshold)

        # Optional dilation so a thin curve isn't easily missed by bbox
        # sampling. Cheap max-filter dilation without scipy.
        if dilate_px > 0:
            m = mask
            for _ in range(dilate_px):
                shifted = np.zeros_like(m)
                shifted[1:, :]  |= m[:-1, :]
                shifted[:-1, :] |= m[1:, :]
                shifted[:, 1:]  |= m[:, :-1]
                shifted[:, :-1] |= m[:, 1:]
                m = m | shifted
            mask = m

        return mask, (x0, by0, x1, by1)
    finally:
        for t, v in hidden:
            t.set_visible(v)
        fig.canvas.draw()


def assert_data_free(fig, ax, *, max_occupied_frac=0.15,
                     tick_labels_too=False):
    """Raise if any in-axes text sits on data-drawn pixels.

    Parameters
    ----------
    max_occupied_frac : float
        Allowed fraction of a text bbox's area that may overlap data
        pixels. Default 0.15 (15 %). Catches labels buried in a KDE peak
        or violin body while tolerating peripheral captions that graze a
        few short bars.
    tick_labels_too : bool
        If True, also check tick labels. Default False — tick labels are
        typically placed outside the plot area by matplotlib and rarely
        collide with data.
    """
    mask, (mx0, my0, mx1, my1) = _compute_data_mask(fig, ax)
    if mask.size == 0:
        return

    fig_h = int(round(fig.bbox.height))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    axes_area = float(mask.shape[0] * mask.shape[1])

    def _check(t, kind):
        if not t.get_visible():
            return None
        s = t.get_text() if hasattr(t, "get_text") else ""
        if not s or not s.strip():
            return None
        # Skip texts with an opaque bbox background — they visually mask
        # the data beneath them, so overlap doesn't hide the label. Common
        # for annotation boxes with `bbox=dict(facecolor="white", ...)`.
        try:
            bp = t.get_bbox_patch()
        except Exception:
            bp = None
        if bp is not None and bp.get_visible():
            fc = bp.get_facecolor()
            if fc is not None and len(fc) == 4 and fc[3] > 0.5:
                return None
        try:
            tb = t.get_window_extent(renderer=renderer)
        except Exception:
            return None
        # Convert text bbox from display (y=0 at bottom) to mask coords
        # (y=0 at top), then subtract mask origin.
        tx0 = max(0, int(math.floor(tb.x0 - mx0)))
        tx1 = min(mask.shape[1], int(math.ceil(tb.x1 - mx0)))
        ty1 = max(0, int(math.floor(fig_h - tb.y1 - my0)))
        ty0 = min(mask.shape[0], int(math.ceil(fig_h - tb.y0 - my0)))
        if tx1 <= tx0 or ty0 <= ty1:
            return None
        sub = mask[ty1:ty0, tx0:tx1]
        if sub.size == 0:
            return None
        # Skip small data markers (median values on violins, count labels
        # on bar tops, etc.). Two heuristics — pass either → skip:
        #   1. bbox < 1.5 % of axes area (small in absolute terms)
        #   2. text ≤ 8 chars (single number / short value label)
        # These catch the common "0.86" / "n=190" / "1,166" style markers
        # that panel code intentionally places on data. Longer floating
        # annotations still get checked.
        if sub.size / axes_area < 0.015:
            return None
        if len(s.strip()) <= 8:
            return None
        occ = float(sub.sum()) / sub.size
        return (s, kind, occ)

    offenders = []
    for t in ax.texts:
        r = _check(t, "axes_text")
        if r is not None and r[2] > max_occupied_frac:
            offenders.append(r)
    if tick_labels_too:
        for t in ax.get_xticklabels() + ax.get_yticklabels():
            r = _check(t, "tick")
            if r is not None and r[2] > max_occupied_frac:
                offenders.append(r)

    if offenders:
        lines = [
            f"  '{s[:40]}' ({k}): {occ*100:.1f}% of bbox area sits on data"
            for (s, k, occ) in offenders
        ]
        raise AssertionError(
            "In-axes text overlaps data-drawn pixels (KDE fill, violin body, "
            "bars, curve). Threshold: "
            f"{max_occupied_frac*100:.0f}%.\n"
            + "\n".join(lines)
            + "\n\nFix: move the text into a data-free region "
              "(call find_open_regions(fig, ax) to locate one), lift the axes "
              "ylim above the peak, or place the text in an ax.legend() "
              "outside the plot area."
        )


def find_open_regions(fig, ax, *, min_run_px=None):
    """Return the un-occupied region of `ax` as an axes-fraction bbox list.

    Rasterizes the axes (via `_compute_data_mask`) and returns the largest
    rectangle of un-occupied pixels in each of nine quadrant-anchor slots
    (TL, TC, TR, ML, MC, MR, BL, BC, BR) — expressed in axes-fraction
    coords `(x0, y0, x1, y1)` with (0,0)=bottom-left, (1,1)=top-right of
    the axes, matching `transform=ax.transAxes`.

    Callers can pick the slot whose free region fits their text's needs:

        free = find_open_regions(fig, ax)
        top_left = free["TL"]        # (x0, y0, x1, y1) in axes fraction
        if top_left is not None and (top_left[3] - top_left[1]) > my_height:
            ax.text(top_left[0] + 0.02, top_left[3] - 0.02,
                    "my label", transform=ax.transAxes,
                    ha="left", va="top")

    Slots with no clear region return None.
    """
    mask, (mx0, my0, mx1, my1) = _compute_data_mask(fig, ax)
    if mask.size == 0:
        return {k: None for k in
                ("TL", "TC", "TR", "ML", "MC", "MR", "BL", "BC", "BR")}
    H, W = mask.shape
    if min_run_px is None:
        min_run_px = max(10, min(H, W) // 20)

    def _largest_free_rect_in(sub):
        """Return the largest all-False axis-aligned rectangle in `sub`.

        Uses the classical maximal-rectangle-in-histogram DP. Returns
        (r0, c0, r1, c1) in local sub-array coords, or None.
        """
        if sub.size == 0 or sub.all():
            return None
        rows, cols = sub.shape
        heights = np.zeros(cols, dtype=np.int32)
        best = (0, 0, 0, 0, 0)  # area, r0, c0, r1, c1
        for r in range(rows):
            # heights[c] = current column's un-occupied run height ending at r
            for c in range(cols):
                heights[c] = 0 if sub[r, c] else heights[c] + 1
            # Largest rectangle in histogram `heights`
            stack = []
            for c in range(cols + 1):
                cur_h = heights[c] if c < cols else 0
                start = c
                while stack and stack[-1][1] > cur_h:
                    idx, h = stack.pop()
                    area = h * (c - idx)
                    if area > best[0] and h >= min_run_px and (c - idx) >= min_run_px:
                        # r0 = r - h + 1, r1 = r
                        best = (area, r - h + 1, idx, r, c - 1)
                    start = idx
                stack.append((start, cur_h))
        if best[0] == 0:
            return None
        return best[1:]  # (r0, c0, r1, c1)

    def _slot_to_frac(rect_local, x_off, y_off_top):
        if rect_local is None:
            return None
        r0, c0, r1, c1 = rect_local
        # Convert to axes-frac coords (0,0)=bottom-left, (1,1)=top-right.
        x0f = (c0 + x_off) / W
        x1f = (c1 + 1 + x_off) / W
        # Rows count from top of image; ax fraction Y counts from bottom.
        y0f = 1.0 - (r1 + 1 + y_off_top) / H
        y1f = 1.0 - (r0 + y_off_top) / H
        return (x0f, y0f, x1f, y1f)

    # Divide the axes into nine quadrants and find largest free rect in each.
    Hs = [(0, H // 3), (H // 3, 2 * H // 3), (2 * H // 3, H)]
    Ws = [(0, W // 3), (W // 3, 2 * W // 3), (2 * W // 3, W)]
    names = [
        ("TL", 0, 0), ("TC", 0, 1), ("TR", 0, 2),
        ("ML", 1, 0), ("MC", 1, 1), ("MR", 1, 2),
        ("BL", 2, 0), ("BC", 2, 1), ("BR", 2, 2),
    ]
    out = {}
    for name, ri, ci in names:
        r0, r1 = Hs[ri]
        c0, c1 = Ws[ci]
        sub = mask[r0:r1, c0:c1]
        rect = _largest_free_rect_in(sub)
        out[name] = _slot_to_frac(rect, c0, r0) if rect else None
    return out


def assert_tick_labels_disjoint(fig, *, tol_px=1.0):
    """Raise if two tick labels of the same axis overlap in bbox.

    Matplotlib does not auto-rotate or truncate colliding tick labels;
    long two-line xticklabels can silently overlap adjacent ticks. This
    check walks adjacent same-axis tick-label pairs and flags overlap.

    Parameters
    ----------
    tol_px : float
        Bbox overlap smaller than this (in display pixels) is treated
        as touching-not-overlapping and ignored.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    def _bbox(t):
        try:
            bb = t.get_window_extent(renderer=renderer)
        except Exception:
            return None
        if bb.width == 0 or bb.height == 0:
            return None
        return (bb.x0, bb.y0, bb.x1, bb.y1)

    def _overlaps(a, b, tol):
        return (a[0] < b[2] - tol and a[2] > b[0] + tol and
                a[1] < b[3] - tol and a[3] > b[1] + tol)

    offenders = []
    for ax_idx, ax in enumerate(fig.axes):
        if not getattr(ax, "axison", True):
            continue
        for side, get_labels in (
            ("xtick", ax.get_xticklabels),
            ("ytick", ax.get_yticklabels),
        ):
            labels = [
                t for t in get_labels()
                if t.get_visible()
                and t.get_text() is not None
                and t.get_text().strip()
            ]
            # Only adjacent pairs — same axis, tick positions are ordered,
            # a label pair that's not adjacent almost never overlaps.
            for i in range(len(labels) - 1):
                bb1 = _bbox(labels[i])
                bb2 = _bbox(labels[i + 1])
                if bb1 is None or bb2 is None:
                    continue
                if _overlaps(bb1, bb2, tol_px):
                    offenders.append(
                        f"  ax{ax_idx} {side}: "
                        f"'{labels[i].get_text()[:24]}' overlaps "
                        f"'{labels[i + 1].get_text()[:24]}'"
                    )

    if offenders:
        raise AssertionError(
            "Tick labels overlap within an axis (matplotlib doesn't "
            "auto-rotate/truncate). Fix: rotate xticks, shorten labels, "
            "widen the axis, or thin the tick set.\n"
            + "\n".join(offenders)
        )


def _dilate_mask(mask, n_iter):
    """4-neighbor max-filter dilation, n_iter iterations. In-place safe.
    Returns a new bool array of the same shape as `mask`."""
    if n_iter <= 0 or mask.size == 0:
        return mask.copy()
    m = mask.copy()
    for _ in range(n_iter):
        shifted = np.zeros_like(m)
        shifted[1:, :]  |= m[:-1, :]
        shifted[:-1, :] |= m[1:, :]
        shifted[:, 1:]  |= m[:, :-1]
        shifted[:, :-1] |= m[:, 1:]
        m = m | shifted
    return m


def assert_legend_clear(fig, ax, *, clearance_px=4,
                        overlap_tolerance_px=0):
    """Raise if the legend bbox on `ax` overlaps (or comes within
    `clearance_px` of) any data-drawn pixel.

    Zero-tolerance by default: a legend that even TOUCHES the boundary
    of a data element (a curve, a KDE fill edge, a bar top) fails. The
    clearance dilation grows the data mask by `clearance_px` pixels
    around every data pixel — the legend must have that much breathing
    room from any visual element, not just avoid pixel intersection.

    Parameters
    ----------
    clearance_px : int
        Pixels of keep-out zone around every data pixel. Default 4
        (≈1 pt at DPI=300).
    overlap_tolerance_px : int
        Number of intersecting pixels that are still tolerated. Default
        0 — any overlap fails. Callers who have reviewed a small overlap
        and decided it's acceptable can raise this per-figure.

    On failure, the error message reports:
      - How many pixels of overlap were detected
      - The axes-fraction bbox of the largest open region on this
        axis (from `find_open_regions`) as a concrete alternative slot
      - A concrete `bbox_to_anchor` suggestion derived from that slot
    """
    leg = ax.get_legend()
    if leg is None or not leg.get_visible():
        return

    mask, (mx0, my0, mx1, my1) = _compute_data_mask(fig, ax)
    if mask.size == 0:
        return

    dilated = _dilate_mask(mask, clearance_px)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    try:
        lb = leg.get_window_extent(renderer=renderer)
    except Exception:
        return
    if lb.width <= 0 or lb.height <= 0:
        return

    fig_h = int(round(fig.bbox.height))
    lx0 = max(0, int(math.floor(lb.x0 - mx0)))
    lx1 = min(mask.shape[1], int(math.ceil(lb.x1 - mx0)))
    ly1 = max(0, int(math.floor(fig_h - lb.y1 - my0)))
    ly0 = min(mask.shape[0], int(math.ceil(fig_h - lb.y0 - my0)))
    if lx1 <= lx0 or ly0 <= ly1:
        # Legend is entirely outside the axes region — nothing on the
        # axes data to conflict with.
        return

    sub = dilated[ly1:ly0, lx0:lx1]
    n_overlap = int(sub.sum())
    if n_overlap <= overlap_tolerance_px:
        return

    # Failure — build a helpful message with an alternative slot suggestion.
    free = find_open_regions(fig, ax)
    slot_priority = ["TR", "TL", "BR", "BL", "TC", "MR", "ML", "BC", "MC"]
    best_slot = None
    best_area = 0.0
    for name in slot_priority:
        rect = free.get(name)
        if rect is None:
            continue
        area = (rect[2] - rect[0]) * (rect[3] - rect[1])
        if area > best_area:
            best_area = area
            best_slot = (name, rect)

    if best_slot is not None:
        name, (fx0, fy0, fx1, fy1) = best_slot
        anchor_map = {
            "TR": ("upper right", fx1, fy1),
            "TL": ("upper left",  fx0, fy1),
            "BR": ("lower right", fx1, fy0),
            "BL": ("lower left",  fx0, fy0),
            "TC": ("upper center", (fx0 + fx1) / 2, fy1),
            "MR": ("right",       fx1, (fy0 + fy1) / 2),
            "ML": ("left",        fx0, (fy0 + fy1) / 2),
            "BC": ("lower center", (fx0 + fx1) / 2, fy0),
            "MC": ("center",      (fx0 + fx1) / 2, (fy0 + fy1) / 2),
        }
        loc, ax_x, ax_y = anchor_map[name]
        alt = (f"  Largest open slot: {name} "
               f"(axes-fraction bbox=({fx0:.2f}, {fy0:.2f}, {fx1:.2f}, {fy1:.2f})).\n"
               f"  Try: ax.legend(..., loc={loc!r}, "
               f"bbox_to_anchor=({ax_x:.2f}, {ax_y:.2f}), "
               f"bbox_transform=ax.transAxes)")
    else:
        alt = ("  No open region found on this axis. Move the legend "
               "OUTSIDE the axes (e.g., bbox_to_anchor=(1.02, 1.0), "
               "loc='upper left') or expand the panel to make room.")

    raise AssertionError(
        f"Legend overlaps data (or comes within {clearance_px} px of it): "
        f"{n_overlap} pixel(s) of clearance-zone intersection "
        f"(tolerance {overlap_tolerance_px}).\n"
        f"{alt}\n\n"
        f"To accept this overlap intentionally, pass "
        f"overlap_tolerance_px=<N> to assert_legend_clear "
        f"(or the corresponding kwarg to render_and_validate)."
    )

"""Geometry & balance helpers for matplotlib grant figures.

Two complementary tools:

  1. `horizontal_layout` / `vertical_layout` — derive named coordinates
     from a sequence of (name, width) items, so symmetric layouts are
     structural (not coincidental). Eliminates a class of bugs where
     the left arrow ends up a different length than the right because
     two free coordinates drifted apart.

  2. `geometry_summary` — given a layout you've already built, print a
     one-screen layout report and check declared symmetry properties
     (equal widths, equal heights, mirror pairs). Run before every
     render — a quick visual confirmation that "left mirrors right"
     and "siblings match" without re-rendering.

Use these in addition to (not instead of) the layout validator
(`validate_figure_layout.py`). The validator catches rendered overlap;
this catches asymmetry of intent.
"""
from __future__ import annotations

from collections import OrderedDict
from typing import Iterable


# ─────────────────────────────────────────────────────────────
# Derivation helpers — make symmetry structural
# ─────────────────────────────────────────────────────────────
def horizontal_layout(start, items):
    """Compute named horizontal coords from a sequence of (name, width) pairs.

    `items` is a list of (name, width) tuples. Names are arbitrary strings;
    widths are summed left-to-right starting from `start`.

    Returns (coords, end_x) where coords[name] = (x0, x1) for each named
    item. Use names like 'gap_*' for spacers — they appear in the dict but
    are usually only referenced when verifying.

    Example
    -------
    >>> coords, end_x = horizontal_layout(0.30, [
    ...     ('SAEC',    4.00),
    ...     ('gap_l1',  0.05),
    ...     ('arrow_l', 0.70),
    ...     ('gap_l2',  0.05),
    ...     ('INS',     4.20),
    ...     ('gap_r1',  0.05),
    ...     ('arrow_r', 0.70),
    ...     ('gap_r2',  0.05),
    ...     ('LTRC',    4.00),
    ... ])
    >>> coords['arrow_l']  # arrow_l_x0, arrow_l_x1
    (4.35, 5.05)
    >>> coords['arrow_r']
    (9.55, 10.25)
    >>> # arrows are guaranteed equal-length because they share the same width.
    """
    coords = OrderedDict()
    x = start
    for name, width in items:
        coords[name] = (x, x + width)
        x += width
    return coords, x


def vertical_layout(top, items):
    """Mirror of horizontal_layout for top-down vertical stacks.

    Items consume vertical space starting at `top` and going down. Returned
    coords are (y_top, y_bot) with y_top > y_bot (matplotlib convention).
    """
    coords = OrderedDict()
    y = top
    for name, height in items:
        coords[name] = (y, y - height)
        y -= height
    return coords, y


# ─────────────────────────────────────────────────────────────
# Summary + balance checks
# ─────────────────────────────────────────────────────────────
def _width(span):
    return abs(span[1] - span[0])


def geometry_summary(
    *,
    page=None,
    columns=None,
    rows=None,
    arrows=None,
    require_equal_widths=None,
    require_equal_heights=None,
    require_arrows_equal=False,
    require_symmetric_gaps=False,
    tol=0.01,
    print_summary=True,
):
    """Print a layout summary and return a list of balance violations.

    Parameters
    ----------
    page : (page_w, page_h) tuple, optional
        Figure dimensions, used only for the printed header.
    columns : list of (name, x0, x1) tuples
        Vertical bands (named columns).
    rows : list of (name, y0, y1) tuples
        Horizontal bands (named rows).
    arrows : list of (name, x0, x1) tuples
        Block arrows whose horizontal span you want listed and (optionally)
        verified equal.
    require_equal_widths : list of [name, name, ...] groups
        Each inner list names columns whose widths should match.
    require_equal_heights : list of [name, name, ...] groups
        Each inner list names rows whose heights should match.
    require_arrows_equal : bool
        If True, all arrows in `arrows` must have the same width.
    require_symmetric_gaps : bool
        If True, gaps between consecutive columns are checked for symmetry
        about the middle (gap[0] == gap[-1], gap[1] == gap[-2], ...).
    tol : float
        Equality tolerance for all checks.

    Returns
    -------
    list of str — human-readable violations (empty list if everything balances).
    """
    issues = []
    lines = []

    if page is not None:
        lines.append(f"Page: {page[0]:.2f} × {page[1]:.2f}")

    if columns:
        col_widths = [(n, _width((x0, x1))) for (n, x0, x1) in columns]
        col_str = "  |  ".join(f"{n} {w:.2f}" for n, w in col_widths)
        lines.append(f"Columns:  {col_str}")
        # Compute gaps between consecutive columns
        sorted_cols = sorted(columns, key=lambda c: c[1])
        gaps = []
        for i in range(len(sorted_cols) - 1):
            g = sorted_cols[i+1][1] - sorted_cols[i][2]
            gaps.append(g)
        if gaps:
            lines.append("Col gaps: " + "  ".join(f"{g:.2f}" for g in gaps))

    if arrows:
        arrow_widths = [(n, _width((x0, x1))) for (n, x0, x1) in arrows]
        a_str = "   ".join(f"{n}={w:.2f}" for n, w in arrow_widths)
        lines.append(f"Arrows:   {a_str}")

    if rows:
        row_heights = [(n, _width((y0, y1))) for (n, y0, y1) in rows]
        r_str = "  ".join(f"{n} {h:.2f}" for n, h in row_heights)
        lines.append(f"Rows:     {r_str}")

    # ── checks ──
    def _by_name(triples):
        return {n: (a, b) for (n, a, b) in (triples or [])}

    cols_by_name = _by_name(columns)
    rows_by_name = _by_name(rows)
    arrows_by_name = _by_name(arrows)

    for group in (require_equal_widths or []):
        widths = [(n, _width(cols_by_name[n])) for n in group if n in cols_by_name]
        if widths and max(w for _, w in widths) - min(w for _, w in widths) > tol:
            detail = ", ".join(f"{n}={w:.2f}" for n, w in widths)
            issues.append(f"unequal column widths in {group}: {detail}")

    for group in (require_equal_heights or []):
        heights = [(n, _width(rows_by_name[n])) for n in group if n in rows_by_name]
        if heights and max(h for _, h in heights) - min(h for _, h in heights) > tol:
            detail = ", ".join(f"{n}={h:.2f}" for n, h in heights)
            issues.append(f"unequal row heights in {group}: {detail}")

    if require_arrows_equal and arrows:
        widths = [_width(arrows_by_name[n]) for n in arrows_by_name]
        if max(widths) - min(widths) > tol:
            detail = ", ".join(f"{n}={_width(arrows_by_name[n]):.2f}"
                               for n in arrows_by_name)
            issues.append(f"unequal arrow widths: {detail}")

    if require_symmetric_gaps and columns:
        sorted_cols = sorted(columns, key=lambda c: c[1])
        gaps = [sorted_cols[i+1][1] - sorted_cols[i][2]
                for i in range(len(sorted_cols) - 1)]
        for i in range(len(gaps) // 2):
            j = len(gaps) - 1 - i
            if abs(gaps[i] - gaps[j]) > tol:
                issues.append(
                    f"gap asymmetry: gap[{i}]={gaps[i]:.2f} vs gap[{j}]={gaps[j]:.2f}")

    if print_summary:
        print("=" * 60)
        for line in lines:
            print(line)
        if issues:
            print("BALANCE ISSUES:")
            for iss in issues:
                print(f"  ✗ {iss}")
        else:
            print("Balance: all checks passed ✓")
        print("=" * 60)

    return issues

"""Reusable matplotlib primitives for grant / scientific figures.

Originated in the NMD 2026 overview figure. Functions are pure draw-on-axes
helpers — no globals, no figure ownership. Pass an `ax` plus geometry; each
function returns None and adds patches/text to the axes.

Categories:
  - Layout cards / panels:    aim_card, ds_box
  - Flow arrows:              block_arrow_right, block_arrow_left, thin_arrow
  - Bullet / timeline content: render_bullets, render_insights_timeline
  - Gene model:               exon_pill, intron_segment
  - Airway biology:           draw_ali_section, draw_cilia, smoke_puff
  - Misc:                     rbox

Companion conventions (typically defined per-figure, not here):
  HEADER_FS = 18         # one header size
  BODY_FS   = 14         # one body size; body text is never bold
  font.sans-serif        # Arial first; matplotlib sees both 400 and 700 faces
"""

import numpy as np
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, Polygon


# ── Default colors (override per-figure) ──
C_HEADER_DEFAULT = '#2C3E50'
C_SMOKE_DEFAULT  = '#7F8C8D'

# Shared layout constants for cards
HEADER_H_DEFAULT = 0.78
STRIPE_W_DEFAULT = 0.15


# ─────────────────────────────────────────────────────────────
# Generic boxes / arrows
# ─────────────────────────────────────────────────────────────
def rbox(ax, x, y, w, h, text, fc, ec='#555', fs=10, fw='normal',
         tc='white', alpha=1.0, lw=1.2, pad=0.06):
    """Rounded box with centered text."""
    box = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad={pad}",
                         facecolor=fc, edgecolor=ec, linewidth=lw,
                         alpha=alpha, zorder=2)
    ax.add_patch(box)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fs, fontweight=fw, color=tc, zorder=3)


def thin_arrow(ax, x1, y1, x2, y2, color='#444', lw=1.6, ms=14):
    """Thin annotation arrow from (x1,y1) to (x2,y2)."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color,
                                lw=lw, mutation_scale=ms), zorder=4)


def block_arrow_right(ax, xL, xR, yC, h=0.65, color=C_HEADER_DEFAULT,
                      alpha=1.0, ec='#222', lw=1.0, head_frac=0.40,
                      zorder=3):
    """Filled block arrow pointing right. xL/xR define horizontal span."""
    body_w = (xR - xL) * (1 - head_frac)
    half_h = h / 2
    head_half = h * 0.85
    body_xL, body_xR = xL, xL + body_w
    head_x_tip = xR
    verts = [
        (body_xL, yC - half_h), (body_xR, yC - half_h),
        (body_xR, yC - head_half), (head_x_tip, yC),
        (body_xR, yC + head_half), (body_xR, yC + half_h),
        (body_xL, yC + half_h),    (body_xL, yC - half_h),
    ]
    ax.add_patch(Polygon(verts, closed=True, facecolor=color, edgecolor=ec,
                         alpha=alpha, lw=lw, zorder=zorder))


def block_arrow_left(ax, xL, xR, yC, h=0.65, color=C_HEADER_DEFAULT,
                     alpha=1.0, ec='#222', lw=1.0, head_frac=0.40,
                     zorder=3):
    """Filled block arrow pointing left. Mirror of block_arrow_right."""
    body_w = (xR - xL) * (1 - head_frac)
    half_h = h / 2
    head_half = h * 0.85
    body_xR = xR
    body_xL = xR - body_w
    head_x_tip = xL
    verts = [
        (body_xR, yC - half_h),  (body_xL, yC - half_h),
        (body_xL, yC - head_half), (head_x_tip, yC),
        (body_xL, yC + head_half), (body_xL, yC + half_h),
        (body_xR, yC + half_h),    (body_xR, yC - half_h),
    ]
    ax.add_patch(Polygon(verts, closed=True, facecolor=color, edgecolor=ec,
                         alpha=alpha, lw=lw, zorder=zorder))


# ─────────────────────────────────────────────────────────────
# Card primitives (aim card, data-source box)
# ─────────────────────────────────────────────────────────────
def aim_card(ax, x, y, w, h, color, header_text=None, header_subtitle=None,
             header_align='left', fs_header=15, fs_sub=11.5,
             header_h=HEADER_H_DEFAULT, stripe_w=STRIPE_W_DEFAULT):
    """White-body card with colored top header band + colored left stripe.

    The stripe and header band share an edge (no overlap), so the layout
    validator's shape-shape overlap check stays clean.
    """
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle="round,pad=0.0,rounding_size=0.10",
                 facecolor='white', edgecolor=color, linewidth=1.6,
                 zorder=2))
    band_x = x + stripe_w
    band_w = w - stripe_w
    ax.add_patch(Rectangle((band_x, y + h - header_h), band_w, header_h,
                 facecolor=color, edgecolor='none', alpha=1.0, zorder=3))
    ax.add_patch(Rectangle((x, y), stripe_w, h,
                 facecolor=color, edgecolor='none', alpha=1.0, zorder=3))
    if header_text is not None:
        if header_align == 'center':
            tx = band_x + band_w/2
            ha = 'center'
        else:
            tx = band_x + 0.30
            ha = 'left'
        if header_subtitle:
            ax.text(tx, y + h - 0.22, header_text,
                    ha=ha, va='top', fontsize=fs_header,
                    fontweight='bold', color='white', zorder=5)
            ax.text(tx, y + h - 0.50, header_subtitle,
                    ha=ha, va='top', fontsize=fs_sub,
                    color='white', alpha=0.97, zorder=5)
        else:
            ax.text(tx, y + h - header_h/2, header_text,
                    ha=ha, va='center', fontsize=fs_header,
                    fontweight='bold', color='white', zorder=5)


def ds_box(ax, x, y, w, h, color, title, sub=None,
           fs_title=12, fs_sub=11, stripe_w=STRIPE_W_DEFAULT,
           text_color=C_HEADER_DEFAULT, sub_italic=True):
    """Data-source box: white body + colored left stripe.

    Title sits at the top in `color` (bold). Optional `sub` line sits at
    the bottom of the box, italic by default, in `text_color`.
    """
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle="round,pad=0.0,rounding_size=0.08",
                 facecolor='white', edgecolor=color, linewidth=1.6,
                 zorder=2))
    ax.add_patch(Rectangle((x, y), stripe_w, h,
                 facecolor=color, edgecolor='none', alpha=1.0, zorder=3))
    ax.text(x + stripe_w + (w - stripe_w)/2, y + h - 0.40, title,
            ha='center', va='center', fontsize=fs_title,
            fontweight='bold', color=color)
    if sub:
        ax.text(x + stripe_w + (w - stripe_w)/2, y + 0.40, sub,
                ha='center', va='center', fontsize=fs_sub,
                color=text_color,
                style='italic' if sub_italic else 'normal')


# ─────────────────────────────────────────────────────────────
# Bullet / timeline content
# ─────────────────────────────────────────────────────────────
def render_bullets(ax, x0, y_top, bullets, color, line_h=0.42, fs=12,
                   label_dx=0.45, body_color=C_HEADER_DEFAULT):
    """Render a list of (label, text) bullets from y_top down.
    Label drawn bold in `color`; text drawn regular in `body_color`.
    Returns the y-position after the last bullet.
    """
    cur_y = y_top
    for lbl, txt in bullets:
        ax.text(x0, cur_y, lbl, ha='left', va='top',
                fontsize=fs, fontweight='bold', color=color)
        ax.text(x0 + label_dx, cur_y, txt, ha='left', va='top',
                fontsize=fs, color=body_color)
        n_lines = txt.count('\n') + 1
        cur_y -= line_h * n_lines + 0.10
    return cur_y


def render_insights_timeline(ax, x_dot, y_centers, bullets, color,
                              dot_r=0.22, fs=13, text_dx=0.45,
                              body_color=C_HEADER_DEFAULT):
    """Numbered timeline: colored dots connected by a vertical line, with
    bullet text floating to the right. No enclosing box.

    `y_centers` is a list of y positions for each bullet's dot center; the
    vertical connector runs from the first to the last dot.
    `bullets` is a list of (label, text) tuples — label goes inside the dot,
    text is rendered to the right.
    """
    for i in range(len(y_centers) - 1):
        y_a = y_centers[i]   - dot_r - 0.02
        y_b = y_centers[i+1] + dot_r + 0.02
        ax.plot([x_dot, x_dot], [y_a, y_b],
                '-', color=color, lw=2.5, zorder=2, alpha=0.55,
                solid_capstyle='round')
    for (yc, (lbl, txt)) in zip(y_centers, bullets):
        ax.add_patch(Circle((x_dot, yc), dot_r,
                            facecolor=color, edgecolor='white', lw=2.0,
                            zorder=3))
        ax.text(x_dot, yc, lbl.rstrip('.'),
                ha='center', va='center',
                fontsize=fs - 2, fontweight='bold', color='white',
                zorder=4)
        ax.text(x_dot + text_dx, yc, txt,
                ha='left', va='center',
                fontsize=fs, color=body_color, zorder=3)


# ─────────────────────────────────────────────────────────────
# Gene-model primitives
# ─────────────────────────────────────────────────────────────
def exon_pill(ax, x, y, w, h, color, ec='#333', label=None, fs=9, tc='white'):
    """Rounded exon pill. Optional centered label."""
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.16",
                         facecolor=color, edgecolor=ec, linewidth=1.0, zorder=3)
    ax.add_patch(box)
    if label:
        ax.text(x + w/2, y + h/2, label, ha='center', va='center',
                fontsize=fs, fontweight='bold', color=tc, zorder=4)


def intron_segment(ax, x1, x2, y, color='#888', lw=1.0):
    """Thin intron line connecting two exons. No-op if span is too small."""
    if x2 - x1 > 0.05:
        ax.plot([x1, x2], [y, y], '-', color=color, lw=lw, zorder=2)


# ─────────────────────────────────────────────────────────────
# Airway biology primitives
# ─────────────────────────────────────────────────────────────
def smoke_puff(ax, cx, cy, scale=1.0, color=C_SMOKE_DEFAULT, alpha=0.55):
    """Cluster of overlapping circles forming a smoke puff."""
    offsets = [(-0.18, 0.0), (0.0, 0.05), (0.18, 0.0),
               (-0.09, 0.10), (0.09, 0.10), (0.0, -0.06)]
    for dx, dy in offsets:
        ax.add_patch(Circle((cx + dx*scale, cy + dy*scale),
                            0.11*scale, facecolor=color, edgecolor='none',
                            alpha=alpha, zorder=3))


def draw_cilia(ax, x_left, x_right, y_top, n=14, h=0.18, color=C_HEADER_DEFAULT):
    """Wavy cilia tufts spread between x_left and x_right, rooted at y_top."""
    xs = np.linspace(x_left + 0.04, x_right - 0.04, n)
    for x in xs:
        t = np.linspace(0, 1, 10)
        wy = y_top + h * t
        wx = x + 0.04 * np.sin(t * np.pi * 2.0)
        ax.plot(wx, wy, '-', color=color, lw=0.9, zorder=4, alpha=0.85)


# ── ALI cross-section: cell-type palette ──
ALI_BASAL_C    = '#4A90E2'   # blue, dome-shaped basal cells
ALI_CLUB_C     = '#E8C8A8'   # tan/peach club cells
ALI_CILIATED_C = '#F2DEB6'   # light cream ciliated cells
ALI_GOBLET_C   = '#7FCDBC'   # teal goblet cells
ALI_NUC_NAVY   = '#0E2B45'   # nucleus for basal
ALI_NUC_BROWN  = '#5D4037'   # nucleus for club / ciliated
ALI_NUC_GREEN  = '#1F5641'   # nucleus / mucus granules for goblet
ALI_CILIA_C    = '#3E3530'   # cilia hairs (dark grey-brown)
ALI_MEDIA_C    = '#FCE0E5'   # pink culture media
ALI_WALL_C     = '#CFD3D7'   # vessel walls


def draw_ali_section(ax, x, y, w, h, stage='full'):
    """Cross-section of an HBEC ALI culture (pseudostratified airway epithelium).

    Anatomy:
      - Vessel walls (light grey) on the left and right.
      - Pink basal media in the lower chamber.
      - Porous membrane spanning the chamber at ~20% of vessel height.
      - Cells anchored on the membrane (every cell connects — no multilayer).

    Stages:
      stage='basal' : primary HBECs, cuboidal cells submerged in apical media.
      stage='inter' : partially differentiated; basal + early club; no cilia,
                      no goblets, no apical media. Heights closer together.
      stage='full'  : mature pseudostratified epithelium (D35) — all four
                      cell types: basal (short dome), club, ciliated (cilia),
                      goblet (mucus granules).
    """
    wall_w = 0.04
    ax.add_patch(Rectangle((x - wall_w, y), wall_w, h,
                facecolor=ALI_WALL_C, edgecolor=None, zorder=2))
    ax.add_patch(Rectangle((x + w, y), wall_w, h,
                facecolor=ALI_WALL_C, edgecolor=None, zorder=2))

    mem_y     = y + h * 0.20
    cell_zone = h * 0.78

    ax.add_patch(Rectangle((x, y), w, mem_y - y,
                facecolor=ALI_MEDIA_C, edgecolor=None, zorder=1))
    ax.plot([x, x + w], [mem_y, mem_y], '-',
            color='#5D4037', lw=1.0, zorder=3)
    for tx in np.linspace(x + 0.05, x + w - 0.05, 7):
        ax.plot([tx, tx], [mem_y - 0.02, mem_y + 0.02], '-',
                color='#5D4037', lw=0.5, zorder=3)

    n = max(7, int(round(w * 7.5)))
    cw = w / n
    cell_w = cw * 0.86

    def _draw_cell(cx, ch, color, nuc_color, top_feat=None):
        top_y = mem_y + ch
        if color == ALI_BASAL_C:
            shoulder = max(cell_w * 0.18, 0.02)
            poly = [
                (cx - cell_w/2, mem_y),
                (cx - cell_w/2, top_y - shoulder),
                (cx,            top_y),
                (cx + cell_w/2, top_y - shoulder),
                (cx + cell_w/2, mem_y),
            ]
            ax.add_patch(Polygon(poly, closed=True,
                        facecolor=color, edgecolor=ALI_NUC_BROWN,
                        lw=0.4, zorder=4))
            nuc_y = mem_y + ch * 0.45
        else:
            ax.add_patch(Rectangle((cx - cell_w/2, mem_y), cell_w, ch,
                        facecolor=color, edgecolor=ALI_NUC_BROWN,
                        lw=0.4, zorder=4))
            nuc_y = mem_y + ch * (0.30 if ch > cell_zone * 0.65 else 0.45)
        ax.add_patch(Circle((cx, nuc_y), min(cell_w * 0.22, 0.030),
                    facecolor=nuc_color, alpha=0.80,
                    edgecolor=None, zorder=5))
        if top_feat == 'cilia':
            n_ci = 5
            ci_h = cell_zone * 0.10
            for ci_i in range(n_ci):
                ci_x = cx + (ci_i / (n_ci - 1) - 0.5) * cell_w * 0.80
                ax.plot([ci_x, ci_x], [top_y, top_y + ci_h],
                        '-', color=ALI_CILIA_C, lw=0.7,
                        zorder=6, alpha=0.9)
        elif top_feat == 'mucus':
            for gx_off, gy_off in [
                (-0.22, -0.10), (0.22, -0.10),
                (0.0,  -0.22), (-0.18, -0.32), (0.18, -0.32),
            ]:
                ax.add_patch(Circle(
                    (cx + gx_off * cell_w, top_y + gy_off * cell_zone * 0.5),
                    cell_w * 0.13,
                    facecolor=ALI_NUC_GREEN, edgecolor=None,
                    alpha=0.55, zorder=5))

    if stage == 'basal':
        ax.add_patch(Rectangle((x, mem_y), w, h - (mem_y - y),
                    facecolor=ALI_MEDIA_C, edgecolor=None,
                    alpha=0.55, zorder=1))
        cell_h = cell_zone * 0.24
        for i in range(n):
            cx = x + (i + 0.5) * cw
            _draw_cell(cx, cell_h, ALI_BASAL_C, ALI_NUC_NAVY)

    elif stage == 'inter':
        pattern = [
            (0.30, ALI_BASAL_C, ALI_NUC_NAVY,  None),
            (0.48, ALI_CLUB_C,  ALI_NUC_BROWN, None),
            (0.32, ALI_BASAL_C, ALI_NUC_NAVY,  None),
            (0.52, ALI_CLUB_C,  ALI_NUC_BROWN, None),
            (0.34, ALI_BASAL_C, ALI_NUC_NAVY,  None),
            (0.49, ALI_CLUB_C,  ALI_NUC_BROWN, None),
            (0.30, ALI_BASAL_C, ALI_NUC_NAVY,  None),
            (0.53, ALI_CLUB_C,  ALI_NUC_BROWN, None),
            (0.34, ALI_BASAL_C, ALI_NUC_NAVY,  None),
        ]
        for i in range(n):
            cx = x + (i + 0.5) * cw
            hf, color, nuc_color, top_feat = pattern[i % len(pattern)]
            _draw_cell(cx, cell_zone * hf, color, nuc_color, top_feat)

    else:  # 'full'
        pattern = [
            (0.61, ALI_CILIATED_C, ALI_NUC_BROWN, 'cilia'),
            (0.32, ALI_BASAL_C,    ALI_NUC_NAVY,  None),
            (0.57, ALI_GOBLET_C,   ALI_NUC_GREEN, 'mucus'),
            (0.61, ALI_CILIATED_C, ALI_NUC_BROWN, 'cilia'),
            (0.53, ALI_CLUB_C,     ALI_NUC_BROWN, None),
            (0.32, ALI_BASAL_C,    ALI_NUC_NAVY,  None),
            (0.61, ALI_CILIATED_C, ALI_NUC_BROWN, 'cilia'),
            (0.57, ALI_GOBLET_C,   ALI_NUC_GREEN, 'mucus'),
            (0.61, ALI_CILIATED_C, ALI_NUC_BROWN, 'cilia'),
        ]
        for i in range(n):
            cx = x + (i + 0.5) * cw
            hf, color, nuc_color, top_feat = pattern[i % len(pattern)]
            _draw_cell(cx, cell_zone * hf, color, nuc_color, top_feat)

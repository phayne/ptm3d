"""
Shared matplotlib style for PTM3D publication figures (docs/figures/v2/).

Built on top of `pretty_plots` (https://github.com/phayne/pretty-plots) to get
IDL-style defaults — bold axis lines, full box-frame, inward ticks, IDL Hershey
serif, IDL color tables — with a few PTM3D-specific color constants and
helpers on top.

Use:

    import matplotlib
    matplotlib.use('Agg')
    from ptm.pubstyle import apply_style, COLORS, COLUMN_W, DOUBLE_W

    apply_style()                     # IDL family
    fig, ax = plt.subplots(figsize=(DOUBLE_W, 0.75 * DOUBLE_W))
    ...
    fig.savefig(out_path)
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

try:
    from pretty_plots import idl as _pp_idl
    HAS_PRETTY = True
except Exception:
    _pp_idl = None
    HAS_PRETTY = False

try:
    from pretty_plots import planetary as _pp_planetary
except Exception:  # pretty_plots is optional: use() then sets the same base
    _pp_planetary = None

# AGU / Icarus single- and double-column widths (inches)
COLUMN_W = 3.5
DOUBLE_W = 7.2

# PTM3D color palette — kept compatible with the IDL `colortables` so curves
# read well against IDL's white background and full-box frame
COLORS = {
    "ptm":      "#1f77b4",   # IDL default blue
    "ref":      "#d62728",   # IDL default red
    "neutral":  "#444444",
    "lat0":     "#1f77b4",
    "lat30":    "#2ca02c",
    "lat45":    "#ff7f0e",
    "lat60":    "#9467bd",
    "lat85":    "#d62728",
    "shadow":   "#cccccc",
    "envelope": "#cfd8e8",
    "good":     "#2ca02c",
    "bad":      "#d62728",
}

LATITUDE_COLORS = {
    0:  COLORS["lat0"],
    30: COLORS["lat30"],
    45: COLORS["lat45"],
    60: COLORS["lat60"],
    85: COLORS["lat85"],
}


def apply_style(variant: str | None = None) -> None:
    """Apply the IDL family of pretty_plots, then a few PTM3D refinements.

    `variant` may be one of pretty_plots' IDL variants: aas / aa / nature /
    icarus / hershey / latex. Defaults to the bare IDL family.
    """
    if HAS_PRETTY:
        _pp_idl.use(variant=variant)
    else:
        # Minimal IDL-ish fallback if pretty_plots isn't installed.
        mpl.rcParams.update({
            "font.family": "sans-serif",
            "lines.linewidth": 1.8,
            "axes.linewidth": 1.4,
            "axes.spines.top":   True,
            "axes.spines.right": True,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top":   True,
            "ytick.right": True,
            "savefig.dpi": 300,
        })

    # PTM3D refinements on top of IDL — large fonts for clean publication plots
    mpl.rcParams.update({
        "font.size": 13,
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "axes.titleweight": "bold",
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11,
        "legend.frameon": False,
        "legend.handlelength": 2.0,
        "lines.linewidth": 2.0,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.06,
        "figure.dpi": 110,
    })


# ---------------------------------------------------------------------------
# Figure style of the tutorials and the documentation
# ---------------------------------------------------------------------------
# One call, use(), gives every figure in docs/ the same look. Its base is the
# `planetary` family of pretty_plots (Hayne 2017; MIT license): serif text and
# 2-point lines. The font sizes are scaled from that family so that a figure
# drawn at COLUMN_W or DOUBLE_W carries 10 to 12 point text at its printed
# size. When pretty_plots is not installed the same base values are set
# directly, so the figure looks the same; pretty_plots is optional.

#: Okabe and Ito (2008) colors, distinguishable with the common forms of color
#: blindness. Yellow is left out because it does not show on white.
BLACK, BLUE, VERMILLION, GREEN, ORANGE, SKY_BLUE, PURPLE = (
    "#000000", "#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7")

#: The order in which curves of one kind take colors and line styles, so a
#: figure also reads in grayscale: latitudes from the equator to the pole,
#: thermal inertias from low to high, and so on.
PALETTE = (BLACK, BLUE, VERMILLION, GREEN, ORANGE, PURPLE)
LINESTYLES = ("-", "--", "-.", ":", (0, (7, 2)), (0, (3, 1, 1, 1, 1, 1)))
MARKERS = ("o", "s", "^", "D", "v", "P")

#: Quantities drawn in more than one tutorial always take these styles.
STYLES = {
    "radiation_on":  dict(color=VERMILLION, linestyle="-"),
    "radiation_off": dict(color=BLUE, linestyle="--"),
    "brightness":    dict(color=PURPLE, linestyle="-"),    # band brightness temperature
    "surface":       dict(color=GREEN, linestyle="--"),    # surface temperature
    "model":         dict(color=BLACK, linestyle="-"),
    "fit":           dict(color=ORANGE, linestyle="-"),
    "coarse":        dict(color="0.35", linestyle=":"),    # coarser time steps
    "observed":      dict(color=BLACK, marker="o", markerfacecolor="white",
                          linestyle="none"),
}

#: Shaded time intervals, gray so they never compete with the data.
SHADES = {
    "occultation": dict(facecolor="0.78", edgecolor="none", zorder=0.2),
    "eclipse":     dict(facecolor="0.91", edgecolor="none", zorder=0.1),
    "transit":     dict(facecolor="none", edgecolor="0.55", hatch="////",
                        linewidth=0.0, zorder=0.1),
}

#: Color maps: temperature, any other quantity of one sign, and a signed
#: difference, always drawn about zero (use a norm symmetric about zero).
CMAP_TEMPERATURE = "inferno"
CMAP_SCALAR = "viridis"
CMAP_DIFFERENCE = "RdBu_r"

#: Resolution of saved raster figures, in dots per inch at the designed size.
SAVE_DPI = 300

_BASE = {  # the pretty_plots `planetary` base, used when it is not installed
    "font.family": "serif",
    "lines.linewidth": 2,
}

_SCALED = {
    "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
    "mathtext.fontset": "dejavuserif",      # math in the same face as the text
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
    "lines.markersize": 6,
    "axes.linewidth": 1.0,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "xtick.major.size": 5, "ytick.major.size": 5,
    "xtick.minor.size": 2.5, "ytick.minor.size": 2.5,
    "xtick.major.width": 1.0, "ytick.major.width": 1.0,
    "xtick.minor.width": 0.8, "ytick.minor.width": 0.8,
    "legend.frameon": False,
    "image.cmap": CMAP_TEMPERATURE,
    "savefig.dpi": SAVE_DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
    "pdf.fonttype": 42, "ps.fonttype": 42,  # text stays text in PDF and SVG
    "svg.fonttype": "none",
}


def use() -> None:
    """Apply the figure style of the PTM3D documentation."""
    from cycler import cycler
    if _pp_planetary is not None:
        _pp_planetary.use()
    else:
        mpl.rcParams.update(_BASE)
    mpl.rcParams.update(_SCALED)
    mpl.rcParams["axes.prop_cycle"] = (cycler(color=list(PALETTE))
                                       + cycler(linestyle=list(LINESTYLES)))


def panel_labels(axes, start="a", titles=None):
    """Label panels (a), (b), ... above their top-left corners, clear of the data.

    ``titles``, one per panel, are written after the label on the same line.
    """
    import numpy as np
    axes = np.ravel(axes)
    titles = list(titles) if titles is not None else [""] * len(axes)
    for i, (ax, title) in enumerate(zip(axes, titles)):
        text = f"({chr(ord(start) + i)})" + (f"  {title}" if title else "")
        ax.text(0.0, 1.015, text, transform=ax.transAxes,
                ha="left", va="bottom", fontsize=mpl.rcParams["axes.labelsize"])


def save(fig, path, dpi=SAVE_DPI):
    """Save a figure; a PNG is written at ``dpi`` and losslessly recompressed.

    Line plots use the default 300 dpi. Shaded maps are saved at 200 dpi, which
    keeps them sharp on a screen at a third of the file size.
    """
    from pathlib import Path
    path = Path(path)
    fig.savefig(path, dpi=dpi, metadata=_metadata(path.suffix))
    if path.suffix.lower() == ".png":
        try:
            from PIL import Image
        except ImportError:
            return path
        with Image.open(path) as im:
            im = im.convert("RGB")
            im.save(path, optimize=True, dpi=(dpi, dpi))
    return path


def _metadata(suffix):
    """No software or creation-date stamps, so a regenerated file differs only if
    the figure does."""
    suffix = suffix.lower()
    if suffix == ".png":
        return {"Software": None}
    if suffix == ".pdf":
        return {"Creator": None, "Producer": None, "CreationDate": None}
    if suffix == ".svg":
        return {"Creator": None, "Date": None}
    return None


def panel_label(ax, label, x=0.02, y=0.97, **kw):
    """Place a panel label like '(a)' in the top-left corner."""
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=10, fontweight="bold",
        ha="left", va="top",
        **kw,
    )


def savefig_pinned(fig, out_stem, formats=(".png", ".pdf"), canvas_px=None, dpi=300):
    """Save `fig` to `out_stem` + each of `formats`.

    If `canvas_px=(w, h)` is given, pin the saved canvas to that EXACT pixel
    size at `dpi` (matching a reference/annotated sibling 1:1) instead of the
    default `savefig.bbox="tight"` rcParam set in `apply_style()`. That rcParam
    crops to rendered-CONTENT extent, so removing a title/stat-box SHRINKS the
    crop even for an otherwise pixel-identical figure -- `bbox_inches=None`
    does NOT disable it; an explicit full-canvas `Bbox` is required. This is
    the pattern validated on the Diviner seasonal comparison /
    the paper scatter figure's `--pubclean-canvas-px`
    (manuscript-v4 review, 2026-07-22): pin every publication-clean variant to
    its annotated sibling's own shipped pixel size so the manuscript swap is
    exact-dims, not just visually close.

    Returns the list of Paths written.
    """
    from pathlib import Path
    out_stem = Path(out_stem)
    written = []
    if canvas_px is None:
        for suf in formats:
            p = out_stem.with_suffix(suf)
            fig.savefig(p, dpi=dpi)
            written.append(p)
        return written
    from matplotlib.transforms import Bbox
    w_px, h_px = canvas_px
    w_in, h_in = w_px / dpi, h_px / dpi
    # Center a fixed w_in x h_in window on the figure's ACTUAL rendered content
    # extent (fig.get_tightbbox(), in inches) rather than anchoring at (0,0) of
    # a resized canvas -- a plain (0,0)-origin Bbox silently clips content that
    # sits left/below the figure's nominal origin once tight_layout has been
    # called with a different (e.g. title-reserving) rect than the pinned
    # target aspect implies (found 2026-07-22 on the indirect-bands map figure:
    # the left panel title got clipped). This is the same fix already applied
    # to the Diviner seasonal comparison's --overview mode for a similar gotcha.
    fig.canvas.draw()
    tb = fig.get_tightbbox(fig.canvas.get_renderer())
    cx, cy = (tb.x0 + tb.x1) / 2.0, (tb.y0 + tb.y1) / 2.0
    bbox = Bbox.from_bounds(cx - w_in / 2.0, cy - h_in / 2.0, w_in, h_in)
    for suf in formats:
        p = out_stem.with_suffix(suf)
        fig.savefig(p, dpi=(dpi if suf == ".png" else None), bbox_inches=bbox)
        written.append(p)
    return written


def load_face_mask(path):
    """Load a boolean per-face mask (uint32 n, then n uint8) written by
    make_psr_mask.py — True = spurious-lit face to reclassify as PSR."""
    import numpy as np
    with open(path, "rb") as fh:
        n = int(np.fromfile(fh, np.uint32, 1)[0])
        return np.fromfile(fh, np.uint8, n).astype(bool)

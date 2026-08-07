"""
GILDAS/GREG colour lookup tables as matplotlib colormaps.

Registers:

- ``rainbow3`` (+ ``rainbow3_r``): the GILDAS ``rainbow3`` LUT verbatim - the
  most popular GREG/CLASS colour table (``lut rainbow3``), 128 entries running
  black→blue→cyan→green→yellow→orange→red→magenta→white.
- ``rainbow3_br`` (+ ``rainbow3_br_r``): a truncated blue→red slice of the same
  LUT - the leading pure-black entries and the trailing magenta/pink/white tail
  are dropped, keeping dark-blue→cyan→green→yellow→orange→red. This matches the
  common "blue background, red peaks" look, because the full LUT's default
  min/max scaling otherwise sends noise to black and peaks to magenta/white.

All available anywhere a colormap name is accepted, e.g. ``--colormap rainbow3_br``.
The RGB values are copied verbatim from GILDAS' ``rainbow3.lut`` and embedded here
so the colormaps work without a GILDAS installation.
"""

from matplotlib.colors import LinearSegmentedColormap

# 128 RGB triplets (0..1), verbatim from GILDAS kernel/etc/lut-desc/rainbow3.lut
_RAINBOW3_RGB = [
    (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0.0784),
    (0, 0, 0.1569), (0, 0, 0.2392), (0, 0, 0.3176), (0, 0, 0.4),
    (0, 0, 0.4784), (0, 0, 0.5569), (0, 0, 0.6392), (0, 0, 0.7176),
    (0, 0, 0.8), (0, 0, 0.8784), (0, 0, 0.9569), (0, 0.0314, 1),
    (0, 0.0941, 1), (0, 0.1569, 1), (0, 0.2196, 1), (0, 0.2863, 1),
    (0, 0.349, 1), (0, 0.4118, 1), (0, 0.4784, 1), (0, 0.5176, 1),
    (0, 0.5569, 1), (0, 0.6, 1), (0, 0.6392, 1), (0, 0.6824, 1),
    (0, 0.7216, 1), (0, 0.7608, 1), (0, 0.7961, 1), (0, 0.8275, 1),
    (0, 0.8549, 1), (0, 0.8824, 1), (0, 0.9137, 1), (0, 0.9412, 1),
    (0, 0.9686, 1), (0, 1, 1), (0, 1, 0.9608), (0, 1, 0.9216),
    (0, 1, 0.8824), (0, 1, 0.8431), (0, 1, 0.8039), (0, 1, 0.7647),
    (0, 1, 0.7255), (0, 1, 0.6588), (0, 1, 0.5647), (0, 1, 0.4706),
    (0, 1, 0.3765), (0, 1, 0.2784), (0, 1, 0.1843), (0, 1, 0.0902),
    (0, 1, 0), (0.0941, 1, 0), (0.1882, 1, 0), (0.2824, 1, 0),
    (0.3765, 1, 0), (0.4706, 1, 0), (0.5647, 1, 0), (0.6588, 1, 0),
    (0.7255, 1, 0), (0.7647, 1, 0), (0.8039, 1, 0), (0.8431, 1, 0),
    (0.8824, 1, 0), (0.9216, 1, 0), (0.9608, 1, 0), (1, 1, 0),
    (0.9961, 0.9608, 0), (0.9961, 0.9255, 0), (0.9922, 0.8902, 0), (0.9922, 0.851, 0),
    (0.9922, 0.8157, 0), (0.9882, 0.7804, 0), (0.9882, 0.7412, 0), (0.9882, 0.7059, 0),
    (0.9882, 0.6745, 0), (0.9922, 0.6392, 0), (0.9922, 0.6039, 0), (0.9922, 0.5725, 0),
    (0.9961, 0.5373, 0), (0.9961, 0.502, 0), (1, 0.4706, 0), (1, 0.4078, 0),
    (1, 0.3451, 0), (1, 0.2824, 0), (1, 0.2157, 0), (1, 0.1529, 0),
    (1, 0.0902, 0), (1, 0.0275, 0), (1, 0, 0), (1, 0, 0),
    (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0),
    (1, 0, 0), (1, 0, 0), (1, 0, 0.0941), (1, 0, 0.1882),
    (1, 0, 0.2824), (1, 0, 0.3765), (1, 0, 0.4706), (1, 0, 0.5647),
    (1, 0, 0.6588), (1, 0, 0.7255), (1, 0, 0.7725), (1, 0, 0.8157),
    (1, 0, 0.8627), (1, 0, 0.9059), (1, 0, 0.9529), (1, 0, 1),
    (1, 0.0706, 1), (1, 0.1412, 1), (1, 0.2157, 1), (1, 0.2863, 1),
    (1, 0.3608, 1), (1, 0.4314, 1), (1, 0.4863, 1), (1, 0.5216, 1),
    (1, 0.5608, 1), (1, 0.5961, 1), (1, 0.6314, 1), (1, 0.6667, 1),
    (1, 0.7059, 1), (1, 0.7882, 1), (1, 0.8706, 1), (1, 0.9569, 1),
]


def _truncate_blue_to_red(rgb):
    """Slice the LUT to dark-blue → ... → first pure red.

    Drops the leading pure-black entries and everything from the first pure red
    onward (the red plateau tail plus the magenta/pink/white top), so the
    resulting colormap spans only the blue→red portion.
    """
    start = 0
    while start < len(rgb) and tuple(rgb[start]) == (0, 0, 0):
        start += 1
    end = len(rgb) - 1
    for i in range(start, len(rgb)):
        if tuple(rgb[i]) == (1, 0, 0):
            end = i
            break
    return rgb[start:end + 1]


def register_gildas_luts():
    """Register GILDAS LUTs as matplotlib colormaps (idempotent).

    Makes ``rainbow3``/``rainbow3_r`` (full LUT) and ``rainbow3_br``/
    ``rainbow3_br_r`` (blue→red truncation) available via ``plt.colormaps[...]``
    and hence the ``--colormap`` option. Safe to call multiple times.
    """
    import matplotlib

    # matplotlib >=3.5: matplotlib.colormaps; older: matplotlib.cm.register_cmap
    try:
        registry = matplotlib.colormaps
        have_registry = True
    except AttributeError:  # very old matplotlib
        registry = None
        have_registry = False

    full = LinearSegmentedColormap.from_list('rainbow3', _RAINBOW3_RGB, N=256)
    br = LinearSegmentedColormap.from_list(
        'rainbow3_br', _truncate_blue_to_red(_RAINBOW3_RGB), N=256)

    for base in (full, br):
        for cm in (base, base.reversed()):
            if have_registry:
                if cm.name not in registry:
                    registry.register(cm)
            else:  # pragma: no cover - legacy path
                import matplotlib.cm as mcm
                try:
                    mcm.get_cmap(cm.name)
                except ValueError:
                    mcm.register_cmap(name=cm.name, cmap=cm)

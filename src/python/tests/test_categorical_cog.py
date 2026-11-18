#!/usr/bin/env python3
"""Class and layer rasters hold only class and layer values (PTM-983).

``ice_retention_class`` and ``ice_table_layer`` are labels. The COG exporter
mapped them like a temperature: barycentric (or inverse-distance) face ->
pixel interpolation and AVERAGE overviews. Measured on a polar product, 1.9%
of the class pixels and 2.9% of the layer pixels were values such as 0.242 or
1.0002, which are no class at all. A categorical field (``categorical`` in
``ptm.results.STATISTICS_FIELDS``) now takes the nearest face, with the same
footprint as the continuous fields, and MODE overviews.

1. ``MeshRegridder.regrid_nearest`` returns only face values, on exactly the
   footprint of ``regrid_linear`` (and of ``regrid`` for ``idw``); the same
   field through ``regrid_linear`` is NOT all labels, so the check can fail.
2. The exported class and layer COGs hold only valid values at the full
   resolution and at every overview level, and say how they were resampled.
3. A continuous field is exported exactly as before.

Arms 2 and 3 need GDAL and rasterio and skip without them. The fixture is
synthetic, at invented coordinates; no site appears here.

    python3 src/python/tests/test_categorical_cog.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))

import numpy as np                                              # noqa: E402

CLASSES = {0.0, 1.0, 2.0, 3.0}
LAYERS = set(float(x) for x in range(-1, 23))

FAILURES, SKIPS = [], []


class Skip(Exception):
    pass


def check(name, fn, *args):
    try:
        detail = fn(*args)
    except Skip as exc:
        SKIPS.append(name)
        print(f"SKIP {name} :: {exc}")
    except AssertionError as exc:
        FAILURES.append(name)
        print(f"FAIL {name} :: {exc}")
    except Exception as exc:                                    # noqa: BLE001
        FAILURES.append(name)
        print(f"FAIL {name} :: {type(exc).__name__}: {exc}")
    else:
        print(f"PASS {name} :: {detail or ''}")


def _fixture():
    """Face centers on a jittered polar patch (invented coordinates) with
    class and layer fields in blocks, so that many pixels sit between faces
    of different labels."""
    try:
        from ptm.results import SCIPY_AVAILABLE, MeshRegridder
    except Exception as exc:                                    # noqa: BLE001
        raise Skip(f"ptm.results not importable ({exc})")
    if not SCIPY_AVAILABLE:
        raise Skip("scipy is not installed")
    rng = np.random.default_rng(7)
    n = 40
    lat0, lon0 = -80.0, 20.0
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    lat = lat0 + 0.002 * ii.ravel() + 0.0006 * rng.standard_normal(n * n)
    lon = lon0 + 0.012 * jj.ravel() + 0.003 * rng.standard_normal(n * n)
    cls = (((ii // 5) + (jj // 7)) % 4).ravel().astype(np.float32)
    layer = (((ii // 3) * 3 + (jj // 4)) % 24 - 1).ravel().astype(np.float32)
    depth = (0.01 + 0.002 * ii + 0.001 * jj).ravel().astype(np.float32)
    reg = MeshRegridder(lat, lon, resolution=0.0006)
    return reg, cls, layer, depth


def _labels_only(grid, allowed):
    v = grid[np.isfinite(grid)]
    bad = v[~np.isin(v, list(allowed))]
    return v.size, bad


def check_nearest_regrid_gives_labels_on_the_same_footprint():
    reg, cls, layer, _ = _fixture()
    for field, vals, allowed in (("class", cls, CLASSES), ("layer", layer, LAYERS)):
        lin = reg.regrid_linear(vals)
        near = reg.regrid_nearest(vals, footprint="linear")
        n, bad = _labels_only(near, allowed)
        assert n > 0 and bad.size == 0, f"{field}: {bad.size} non-label pixels of {n}"
        assert np.array_equal(np.isfinite(near), np.isfinite(lin)), f"{field}: footprint differs"
        _, bad_lin = _labels_only(lin, allowed)
        assert bad_lin.size > 0, f"{field}: interpolation gave only labels; the fixture cannot fail"
        idw = reg.regrid(vals)
        near_i = reg.regrid_nearest(vals, footprint="idw")
        assert np.array_equal(np.isfinite(near_i), np.isfinite(idw)), f"{field}: idw footprint differs"
        assert _labels_only(near_i, allowed)[1].size == 0
    return (f"{n} pixels, all labels; interpolation would give "
            f"{bad_lin.size} non-label pixels on the same grid")


def _raster():
    try:
        import ptm_cog
        import rasterio                                         # noqa: F401
    except Exception as exc:                                    # noqa: BLE001
        raise Skip(f"raster stack not importable ({exc})")
    if not ptm_cog.GDAL_AVAILABLE:
        raise Skip("osgeo.gdal is not importable")
    try:
        ptm_cog._find_gdal_tool("gdaladdo")
    except RuntimeError as exc:
        raise Skip(str(exc))
    return ptm_cog


def _all_levels(path):
    import rasterio
    with rasterio.open(path) as ds:
        nd, out = ds.nodata, []
        for i, f in enumerate([1] + ds.overviews(1)):
            a = ds.read(1, out_shape=(-(-ds.height // f), -(-ds.width // f))) if i else ds.read(1)
            out.append(a[(a != nd) & np.isfinite(a)])
        return out, ds.tags(1), len(ds.overviews(1))


def check_class_and_layer_cogs_hold_only_labels_at_every_level(tmp):
    ptm_cog = _raster()
    reg, cls, layer, _ = _fixture()
    msg = []
    for field, vals, allowed in (("ice_retention_class", cls, CLASSES),
                                 ("ice_table_layer", layer, LAYERS)):
        out = tmp / f"{field}.tif"
        ptm_cog.export_field(out, vals, reg, field, blocksize=64)
        levels, tags, n_ovr = _all_levels(out)
        assert n_ovr >= 1, f"{field}: no overview levels; the check would not see them"
        for k, v in enumerate(levels):
            bad = v[~np.isin(v, list(allowed))]
            assert v.size > 0 and bad.size == 0, f"{field} level {k}: {bad.size} of {v.size} not labels"
        assert "nearest face" in tags.get("PTM_RESAMPLING", ""), tags
        msg.append(f"{field} {len(levels)} levels")
    return ", ".join(msg)


def check_a_continuous_field_is_exported_as_before(tmp):
    ptm_cog = _raster()
    reg, _, _, depth = _fixture()
    out = tmp / "ice_depth.tif"
    grid, tags = ptm_cog.export_field(out, depth, reg, "ice_depth", blocksize=64)
    expect = reg.regrid_linear(depth)[::-1, :]
    fin = np.isfinite(expect)
    assert np.array_equal(grid[fin], expect[fin].astype(np.float32)), "continuous pixels changed"
    assert "PTM_RESAMPLING" not in tags
    return "barycentric pixels unchanged"


def main():
    with tempfile.TemporaryDirectory(prefix="ptm_cat_cog_") as d:
        tmp = Path(d)
        check("nearest regrid: labels only, same footprint",
              check_nearest_regrid_gives_labels_on_the_same_footprint)
        check("class and layer COGs: labels at every level",
              check_class_and_layer_cogs_hold_only_labels_at_every_level, tmp)
        check("continuous field unchanged", check_a_continuous_field_is_exported_as_before, tmp)
    print(f"\n[test_categorical_cog] {3 - len(FAILURES) - len(SKIPS)} passed, "
          f"{len(FAILURES)} failed, {len(SKIPS)} skipped")
    if FAILURES:
        print("[test_categorical_cog] FAIL")
        return 1
    print("[test_categorical_cog] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

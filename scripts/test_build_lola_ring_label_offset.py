#!/usr/bin/env python3
"""Regression test for the elevation scale/offset handling in scripts/build_lola_ring.py.

A LOLA polar product opened through its PDS label (e.g. LDEM_85S_40M.LBL) reports
Scale 0.5 and Offset 1737400 through GDAL. The label defines

    HEIGHT           = DN * SCALING_FACTOR
    PLANETARY_RADIUS = DN * SCALING_FACTOR + OFFSET   (OFFSET = reference sphere radius)

The builder places vertices at XYZ = (R + height) * spherical(lat, lon), so adding the
OFFSET to the height put every ring vertex near 2 R. This test runs the real main() on
small synthetic polar-stereographic rasters (an in-memory stand-in for osgeo.gdal, so
no GDAL install is needed) and checks the radial distance of every written vertex:

  1. Int16, scale 0.5, offset 1737400 (label style), no options: |XYZ| = R + height.
  2. The same raster with --elev-scale 0.5 --elev-offset 0: byte-identical to case 1.
  3. Float32 meters, no scale/offset metadata: heights unchanged.
  4. Float32 with a 100 m metadata offset (a height offset): applied as a height.
  5. Int16, scale 0.5, offset 1738000 (another lunar reference sphere): +600 m.
  6. Offset 3396190 m (not a lunar radius): the builder refuses.
  7. Explicit --elev-offset 1737400: still added to every height (override meaning kept).

Pure-python, no GPU, no sim binary -- runnable anytime:
    python3 scripts/test_build_lola_ring_label_offset.py
"""
import contextlib
import io
import sys
import tempfile
import types
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R = 1737400.0

# 400 x 400 px at 500 m, centered on the pole (+/-100 km).
DX = 500.0
NPX = 400
ULX, ULY = -NPX * DX / 2, NPX * DX / 2
SITE_X, SITE_Y = 20000.0, 10000.0
R_EXT_KM = 60.0


def _height_m(rows, cols):
    """Synthetic terrain in meters, a multiple of 0.5 m so DN = 2 h is exact."""
    h = 3000.0 * np.sin(cols / 37.0) - 2500.0 * np.cos(rows / 53.0)
    return np.round(2.0 * h) / 2.0


class _Raster:
    """Configuration of the fake band: data type, metadata scale/offset, DN encoding."""
    dtype = "Float32"
    scale = None
    offset = None
    dn_per_m = 1.0


def _fake_gdal(cfg):
    g = types.ModuleType("osgeo.gdal")
    g.GDT_Byte, g.GDT_UInt16, g.GDT_Int16, g.GDT_UInt32, g.GDT_Int32, g.GDT_Float32 = 1, 2, 3, 4, 5, 6
    g.UseExceptions = lambda: None
    g.GetDataTypeName = lambda t: str(t)

    class Band:
        @property
        def DataType(self):
            return g.GDT_Int16 if cfg.dtype == "Int16" else g.GDT_Float32
        def GetScale(self): return cfg.scale
        def GetOffset(self): return cfg.offset
        def ReadAsArray(self, xoff, yoff, xsize, ysize):
            rows, cols = np.meshgrid(np.arange(yoff, yoff + ysize),
                                     np.arange(xoff, xoff + xsize), indexing="ij")
            dn = _height_m(rows, cols) * cfg.dn_per_m
            return dn.astype(np.int16 if cfg.dtype == "Int16" else np.float32)

    class DS:
        RasterXSize = NPX
        RasterYSize = NPX
        def GetGeoTransform(self): return (ULX, DX, 0.0, ULY, 0.0, -DX)
        def GetRasterBand(self, i): return Band()

    g.Open = lambda path: DS()
    pkg = types.ModuleType("osgeo")
    pkg.gdal = g
    return pkg, g


def _read_ply(path):
    b = Path(path).read_bytes()
    h = b.index(b"end_header\n") + 11
    nv = int(b[:h].decode().split("element vertex ")[1].split()[0])
    vd = np.dtype([('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('nx', '<f4'), ('ny', '<f4'),
                   ('nz', '<f4'), ('lat', '<f8'), ('lon', '<f8'), ('el', '<f8')])
    return np.frombuffer(b, vd, nv, h), b


def _true_height(lat, lon):
    """Height of the synthetic terrain at each written vertex (its source pixel center)."""
    colat = np.radians(90.0 + lat)
    rho = 2.0 * R * np.tan(colat / 2.0)
    px = rho * np.sin(np.radians(lon))
    py = rho * np.cos(np.radians(lon))
    cols = np.round((px - ULX) / DX - 0.5)
    rows = np.round((ULY - py) / DX - 0.5)
    return _height_m(rows, cols)


def _run(blr, cfg, extra, out):
    """Run the builder's main() on cfg; return (exit message or None, stdout)."""
    rho = np.hypot(SITE_X, SITE_Y)
    site_lat = float(-(90.0 - np.degrees(2.0 * np.arctan(rho / (2.0 * R)))))
    site_lon = float(np.degrees(np.arctan2(SITE_X, SITE_Y)))
    pkg, g = _fake_gdal(cfg)
    blr.gdal = g
    sys.argv = ["build_lola_ring.py", "--lola", "synthetic.LBL",
                "--site-lat", repr(site_lat), "--site-lon", repr(site_lon),
                "--r-core-km", "5", "--r-ext-km", str(R_EXT_KM),
                "--max-points", "8000", "--allow-unregistered", "--out", str(out)] + extra
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            blr.main()
        except SystemExit as e:
            return str(e), buf.getvalue()
    return None, buf.getvalue()


def main():
    pkg, g = _fake_gdal(_Raster())
    sys.modules["osgeo"] = pkg
    sys.modules["osgeo.gdal"] = g
    sys.path.insert(0, str(HERE))
    import build_lola_ring as blr

    fails = 0

    def check(ok, msg):
        nonlocal fails
        print(("ok:   " if ok else "FAIL: ") + msg)
        fails += 0 if ok else 1

    def radial_case(name, cfg, extra, shift_m, td):
        """Assert |XYZ| = R + true height + shift_m and el = true height + shift_m."""
        out = Path(td) / f"{name}.ply"
        err, log = _run(blr, cfg, extra, out)
        if err is not None:
            check(False, f"{name}: builder exited: {err}")
            return None
        v, raw = _read_ply(out)
        r = np.sqrt(v['x']**2 + v['y']**2 + v['z']**2)
        h = _true_height(v['lat'], v['lon']) + shift_m
        dr = np.abs(r - (R + h)).max()
        de = np.abs(v['el'] - h).max()
        check(dr < 1e-6 and de < 1e-6,
              f"{name}: |XYZ|/R in [{r.min()/R:.6f}, {r.max()/R:.6f}], "
              f"max | |XYZ| - (R + h) | = {dr:.2e} m, max |el - h| = {de:.2e} m")
        return raw, log

    with tempfile.TemporaryDirectory() as td:
        label = _Raster(); label.dtype = "Int16"; label.scale = 0.5; label.offset = R
        label.dn_per_m = 2.0
        res1 = radial_case("1 label scale 0.5, offset 1737400, no options", label, [], 0.0, td)
        res2 = radial_case("2 same, --elev-scale 0.5 --elev-offset 0", label,
                           ["--elev-scale", "0.5", "--elev-offset", "0"], 0.0, td)
        if res1 and res2:
            check(res1[0] == res2[0], "2 is byte-identical to 1")
            check("radius of its reference sphere" in res1[1],
                  "1 reports that the label offset is a reference radius")

        tif = _Raster()
        radial_case("3 Float32 meters, no metadata", tif, [], 0.0, td)

        small = _Raster(); small.offset = 100.0
        radial_case("4 Float32, metadata offset 100 m", small, [], 100.0, td)

        other = _Raster(); other.dtype = "Int16"; other.scale = 0.5; other.offset = 1738000.0
        other.dn_per_m = 2.0
        radial_case("5 label offset 1738000 (another lunar sphere)", other, [], 600.0, td)

        mars = _Raster(); mars.dtype = "Int16"; mars.scale = 0.5; mars.offset = 3396190.0
        mars.dn_per_m = 2.0
        err, _ = _run(blr, mars, [], Path(td) / "6.ply")
        check(err is not None and "REFUSING" in err and not (Path(td) / "6.ply").exists(),
              f"6 non-lunar reference radius 3396190 m is refused (exit: {err!r:.60})")

        raw7 = radial_case("7 explicit --elev-offset 1737400 (override meaning kept)", label,
                           ["--elev-offset", str(R)], R, td)
        if raw7:
            check("WARNING: --elev-offset" in raw7[1], "7 warns that the offset is a radius")

    print("PASS" if fails == 0 else f"{fails} FAILURE(S)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

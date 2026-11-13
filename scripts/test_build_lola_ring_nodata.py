#!/usr/bin/env python3
"""Regression test: scripts/build_lola_ring.py never builds terrain from missing data.

Synthetic polar-stereographic rasters (an in-memory stand-in for osgeo.gdal, so no
GDAL install is needed) with a 10 km x 10 km block of missing pixels inside the
ring:

  1. Int16, scale 0.5, reference radius 1737400 (a LOLA label), NoData -32768:
     no vertex may come from a fill pixel (it would sit at -16384 m), and no
     triangle may lie across the block (its heights would be invented).
  2. Float32 meters with NaN in the block: the same.
  3. The same Int16 raster without the block: every triangle is kept, and the
     ring is byte-identical with and without a declared NoData value.

    python3 scripts/test_build_lola_ring_nodata.py
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
DX = 500.0
NPX = 400
ULX, ULY = -NPX * DX / 2, NPX * DX / 2
SITE_X, SITE_Y = 20000.0, 10000.0
R_EXT_KM = 60.0
# Missing block: columns/rows of the raster (10 km x 10 km), about 25 km from the site.
BLOCK_C = (270, 290)
BLOCK_R = (150, 170)


def _height_m(rows, cols):
    h = 3000.0 * np.sin(cols / 37.0) - 2500.0 * np.cos(rows / 53.0)
    return np.round(2.0 * h) / 2.0


def _fake_gdal(dtype, scale, offset, nodata, block):
    g = types.ModuleType("osgeo.gdal")
    g.GDT_Byte, g.GDT_UInt16, g.GDT_Int16, g.GDT_UInt32, g.GDT_Int32, g.GDT_Float32 = 1, 2, 3, 4, 5, 6
    g.UseExceptions = lambda: None
    g.GetDataTypeName = lambda t: str(t)

    class Band:
        DataType = g.GDT_Int16 if dtype == "Int16" else g.GDT_Float32
        def GetScale(self): return scale
        def GetOffset(self): return offset
        def GetUnitType(self): return ""
        def GetNoDataValue(self): return nodata
        def ReadAsArray(self, xoff, yoff, xsize, ysize):
            rows, cols = np.meshgrid(np.arange(yoff, yoff + ysize),
                                     np.arange(xoff, xoff + xsize), indexing="ij")
            h = _height_m(rows, cols)
            inb = ((cols >= BLOCK_C[0]) & (cols < BLOCK_C[1]) &
                   (rows >= BLOCK_R[0]) & (rows < BLOCK_R[1]))
            if dtype == "Int16":
                dn = (2.0 * h).astype(np.int16)
                if block:
                    dn[inb] = np.int16(-32768)
                return dn
            f = h.astype(np.float32)
            if block:
                f[inb] = np.nan
            return f

    class DS:
        RasterXSize = NPX
        RasterYSize = NPX
        def GetGeoTransform(self): return (ULX, DX, 0.0, ULY, 0.0, -DX)
        def GetMetadataDomainList(self): return [""]
        def GetRasterBand(self, i): return Band()

    g.Open = lambda path: DS()
    pkg = types.ModuleType("osgeo")
    pkg.gdal = g
    return pkg, g


def _read_ply(path):
    b = Path(path).read_bytes()
    h = b.index(b"end_header\n") + 11
    head = b[:h].decode()
    nv = int(head.split("element vertex ")[1].split()[0])
    nf = int(head.split("element face ")[1].split()[0])
    vd = np.dtype([('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('nx', '<f4'), ('ny', '<f4'),
                   ('nz', '<f4'), ('lat', '<f8'), ('lon', '<f8'), ('el', '<f8')])
    v = np.frombuffer(b, vd, nv, h)
    f = np.frombuffer(b, np.dtype([('n', 'u1'), ('i', '<i4', 3)]), nf, h + nv * vd.itemsize)
    return v, f['i'], b


def _xy(lat, lon):
    rho = 2.0 * R * np.tan(np.radians(90.0 + lat) / 2.0)
    return rho * np.sin(np.radians(lon)), rho * np.cos(np.radians(lon))


def _run(blr, pkg_g, out):
    rho = np.hypot(SITE_X, SITE_Y)
    site_lat = float(-(90.0 - np.degrees(2.0 * np.arctan(rho / (2.0 * R)))))
    site_lon = float(np.degrees(np.arctan2(SITE_X, SITE_Y)))
    blr.gdal = pkg_g[1]
    sys.argv = ["build_lola_ring.py", "--lola", "synthetic.LBL",
                "--site-lat", repr(site_lat), "--site-lon", repr(site_lon),
                "--r-core-km", "5", "--r-ext-km", str(R_EXT_KM),
                "--max-points", "8000", "--allow-unregistered", "--out", str(out)]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        blr.main()
    return buf.getvalue()


def main():
    pkg, g = _fake_gdal("Int16", 0.5, R, -32768.0, True)
    sys.modules["osgeo"] = pkg
    sys.modules["osgeo.gdal"] = g
    sys.path.insert(0, str(HERE))
    import build_lola_ring as blr

    # block bounds in projected meters
    bx0, bx1 = ULX + BLOCK_C[0] * DX, ULX + BLOCK_C[1] * DX
    by1, by0 = ULY - BLOCK_R[0] * DX, ULY - BLOCK_R[1] * DX
    fails = 0

    def check(ok, msg):
        nonlocal fails
        print(("ok:   " if ok else "FAIL: ") + msg)
        fails += 0 if ok else 1

    with tempfile.TemporaryDirectory() as td:
        for tag, cfg in (("1 Int16 label, NoData -32768", ("Int16", 0.5, R, -32768.0, True)),
                         ("2 Float32 meters, NaN", ("Float32", None, None, None, True))):
            out = Path(td) / "ring.ply"
            _run(blr, _fake_gdal(*cfg), out)
            v, f, _ = _read_ply(out)
            check(v['el'].min() > -6000.0,
                  f"{tag}: lowest vertex {v['el'].min():.1f} m (a fill pixel would give -16384 m)")
            x, y = _xy(v['lat'], v['lon'])
            inside_v = (x > bx0) & (x < bx1) & (y > by0) & (y < by1)
            check(not inside_v.any(), f"{tag}: {int(inside_v.sum())} vertices inside the missing block")
            cx, cy = x[f].mean(axis=1), y[f].mean(axis=1)
            inside_f = (cx > bx0 + DX) & (cx < bx1 - DX) & (cy > by0 + DX) & (cy < by1 - DX)
            check(not inside_f.any(),
                  f"{tag}: {int(inside_f.sum())} triangles laid across the missing block "
                  f"({f.shape[0]} triangles in the ring)")

        a, b = Path(td) / "a.ply", Path(td) / "b.ply"
        log_a = _run(blr, _fake_gdal("Int16", 0.5, R, -32768.0, False), a)
        _run(blr, _fake_gdal("Int16", 0.5, R, None, False), b)
        check(a.read_bytes() == b.read_bytes(),
              "3 without missing pixels, the ring is byte-identical with and without a "
              "declared NoData value")
        check("missing pixels" not in log_a, "3 no triangles dropped when nothing is missing")

    print("PASS" if fails == 0 else f"{fails} FAILURE(S)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

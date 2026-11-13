#!/usr/bin/env python3
"""Build a LOLA occluder-ring PLY WITHOUT dtm_to_ply_delaunay (which is gone
from the evicted dtm/ tree). Replicates the polar-site campaign ring builder's output:
a coarse LOLA annulus [R_core, R_ext] around the site, Delaunay-triangulated
in the south-polar-stereo plane, written in the exact dtm PLY schema
(binary_le; dbl x,y,z / flt nx,ny,nz / dbl lat,lon,elev; faces uchar+int*3),
in the selenocentric frame XYZ = (R+elev)*spherical(lat,lon) — verified to
0.000 m against site A's dtm-built ring.

  python3 build_lola_ring.py --lola <LOLA.tif> --site-lat -85.0 \
      --site-lon 30.0 --r-core-km 3 --r-ext-km 81.6712 \
      --max-points 80000 --core-ply core.ply --out ring.ply

The inner hole is registered to the core mesh's ACTUAL outer boundary
polygon (--core-ply). The legacy nominal-circle centroid clip kept every
triangle whose centroid sat outside R_core, so its inward-most vertex could
protrude ~130 m inside — laying coarse ring surface on top of the core's own
terrain, which at grazing polar sun casts a spurious hard shadow on core edge
faces. See scripts/ring_registration.py.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
from osgeo import gdal
gdal.UseExceptions()

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ring_registration as rreg  # noqa: E402
import dem_elevation  # noqa: E402

R = 1737400.0  # Moon 2015 sphere radius (m), matches the LOLA base CRS

def geodesic_m(slat, slon, lat, lon):
    """Great-circle distance (m) on the R sphere; all args radians."""
    ca = np.sin(slat)*np.sin(lat) + np.cos(slat)*np.cos(lat)*np.cos(lon-slon)
    return R * np.arccos(np.clip(ca, -1.0, 1.0))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lola", required=True)
    ap.add_argument("--site-lat", type=float, required=True)  # deg, south negative
    ap.add_argument("--site-lon", type=float, required=True)  # deg east
    ap.add_argument("--r-core-km", type=float, default=3.0)
    ap.add_argument("--r-ext-km", type=float, required=True)
    ap.add_argument("--max-points", type=int, default=80000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--elev-scale", type=float, default=None,
                    help="Meters per stored DEM unit. Default: the file's scale metadata "
                         "times its unit (a LOLA polar product opened through its PDS .LBL "
                         "gives 0.5); 1 for a float raster without it. Required for an "
                         "integer raster without scale metadata. See scripts/dem_elevation.py.")
    ap.add_argument("--elev-offset", type=float, default=None,
                    help="Meters added to every height. Default: read band.GetOffset() "
                         "(or 0). A metadata offset larger than 100 km is the radius of "
                         "the DEM's reference sphere (PDS OFFSET: PLANETARY_RADIUS = "
                         "DN*SCALING_FACTOR + OFFSET), so heights are taken relative to "
                         "the 1737.4 km sphere instead. See scripts/dem_elevation.py.")
    ap.add_argument("--core-ply", default=None,
                    help="Core mesh PLY. Its actual outer boundary polygon "
                         "defines the ring's inner hole. Required unless "
                         "--allow-unregistered.")
    ap.add_argument("--seam-margin-m", type=float, default=0.0,
                    help="Extra clearance beyond the core boundary (m).")
    ap.add_argument("--allow-unregistered", action="store_true",
                    help="Use the legacy nominal-circle centroid clip. The "
                         "resulting ring overlaps the core mesh and casts "
                         "spurious hard shadows on its edge faces.")
    a = ap.parse_args()
    if a.core_ply is None and not a.allow_unregistered:
        raise SystemExit(
            "[build_lola_ring] REFUSING to proceed: no --core-ply. The ring's "
            "inner hole must be registered to the core mesh's actual outer "
            "boundary; the legacy nominal-circle centroid clip leaves up to "
            "~130 m of ring surface lying on top of the core and casts a "
            "spurious hard shadow on its edge faces. Pass --core-ply, or "
            "--allow-unregistered to reproduce the old geometry deliberately.")
    rng = np.random.default_rng(a.seed)

    slat, slon = np.radians(a.site_lat), np.radians(a.site_lon)
    # site projected coords (south-polar stereo, lat_0=-90 lon_0=0):
    #   rho = 2R tan(colat/2), x = rho sin(lon), y = rho cos(lon)
    colat = np.radians(90.0 - abs(a.site_lat))
    rho = 2.0*R*np.tan(colat/2.0)
    x0, y0 = rho*np.sin(slon), rho*np.cos(slon)

    ds = gdal.Open(a.lola)
    gt = ds.GetGeoTransform()      # (ulx, dx, 0, uly, 0, dy<0)
    ulx, dx, _, uly, _, dyg = gt
    dy = abs(dyg)
    band = ds.GetRasterBand(1)

    # read a window of half-size R_ext*1.05 around the site
    half = int(a.r_ext_km*1.05*1000.0/dx) + 2
    col0 = max(0, int((x0-ulx)/dx) - half)
    row0 = max(0, int((uly-y0)/dy) - half)
    # Clip each axis to the raster on its own. A single square side clipped by
    # min(x-room, y-room) let a site near ONE raster edge silently lose ring
    # coverage on the OTHER axis too (scripts/test_build_lola_ring_window.py).
    ncol = min(2*half, ds.RasterXSize-col0)
    nrow = min(2*half, ds.RasterYSize-row0)
    want0c, want0r = int((x0-ulx)/dx) - half, int((uly-y0)/dy) - half
    if (want0c < 0 or want0r < 0 or want0c + 2*half > ds.RasterXSize
            or want0r + 2*half > ds.RasterYSize):
        print(f"[build_lola_ring] WARNING: the requested {a.r_ext_km} km ring "
              f"extends past the DEM edge; ring coverage is cut there (window "
              f"cols {col0}..{col0+ncol}, rows {row0}..{row0+nrow} of "
              f"{ds.RasterXSize} x {ds.RasterYSize}).")
    arr = band.ReadAsArray(col0, row0, ncol, nrow).astype(np.float64)

    # Stored value -> height in meters, by the rule shared with the C++ DEM tools
    # (scripts/dem_elevation.py). GDAL's ReadAsArray returns the stored values;
    # a PDS label's SCALING_FACTOR/OFFSET/UNIT are metadata only. Using raw LOLA
    # DN as meters once DOUBLED a ring; adding the label's OFFSET (the reference
    # radius) once put one at 2 R.
    m_scale, m_offset, m_unit, m_int, nodata = dem_elevation.gdal_band_metadata(ds, band, gdal)
    try:
        tf = dem_elevation.resolve(m_scale, m_offset, m_unit, m_int,
                                   a.elev_scale, a.elev_offset, R)
    except dem_elevation.DemElevationError as e:
        raise SystemExit(f"[build_lola_ring] REFUSING to proceed: {a.lola}: {e}")
    for w in tf.warnings:
        print(f"[build_lola_ring] WARNING: {w}")
    print(f"[build_lola_ring] elevation: {tf.describe()}")
    scale, offset = tf.scale, tf.offset
    # Missing data (declared NoData, NaN, float fill) never becomes a vertex,
    # and no triangle is laid across it (see below).
    missing = dem_elevation.missing_mask(arr, nodata)
    arr = arr * scale + offset

    cols = np.arange(col0, col0+ncol)
    rows = np.arange(row0, row0+nrow)
    PX = ulx + (cols+0.5)*dx
    PY = uly - (rows+0.5)*dy
    PXg, PYg = np.meshgrid(PX, PY)
    px = PXg.ravel(); py = PYg.ravel(); el = arr.ravel()
    good = ~missing.ravel()
    px, py, el = px[good], py[good], el[good]

    # pixel -> lat/lon (inverse polar stereo)
    rho_p = np.sqrt(px*px + py*py)
    colat_p = 2.0*np.arctan2(rho_p, 2.0*R)
    lat = -(90.0 - np.degrees(colat_p))            # south negative, deg
    lon = np.degrees(np.arctan2(px, py))           # east, deg
    d = geodesic_m(slat, slon, np.radians(lat), np.radians(lon))

    disk = d <= a.r_ext_km*1000.0                  # mesh the full disk, drop core later
    px, py, el, lat, lon = px[disk], py[disk], el[disk], lat[disk], lon[disk]
    d = d[disk]

    # subsample to the point budget. Weight prob ~ 1/dist so the near-field
    # (which sets the shadowing horizon) stays dense, ~matching dtm's
    # density-exp=1.0 linear falloff. (Uniform-area sampling under-resolves
    # the near horizon and lets grazing sun leak -> validated-failing.)
    if px.size > a.max_points:
        w = 1.0/np.maximum(d, 1000.0); w /= w.sum()
        sel = rng.choice(px.size, a.max_points, replace=False, p=w)
        px, py, el, lat, lon, d = px[sel], py[sel], el[sel], lat[sel], lon[sel], d[sel]

    # selenocentric XYZ = (R+elev)*spherical(lat,lon)  [verified vs dtm: 0.000 m]
    la, lo = np.radians(lat), np.radians(lon)
    rr = R + el
    X = rr*np.cos(la)*np.cos(lo); Y = rr*np.cos(la)*np.sin(lo); Z = rr*np.sin(la)

    # Delaunay in the projection plane (same as dtm tool)
    from scipy.spatial import Delaunay
    tri = Delaunay(np.column_stack([px, py]))
    faces = tri.simplices.astype(np.int64)

    # Where data are missing the ring has a hole: drop every triangle whose
    # centroid or an edge midpoint falls on a missing pixel, instead of letting
    # Delaunay bridge the gap with heights the DEM does not have.
    if missing.any():
        fx, fy = px[faces], py[faces]
        probes_x = np.column_stack([fx.mean(axis=1), (fx[:, 0]+fx[:, 1])/2,
                                    (fx[:, 1]+fx[:, 2])/2, (fx[:, 2]+fx[:, 0])/2])
        probes_y = np.column_stack([fy.mean(axis=1), (fy[:, 0]+fy[:, 1])/2,
                                    (fy[:, 1]+fy[:, 2])/2, (fy[:, 2]+fy[:, 0])/2])
        pc = np.clip(np.floor((probes_x - ulx)/dx).astype(np.int64) - col0, 0, ncol-1)
        pr = np.clip(np.floor((uly - probes_y)/dy).astype(np.int64) - row0, 0, nrow-1)
        over_gap = missing[pr, pc].any(axis=1)
        print(f"[build_lola_ring] {int(missing.sum())} missing pixels in the DEM window: "
              f"dropped {int(over_gap.sum())} of {faces.shape[0]} triangles that cross them")
        faces = faces[~over_gap]

    # punch the inner hole
    core = None
    if a.core_ply is not None:
        core = rreg.core_boundary_profile(Path(a.core_ply), a.site_lat,
                                          a.site_lon, body_r_m=R)
        keep, _, _ = rreg.clip_ring_to_core(lat, lon, faces, core,
                                            a.site_lat, a.site_lon,
                                            margin_m=a.seam_margin_m,
                                            body_r_m=R)
        print(f"[build_lola_ring] inner hole registered to "
              f"{Path(a.core_ply).name} (boundary r = {core.r.min():.1f} .. "
              f"{core.r.max():.1f} m): kept {int(keep.sum())} of "
              f"{faces.shape[0]} disk triangles")
        faces = faces[keep]
    else:
        print("[build_lola_ring] WARNING: UNREGISTERED ring (legacy "
              "nominal-circle centroid clip). Ring geometry will overlap the "
              "core mesh and cast spurious hard shadows on its edge faces. "
              "Do not use this for science.")
        # drop triangles whose centroid is inside R_core (geodesic) -> annulus
        clat = la[faces].mean(axis=1); clon = lo[faces].mean(axis=1)
        cd = geodesic_m(slat, slon, clat, clon)
        faces = faces[cd >= a.r_core_km*1000.0]

    # renumber used vertices
    used = np.unique(faces)
    remap = np.full(px.size, -1, dtype=np.int64); remap[used] = np.arange(used.size)
    faces = remap[faces]
    X, Y, Z, lat, lon, el = X[used], Y[used], Z[used], lat[used], lon[used], el[used]

    # outward radial normals (irrelevant to shadow BVH, but fill the schema)
    nrm = np.column_stack([X, Y, Z]); nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)

    nv, nf = X.size, faces.shape[0]
    vdt = np.dtype([('x','<f8'),('y','<f8'),('z','<f8'),('nx','<f4'),('ny','<f4'),
                    ('nz','<f4'),('lat','<f8'),('lon','<f8'),('el','<f8')])
    varr = np.empty(nv, dtype=vdt)
    varr['x'],varr['y'],varr['z'] = X,Y,Z
    varr['nx'],varr['ny'],varr['nz'] = nrm[:,0],nrm[:,1],nrm[:,2]
    varr['lat'],varr['lon'],varr['el'] = lat,lon,el
    frows = np.zeros((nf,13), dtype=np.uint8)
    frows[:,0] = 3
    frows[:,1:13] = faces.astype('<i4').view(np.uint8).reshape(nf,12)

    hdr = (b"ply\nformat binary_little_endian 1.0\n"
           b"comment DTM to PLY with Delaunay triangulation for PTM3D\n"
           b"element vertex %d\n"
           b"property double x\nproperty double y\nproperty double z\n"
           b"property float nx\nproperty float ny\nproperty float nz\n"
           b"property double latitude\nproperty double longitude\nproperty double elevation\n"
           b"element face %d\nproperty list uchar int vertex_indices\nend_header\n"
           % (nv, nf))
    with open(a.out, "wb") as f:
        f.write(hdr); varr.tofile(f); frows.tofile(f)
    print(f"wrote {a.out}: {nv} verts, {nf} faces, "
          f"elev[{el.min():.0f},{el.max():.0f}] dist_km[{(geodesic_m(slat,slon,la[used] if False else np.radians(lat),np.radians(lon)).min())/1e3:.2f},"
          f"{geodesic_m(slat,slon,np.radians(lat),np.radians(lon)).max()/1e3:.2f}]")

    # registration gate: re-measure what was actually written
    if core is not None:
        rlat, rlon, _, rtris = rreg.read_ply(Path(a.out))
        rep = rreg.registration_report(rlat, rlon, rtris, core,
                                       a.site_lat, a.site_lon, body_r_m=R)
        rreg.print_report(rep, name=Path(a.out).name)
        if not rreg.gate(rep, name=Path(a.out).name):
            # Do not leave a rejected ring at the requested path: the caller
            # sees exit 1, but anything that later globs for the file finds a
            # ring that failed its own gate.
            try:
                Path(a.out).unlink()
                print(f"[build_lola_ring] removed the rejected ring {a.out}")
            except OSError as e:
                print(f"[build_lola_ring] WARNING: could not remove the "
                      f"rejected ring {a.out}: {e}. DELETE IT BY HAND.")
            raise SystemExit(1)

if __name__ == "__main__":
    main()

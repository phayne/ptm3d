# Building a mesh from a digital terrain model

This page describes how a digital terrain model (DTM) raster in GeoTIFF format is converted to
a PLY mesh for `ptm_thermal_sim`, from the command line and from the DEM module of the GUI
dashboard.

The converter triangulates the terrain and can grade the mesh density around a region of
interest, so that the mesh is fine at the site under study and coarse far from it (section 4).
There is no separate mesh-refinement step. The
[Status and known gaps](#status-and-known-gaps) section lists what the tools do not do.

---

## 1. The converter

The conversion is done by one program, which is built into `build/bin/`.

| Program | Source | Function |
|---|---|---|
| `dtm_to_ply_delaunay` | `src/cpp/tools/dtm_to_ply_delaunay.cpp` | GeoTIFF DTM to PLY. Delaunay triangulation of the projected point cloud, region-of-interest extraction, radial density grading, edge filtering, mesh quality control. |

A second program, `dem_mesh_generator` (`src/cpp/tools/dem_mesh_generator.cpp`), is used by
`data/dem/copernicus/README.md`. The GUI does not use it and this page does not describe it,
except that it reads heights in the same way (section 3, "Heights"). Its sixth and seventh
arguments set the scale and the offset.

Both programs need GDAL when they are built. Without GDAL, CMake prints

```
CMake Warning: GDAL not properly configured. SKIPPING these targets:
dem_mesh_generator, dtm_to_ply_delaunay. ...
```

and configures the rest of the model without them. To build them, install GDAL and run
`cmake -S . -B build` again.

### Build

```bash
# Prerequisite (macOS/Linux, recommended):
conda install -c conda-forge gdal

cd <repo>
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target dtm_to_ply_delaunay -j 8
```

`dtm_to_ply_delaunay` also uses [CDT](https://github.com/artem-ogre/CDT) 1.4.4, a header-only
Delaunay library included in `third_party/CDT/` (MPL-2.0), so there is nothing further to
install. The converter is compiled with `-fno-fast-math -ffp-contract=off` in place of the
`-ffast-math` used for the rest of the model. The geometric predicates of the triangulation
require exact IEEE-754 arithmetic, and without it the triangulation can be wrong for nearly
degenerate triangles, which thinning around a region of interest produces.

---

## 2. Examples

The first example converts a whole DTM at half resolution.

```bash
./build/bin/dtm_to_ply_delaunay input_dtm.tif out_mesh.ply \
    --downsample 2 \
    --body-radius 1737400
```

The second extracts a region of interest of 2 km radius around a site, with a limit on the
number of points and a core at full resolution.

```bash
./build/bin/dtm_to_ply_delaunay input_dtm.tif site_mesh.ply \
    --body-radius 1737400 \
    --roi-center -89.95 0.0 \
    --roi-radius 2.0 \
    --roi-max-points 8000 \
    --roi-core-frac 0.3
```

The second command is the typical use. The mesh is fine near the site, where shadows and
temperature gradients must be resolved, and coarse elsewhere.

---

## 3. `dtm_to_ply_delaunay` options

```
Usage: dtm_to_ply_delaunay input.tif output.ply [options]
```

### Sampling

| Flag | Default | Meaning |
|---|---|---|
| `--downsample N` | 1 | Take every Nth pixel in each direction before anything else. Applied to the whole raster, ROI or not. |
| `--body-radius METERS` | auto | Body radius used to place vertices in Cartesian coordinates and to measure geodesic distances. Auto-detected from the GeoTIFF's CRS when omitted; pass it explicitly if the CRS is missing or wrong. Moon = 1737400. |
| `--elev-scale S` | from the file | Meters per stored value. Needed for an integer file without scale metadata. |
| `--elev-offset METERS` | from the file | Meters added to every height. |

### Heights

A DEM stores numbers, and the file or its label gives the conversion to heights. Open a LOLA
product from the PDS through its label. The integer products store heights in steps of 0.5 m, the
floating-point products store kilometers, and every label gives the radius of the reference
sphere. A JPEG 2000 file is opened through its own `_jp2.lbl` label, because the image alone does
not record its scale. The LOLA GeoTIFFs from the Planetary Geodesy Data Archive and the SfS
GeoTIFFs store heights in meters and carry no such information. The converter reads the scale,
the offset and the unit from the file, and the log reports the values it used and where they
came from.

```
Elevation: height (m) = stored value * 0.5 [file metadata, unit METER] + 0 [file metadata: reference radius 1.7374e+06 m, unit METER]
```

An offset larger than 100 km is the radius of the reference sphere. Heights are then taken
relative to that sphere and shifted by its difference from the body radius. A file in kilometers
is converted to meters when its label gives the unit. A floating-point file without a scale is taken
to be in meters. The run stops with an error for an integer file without a scale, for a reference
radius more than 20 km from the body radius, and for an offset equal to the body radius in
kilometers when the file gives no unit. `--elev-scale` and `--elev-offset` set the values
directly. Pixels marked with the file's NoData value never become vertices.

`dem_mesh_generator` and `scripts/build_lola_ring.py` read heights the same way.

### Region of interest and density grading

| Flag | Default | Meaning |
|---|---|---|
| `--roi-center LAT LON` | (off) | Enables ROI mode and sets its center, in degrees. |
| `--roi-radius KM` | 0 | Hard cutoff: points beyond this geodesic distance from the center are dropped entirely. |
| `--roi-max-points N` | 0 (no budget) | Target number of vertices. The falloff strength `beta` is found by bisection so that the expected number of retained points is N or fewer. |
| `--roi-core-frac F` | 0.15 | Fraction of `--roi-radius` kept at full resolution. `0.3` with a 2 km radius means the inner 0.6 km is untouched and grading starts outside it. Clamped to [0, 1]. If the core alone exceeds the budget, a message reports this and the entire ROI is thinned. |
| `--density-exp P` | 2.0 | Falloff exponent. Keep probability outside the core is `1 / (1 + beta·d^P)²` where `d` is the normalized distance from the core edge to the rim. 1.0 = linear, 2.0 = quadratic, 0.5 = square root. Lower = gentler falloff. Floored at 0.1. |

Points are thinned at random with the keep probability above, and not at a regular stride, so
the density changes smoothly with distance. The random number generator starts from a fixed
seed (42), so the same command on the same input gives the same mesh.

### Triangulation

| Flag | Default | Meaning |
|---|---|---|
| `--max-edge METERS` | auto | Drop triangles with any projected edge longer than this. `0` disables filtering entirely (no holes, but long thin triangles across ROI gaps survive). Omit for auto: `pixel_size × downsample × (1 + beta) × 3`. |
| `--swap-axes` | off | Swap X/Y when inverting a polar stereographic projection. Use when a polar DTM appears mirrored or rotated 90°. |

### Mesh quality control

Quality control follows the triangulation. Slivers are removed first, the face normals are then
given a consistent orientation by a breadth-first flood fill, and small holes are filled last.

| Flag | Default | Meaning |
|---|---|---|
| `--no-qc` | (QC on) | Skip the entire quality-control sweep. |
| `--qc-sliver-angle DEG` | 5.0 | Remove boundary triangles whose minimum angle is below this. Interior slivers are kept, because removing them would open holes. |
| `--qc-max-hole N` | 20 | Fill boundary loops of at most N edges (triangle fill for 3, centroid fan otherwise). `0` disables hole filling. |

### Not implemented

| Flag | Status |
|---|---|
| `--south-pole-coords` | Accepted, with no effect. It applies no coordinate correction, and meshes built with and without it are identical. A warning is printed when it is given. To swap the axes of a polar stereographic projection, use `--swap-axes`. |

An unrecognized or incomplete argument produces a warning on stderr, and the conversion
continues.

### Exit codes

| Code | Meaning |
|---|---|
| 0 | A mesh with at least one face was written. |
| 1 | The input could not be opened, its heights cannot be read without `--elev-scale` or `--elev-offset` (see "Heights"), or the PLY could not be written. |
| 2 | The mesh came out empty, and no file was written. |

A mesh without faces means that the request could not be met. The usual causes are an
`--roi-center` outside the raster (often a wrong sign of latitude on a polar DEM), an
`--roi-radius` smaller than the pixel size, a `--max-edge` shorter than the point spacing, which
removes every triangle, or a region that contains only NoData. In this case the converter
writes no file, lists these causes, and reports the number of points skipped and triangles
filtered. The GUI uses the exit code to decide whether a conversion succeeded.

### Output format

The output is a binary little-endian PLY file.

- Vertices are `double x y z`, `float nx ny nz`, `double latitude longitude elevation`.
  Latitude, longitude and elevation are stored in double precision with the mesh, so later
  steps do not have to derive them again.
- Faces are a `uchar` count followed by three `int` vertex indices.

The Cartesian coordinates are body-fixed, with +X through (0°N, 0°E), +Y through (0°N, 90°E)
and +Z through the north pole. The radius of a vertex is `body_radius + elevation`.

---

## 4. Grading mesh resolution with the ROI controls

The region-of-interest options are the way to obtain a mesh that is fine at a site and coarse
at its rim within a given number of points. The table relates each property of the mesh to
the option that sets it.

| Property of the mesh | `dtm_to_ply_delaunay` flag |
|---|---|
| Center the fine region on a site | `--roi-center LAT LON` |
| Radius kept at full resolution | `--roi-core-frac` × `--roi-radius` |
| Radius beyond which nothing is kept | `--roi-radius KM` |
| A budget on mesh size | `--roi-max-points N` (vertices; faces ≈ 2×) |
| The shape of the falloff between them | `--density-exp P` |
| An upper bound on triangle edge length | `--max-edge METERS` |

As an example, a synthetic south-polar lunar DTM of 400×400 pixels at 20 m/pixel was converted
with `--roi-center -89.95 0.0 --roi-radius 2.0 --roi-max-points 8000 --roi-core-frac 0.3`.
The vertex density in rings around the center is given below.

| Ring from ROI center | Vertices | Vertex density |
|---|---|---|
| 0.0 – 0.6 km (the protected core) | 2822 | 2495 verts/km² |
| 0.6 – 1.2 km | 4280 | 1262 verts/km² |
| 1.2 – 1.8 km | 788 | 139 verts/km² |
| 1.8 – 2.4 km | 95 | 12 verts/km² |

The mesh has 7985 vertices and 15916 faces, within the budget of 8000 vertices, and its
farthest vertex is 1.9996 km from the center for the requested cutoff of 2.0 km. The vertex
density at the core is 208 times that at the rim.

---

## 5. GUI dashboard

```bash
cd gui/
python run_dashboard.py
```

The dashboard is an optional launcher for several modules. The thermal application alone is
started with `run_gui.py`, and the DEM module is available only in the dashboard. Select the
**DEM** module in the sidebar. It has three pages, named Preview, Convert and Inspect.

**Preview page** (`dem_panels/dem_preview_panel.py`)

1. Click **Browse…** and select a GeoTIFF DTM. The panel reads it with the GDAL Python
   bindings (`osgeo`) and draws an elevation map. Without `osgeo` the panel accepts the file
   and shows a notice in place of the map. The conversion does not depend on `osgeo`, because
   the converter is linked with GDAL.
2. Check **Enable ROI extraction** and set the center latitude, longitude and radius, or click
   on the preview to place the center.

**Convert page** (`dem_panels/dem_config_panel.py`)

- *Grid & Body*: downsample factor, body (sets the radius), body radius.
- *ROI*: core fraction, density exponent. The center and radius come from the
  Preview page.
- *Triangulation*: max-edge mode (auto / disabled / explicit) and its value.
- *Coordinate corrections*: swap axes, and south pole coords (this checkbox has no effect,
  see section 3).
- *Quality control*: enable, sliver angle, max hole edges.
- *Output*: the destination `.ply`.

Press **Convert**. The status line shows the converter's output as it runs, including the
percentage completed.

**Inspect page** (`dem_panels/dem_mesh_panel.py`)

This page gives statistics of the converted mesh, which is loaded automatically. The 3D viewer
is the right pane of the dashboard.

There is no page for mesh refinement. Use the ROI controls of the Convert page (section 4).

### Converter search path of the GUI

`gui/dem_core/dem_runner.py` searches the following locations in order.

1. `<repo>/build/bin/<name>`, which `cmake --build build` produces.
2. `<repo>/../dtm/<name>`, a standalone `dtm/` directory beside the repository. A converter
   found there may have been built for a different CPU architecture.
3. `PATH`.

If the converter is not found, the module prints the `cmake` commands that build it and notes
that GDAL is required.

### Platform support

The Convert action of the DEM module works on macOS and Linux. It reads the converter's output
with `select.select()` on pipes, which Windows does not support. On Windows the module reports
that the action is unavailable, and the conversion is run from the command line.

---

## 6. Using the mesh in a simulation

Give the PLY file as `mesh.import_file`.

```yaml
mesh:
  import_file: site_mesh.ply
```

Check two properties of a mesh derived from a DEM before using it.

- Topology. A mesh can have the expected number of faces, radial extent and elevations and
  still be non-manifold or contain interior holes. Count the boundary edges, which should all
  lie on the outer rim, and the non-manifold edges, of which there should be none. To obtain a
  smaller mesh, convert the raster again with fewer points in preference to decimating an
  existing mesh, because a new triangulation has no such defects.
- Units. Meshes derived from DEMs are in meters, and asteroid `.obj` shape models are in
  kilometers. A mesh read in the wrong unit is wrong in scale by a factor of 1000, and so are
  the offsets of the shadow rays.

---

## Status and known gaps

| Item | Status |
|---|---|
| DTM to PLY conversion | Available. |
| ROI extraction + radial density grading | Available (`dtm_to_ply_delaunay --roi-*`). |
| Mesh quality control (slivers, normals, holes) | Available. |
| Empty-mesh handling | The converter writes no file and exits with code 2 when the mesh has no faces. |
| Adaptive mesh refinement | Not available. Use ROI density grading (section 4). |
| `--south-pole-coords` | No effect. A warning is printed when it is given. |
| Windows support for the GUI Convert action | Not supported (POSIX pipe I/O). The command line is available wherever the converter can be built. |

## Tests

The test for the converter is `src/python/tests/test_dem_tools.py`. It generates a synthetic
lunar DTM in south polar stereographic projection, so it needs no data files. The test confirms
that the converter reports `Detected POLAR STEREOGRAPHIC projection`, so that the projected
coordinates are the ones tested. When a program has not been built, the tests that need it are
skipped with a message that names the program and the `cmake --build` command.

The converter must be compiled with `-fno-fast-math -ffp-contract=off` for its output to be
reproducible. With `-ffast-math` the single-precision vertex normals change in the last one or
two bits (at most 2.98e-7), while the double-precision coordinates and the faces are unchanged.

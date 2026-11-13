# Core mesh plus occluder ring

## Overview

High-resolution site studies need distant topography (massifs and crater rims
out to the physics horizon, typically 20 km) to shadow the region of interest,
but the temperature of that terrain does not need to be computed. Meshing it at
the resolution of the core multiplies the face count several-fold and makes the
thermal solution, the view-factor matrix and the output volume impractically
large. A core plus ring of 2.8 M faces requires roughly 60 GB and 605 M stored
view-factor entries, compared with about 1.8 M faces for the core alone.

The occluder ring separates the two roles. The simulation runs on the core mesh
exactly as it would without a ring, and a second, coarse mesh is given to the
illumination ray tracer as shadow-casting geometry only.

| Mesh | Source | Role | Typical size |
|---|---|---|---|
| core | site DEM at full resolution | thermal, radiation, statistics, output | ~10<sup>6</sup> faces |
| ring | coarse polar DEM annulus, inner radius to physics horizon | occluder triangles in the shadow-ray BVH only | ~10<sup>5</sup> faces |

Ring triangles are appended to the ray-tracing acceleration structure after the
core faces, so their primitive indices are at or above the core face count and
never index a flux, temperature or output array. The shadow test requires only an
any-hit query: a ray from a core face that intersects a ring triangle is shadowed.

Everything after the illumination calculation is unchanged. `n_faces` equals the
face count of the core mesh in the thermal solver, in the N×N view-factor
operator, in the statistics and in the VTK and binary output.

## Configuration

```yaml
mesh:
  import_file:   core.ply     # solved
  occluder_file: ring.ply     # shadow-casting only
```

With an empty or absent `occluder_file` the ring is not used, and the results are
identical to those of a run on the core alone. When the file is set but cannot be
loaded, a warning that names the file is issued and the run continues without
shadowing by distant terrain.

A successful load reports the geometry that was registered:

```
Loading occluder mesh (shadow-only): ring.ply
Occluder mesh registered: 61490 verts, 121872 tris (shadow BVH only)
```

`mesh.scale_factor` is applied to the occluder as well as the core, so both
meshes must be in the same length unit.

The occluder ring is available with both the CPU (Embree) and the GPU ray tracers.

## Building the ring

`scripts/build_lola_ring.py` builds a coarse annulus from a polar DEM. Open a polar
PDS product through its label (`.LBL`). The label gives the height units and the
radius of the reference sphere, and the script converts the stored integers to
heights above the 1,737,400 m sphere. It reads heights by the same criterion as the
mesh converter ([dem-to-mesh-workflow.md](dem-to-mesh-workflow.md), "Heights"), so an
integer raster without scale metadata needs `--elev-scale`. Where the DEM has no data,
the ring has a hole. No vertex is taken from a NoData pixel, and no triangle spans one.

Coarse facets are adequate. Occluders only have to block rays, so the native DEM
resolution of 80-160 m over the annulus is sufficient, and the added cost of ray
tracing is small.

The script first builds a coarse full disk from the clipped DEM of the site and then
removes the inner hole.

### Registration of the inner hole to the core mesh

The inner hole of the ring is not a nominal circle of radius `R_core`. It is
clipped against the actual outer boundary polygon of the core mesh, which
`--core-ply` supplies. `scripts/ring_registration.py` extracts that polygon,
performs the clip, and checks the result. A ring is not written without
registration unless `--allow-unregistered` is passed.

Registration is necessary because a ring triangle retained on the basis of its
centroid, lying outside a nominal circle of radius `R_core`, can have its innermost
vertex inside `R_core` by up to its centroid-to-vertex distance (100-130 m at LOLA
ring resolution). The coarse ring surface then overlaps the terrain of the core in
plan view. At lunar polar sun elevations of 0-5 deg, the overlapping edge casts a
hard, nearly permanent spurious shadow on the edge faces of the core. At lunar polar
site A an unregistered ring covered 2,730 core faces. The fraction of severely cold
faces was 0.165 under the overlap against 0.031 on clear edge faces, and no severely
cold edge face was more than 200 m from the overlap. Unregistered rings at other
polar sites had overlaps of similar size (97 m at lunar polar site B, 88 m at lunar
polar site C, and 78 m for a wider ring at site A).

Registration replaces the overlap with a small gap. The clip removes whole
triangles, so the inner boundary of the ring lies a fraction of one ring facet
outside the rim of the core (at lunar polar site A, 0.4 m minimum / 46 m median /
167 m maximum). This band contains no occluder geometry, so an edge face of the core
can be slightly under-shadowed at grazing sun. The gap is much the smaller error.
An overlap produces a hard permanent shadow with an effect of 50-300 K, while the
gap omits ~100 m of nearly flat LOLA terrain. The check prints the gap as a function
of azimuth so that it can be inspected.

Clipping a ring again does not require the view factors to be rebuilt, because the
ring is shadow-only geometry and does not enter the view-factor matrix.

## Limitations

The ring only casts shadows. Thermal emission and reflected sunlight from it are not
modeled. This omission is not small. On the floor of Haworth crater nearly all of the
terrain infrared flux comes from sunlit terrain 10–60 km away (Hayne et al., in
revision). The core mesh must therefore include all warm terrain that the region of
interest can see. The ring should represent only terrain whose effect is shadowing.

A study that requires IR from distant terrain would need the ring to be added to the
radiation exchange as a one-way source, which is not implemented.

## Verification of a ring

1. Occlusion is active. A run with a ring that contains a known massif shadows core
   faces that a run without the ring leaves lit. Check the shadow geometry against
   the Sun vector at a known epoch.
2. No change without a ring. A run with no `occluder_file` gives results identical
   to those of a run on the core alone.
3. Ring elevations are correct. `scripts/audit_ring_horizon.py` compares the
   horizon profile of the ring against the source DEM. A ring built with the wrong
   elevation scale has a horizon that is systematically too high and over-shadows
   the core at grazing sun.
4. The inner hole is registered. Run the registration check on any ring before using it:

   ```
   python3 scripts/ring_registration.py --core core.ply --ring ring.ply \
       --site-lat <lat> --site-lon <lon>
   ```

   It exits with a non-zero status if any ring geometry lies inside the outer
   boundary of the core, and prints the gap as a function of azimuth in either
   case. A ring that fails this test produces a discontinuous cold rim on the edge
   faces of the core, which checks of mesh quality, coverage or face count do not
   detect.

## References

Hayne, P. O., et al. (in revision). Advanced 3D illumination and thermal modeling for planetary surfaces: Application to the Moon3D test dataset. *Journal of Geophysical Research: Planets*.

# Flux cubes

A flux cube is a virtual sensor placed on the terrain of a site. At every output step it reports
the irradiance incident on six oriented faces, separated into components, together with the
temperature of the mesh face beneath it. The cube is a sensor and not a body, which has the
following consequences.

- It does not occlude, emit, or scatter radiation.
- It is not a row or a column of the view-factor matrix and does not enter the radiosity
  solution.
- Cubes do not see each other, so cube-to-cube coupling is zero.
- Every reported quantity is an incident irradiance in W/m^2. No albedo, emissivity, or cube
  material is applied. Users apply their own optical model to the incident fluxes.

Enabling flux cubes does not change the temperature or flux fields of the mesh.

---

## 1. Cube geometry

By default the height of the cube center and the cube edge are defined relative to the local
facet size of the mesh.

```
  local_median_edge = median over CORE faces whose centroid lies inside the grid disk
                      of the equivalent edge  e = sqrt(4A/sqrt(3))
                      (the edge of an equilateral triangle of the same area A)

  height_m = height_factor * local_median_edge     (default height_factor = 0.1)
  edge_m   = edge_factor   * local_median_edge     (default edge_factor   = 0.1)
```

`height_m` is a plumb-line height along the site vertical `U`, as given by a leveled lander leg.
It is not a height along the local face normal. The two differ by `h(1 - cos slope)`, which is
less than a millimeter at the 4 degree median inner slope of lunar polar site A.

The geometry tool computes `local_median_edge` from the mesh, prints it, and writes it into the
file header, so that the height and edge of the cubes can be derived again from the file.

Absolute values take precedence over the factors. When `height_m` or `edge_m` is given explicitly
(CLI `--height-m` / `--edge-m`), the file header records `height_factor = 0`.

> Coarse meshes need the absolute values. The factor is a fraction of the local facet size, so
> the cube is close to the ground only where the mesh is fine. Lunar polar site A has a
> high-resolution SfS core (median equivalent edge 1.083 m in the grid region) and resolves to
> h = 0.1083 m. Lunar polar site C has only the 80 m/px LOLA base and an 8,746-face core, where
> a height of 0.1 times the facet edge is meters above the ground. Set `height_m` explicitly at
> such sites.

The values at lunar polar site A (2,397,890-face core, 724,593 faces inside the 500 m grid disk)
are as follows.

| quantity | value |
|---|---|
| local median equivalent edge | 1.082969 m |
| resolved height | 0.108297 m |
| resolved edge | 0.108297 m |
| cubes placed (25 m spacing, 500 m radius) | 1,257 |
| frame max tilt across the grid | 0.016438 deg |

---

## 2. The largest view factor, `f_max`

`f_max` is the largest single view-factor entry in the row of a patch, i.e., the fraction of the
terrain flux incident on that cube face that comes from one mesh facet. `n_faces_50pct` is the
number of partner facets needed to reach half of the terrain weight.

At a height of 0.1 times the local facet edge, a large `f_max` on the nadir face is expected. A
cube ~11 cm above a ~1.1 m facet has a nadir face whose view factor to the facet beneath it is
close to unity.

The distribution of `f_max` at lunar polar site A, with `h = 0.1083 m` and 65,536 rays per patch,
is given below.

| face | min | p05 | p50 | p95 | max |
|---|---|---|---|---|---|
| zenith | 0.0000 | 0.0000 | 0.0001 | 0.0021 | 0.0350 |
| **nadir** | 0.2107 | 0.5013 | **0.8521** | 0.9765 | 0.9949 |
| north | 0.0713 | 0.1841 | 0.3207 | 0.4350 | 0.4895 |
| south | 0.0885 | 0.1980 | 0.3486 | 0.4673 | 0.5600 |
| east | 0.0814 | 0.2019 | 0.3474 | 0.4893 | 0.6465 |
| west | 0.0837 | 0.1785 | 0.3146 | 0.4326 | 0.5046 |

The nadir `n_faces_50pct` is 1 at both the median and the minimum. A warning threshold of
`f_max > 0.25` would be exceeded by 99.6% of nadir patches and 69.4% of all patches, so it would
carry no information. The tool therefore reports `f_max` and `n_faces_50pct` for every cube and
issues no warning when the height is set by `height_factor`. The values should be inspected when
an absolute `height_m` places a cube close to the ground.

### Near-field error bound

The error of the nadir irradiance depends on the cube height. At `h = 1 m` the nadir face
averages over ~6 facets, whereas at `h = 0.1083 m` it sees mostly one. Lowering the cube reduces
the numerical error and increases the sensitivity to the temperature of a single facet.

The ray-traced row integrates the facet geometry exactly, so the number of facets in view
introduces no error. The remaining error comes from the variation of radiosity within a facet,
which the model does not resolve. At the median `f_max = 0.852` the neighboring facets carry the
remaining 0.148, and the relative facet-to-facet spread of `B_ir` within 1 m is 2.7% at this
site. The reported nadir irradiance therefore differs from the radiosity of the underlying facet
by

```
  0.148 x 2.7%  ~  0.40%
```

Because ~85% of the nadir irradiance comes from one facet, the nadir face carries the full
temperature error of that facet. At a height of 1 m the same error would be averaged over ~6
facets. The facet-to-facet temperature spread within 1 m is 1.51 K at this site. The nadir
reading at 11 cm is thus a nearly exact sample of one facet, and it retains the facet-scale
scatter of the thermal solution.

We do not add a near-field model below the facet scale, because the thermal solution contains no
temperature information at that scale.

---

## 3. Sampling and the ray budget

Rays are drawn with a cosine weighting about each patch normal, using a stratified Shirley-Chiu
concentric mapping with a Cranley-Patterson rotation for each patch. Every ray then carries
equal weight and

```
  F[c][p][j] = (rays hitting core face j) / N
  f_sky      = (rays escaping)            / N
  f_ring     = (rays hitting the ring)    / N
```

so that `sum_j F + f_sky + f_ring = 1` exactly. For Lambertian emitters of radiosity `B_j` the
irradiance on the patch is `E = sum_j F[j] B_j`.

The rotation is applied independently to each patch. Without it every patch would reuse the same
sample pattern, whose largest error lies along the 45/135 degree diagonal. With 7,542 patches on
a regular 25 m lattice, a shared pattern would produce coherent stripes in the flux map instead
of uncorrelated noise. The rotation adds no cost.

### Accuracy metric

`closure_residual` is not a measure of sampling accuracy. The estimator sums to 1 by
construction, so the residual is zero on every patch at every ray budget. It is retained in the
file as a bookkeeping check, labeled as such, which detects a ray channel that has been dropped
or counted twice.

Sampling accuracy is given by the split-half total variation between two independently rotated
half-rows,

```
  splithalf_tv = 0.5 * sum_j |F_a[j] - F_b[j]|
```

which bounds the disagreement in flux directly, `|E_a - E_b| <= splithalf_tv * max_j B_j`.

### Ray budget at lunar polar site A

The tables below give the convergence and cost for the full grid at lunar polar site A on eight
threads. The first lists the median split-half total variation for each face orientation at
three ray budgets.

| face | 4,096 | 16,384 | 65,536 |
|---|---|---|---|
| zenith | 1.95e-3 | 1.65e-3 | 1.30e-3 |
| nadir | 8.30e-3 | 4.70e-3 | 2.78e-3 |
| north | 3.39e-2 | 2.47e-2 | 1.77e-2 |
| south | 2.93e-2 | 1.94e-2 | 1.27e-2 |
| east | 3.69e-2 | 2.61e-2 | 1.86e-2 |
| west | 2.81e-2 | 1.89e-2 | 1.21e-2 |

The second lists the cost. Peak memory was 0.62 GB in all three cases.

| n_rays | rays | wall | cache | nnz |
|---|---|---|---|---|
| 4,096 | 30.9 M | 6 s | 6.0 MB | 737,488 |
| 16,384 | 123.6 M | 14 s | 16.5 MB | 2,113,257 |
| 65,536 | 494.5 M | 46 s | 47.2 MB | 6,131,014 |

The default is `n_rays = 65536`. The throughput is 11.9-20 M rays/s on eight threads. The rays
are short and all originate inside one 1 km patch, so the acceleration structure remains in the
processor cache, and the throughput is higher than for whole-body illumination. Because the
whole calculation takes less than a minute, the factor of ~1.4 in accuracy on the side faces,
which converge most slowly, is obtained at little cost.

### Truncation

`truncate_rel` defaults to 0, which keeps every partner face. Truncation can introduce a large
flux error. Dropping a weight `w_res` changes the irradiance by `|dE| <= w_res * max_j B_j`, and
`B_ir` reaches `sigma*(363.5 K)^4 ~ 990 W/m^2` at this site. A dropped weight of `w_res = 0.002`
then corresponds to 2.0 W/m^2, which is larger than the entire Earthshine flux. When
`truncate_rel` is set above zero, the tool prints the dropped weight and the resulting bound for
each face orientation.

---

## 4. The occluder ring

Rays that hit the occluder ring are counted in a separate channel, `f_ring`. They are not added
to `f_sky`. The ring has no thermal solution and therefore no radiosity, and counting it as sky
would remove a source of infrared flux from the budget without a record.

For a 500 m grid inside a 3 km core, `f_ring` is small but not zero. The values at lunar polar
site A with the 27.3 km ring are given below.

| face | median f_ring | max over all patches |
|---|---|---|
| zenith | 0.0002 | |
| nadir | 0.0000 | |
| north / south / east / west | 0.0008 - 0.0064 | **0.0397** |

The ring is visible because the cubes are low. The cosine lobe of a side patch points at the
horizon, so rays that clear the local terrain travel ~2.5 km to the rim of the core and reach the
ring beyond it. Up to ~4% of the terrain infrared flux on a side face is therefore absent from
the product. For ring terrain at ~100 K this is of order 0.1-0.2 W/m^2, which is larger than the
infrared flux from the primary body at this site (0.010-0.047 W/m^2).

The geometry tool prints `f_ring` each time it runs. The missing flux can be accepted and
documented, or removed by extending the thermal core or by assigning a radiosity to the ring.

---

## 5. Tooling

```
  ptm_flux_cube_geom     one-time, CPU-only geometry build  -> <site>.fcg
      no GPU lock, no OpenCL context; safe alongside production GPU windows.
      DO NOT overlap a view-factor assemble step (those peak above 300 GB RSS).
```

The simulation loads a `.fcg` file and stops with an error if any key differs from the run (core
geometry hash, core face count, ring hash, grid parameters, ray count, sampler, seed). It does
not rebuild the geometry itself, because a hemisphere sweep over 7,542 patches would add an
unplanned cost to a long run.

The geometry tool stops with an error in the following cases.

- A lattice point has no terrain below it. Dropping the cube would change the cube count between
  runs of the same configuration.
- A lattice point falls on the occluder ring, where there is no `T_face` and no radiosity.
- The site is within 0.1 degrees of the spin pole, where local north is undefined and the face
  names would be arbitrary.
- `--ring` is given but the file cannot be read. The ring changes every sky fraction in the file.

### Validation

`scripts/flux_cube_closed_form_check.py` builds analytic meshes and compares the geometry with
closed-form solutions, using the binomial error of the estimator as the tolerance.

- Infinite flat plane, `h = 0.1 m`. The exact `f_sky` is 1 (zenith), 0 (nadir), and 1/2 (each
  side). The computed values are 1.00000 / 0.00000 / 0.50013-0.50015, against a 3-sigma tolerance
  of 0.01172.
- Finite square plane, `L = 10 m`, `h = 5 m`. The cosine-weighted nadir row sum must fall between
  the values for the inscribed and circumscribed disks, `0.8000 <= sum_j F <= 0.8889`. The
  computed value is 0.83130. A sampler uniform over the hemisphere would give 0.5528, so this
  case tests the cosine weighting.

---

## 6. The flux-cube pass

Flux cubes are evaluated once, after the radiation coupling loop has converged and immediately
before the location export. They are not evaluated inside the loop, for three reasons. First,
only the converged solution is a product, and the loop is the most expensive part of the
calculation. Second, inside the loop the band series are Picard-blended and `B_ir` carries an
emission term that changes with `T` between iterations, so a blended `B` does not correspond to
any single state. Third, after the loop the CPU and multi-GPU solvers have written the same
series, so one evaluation serves both.

### Reconstruction of the radiosity

The converged radiosity vectors are not retained after the radiosity solution. The pass
reconstructs them from the stored fluxes by inverting the two closing identities of the solver.

```
  B_ir [j] = eps_j*sigma*T_j^4 + (1 - eps_j) * ( Q_ind_ir [j] / eps_j    )
  B_vis[j] = A_src_j * Q_dir_j +      A_j    * ( Q_ind_vis[j] / (1 - A_j) )
```

The inversion is exact. It is implemented once, in
`src/cpp/radiation/radiosity_reconstruction.h`, and used by both the pass and its test. The
reconstructed `B` agrees with the solver's own values to 7.7e-17 relative.

`A_src(t)` is not stored. The pass recomputes it with `applyIncidenceAlbedo()` from quantities
that remain available after the solution, which avoids storing `n_faces x n_steps` values.

### Error conditions

Each of the following conditions stops the run with an error message that gives the reason.

| condition | reason |
|---|---|
| cache built for a different mesh | a cache from another mesh would give wrong view factors. The message prints the command that regenerates it |
| `equilibration_steps + n_steps != total_time_steps` | cube rows would not correspond to the lines of `time_steps.txt` |
| gray path (no band-separated series) | a single mixed radiosity cannot be separated into visible and infrared components |
| no ray-occlusion function | without it every cube would be treated as unshadowed, and shadowed cubes would receive direct solar flux |
| `eps` or `(1-A)` collapsing | the inversion is ill-posed when the emissivity or the absorptivity approaches zero |

A failure of the pass stops the run, so a run with flux cubes enabled does not finish without
the product.

### Units

`sun_positions` and `primary_positions` are in kilometers. The mesh is in meters for lunar and
Europa DEMs and in kilometers for asteroid shape models. `PassInputs::mesh_to_km` carries the
conversion, which is set by the same magnitude test that the illumination model uses. With the
wrong conversion no error is raised. The direction to the primary points into the ground, every
disk sample is occluded, and the flux from the primary body is exactly zero. The pass prints a
warning whenever the primary is in the ephemeris and no cube sees it at any epoch.

### Open meshes and the body sphere

An open mesh has no geometry below its own terrain, so a cube ray aimed below the limb that
misses the mesh would be counted as clear sky. On open meshes every cube ray is therefore tested
first against the body sphere of radius R_occ, with the same test that the model applies to
direct sunlight (`IlluminationModel::probeBodySphereBlocks`), and then against the mesh. This
applies to the Sun probes, the solar disk samples and the primary-disk samples. Closed meshes
skip the sphere test. Where the mesh already blocks every ray that the sphere blocks, the
product is unchanged.

Each pass also runs a check for every site (each named location, or `lattice`). From every
cube it casts a mesh-only ray toward the edge of the spherical cap, `acos(R_occ/|x|)` below the
local horizontal, at 72 azimuths. It then finds by bisection the mesh horizon of the cube with
the lowest horizon. Both results are written to the log and to `provenance.txt`, together with
the number of cube rays blocked by the sphere.

```
flux_cube_body_cull[<stem>]=on (R_occ = ...)
flux_cube_horizon_vs_cap[<stem>][<site>]=N cubes; cap depression a..b deg; cubes with a
    mesh-open azimuth at the cap edge (72 azimuths): k; worst cube min mesh horizon h deg ...
flux_cube_body_cull_effect[<stem>]=n cube ray(s) blocked by the body sphere
```

`k = 0` means that the mesh horizon is above the cap in every scanned direction, so that no ray
can pass between the mesh and the sphere. The count of blocked rays includes rays that the mesh
would also have blocked. These records are not written to the `.fcb` file.

## 7. Configuration

```yaml
flux_cubes:
  enabled: false            # default OFF; ON must be byte-neutral to every mesh field
  cache_file: ''            # REQUIRED when enabled; built by ptm_flux_cube_geom
  write_binary: true        # flux_cubes/flux_cubes.fcb
  write_csv: true           # flux_cubes/csv/cube_XXXXX.csv
  write_csv_max_cubes: 0    # 0 = every cube; cap it on long production windows
  shadow_probe_rays: 5      # 1 = center only; 5 = center + four rim points
```

The block is at the top level of the configuration and is not nested under `output:`. The log
names any key in the block that the parser does not use. The geometry keys (`spacing_m`,
`height_factor`, `n_rays`, ...) are accepted but have no effect, and the log reports this. Their
values are taken from the header of the cache file.

## 8. Outputs

| file | producer |
|---|---|
| `flux_cubes/flux_cubes.fcb` | engine, float64 flux `[cube][face][component][step]` |
| `flux_cubes/csv/cube_NNNNN.csv` | engine, one row per output step, 36 flux columns |
| `<site>_flux_cubes.h5` | `scripts/ptm_flux_cubes_to_hdf5.py` (h5py, `europa_viz`) |

The HDF5 file is written by a Python script, because the C++ code is not linked against HDF5 by
default.

The `--qa` mode of the packaging script applies four checks. It requires
`total == sum(components)`, no NaN or negative irradiance, closure of each patch, and
`nadir terrain_ir >= eps*sigma*T_face^4 * f_max`. The last check detects an incorrect mapping
between cubes and mesh faces, which the first three cannot.

### Compatibility with `skip_time_series`

The chunked runner (`scripts/run_chunked_simulation.py`) always sets
`output.skip_time_series = True`. Flux cubes read the band series held in memory and not the
per-step frames on disk, so the two are compatible. For the crater80N configuration the `.fcb`
file is identical (md5 `0133a0e073a3162ee21ce2dc38066b59`) with `skip_time_series` true and
false, while the number of per-step frames changes from 118 to 0.

The pass does not combine chunks. Under the chunked runner each chunk runs its own pass over its
own part of the window and writes its own `flux_cubes.fcb`, and `ptm_merge_chunks` does not merge
these files. Run a window that needs cubes as a single simulation, or treat the files of the
individual chunks as separate products.

## 9. Validation status

| check | result |
|---|---|
| radiosity reconstruction vs the solver's own `B` | 7.7e-17 / 9.4e-17 relative (bar 1e-12) |
| closed forms, infinite plane | `f_sky` 1.00000 / 0.00000 / 0.50014, 3-sigma bar 0.00586 |
| cosine weighting, analytic bracket | `sum_j F` = 0.83080 in [0.8000, 0.8889]; uniform would give 0.5528 |
| flush limit, h = 1.0 / 0.1 / 0.01 m | `f_max` 0.1307 -> 0.9095 -> 0.9990, monotone |
| pass unit test (3 cases, 5 negative controls) | 0 failures; nadir `terrain_ir` = 904.784902 W/m^2 at 360 K |
| Ingersoll bowl crater | PASSED, 0.492713% |
| mesh fields with cubes on and off | identical |
| body sphere, unit test | Sun 5 deg below the horizontal, past a cap 2.56 deg down: east face 0 (S cos 5 deg without the sphere); Sun 1 deg down: identical with and without the sphere |
| body sphere, crater80N | cubes on vs off: 855 surface files identical; with vs without the sphere test: surface files identical and all three `.fcb` identical apart from the header; 0 open azimuths, 37,715-37,720 sphere-blocked cube rays per height, all of which the mesh also blocks |

The most important check is that the cubes do not perturb the simulation. For the 12,800-face
crater80N configuration, the 862 physics output files have the same aggregate md5,
`cab1be996c71c52543c006822e6a9deb`, in four cases. These are the model compiled without the flux-cube code, the model with
no `flux_cubes` block, with `enabled: false`, and with `enabled: true`.

This test covers one configuration, the Moon as a built-in body with the CPU `SparseVFOperator`
solver, two-band radiosity, and a single run without chunks. It has not been repeated for
Europa, generic bodies, the gray path, the multi-GPU solver, the H-matrix solver, or chunked and
warm-started runs.

---

## 10. Shadowed cubes and permanently shadowed regions

A permanently shadowed region (PSR) is dark in all seasons. A cube that receives no direct solar
flux during one evaluated window is called window-dark. This is a weaker condition, because a
cube that is dark in a winter window may be sunlit in summer.

The two categories answer different questions.

- A window-dark cube receives only indirect flux. The direct solar flux is zero, the terrain
  infrared flux falls to tens of mW/m², and the infrared flux from the primary body is not
  negligible. A cube in a PSR is in the same physical regime, so window-dark cubes test the
  calculation of indirect flux in shadow.
- A PSR is a property of the site. It is determined by the site's `psr_mask` across all seasons
  and not by any single window.

`scripts/flux_cube_strata.py` reports the first category as `window_dark` and does not label it
PSR.

### Reachability of a PSR

Before placing a lattice to sample a PSR, determine whether any lattice point can fall on one.
This requires only the site's `psr_mask` raster, which `scripts/ptm3d_stats_to_cog.py` derives
from `direct_solar_max`.

| quantity | how to measure it |
|---|---|
| PSR pixels within the lattice radius | count of mask pixels inside the radius |
| nearest PSR pixel to the site | distance from the lattice origin |
| distinct PSR patches inside the grid | connected-component labeling of the mask |
| largest patch | its area, its extent, and its centroid in site ENU |
| PSR areal fraction of the grid | PSR area over the lattice disk area |
| lattice points landing on PSR | sample the mask at every lattice point, and compare with the count expected by area |

A COG mask is 0 or 1 almost everywhere, but linear interpolation gives fractional values along
its boundary. Report both an inclusive (`>= 0.5`) and a strict (`== 1`) count.

The local east-north-up frame of the site is rotated from the axes of the map projection by the
site longitude. In a south-polar stereographic projection with `lon_0 = 0`, local north is the
radial outward direction and not projected +Y. Omitting the rotation displaces a patch centroid
by tens of meters and changes which lattice points appear to be on a PSR.

When the largest patch is smaller than the lattice spacing, the lattice may place no point on a
PSR even though PSR pixels lie inside its radius, and a finer grid over the whole site changes
this little. Reaching such a patch requires a dense array centered on it (§11).

## 11. Targeted PSR array

A dense array centered on a PSR patch can be built with `ptm_flux_cube_geom` as it is, at the
cost of one additional run. An array of 1.0 m spacing and 10 m radius centered on the largest
patch inside the grid is an example.

Such arrays are not part of the standard flux-cube product. Landed hardware and crew are not
expected to remain inside a permanently shadowed region for long, so the fluxes there are of
limited engineering interest. The reachability calculation of §10 shows whether the standard
lattice samples a PSR at a given site.

---

## 12. Named locations

Cubes can be placed at named points of given latitude and longitude, in addition to the grid or
in place of it. Outputs are identified by name, as for the `locations:` block of the thermal
model.

### Building the geometry

```
ptm_flux_cube_geom --core data/meshes/lunar_crater_85N.ply \
    --center-lat 85.0 --center-lon 0.0 \
    --spacing 25 --radius 400 \
    --locations pads.csv --out site.fcg
```

The example uses the synthetic bowl crater included in the repository (1 km across, 200 m deep,
centered at 85°N, 0°E). On real terrain, add `--ring <ring.ply>` for far-field shadowing.

`pads.csv` has one entry per line, and the header is optional.

```
name,lat_deg,lon_deg,height_m,edge_m
pad_A,85.0000,0.0000
pad_B,85.0110,0.0000,2.0,1.0
```

`pad_A` is at the center of the crater floor, and `pad_B` is on the north wall, about 330 m away.

`height_m` and `edge_m` are optional for each entry. When omitted, each takes the site-relative
value (§1). `--radius <0` disables the lattice, so that only the named locations are built. `--radius 0` gives a single cube at the origin.

Named cubes follow the lattice cubes in the same cache, so a single run covers both. They carry
the grid index `-2147483648`, which distinguishes them from lattice cubes, and a block of names
is appended after the JSON in both `.fcg` and `.fcb`. A reader that does not use the names can
parse everything before the block at the same offsets and ignore the rest.

### Declaring them in the configuration

```yaml
flux_cubes:
  enabled: true
  cache_file: site.fcg
  locations:
    - {name: pad_A, lat_deg: 85.0000, lon_deg: 0.0000}
    - {name: pad_B, lat_deg: 85.0110, lon_deg: 0.0000, height_m: 2.0, edge_m: 1.0}
```

The list in the configuration does not place any cube. Placement is done only by
`ptm_flux_cube_geom`. The simulation checks that the cache contains every name in the list and
stops with an error otherwise. This prevents a run from using a cache that was built for a
different set of locations after the YAML has been edited.

Named locations are written to `flux_cubes/csv/cube_<name>.csv`, with `name=` and
`placement=named_location` in the header. Lattice cubes are identified by index.
`write_csv_max_cubes` limits the number of lattice files but does not remove a named location.

### Consistency check

`scripts/flux_cube_location_check.py` requires that a cube placed by name at the latitude and
longitude of the site is identical to the lattice cube at the origin. The comparison covers every
`CubeRecord` field apart from the grid indices, all six `PatchRecord`s, and all view-factor
columns and values, and it is exact. The named and lattice placements use the same routines for
the direction and the site frame, so that the same coordinates give the same geometry to the
last bit.

---

## 13. Known limitations

### Zenith terrain infrared
`scripts/flux_cube_vf_parity.py` compares the upward view-factor row of a cube, obtained by ray
tracing, with the row of the model's view-factor matrix for the facet beneath it. Over 40 level
facets the median ratio is 0.978–0.99. The ~2% difference is within the sampling scatter.

The comparison requires a view-factor matrix with correctly normalized entries. A matrix from an
offline assembler that divides entries by the hit count instead of by 2 has far-field view
factors too large by up to a factor of 2, and gives a median ratio of 0.513. Flux-cube products
computed from a thermal solution that used such a matrix inherit the error in their absolute
terrain fluxes, and must be recomputed with a corrected matrix.

### Sampling budget
A single `n_rays` applies to every patch. It is chosen from the split-half total variation,
which does not involve the radiosity and converges as ~N^-0.23. Referenced to the flux, the
median sampling error is 0.01–0.04% for the nadir and side faces. For the zenith face it is
3.98% at the median and 18.97% at the 95th percentile. The zenith face sees mostly sky, and the
estimator must resolve its small terrain row, so it is under-sampled at the default budget.

Setting a separate budget for each face orientation from a flux-referenced split-half metric
would remove this limitation. It requires the geometry tool to retain both half-rows, whereas
the cache stores only their sum, and therefore a change of file format.

### Shadow-probe spread
`edge_m` defaults to 0.1 × the local facet scale, so the spatial spread of the shadow probes
depends on the mesh resolution and not on the size of the object being modeled. A fixed scale
such as 0.30 m would represent an avionics box. At lunar polar site A the default gives an edge
of 0.108 m, and the nadir face is then h − edge/2 = 0.054 m above the surface, not h. Use the
`edge_m` value of a named location (§12) when the hardware scale is known.

### Direct solar
Files whose header gives `solar_model = finite_disk` resolve the partial solar disk as the
illumination model does. Files with `point_sun` use a single shadow test toward the center of
the Sun, which can be in error by the full direct beam at shadow edges.

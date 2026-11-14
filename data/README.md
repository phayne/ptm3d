# PTM3D data directory

Small data files (locations, cached small-body ephemerides, tutorial crater
meshes, mission SPK files) are included in the repository. Large science
datasets are not stored in git. They are distributed as versioned archives and
downloaded with `scripts/fetch_data.py`:

```bash
python3 scripts/fetch_data.py --list          # show available datasets
python3 scripts/fetch_data.py --core          # everything needed by the shipped examples (about 1.4 GB)
python3 scripts/fetch_data.py --dataset moon-topo
```

The catalog, `data/data_manifest.json`, records the files, download URLs,
sizes and SHA-256 checksums of every dataset. Downloads can be resumed and are
verified against the checksums, as with `scripts/fetch_kernels.py` for SPICE
kernels.

The datasets are archived on Zenodo as *Supporting data for Planetary
Thermophysical Model 3D (PTM3D)*, concept DOI
[10.5281/zenodo.21693403](https://doi.org/10.5281/zenodo.21693403)
(v1: [10.5281/zenodo.21693404](https://doi.org/10.5281/zenodo.21693404)),
license CC-BY-4.0. Please cite the data DOI, and the source citations given in
the README of each subdirectory, together with the software.

| Subdirectory | Contents | Source |
|---|---|---|
| `kernels/` | SPICE kernels | `scripts/fetch_kernels.py` |
| `topo/moon/` | global lunar topography grids (ico3–ico9) | `fetch_data.py --dataset moon-topo` |
| `property_maps/moon/` | Diviner-derived thermophysical property maps | `fetch_data.py --dataset moon-property-maps` |
| `albedo_maps/` | Moon + Europa albedo maps | `fetch_data.py --dataset albedo-maps` |
| `dem/validation/` | Diviner-validation crater meshes + view factors | `fetch_data.py --dataset validation-meshes` |
| `dem/copernicus/` | Copernicus crater DEM crop (SLDEM2015) | `fetch_data.py --dataset copernicus-dem` |
| `meshes/` | tutorial crater meshes (in the repository) and the Moon3D testbed scan (downloaded) | repository / `fetch_data.py --dataset moon3d-mesh` |
| `spice/` | mission trajectory SPKs, generated when needed | `scripts/fetch_horizons_spk.py` |
| `horizons_cache/` | small-body ephemerides from JPL Horizons | written on first use. The three bodies stored in the repository are kept as a record (see its README) |
| `locations/` | small analysis inputs | repository |
| `asteroids/` | small-body ephemerides, shapes and spin states for the asteroid tutorial and the front-page animation (not tracked by git). The Eros shape of the animation is Gaskell, R., *Gaskell Eros Shape Model V1.1*, NASA Planetary Data System, [doi:10.26033/d0gq-9427](https://doi.org/10.26033/d0gq-9427) (CC0), derived from NEAR images | `examples/python/fetch_asteroid.py`, `examples/python/fetch_eros_shape.py` |
| `external/sutton_DEMs/` | Europa terrain meshes + view factors (Galileo SSI SOCET SET DTMs) for `configs/europa_*` | `fetch_data.py --dataset europa-terrain` |

**Data location.** Set `PTM_DATA_ROOT=/path/to/data` to keep the datasets
outside the repository. Every relative `data/...` path in a configuration file
is then resolved against it (see [docs/ENVIRONMENT.md](../docs/ENVIRONMENT.md)).
Each dataset directory has its own README with the source, citation, processing
and license of the data.

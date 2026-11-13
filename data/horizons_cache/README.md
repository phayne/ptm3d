# JPL Horizons small-body cache (`data/horizons_cache/`)

When an asteroid or comet is defined as a `generic_body` through JPL Horizons, the results of
the query are stored in this directory. Each body has a subdirectory that holds the downloaded
`ephemeris.spk`, a generated `metakernel.tm`, and the `body_config.yaml` derived from the
Horizons orbital elements. A master `index.json` lists the bodies.

The directory is written by `src/python/data_manager.py` and is not edited by hand. A body
that is not in the cache is downloaded when it is first used.

## Bodies stored in the repository

Four bodies are stored in the repository. They are 269 Justitia, 623 Chimaera, 13294 Rockox
and 101955 Bennu. They are kept as a record of the orbit and spin elements used in published
small-body simulations.

269 Justitia, 623 Chimaera and 101955 Bennu can be simulated from a new clone. Their shape
models are in the repository (see the policy below), and their `metakernel.tm` and `index.json`
entries use paths relative to the repository root.

The Justitia kernel covers 2033-12-14T23:58:51 to 2036-05-31T23:58:51 UTC (899.0 days, one
interval, body 20000269). A simulation must lie within this interval, including the
equilibration that precedes its output window.

A new clone must still provide two things.

1. The standard kernels. `$KERNELS/lsk/naif0012.tls`, `$KERNELS/pck/pck00011.tpc`
   and `$KERNELS/spk/de440s.bsp` are not part of the cache. Download them with
   `python3 scripts/fetch_kernels.py`.
2. The working directory. SPICE (`furnsh_c`) interprets the `PATH_VALUES` of a metakernel
   relative to the current working directory, so these simulations must be started from the
   repository root. The model finds the file named by `generic_body.spice_metakernel` from
   any directory, but the paths inside the metakernel are relative to the working directory.

13294 Rockox is the exception. Its `index.json` entry contains absolute paths from the
machine that generated it, and it has no shape model in the repository. It is a record only.

Do not regenerate a `metakernel.tm` that is stored in the repository.
`SPKGenerator.create_metakernel()` writes absolute paths, to avoid the 78-character limit of
SPICE on long paths, so a regenerated file is valid only on the machine that wrote it. The
metakernels in the repository were written by hand with relative paths, as the header of each
states.

## Shape-model policy

A shape model that is loaded by a configuration in the repository is kept in the repository.
Each such file has the following four parts.

1. The file itself, tracked in git.
2. A named `!` exception in `.gitignore`, with the exact file name. The `*.obj` entry of
   `.gitignore` is intended for C/C++ object files and also matches every Wavefront shape
   model. The exception is written explicitly, in preference to adding the file with
   `git add -f`, so that it remains visible and continues to apply if the file is added again.
3. A machine-readable provenance file beside it, `<stem>_provenance.json`, with the source,
   archive, version or model id, URL, md5, sha256, byte count, vertex and face counts, units,
   and the date retrieved. `gui_core.body_fetch.shape_provenance()` reads this file name.
4. A `PROVENANCE.md` in the directory of the body, which explains what the file is, what it
   supersedes, and what has and has not been verified.

The shape models in the repository are listed below.

| body | file | what it is |
|---|---|---|
| (269) Justitia | `shape.obj` | DAMIT model 4966 (asteroid id 2858), the original |
| (269) Justitia | `shape_reconstructed_20260704.obj` | the superseded 2026-07-04 hand reconstruction, geometrically identical to the original — max vertex deviation 0.000000000 km |
| (623) Chimaera | `shape.obj` | DAMIT model 9110, the original (byte-identical to a 2026-08-31 re-download) |
| (623) Chimaera | `shape_refined.obj` | 2-pass Loop subdivision of the above; **this** is what `body_config.yaml` loads, and its generator is not in this repository |
| (101955) Bennu | `g_06370mm_alt_obj_0000n00000_v021.obj` | OSIRIS-REx OLA v21 global model |
| (101955) Bennu | `bennu_radar_shape.obj` | the superseded Nolan et al. (2013) radar model |

Superseded models are kept because they are the geometry of results that have been published
from this repository, and those results can then be checked again.

All four parts are required when the shape model of a body is added. Without the provenance
files the origin of a mesh cannot be traced, and a mesh that is not tracked exists only on one
machine.

Bennu differs from the other bodies. It is not a Horizons entry and has no `body_config.yaml`,
`metakernel.tm` or `index.json` record, because Horizons does not generate an SPK for it. Its
directory holds three files, of which a simulation loads only one.

| file | status |
|---|---|
| `g_06370mm_alt_obj_0000n00000_v021.obj` | **live** — the OSIRIS-REx OLA v21 global shape model, 6.37 m sampling, 49,152 facets, km |
| `bennu_radar_shape.obj` | archived — the pre-encounter Nolan et al. (2013) radar model, superseded 2026-08-28 |
| `ephemeris.spk` | archived — JPL solution 76, superseded 2026-08-27 |

A simulation does not load either archived file. Both are kept because previous Bennu
results from this repository were computed with them, and
`docs/tutorials/notebooks/04_asteroid_spice.ipynb` names both paths.
`data/kernels/bennu.tm` loads `spk/sb-101955-118_long.bsp`, which
`scripts/fetch_kernels.py` downloads, and `examples/quickstart_asteroid_bennu.yaml` loads the
OLA shape model. `bodies/101955_Bennu_(1999_RQ36)/PROVENANCE.md` gives the checksums, the
citations, the validation of the mesh against closed-form results, and the differences between
the models.

The Bennu shape models follow the same policy as the others. They are tracked, named as
exceptions in `.gitignore`, and accompanied by the provenance file and `PROVENANCE.md`.

## Regenerating a body locally

Download the body in the usual way for a `generic_body`, with the Horizons panel of the GUI or
with `DataManager` directly. This rewrites `bodies/<body>/` and the corresponding `index.json`
entry with paths in your own checkout. `data/horizons_cache/` is located relative to the
directory that contains `CMakeLists.txt`, so no configuration is needed. Give the `cache_dir`
argument to use another location.

## `surface.albedo` is a Bond albedo

The model reads `generic_body.surface.albedo` as a Bond albedo. The absorbed solar flux is
proportional to `1 - albedo`. The JPL Small-Body Database publishes the geometric albedo p_V,
which is larger. For (101955) Bennu, p_V is 0.044 and the Bond albedo is 0.0118.

Every cached body therefore carries the converted value and the convention it follows.

```yaml
  surface:
    albedo: 0.014526                    # the BOND albedo the engine consumes
    albedo_type: bond                   # so no reader has to guess
    albedo_geometric_pv: 0.037          # what SBDB published
    albedo_slope_parameter_g: 0.15      # the G used, and whether it was assumed
    albedo_provenance: "Bond albedo A_B = p_V x q, q = 0.290 + 0.684 G ..."
```

The conversion is `A_B = p_V q` with `q = 0.290 + 0.684 G` (Bowell et al. 1989, *Asteroids II*;
Harris 1998, *Icarus* 131, 291). `scripts/fix_cached_body_albedo.py` converts an existing
cache. It leaves unchanged a file that already contains `albedo_type: bond`, so it can be run
more than once. A file with no `albedo_type` holds p_V, and the GUI converts it when the file
is read and reports the conversion.

`fix_cached_body_albedo.py` converts the cache of this repository by default and requires
`--cache-dir` for any other. Software outside PTM3D that reads `surface.albedo` as a geometric
albedo should keep its own cache, which this script does not change unless it is named.

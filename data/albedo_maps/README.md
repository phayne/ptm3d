# Albedo maps (`albedo_maps/`)

These are Bond albedo maps for the Moon and Europa. For each icosphere level there is a
per-face binary file (`Moon_Albedo_res3..9.bin`, `Europa_Albedo_res3..9.bin`), and the source
FITS maps from which they were sampled are included. The maps are used when
`illumination.use_albedo_map: true`, and their directory is set by `albedo_maps_directory`
(default `data/albedo_maps`).

The maps are not stored in the repository. Download them with the following command.

```bash
python3 scripts/fetch_data.py --dataset albedo-maps   # ~150 MB
```

## Source

- Moon: normal-albedo mapping from LRO instruments (Diviner, LROC and LOLA), converted to Bond
  albedo with the incidence-angle parametrization of Hayne et al. (2017, Appendix).
- Europa: Bond albedo maps derived from Voyager and Galileo.
- Original mission data: NASA PDS (public domain).

## Processing

`albedo_map_generator` regrids the source maps to the face centers of each icosphere. It samples
the source FITS at the latitude (north positive) and IAU longitude of each face, which is the
convention used wherever the model samples a map (`ptm_thermal_sim.cpp`, `thermal_model.cpp`,
`albedo_manager.cpp`, the DEM sampler) and in the Moon camera file, which was compared with
photographs. Each per-face binary is an 8-byte little-endian face count followed by one float64
(Bond albedo) per face.

## Registration of the lunar maps

The lunar FITS map `Moon-BondAlbedo.fits` is registered to its header and to IAU_MOON. Column 0
is at IAU longitude 0, and longitude increases eastward from 0 to 360°. The nearside maria have
albedos of ≈ 0.07–0.09, and the farside highlands and the Tycho region ≈ 0.17. The per-face
files generated from it agree with `Moon_Albedo_res9_camera.bin` to RMS ≈ 0.005.

Lunar maps distributed with previous releases had two registration errors. The FITS raster was
offset by 180° in longitude from its header, which placed the dark nearside maria on the bright
farside highlands, and the per-face files were sampled at −latitude, which exchanged the
hemispheres (RMS difference ≈ 0.031 from the camera file). If you have a local copy from a
previous release, download the maps again. The Europa maps were not affected. The topographic
grids in `data/topo/` were also not affected.

`scripts/fix_moon_albedo_lon_registration.py` applies the longitude correction to the lunar
FITS file, and `scripts/verify_moon_albedo_registration.py` verifies the registration. Both
scripts need `numpy` and `astropy` (`python3 -m pip install -r scripts/requirements.txt`),
which the model itself does not need. The per-face files are generated from the FITS map as
follows.

```bash
build/bin/albedo_map_generator data/albedo_maps/Moon-BondAlbedo.fits \
    data/albedo_maps 3 9
```

See also `docs/KNOWN_LIMITATIONS.md`.

### Registration check when a map is loaded

A map with a longitude offset or exchanged hemispheres has the same distribution of albedo
values as the correct map. Its range, its mean and the reported effective albedo `A_eff`
therefore do not reveal the error. On `moon_global_res7`, the 180° offset changed the global
mean `A_eff` by 0.0001 and the global mean `T_max` by 0.002 K, but changed individual faces by
up to 15 K.

The model therefore checks the registration spatially wherever it loads a per-face albedo
map. This is done for the absorbed flux on the CPU, for the GPU thermal solver, and for the
source albedo of the two-band radiosity calculation. The log reports the result.

```
[albedo] Lunar map registration check: maria 0.0807 < highlands 0.1597 (contrast 0.0790, ...) -- OK
```

If the maria are brighter than the highlands, the model prints a warning that names the file
and continues. It does not correct the map, because the correction belongs in the data file,
and a correction applied to every map would be wrong for Europa. The check applies only to
bodies with a known pattern of dark and bright terrain (at present the Moon), and only when
the mesh resolves the reference sites. It is not made for regional DEM meshes or other bodies.
The two values it reports cannot distinguish a longitude offset from exchanged hemispheres.
For that, run `scripts/verify_moon_albedo_registration.py`, which compares the map with the
source FITS under all four candidate transformations and names the axis. The test of the check
is `build/bin/test_albedo_registration_check`.

Each of the three loads is checked independently, because a configuration can reach the
radiosity calculation without the absorbed-flux check having run. The per-location depth
profiles use the same absorbed flux (`AbsorbedFluxComposer`), so the check covers them as well.

### Verifying a local copy

The check above is made only when a run loads a map. To verify a local copy beforehand, use
the following commands.

```bash
python3 scripts/fetch_data.py --dataset albedo-maps --verify   # exit 1 if stale
python3 scripts/verify_moon_albedo_registration.py \
    --bins-dir data/albedo_maps --topo-dir data/topo/moon --geographic
```

`--verify` computes the checksum of the downloaded archive and of every file unpacked from
it. The unpacked files are included because the `.unpacked` marker records only that the
archive was unpacked once, so files from a previous release could otherwise remain beside a
current archive.

The `--geographic` mode of the verification script needs no FITS file. The nearside is
dominated by maria and must be darker than the farside. Correct maps give ≈ 0.1165 for the
nearside and ≈ 0.1511 for the farside (res6), and maps with the 180° offset give the inverse.

## License

The maps are derived from NASA PDS public-domain data and distributed under CC-BY-4.0.

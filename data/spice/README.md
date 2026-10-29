# Mission trajectory kernels (`spice/`)

The place for small SPICE SPK (trajectory) files that a camera or flyby
configuration names, for example a spacecraft's Moon-centered trajectory. None
ship with the repository; the large generic kernels live in `data/kernels/` and
are fetched with `scripts/fetch_kernels.py`.

To make one from JPL HORIZONS:

```bash
python3 scripts/fetch_horizons_spk.py --spacecraft <NAIF id> \
    --start <UTC> --stop <UTC> --output data/spice/<mission>_moon.bsp
```

For authoritative, current mission kernels use the NAIF archives at
https://naif.jpl.nasa.gov/naif/data.html.

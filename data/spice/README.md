# Mission trajectory kernels (`spice/`)

This directory is for small SPICE SPK (trajectory) files named by a camera or flyby
configuration, for example the Moon-centered trajectory of a spacecraft. The repository
contains none. The large generic kernels are kept in `data/kernels/` and are downloaded with
`scripts/fetch_kernels.py`.

A trajectory kernel can be made from JPL HORIZONS.

```bash
python3 scripts/fetch_horizons_spk.py --spacecraft <NAIF id> \
    --start <UTC> --stop <UTC> --output data/spice/<mission>_moon.bsp
```

For the official and current kernels of a mission, use the NAIF archives at
https://naif.jpl.nasa.gov/naif/data.html.

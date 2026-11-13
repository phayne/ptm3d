# Lunar thermophysical property maps (`property_maps/moon/`)

These are the per-face maps used by the default lunar property model, for icosphere levels
3–9. They are the H-parameter derived from Diviner (the scale depth of regolith density)
`hpar_moon_icoN.bin`, the maria and highlands masks `mask_moon_icoN.bin`, and the maps of
thermal inertia at 273 K `ti273k_moon_icoN.bin`.

The maps are not stored in the repository. Download them with the following command.

```bash
python3 scripts/fetch_data.py --dataset moon-property-maps   # ~60 MB
```

## Source

- Origin: H-parameter and thermal-inertia mapping from the Diviner Lunar
  Radiometer Experiment.
- Citations: Hayne, P. O., et al. (2017). Global regolith thermophysical
  properties of the Moon from the Diviner Lunar Radiometer Experiment.
  *J. Geophys. Res. Planets*, 122, 2371–2400.
  https://doi.org/10.1002/2017JE005387 ·
  Paige, D. A., et al. (2010). *Science*, 330, 479–482.
- Original data: NASA PDS (public domain).

## Processing

The global maps derived from Diviner are sampled at each icosphere face center. Each file holds
one single-precision value per face in mesh face order, in the binary layout with a 4-byte
header described in `docs/tutorials/03_outputs.md`.

## License

The maps are derived from NASA PDS public-domain data and distributed under CC-BY-4.0.
Please cite Hayne et al. (2017) when you use these maps.

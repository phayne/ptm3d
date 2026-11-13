# Lunar global topography grids (`topo/moon/`)

These are per-face elevation grids for the global Moon icosphere meshes, used when
`mesh.use_topographic: true`. There is one binary file for each icosphere subdivision level,
from `topo_moon_ico3.bin` to `topo_moon_ico9.bin`.

The grids are not stored in the repository. Download them with the following command.

```bash
python3 scripts/fetch_data.py --dataset moon-topo          # ico3-ico9, ~613 MB
```

Grids for `ico10` and `ico11` (1.8 and 7.4 GB) exist for global simulations at very high
resolution and are distributed separately.

## Source

- Product: LRO LOLA global lunar DEM (LDEM), resampled onto the PTM3D
  icosphere face centers.
- Citation: Smith, D. E., et al. (2010). The Lunar Orbiter Laser Altimeter
  investigation on the Lunar Reconnaissance Orbiter mission. *Space Sci. Rev.*,
  150, 209–241. https://doi.org/10.1007/s11214-009-9512-y
- Original data: NASA PDS Geosciences Node (public domain).

## Processing

The LDEM elevations are sampled at each icosphere face center and stored as one
single-precision value per face in mesh face order. The binary layout, with a 4-byte header, is
that of the output statistics files (see `docs/tutorials/03_outputs.md`).

## License

The grids are derived from U.S. government (NASA PDS) public-domain data and are distributed
under CC-BY-4.0. Please cite Smith et al. (2010) when you use the lunar topography.

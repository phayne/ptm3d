# Meshes (`meshes/`)

- `lunar_crater_80N.ply` (12,800 faces) and `lunar_crater_85N.ply` (20,000 faces): synthetic
  bowl-crater meshes, 1 km across and 200 m deep, generated with `crater_mesh_generator` from
  this repository. They are stored in the repository and used by the tutorials
  ([`docs/tutorials/02_crater_radiosity.md`](../../docs/tutorials/02_crater_radiosity.md)).
- `moon3d_cartesian_double.ply`: the Moon3D testbed crater scan, a 1,698,304-face Delaunay
  triangulation of an 849,283-point Moon3D instrument point cloud. The Moon3D paper used a
  509,337-face decimated version of it, which is not distributed. This file (32 MB) is not
  stored in the repository. Download it with the following command.

```bash
python3 scripts/fetch_data.py --dataset moon3d-mesh
```

## Source

The synthetic craters have the spherical-bowl geometry of Ingersoll et al. (1992). The Moon3D
scan is from the testbed measurements of the Moon3D project (Haddad et al. 2026).

## License

The Moon3D scan is distributed under CC-BY-4.0, as archived on Zenodo. The synthetic meshes are
part of the PTM3D distribution (MIT).

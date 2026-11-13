# External datasets (`external/`)

This directory holds the datasets that the configurations in `configs/` read from
`data/external/...`. They are not stored in the repository. The examples in `examples/` do not
need any file in `data/external/`. Their data are downloaded with `scripts/fetch_kernels.py`
and `scripts/fetch_data.py --core`.

## Europa terrain meshes and view factors (`sutton_DEMs/`)

These are the five Europa terrains used by `configs/europa_*_<terrain>.yaml` and by the Europa
terrain card of the GUI. They are Agenor Linea, Cilix crater, Pwyll crater, the pits and chaos
near Rhadamanthys Linea, and the ridged terrain near Yelland Linea. Each is a triangulated mesh
(`<name>.ply`) with its view-factor matrix (`<name>.vf`).

The files are not stored in the repository. Download them with the following command.

```bash
python3 scripts/fetch_data.py --dataset europa-terrain   # ~165 MB
```

The dataset is a separate Zenodo record, listed as `europa-terrain` in
`data/data_manifest.json`, where its DOI is recorded once the record is published. The files
`DATASET_README.md` and `PROVENANCE.json`, unpacked beside the data, describe how each file was
made.

### Source and credit

The meshes are triangulations of the five SOCET SET digital terrain models
of Europa that Michael T. Bland and Donna M. Galuszka (U.S. Geological Survey,
Astrogeology Science Center) made from Galileo Solid State Imager (SSI) images,
under NASA PDART award NNH16AC37I. They are public under CC BY 4.0. Please cite the
dataset

> Bland, Michael (2022). *Europa DEMs and Orthoimages in Support of Bland et al.
> 2021 (Remote Sensing)*. figshare. Dataset.
> <https://doi.org/10.6084/m9.figshare.16904632>

and the paper it supports:

> Bland, M. T., Kirk, R. L., Galuszka, D. M., Mayer, D. P., Beyer, R. A., and
> Fergason, R. L. (2021). "How Well Do We Know Europa's Topography? An
> Evaluation of the Variability in Digital Terrain Models of Europa."
> *Remote Sensing* **13**(24), 5097. <https://doi.org/10.3390/rs13245097>

The meshes were made with `dtm_to_ply_delaunay` at full resolution, on the reference sphere of
the rasters, with the command `dtm_to_ply_delaunay <name>_SOCET.tif <name>_SOCET.ply
--body-radius 1562090`. The view factors were computed by PTM3D with the
radiation settings of `configs/europa_lat00_<terrain>.yaml`. One matrix serves
every latitude, because `mesh.latitude_override` is a rigid rotation. Keep each
mesh and its view factors together.

Read the usage note of the source dataset before publishing results from these meshes. It
states that the products "should only be used to assess and or verify the conclusions
described in Bland et al. 2021 ... as the DTMs themselves have been superseded
by updated versions that are registered to Europa's cartographic coordinate
system." A thermal simulation is not the use stated there. PTM3D uses the relief, which is what
that evaluation examined, and the absolute cartographic position is what has been superseded.
A result from these meshes should therefore be treated as a study of terrain of that type and
not of a mapped location. Use the registered versions if the position must be current.

The directory name `sutton_DEMs/` is not an attribution. It records the colleague who shared
the files, not their authors, and it is kept because the configurations in the repository use
that path.

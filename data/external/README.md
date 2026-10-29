# External datasets (`external/`)

Datasets that configurations in `configs/` read from under `data/external/...`.
They are not stored in git. Everything needed by the runnable examples in
`examples/` is covered by `scripts/fetch_kernels.py` and
`scripts/fetch_data.py --core`; no file in `data/external/` is required for them.

## Europa terrain meshes and view factors (`sutton_DEMs/`)

The five Europa terrains behind `configs/europa_*_<terrain>.yaml` and the GUI's
Europa terrain card: Agenor Linea, Cilix crater, Pwyll crater, the pits and chaos
near Rhadamanthys Linea, and the ridged terrain near Yelland Linea. Each is a
triangulated mesh (`<name>.ply`) with its view-factor matrix (`<name>.vf`).

**Not in git.** Fetch with:

```bash
python3 scripts/fetch_data.py --dataset europa-terrain   # ~165 MB
```

The dataset is its own Zenodo record (`europa-terrain` in `data/data_manifest.json`;
its DOI is recorded there once the record is published). Its `DATASET_README.md`
and `PROVENANCE.json`, unpacked beside the files, record how each file was made.

### Source and credit

The meshes are triangulations of the five SOCET SET digital terrain models
of Europa that Michael T. Bland and Donna M. Galuszka (U.S. Geological Survey,
Astrogeology Science Center) made from Galileo Solid State Imager (SSI) images,
under NASA PDART award NNH16AC37I. They are public under CC BY 4.0. Cite the
dataset:

> Bland, Michael (2022). *Europa DEMs and Orthoimages in Support of Bland et al.
> 2021 (Remote Sensing)*. figshare. Dataset.
> <https://doi.org/10.6084/m9.figshare.16904632>

and the paper it supports:

> Bland, M. T., Kirk, R. L., Galuszka, D. M., Mayer, D. P., Beyer, R. A., and
> Fergason, R. L. (2021). "How Well Do We Know Europa's Topography? An
> Evaluation of the Variability in Digital Terrain Models of Europa."
> *Remote Sensing* **13**(24), 5097. <https://doi.org/10.3390/rs13245097>

The meshes were made with PTM3D's `dtm_to_ply_delaunay` at full resolution, on
the rasters' own sphere: `dtm_to_ply_delaunay <name>_SOCET.tif <name>_SOCET.ply
--body-radius 1562090`. The view factors were computed by the engine from the
radiation settings of `configs/europa_lat00_<terrain>.yaml`. One matrix serves
every latitude, because `mesh.latitude_override` is a rigid rotation. Keep each
mesh and its view factors together.

**Read the source deposit's usage note before you publish on these.** It states
that the products "should only be used to assess and or verify the conclusions
described in Bland et al. 2021 ... as the DTMs themselves have been superseded
by updated versions that are registered to Europa's cartographic coordinate
system." A PTM3D thermal run is not that stated use. What PTM3D consumes is the
relief, which is the part the evaluation was about; the absolute cartographic
placement is the part that has been superseded. So treat a result from these
meshes as a study of a terrain of that character, not of a mapped location, and
prefer the registered versions if you need the placement to be current.

The directory name `sutton_DEMs/` is historical and is **not** an attribution:
it records the colleague who shared the files, not who made them. It is kept
because the shipped configurations reference that path.

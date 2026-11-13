# PTM3D tutorials

Work through the tutorials in order. Each runs on a CPU. Run every command from the repository
root after building PTM3D ([INSTALL.md](../../INSTALL.md)) and downloading the Moon SPICE
kernels.

```bash
python3 scripts/fetch_kernels.py --bodies moon
```

| # | Tutorial | What it covers | Time on four cores | Needs |
|---|---|---|---|---|
| 1 | [A global Moon](01_quickstart.md) | One lunar day on a 1,280-face sphere, the run log, the spin-up estimate, and the temperatures to expect | under a minute | Moon kernels (44 MB) |
| 2 | [Radiation exchange in a crater](02_crater_radiosity.md) | A polar crater with and without scattered sunlight and thermal emission between faces | about a minute | the crater mesh in the repository |
| 3 | [The output files](03_outputs.md) | Every file in an output directory and how to read it | reading only | the output of tutorials 1 and 2 |
| 4 | [The GUI](04_gui.md) | The same simulations through the desktop application | a few minutes | `pip install -r gui/requirements.txt` |
| 5 | [Lunar polar ice from a DEM](05_lunar_polar_ice.md) | A shape-from-shading DEM of the de Gerlache crater rim turned into a mesh with distant terrain, one year of temperatures, water-ice stability depth, and H<sub>2</sub>O and CO<sub>2</sub> cold traps | about 15 minutes | the DEM (1.3 GB download) and a LOLA polar DEM (115 MB) |
| 6 | [Asteroids](06_asteroids.md) | Spheres and an ellipsoid on defined and Horizons orbits, the Bennu shape model, the DAMIT shape of Eros over a flyby window and one orbit, temperatures at sites and with depth, and a thermal inertia sweep | about 25 minutes | Moon kernels (44 MB), the Bennu kernel (81 MB), and 0.3 MB from Horizons and DAMIT |
| 7 | [Europa seen from Earth in the thermal infrared](07_observing_thermal_infrared.md) | Europa through one orbit and a Jupiter eclipse, maps of surface temperature and 8 to 14 μm brightness temperature on the disk seen from Earth, sunlight reflected in the band, the disk brightness temperature through the eclipse, and a fit of thermal inertia and emissivity to synthetic observations | under two minutes for the runs and figures, about 18 minutes for the fit (eight cores) | `fetch_kernels.py --bodies europa` (1.1 GB) |

## Jupyter notebooks

The notebooks in [`notebooks/`](notebooks/) run the model, read its output and make figures.
The copies in the repository include their output. Start with 01, which introduces the reader
used by the others. Notebooks 05, 06 and 07 continue tutorials 7, 5 and 6, on Europa seen from
Earth in the thermal infrared, lunar polar ice and asteroids, and each writes its results to a
table.

```bash
python3 -m venv .venv && source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
python3 -m pip install -e .
jupyter lab docs/tutorials/notebooks/
```

| # | Notebook | What it covers | Time on four cores | Needs |
|---|---|---|---|---|
| 01 | [`01_moon_diurnal`](notebooks/01_moon_diurnal.ipynb) | The global Moon as maps and curves, the same run on LOLA topography, and changing a thermophysical property | a few minutes | Moon kernels. Section 9 also needs the `moon-topo`, `moon-property-maps` and `albedo-maps` datasets (about 700 MB) |
| 02 | [`02_europa_thermal`](notebooks/02_europa_thermal.ipynb) | Europa through a Jupiter eclipse, spin-up length, and the observed albedo map | about ten minutes | `fetch_kernels.py --bodies europa` (1.1 GB), `fetch_data.py --dataset albedo-maps` |
| 03 | [`03_asteroid_thermal`](notebooks/03_asteroid_thermal.ipynb) | A synthetic asteroid on a Keplerian orbit at two thermal inertias | about a minute | the leapseconds kernel only |
| 04 | [`04_asteroid_spice`](notebooks/04_asteroid_spice.ipynb) | (101955) Bennu on its SPICE ephemeris with the OSIRIS-REx shape model | about half an hour | `fetch_kernels.py --bodies bennu` (81 MB) |
| 05 | [`05_europa_thermal_infrared`](notebooks/05_europa_thermal_infrared.ipynb) | Europa seen from Earth at 8 to 14 μm, disk maps of surface and band brightness temperature, the disk brightness temperature through a Jupiter eclipse, and the fit of thermal inertia and emissivity of tutorial 7 | under a minute, after the runs of tutorial 7 | the outputs of tutorial 7 |
| 06 | [`06_lunar_polar_ice`](notebooks/06_lunar_polar_ice.ipynb) | Temperature maps, water and carbon dioxide cold traps, and the depth to stable water ice at the de Gerlache crater rim | under a minute, after the runs of tutorial 5 | the outputs of tutorial 5 |
| 07 | [`07_asteroid_flyby_and_orbit`](notebooks/07_asteroid_flyby_and_orbit.ipynb) | A sphere defined in the configuration and the DAMIT shape of Eros, six hours around a flyby and one whole orbit | a few minutes, plus about 10 minutes for the orbit | the leapseconds kernel and 0.3 MB from DAMIT and Horizons |

Run the notebooks from inside the repository, or set `PTM_ROOT` to the repository directory.
Science data are found under `$PTM_DATA_ROOT` if it is set, and otherwise under the
repository's `data/` directory ([data/README.md](../../data/README.md)).

## Figure style

Every figure in the tutorials and notebooks is drawn in one style, applied by a single call,
`plt = plots.styled()` (from `ptm import plots`), before the figure is made. Text is set in a
serif face at 10 to 12 points, with mathematics in the same face, and figures are drawn at the
width of a journal column (3.5 inches, `plots.COLUMN_W`) or of a full page (7.2 inches,
`plots.DOUBLE_W`), so that the text has the same printed size in every figure. Curves of one
kind, such as latitudes from the equator to the pole, take their colors in a fixed order from a
palette that remains distinct with the common forms of color blindness (Okabe and Ito, 2008),
and each also takes its own line style, so that a figure reads in grayscale. Quantities that
appear in several tutorials keep their style everywhere (`plots.STYLES`): with radiation
exchange a solid vermillion line and without it a dashed blue line, the band brightness
temperature a solid purple line and the surface temperature a dashed green line, a model a
black line, observations open circles, and a fitted model orange. Temperature maps use the
`inferno` color map, other quantities of one sign `viridis`, and differences the diverging
`RdBu_r` map centered on zero. Panels are labeled (a), (b), ... above their upper left corners,
and legends are placed beside or above the data. `plots.save` writes line figures at 300 dots
per inch and shaded maps at 200.

The style is built on the `planetary` style of
[pretty-plots](https://github.com/phayne/pretty-plots) when that package is installed, and
otherwise sets the same values itself, so the figures look the same either way. pretty-plots
needs Python 3.10 or newer and is installed from its repository; the package of the same name
on PyPI is a different one.

```bash
python3 -m pip install "pretty-plots @ git+https://github.com/phayne/pretty-plots"
```

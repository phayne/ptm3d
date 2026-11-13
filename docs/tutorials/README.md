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
| 6 | Asteroids (`06_asteroids.md`) | Asteroids with and without shape models, over a flyby window and over an orbit | (in preparation) | |
| 7 | Observing Europa (`07_observing_europa.md`) | Europa's surface temperatures at an observation time and the disk-integrated emission seen by a telescope | (in preparation) | |

## Jupyter notebooks

The notebooks in [`notebooks/`](notebooks/) run the model, read its output and make figures.
The copies in the repository include their output. Start with 01, which introduces the reader
used by the others.

```bash
python3 -m venv .venv && source .venv/bin/activate
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

Run the notebooks from inside the repository, or set `PTM_ROOT` to the repository directory.
Science data are found under `$PTM_DATA_ROOT` if it is set, and otherwise under the
repository's `data/` directory ([data/README.md](../../data/README.md)).

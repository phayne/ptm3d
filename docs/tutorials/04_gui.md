# Tutorial 4: the PTM3D desktop application

**Application:** `gui/run_gui.py`
**Time:** a few minutes, CPU only

The desktop application builds a configuration, runs the model and displays the results. It
writes the same YAML files the command line reads, so a configuration made in the application
can be run from the command line and the other way round. Build PTM3D first
([INSTALL.md](../../INSTALL.md)).

## Starting the application

```bash
python3 -m pip install -r gui/requirements.txt
cd gui
python3 run_gui.py
```

`python3 run_gui.py --project <file.yaml>` opens a configuration, and
`python3 run_gui.py --results <output directory>` opens a finished run.

The application opens on a choice of body.

![Body selection](images/gui/01_landing_bodies.png)

The Moon, Europa and Mercury have built-in properties. **Asteroid** and **Custom body** describe
any other body. **Open a configuration file…** loads an existing YAML file, and
**All settings (advanced)…** shows every option the model has.

## A global Moon, one lunar day

Choose **Moon**. The run configuration for the Moon opens.

![Moon run configuration](images/gui/02_moon_card.png)

1. Under **Model domain**, choose **Global**, a sphere of triangular faces.
2. Under **Run type**, choose **Time series**, which writes a temperature map at every output time
   as well as the statistics.
3. Set **Subdivision level** to 3 (1,280 faces), the size of [tutorial 1](01_quickstart.md).
4. Under the time window, **Start (UTC)**, **Duration** and **Output samples** set the output
   times. One lunar day with 144 samples is about 5 hours per sample.
5. Under **How long to settle**, **Spin-up (diurnal cycles)** sets the number of lunar days
   integrated before the output window.
6. Set the **Output directory**.

The **YAML Preview** below the form is the configuration the run will use. Change a field and the
corresponding key changes in the preview.

Press **Run** or F5. **Simulation > Check Configuration** (Ctrl+T) checks the configuration
without running it.

## The run log

The console at the bottom shows the model's log. Warnings are yellow and errors red.

![Console at the end of a run](images/gui/03_console_completed.png)

When the run finishes, the console reports `Simulation completed successfully` and the application
opens the output directory in the **Results Browser**. A run that ends with an error does not
open it. Before using the results, check the end of the log for the number of faces and time
steps and for `Statistics export completed successfully`.

## Viewing the results

![Results Browser](images/gui/04_results_browser.png)

The **Results Browser** lists the fields the run wrote and the time of every output step.
Clicking a step shows it in the 3D view.

![3D view](images/gui/05_3d_view.png)

The **3D Visualization** tab colors the mesh by the field chosen under **Field**, either a
statistic such as the maximum temperature or the value at the chosen **Timestep**. **Play** steps
through the output times. With **Select Face Mode**, clicking a face reports its values.

![Time series of one face](images/gui/06_2d_time_series.png)

The **2D Plots** tab shows histograms, temperature against latitude, and the time series of the
face you selected. **Export Data** writes the plotted values to a text file.

## Radiation exchange in a crater

Use **Open a configuration file…** on the body screen, or **File > Open**, to load
`examples/tutorial04_crater80N_coupled.yaml` from [tutorial 2](02_crater_radiosity.md). Then choose
**All settings (advanced)…** and the **Radiation** page.

![Radiation page](images/gui/07_radiation_page.png)

**Enable Radiation Exchange** turns on the exchange of scattered sunlight and thermal emission
between faces. The page also holds the view-factor file and the coupling settings. A box at the
bottom estimates the size of the view-factor matrix for the current mesh, and the application
repeats the estimate before a run that has to compute view factors.

## An asteroid

Choose **Asteroid** on the body screen.

![Asteroid configuration](images/gui/08_asteroid_card.png)

Type a name or number under **Body identity** and press **Look up** to find the body at JPL
Horizons. **Fill in properties (SBDB)** takes its rotation period, spin pole and albedo from the JPL
Small-Body Database, **Fetch kernel from HORIZONS** downloads its ephemeris, and
**DAMIT shape model…** loads a shape model from the DAMIT database. The surface properties of an
asteroid are set directly on this page.

## Saving

**File > Save** (Ctrl+S) writes the configuration as YAML. Saving changes only the keys you
edited and keeps the comments in the file. Each run also writes the complete configuration it
used as `config.yaml` in its output directory. Keep that file with the results.

## Troubleshooting

- **The engine is not found.** The application looks for `build/bin/ptm_thermal_sim` in the
  repository. Build it as described in [INSTALL.md](../../INSTALL.md).
- **The run stops at once with a SPICE error.** Download the kernels, for example
  `python3 scripts/fetch_kernels.py --bodies moon`.
- **The 3D view is empty.** Install `pyvista`, `pyvistaqt` and `vtk`. The rest of the application
  works without them.
- **GPU options.** They apply only to a build with GPU support. On a CPU build, leave them off.

# PTM3D GUI

The PTM3D graphical interface is a PyQt6 application for setting up simulations, running the
model and displaying the results.

## Starting the GUI

Start the GUI with `python run_gui.py` from the `gui/` directory. This starts the thermal
simulation application (`ptm3d_gui.py`).

```bash
cd gui/
python run_gui.py
```

An optional launcher, `run_dashboard.py`, opens the same application together with the DEM and
mesh tools. It has been tested less than `run_gui.py`.

Before a simulation starts, the GUI checks the input files (`gui_core/preflight.py`). Every
path must exist and be readable. A path in a cloud-synced folder (`/Library/CloudStorage/...`)
requires an explicit override, and a file that the cloud service has removed from the local
disk blocks the simulation, because the model would wait indefinitely to read it. Keep the
input data of a simulation on a local disk.

## The interface

### 1. Simulation Configuration

A simulation is set up in a sequence of choices, which are the body, the run type, a card of
settings, and the run. A second view lists every setting. It is hidden until it is opened, and
both views edit the same configuration.

- Start page. The start page has five buttons, for the Moon, Europa, Mercury, Asteroid and
  Custom body. The three built-in bodies are those the model implements, and the list is
  generated from `planetary_properties.h`. The tooltip of each button describes what it
  opens. Where little data is available for a body, its card states this. The Mercury card
  opens with a banner that names the datasets not provided for Mercury. "Open a configuration
  file" and "Reproduce a published run" are the two other starting points.
- Model domain and run type. The domain is Global (an icosphere) or Local mesh (a terrain mesh
  you supply). The run type is Time series, Integrated statistics or, for the Europa terrain
  cases, Parameter grid.
- Card. One card then shows the settings for that body and run type. It has eleven controls
  for a global Moon run, nine for a Europa terrain site and eight for Mercury. Its sections
  are Time window, the mesh, the surface, Radiation exchange, Thermal and Output. The solver
  settings, `thermal.use_gpu` among them, are in Thermal, and Output holds only the location
  of the output files. Each card lists the keys it has set and every choice made
  automatically, together with any keys of an opened file that the card does not show.
- All settings (advanced). This is the complete reference view, with the sections Simulation,
  Mesh, Properties, Thermal, Output, Body, Illumination, Radiation, GPU and Engine Flags. Open
  it from the button on any card, from the start page, or from
  *View > All Settings (Advanced)* (Ctrl+Shift+A). Leave it with the Back button above the
  section list, or with *View > Back to the Body Cards* (Ctrl+Shift+B), which returns to the
  card from which it was opened. The card reads the configuration again on return, so the two
  views always agree. Every section is always listed, and "Advanced options" shows the expert
  groups inside each section.

  A saved configuration does not depend on whether this view is open. Every section is loaded
  and saved in either case.
- YAML Preview. The preview shows the configuration that a run passes to the model. Save does
  not write this file. It edits your project file in place, changes only the keys you changed
  and keeps its comments (this needs `ruamel.yaml`, which is listed in `requirements.txt`).
- Validation. The GUI makes its own checks. In addition, *Simulation > Check Configuration*
  passes the configuration to `ptm_thermal_sim --validate-config` and shows the result.

The cards are organized by body because most simulations differ from the defaults in only a
few keys, and a card shows the controls for those keys.

A file can be opened if it is YAML and names a body under `simulation.body`. Everything else
in it is loaded unchanged, listed above the card, and can be edited under All settings.
*File > Open* always opens the card of the body named in the file, and closes All settings if
it was open.

Each button has a label of one short line, and the detail is in its tooltip. The Mercury card
states that no albedo map, topography or worked example is provided for Mercury.
`src/python/tests/test_gui_tooltip_coverage.py` checks every text string of the application
for internal abbreviations.

### 2. 3D Visualization (PyVista)
- Interactive 3D mesh rendering with zoom, pan, and rotate
- Multiple data fields: Temperature, Illumination, Solar flux
- Face picking with detailed information
- Timestep animation with playback controls
- Export high-resolution images
- Coordinate readout under the canvas, following the pointer

#### The coordinate readout

The line under the canvas names the facet under the pointer and gives its planetocentric
latitude and east-positive longitude. These are the coordinates written to `mesh.ply` by the
simulation, so they are in the body-fixed frame of the run.

The readout gives no number when a coordinate would have no meaning.

| What is under the pointer | What the readout says |
|---|---|
| A facet of a body-centered mesh whose radius matches a known body | `lat … ° S   lon … ° E · Moon, planetocentric, east-positive longitude · facet N` |
| A facet of a body-centered mesh of an unrecognized body (a shape model) | the same, labeled `body-fixed model frame`; latitude is real, and longitude is east-positive from whatever the shape model calls zero, which is not an IAU-ratified prime meridian |
| A local tangent-plane patch, built about its own center | `lat / lon: unavailable — this mesh has no georeferenced frame` |
| Empty space | `lat / lon: — (pointer is off the mesh)` |
| A mesh just opened, pointer not yet moved over the canvas | `lat / lon: — (waiting for the pointer)` |
| The pointer handler itself failed | `lat / lon: unavailable — the pointer handler failed (see the terminal)` |

#### Diagnosing a readout or Select Face Mode that does not follow the pointer

Run the viewer with the pointer record turned on, open the results as usual, move the pointer
across the mesh, then quit.

```bash
PTM_GUI_TRACE=1 python run_gui.py
```

The record is written to `~/ptm3d_pointer_trace.log` and to the terminal. `PTM_GUI_TRACE` also
turns on the flux-cube record. Each `pointer-move` line gives the position received from the
interactor, the display position used for the pick, the facet found, and the facet that the
uncorrected position would have found. Together these distinguish an error in the coordinate
transformation from an error in the pick. The `environment-*` lines name the interactor to
which the observers are attached, the size of the render window in VTK display pixels, the
logical size and device pixel ratio of the canvas, and whether the mesh carries the per-face
coordinate arrays.

### 3. 2D Analysis (Matplotlib)
- Temperature distribution histograms
- Temperature vs. latitude scatter plots
- Time series analysis
- Min/max/mean statistics over time
- Individual face time series on selection

### 4. Results Browser
- Browse simulation output directories
- View available data fields and timesteps
- Load and compare multiple results
- Export to various formats (VTK, NumPy, HDF5, NetCDF, GeoTIFF, COG)

#### Cloud Optimized GeoTIFF export

"Export Data" offers COG beside GeoTIFF, HDF5 and VTU/PVD. It writes Cloud Optimized GeoTIFFs
with `src/python/ptm_cog.py`, which `scripts/ptm3d_stats_to_cog.py` also uses. The files are
single-band Float32 in the GDAL COG layout, with internal tiles of 512 px, LZW compression and
overviews by averaging. One file is written for each statistic under `cog/statistics/`. For a
time series one file is written for each time step under `cog/<field>/`, because a COG holds
one band.

The projection is selected from the mesh. It is polar stereographic when the whole mesh is
beyond 70° latitude, and equirectangular otherwise. The body radius is the median vertex radius
of the mesh. For the polar projection there is one option, "Shared CRS: central meridian 0°".
With a shared central meridian, rasters of different sites overlay one another in a GIS. The
dialog states which meridian will be used.

The resolution is taken from the median face spacing and rounded to a fixed table whose finest
entry is 0.001° (about 30 m on the Moon). On a finer mesh the dialog reports the spacing of the
mesh. Switch the unit selector to m/px to export at that scale.

A raster has one nodata value, and PTM3D writes two special values with different meanings.
`NaN` means that the facet was not evaluated. `+inf` means that the Sun was still not visible
at the top of the shadow-depth search, which is a result. Both become nodata in the main band.
The `+inf` facets are also written beside it as `<field>_unresolved.tif`, the fraction of the
facets behind each pixel that were unresolved, and they are excluded from the main band so that
an unresolved area is not read as a finite depth. Every band carries `PTM_FACES_TOTAL`,
`PTM_FACES_FINITE`, `PTM_FACES_NAN`, `PTM_FACES_POS_INF`, `PTM_FACES_NEG_INF` and
`PTM_SENTINEL_POLICY`, and the README written with the export describes the same convention.

### 5. Console Output
- Real-time simulation progress
- Color-coded messages (info, success, warning, error)
- Auto-scroll with toggle
- Exportable log files

## Installation

### Prerequisites
- Python 3.9 or later
- PTM3D C++ simulation built and available

### Install Dependencies

```bash
cd gui/
pip install -r requirements.txt
```

### Required packages
- `PyQt6` - GUI framework
- `numpy` - Numerical operations
- `PyYAML` - Configuration file handling
- `ruamel.yaml` - Saving a project without losing its comments or layout
  (without it the GUI prints a warning and the saved file loses them)

### Recommended packages

These are listed in `requirements.txt`. The GUI starts without them, with the corresponding view
disabled.

- `pyvista` + `pyvistaqt` - 3D visualization
- `matplotlib` - 2D plotting

### Optional packages
- `pyopencl` - GPU detection
- `h5py` - HDF5 export
- `netCDF4` - NetCDF export
- `rasterio` - GeoTIFF export of results
- `gdal` (`osgeo`) and a system GDAL on PATH (`gdaladdo`, `gdal_translate`) - COG export of
  results. Both are needed. The band is written with the Python bindings and the overviews
  with the command-line tools, because only the command-line tools of GDAL >= 3.12 record
  `OVERVIEW_RESAMPLING=AVERAGE` in the `IMAGE_STRUCTURE` metadata of the COG. Without them
  the COG option is disabled, and the dialog gives the reason.
- `gdal` (`osgeo`) - GeoTIFF DEM preview in the DEM module of the dashboard. It cannot be
  installed with pip alone. Use `conda install -c conda-forge gdal`. The same GDAL
  installation is needed to build the DTM tools that the DEM module runs (see INSTALL.md,
  "DEM / DTM tools").

## Usage

### Launch the GUI

```bash
# From the gui/ directory
python run_gui.py

# Or with a project file
python run_gui.py --project /path/to/config.ptm3d

# Or with results to visualize
python run_gui.py --results /path/to/output/
```

### A first simulation

1. Configure Simulation
   - On the "Configuration" tab, select a body on the start page, then a
     run type. The card that opens holds the settings below
   - Set simulation parameters (duration, timesteps)
   - Configure mesh resolution
   - Set thermal and illumination options
   - Choose output directory

2. Save and Run
   - Save project (Ctrl+S)
   - Click "Run" button or press F5
   - Monitor progress in Console

3. Visualize Results
   - Results load automatically after simulation
   - Switch between 3D view and 2D plots
   - Use timestep slider or Play button for animation
   - Click faces in 3D view for details

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| Ctrl+N | New project |
| Ctrl+O | Open project |
| Ctrl+S | Save project |
| Ctrl+L | Load results |
| F5 | Run simulation |
| Shift+F5 | Stop simulation |
| F6 | Refresh visualization |
| Home | Reset camera |
| Ctrl+` | Toggle console |
| Ctrl+E | Export image |

## Project Structure

```
gui/
├── ptm3d_gui.py           # Main application (thermal simulation)
├── run_gui.py             # Launcher script — the supported entry point
├── dashboard.py           # Multi-module dashboard shell
├── run_dashboard.py       # Dashboard launcher (optional, less tested)
├── requirements.txt       # Python dependencies
├── README.md              # This file
│
├── gui_core/              # Core modules
│   ├── __init__.py
│   ├── config_manager.py  # Configuration handling
│   ├── styles.py          # Modern dark theme
│   └── simulation_runner.py # Background execution
│
├── gui_panels/            # UI panels
│   ├── __init__.py
│   ├── config_panel.py    # Configuration tabs
│   ├── visualization_panel.py  # 3D/2D views
│   ├── results_panel.py   # Results browser
│   └── console_panel.py   # Console output
│
├── dashboard_core/        # Dashboard shell: module registry + chaining
│   ├── module_interface.py  # Base class every dashboard module implements
│   ├── module_registry.py   # Discovery/registration of modules
│   ├── chaining_manager.py  # Hands one module's output to the next
│   └── sidebar.py           # Module navigation
│
├── module_adapters/       # Adapters wrapping each app as a dashboard module
│
├── dem_core/              # DEM module: logic (no Qt widgets)
│   ├── dem_config_manager.py # DEM tool settings
│   ├── dem_runner.py         # Runs dtm_to_ply_delaunay
│   └── dem_tooltips.py       # Flag-level help text
│
└── dem_panels/            # DEM module: UI panels (Preview / Convert / Inspect)
    ├── dem_preview_panel.py  # GeoTIFF preview + ROI picking (needs osgeo/GDAL)
    ├── dem_config_panel.py   # Nav rail + Conversion / ROI / QC options
    ├── dem_mesh_panel.py     # Mesh statistics and quality analysis
    └── dem_results_panel.py  # Generated-mesh summary
```

### DEM module prerequisites

`dem_runner.py` runs one compiled program built from this repository, `dtm_to_ply_delaunay`,
which converts a DTM to a PLY mesh, extracts a region of interest and grades the point density
radially. It is searched for in `build/bin/`, then in a `../dtm/` directory, then on `PATH`.
The program requires GDAL when it is built. Without GDAL, CMake does not build it, and the
Convert action reports that the program is missing and gives the commands to build it.

Two limitations apply to this module.

- There is no mesh-refinement step. Use the Convert page with ROI enabled, which grades the
  point density radially (Core Fraction and Density Exponent).
- The Convert action works on macOS and Linux only. It reads the output of the converter with
  `select.select()` on pipes, which Windows does not support, and on Windows the module
  reports that the action is unavailable.

The conversion, its options and its limitations are described in
[`docs/dem-to-mesh-workflow.md`](../docs/dem-to-mesh-workflow.md).

## Configuration File Format

Projects are saved as YAML files (`.ptm3d` or `.yaml`).

```yaml
simulation:
  start_time: "2030-04-10T00:00:00"
  duration_days: 3.55
  timesteps: 144
  body: Moon

mesh:
  resolution: 5
  import_file: ""
  scale_factor: 1.0

thermal:
  run_thermal_mesh: true
  threads: 0  # auto
  equilibration_orbits: 2
  use_gpu: false

radiation:
  enabled: false
  use_hmatrix: false

output:
  directory: output/my_simulation
  format: binary
  vtk_mode: both
```

## Troubleshooting

### PyVista not showing 3D view
- Ensure VTK is properly installed: `pip install vtk`
- Try: `pip install --upgrade pyvista pyvistaqt`
- On macOS, you may need: `pip install pyobjc`

### High DPI scaling issues
- Set environment variable: `QT_ENABLE_HIGHDPI_SCALING=1`
- Or: `QT_SCALE_FACTOR=1.5`

### Simulation executable not found
- Ensure PTM3D is built: `cd build && make`
- Check the build path in simulation_runner.py

### Memory issues with large meshes
- Reduce mesh resolution
- Disable statistics computation
- Use H-matrix compression for radiation

## Development

### Adding New Configuration Options

1. Add default value in `gui_core/config_manager.py`
2. Add widget in appropriate tab in `gui_panels/config_panel.py`
3. Map to C++ command-line flag in `gui_core/simulation_runner.py`

### Adding New Visualization Features

1. Extend `PyVistaWidget` in `gui_panels/visualization_panel.py`
2. Add controls to the visualization panel
3. Connect signals to main window

## Building a distribution package

To package the model for distribution, use `git archive`, which applies the `export-ignore`
entries of the repository's `.gitattributes`.

```bash
git archive --format=tar.gz -o ptm3d.tar.gz HEAD
```

This leaves out the per-developer settings (`.claude/`, `gui/.claude/`). Copying the working
tree with `tar` or `cp -r` does not apply these entries and would include those directories.

Some input datasets (SPICE kernels, DEMs, view factors) are large and are not included in the
package. The top-level documentation describes how to download them. Keep the input data on a
local disk. If the files of a cloud-synced folder have been removed from the local disk, the
check made before a simulation starts blocks the run.

## License

Part of the PTM3D simulation package.

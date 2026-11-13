# Installing PTM3D

PTM3D is a planetary thermophysical model written in C++. The model runs on CPUs and uses all
available cores. A GPU with OpenCL can be used to accelerate it where suitable hardware exists,
but is not required. This guide describes how to build the model from source on Linux and
macOS. On Windows, use WSL2 (Ubuntu) and follow the Linux instructions.

> The minimal CPU-only build on Linux and macOS is
> ```bash
> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPTM_ENABLE_GPU=OFF \
>       -DPTM_ARCH_NATIVE=OFF -DCSPICE_ROOT=/path/to/cspice
> cmake --build build --target ptm_thermal_sim -j
> ./build/bin/ptm_thermal_sim --help
> ```

## 1. Requirements

| Component | Minimum | Notes |
|-----------|---------|-------|
| CMake | 3.10+ | |
| C++ compiler | C++17 (GCC 10+, Clang 12+, AppleClang) | native Windows (MSVC) is not supported; use WSL2 |
| Eigen3 | 3.3+ | header-only |
| Embree | 4.x | ray tracing for illumination/view factors; stock packages are SSE2-only (see §3 "Embree ISA") |
| CSPICE | N0067+ | NASA NAIF SPICE toolkit (C); not in apt/brew, see below |
| CFITSIO | any | camera/exosphere I/O (can be made optional) |
| yaml-cpp | any | config parsing |
| OpenMP | any | threaded CPU core (`libomp` on macOS) |
| OpenCL | optional | only for `PTM_ENABLE_GPU=ON` |

A laptop with 16 GB of memory is sufficient for meshes of 1,000 to 100,000 faces, such as
global models at low resolution and asteroids. Meshes of 500,000 to 1 million faces, such as
small craters and sites from high-resolution DEMs, need a workstation or a cluster. The memory
a simulation needs can be estimated before it is run.
```bash
python3 scripts/ptm_mem_estimate.py --faces 100000 --steps 360   # -> ~0.80 GiB, OK on 16 GB
```

## 2. Dependencies per platform

### Linux (Debian/Ubuntu)
```bash
sudo apt-get update
sudo apt-get install -y cmake g++ libomp-dev libeigen3-dev libcfitsio-dev \
                        libyaml-cpp-dev libembree-dev pkg-config
```
If `libembree-dev` provides Embree 3 on your system, install Embree 4 from the official release
archive and give its location to CMake with `-Dembree_DIR=` or with
`EMBREE_INCLUDE_DIRS`/`EMBREE_LIBRARIES`.

### macOS (Homebrew)
```bash
brew install cmake embree eigen cfitsio yaml-cpp libomp pkg-config
```
The Homebrew `embree` package is built for SSE2 only, which makes the illumination calculation
somewhat slower on a CPU with AVX or AVX-512. Section 3, "Embree ISA (performance)", describes
how to check which version is installed and how to build one for several instruction sets
beside it.

### DEM / DTM tools (optional, requires GDAL)

The thermal model does not need GDAL. GDAL is needed only to build the two tools that convert
a GeoTIFF digital terrain model into a PTM3D mesh (`dtm_to_ply_delaunay` and
`dem_mesh_generator`), and for the `osgeo` bindings used by the DEM preview of the GUI
dashboard. This section can be skipped if you do not start from a DEM.

```bash
# macOS / Linux, recommended (installs the C++ library AND the Python bindings):
conda install -c conda-forge gdal

# macOS, Homebrew alternative (C++ library only; add the bindings separately):
brew install gdal
python3 -m pip install "GDAL==$(gdal-config --version)"

# Debian/Ubuntu:
sudo apt-get install -y libgdal-dev python3-gdal
```

CMake finds GDAL through its config package or through `gdal-config`. A conda environment is
found automatically when `CONDA_PREFIX` is set. Without GDAL the configuration still succeeds
and prints

```
CMake Warning: GDAL not properly configured. SKIPPING these targets: dem_mesh_generator,
dtm_to_ply_delaunay. ...
```

and everything else is built. If you install GDAL later, run `cmake -S . -B build` again
before building, because the targets are selected when CMake configures the build.

Then build the DEM tool (see section 3).

```bash
cmake --build build --target dtm_to_ply_delaunay -j
```

`dtm_to_ply_delaunay` also uses the header-only [CDT](https://github.com/artem-ogre/CDT)
triangulation library, included in `third_party/CDT/` (MPL-2.0), so there is nothing to install.

The conversion of a DEM to a mesh, with the meaning of each command-line option and the
functions that are not available, is described in
[docs/dem-to-mesh-workflow.md](docs/dem-to-mesh-workflow.md).

### CSPICE (both platforms, manual install)

CSPICE is distributed by NASA NAIF and not by a package manager. It is provided as a prebuilt
static library, with one package for each CPU architecture, so the package must match the
machine on which you build.

| Build host | NAIF package directory | Download |
|---|---|---|
| Linux, x86_64 | `PC_Linux_GCC_64bit` | ~41 MB compressed |
| macOS, Apple Silicon (`uname -m` reports `arm64`) | `MacM1_OSX_clang_64bit` | ~40 MB compressed |
| macOS, Intel (`uname -m` reports `x86_64`) | `MacIntel_OSX_AppleC_64bit` | ~40 MB compressed |

The following commands select and unpack the matching package. Run them in a directory outside
the PTM3D repository, for example `~/software`. They take a few seconds on a fast connection
and leave about 120 MB in `cspice/`.

```bash
case "$(uname -s)/$(uname -m)" in
  Linux/*)       NAIF_PKG=PC_Linux_GCC_64bit        ;;
  Darwin/arm64)  NAIF_PKG=MacM1_OSX_clang_64bit     ;;   # Apple Silicon (M1/M2/M3/M4)
  Darwin/x86_64) NAIF_PKG=MacIntel_OSX_AppleC_64bit ;;   # Intel Mac
  *) echo "No NAIF package mapping for $(uname -s)/$(uname -m)" ;;
esac
URL="https://naif.jpl.nasa.gov/pub/naif/toolkit//C/${NAIF_PKG}/packages/cspice.tar.Z"
curl -sSL "$URL" -o cspice.tar.Z && gzip -d cspice.tar.Z && tar xf cspice.tar && rm cspice.tar
# CMake's find_library wants libcspice.a; NAIF ships lib/cspice.a:
cp cspice/lib/cspice.a cspice/lib/libcspice.a
export CSPICE_ROOT="$PWD/cspice"
```

The configure commands below pass `$CSPICE_ROOT` to CMake. In a new shell, set it again before
running CMake, or add the `export` line with the full path to your shell profile.

On Apple Silicon, verify the architecture of the library before building. CMake's
`find_library` checks only that `libcspice.a` exists, and not the architecture for which it was
compiled. With the Intel package (`MacIntel_OSX_AppleC_64bit`) an arm64 build therefore
configures without error and then fails at the link step. The linker ignores an archive built
for another architecture, and every SPICE symbol (`_furnsh_c`, `_spkezr_c`, …) is reported as
undefined for `arm64`. The following command checks the archive.

```bash
lipo -info "$CSPICE_ROOT/lib/libcspice.a"
# -> ... is architecture: arm64     correct on an Apple Silicon Mac
# -> ... is architecture: x86_64    wrong package — re-download MacM1_OSX_clang_64bit
```

On an Intel Mac the correct architecture is `x86_64`. Running an arm64 build under Rosetta is
not supported, so build for the native architecture.

CMake looks for CSPICE in `-DCSPICE_ROOT=`, then in `$CSPICE_ROOT`/`$CSPICE_DIR`, then in
common installation prefixes. If `libcspice.a` is not found, the configuration stops with an
error message that names the variable to set.

## 3. Configure and build

### Portable CPU core (recommended default)
```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPTM_ENABLE_GPU=OFF \
  -DPTM_ARCH_NATIVE=OFF \
  -DCSPICE_ROOT="$CSPICE_ROOT"
cmake --build build --target ptm_thermal_sim -j
```
`-DPTM_ARCH_NATIVE=OFF` builds portable code without `-march=native`. Omit it (the default is
ON) to optimize for the CPU of the build machine. `PTM_ENABLE_GPU=OFF` is the default.

### GPU-accelerated build (OpenCL-capable hosts)
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPTM_ENABLE_GPU=ON -DCSPICE_ROOT="$CSPICE_ROOT"
cmake --build build --target ptm_thermal_sim ptm_gpu_selfcheck -j
# One-time per host: profile the safe GPU dispatch limits:
./build/bin/ptm_thermal_sim --calibrate-gpu     # writes ~/.ptm3d/machine_profile.json
```
On macOS, `-DENABLE_GPU=ON` is accepted as an alias for `PTM_ENABLE_GPU`. On Apple Silicon use
the CPU build.

#### Checking GPU illumination on your hardware

GPU illumination (`illumination.use_gpu: true`) was validated against CPU Embree on an
AMD Radeon Pro W6800X. Before using it on any other GPU or driver, run the self-check.
```bash
ptm gpu-selfcheck                 # or: ./build/bin/ptm_gpu_selfcheck
ptm gpu-selfcheck --device 1      # one GPU, by its --list-gpu-devices index
```
The self-check needs no data files. It builds three synthetic crater meshes in memory, two with
facets of about 20 m and 2 mm and one of 259,200 triangles with a deep bounding-volume
hierarchy. On each GPU it traces about 41,600 shadow rays of the kinds used in a simulation,
with the same ray-tracing routine, and compares every result with a brute-force calculation in
float64. Rays within a small tolerance of a tie, in which a ray passes exactly through an edge,
are skipped and counted, since either result is correct for them. The check takes a few
seconds per GPU and pauses 500 ms between dispatches (`--pause-ms`), as in the validation runs.

The exit status is `0` for PASS on every GPU and `1` for FAIL, which means that a GPU gave a
wrong result and `use_gpu: false` should be kept on that machine. It is `2` for SKIP, when
there is no GPU build or no device and nothing was checked, and `3` for INVALID, when the
check itself could not be relied on, for example because the reference compiled on the host
disagreed with the brute-force calculation.

#### GPU checklist before a long run (macOS)

A long GPU run on macOS with AMD hardware can hang the GPU, and the machine may then restart.
Three precautions apply on each machine.

1. Disable sleep and display sleep. A transition to display or system sleep resets the GPU
   display mode, which is known to trigger the hang.
   ```bash
   sudo pmset -a disablesleep 1 displaysleep 0     # before a campaign
   sudo pmset -a disablesleep 0 displaysleep 30    # after
   ```
   At startup the model reads `pmset -g custom` and prints a warning if sleep is enabled. It
   does not change any setting.
2. Calibrate the GPU dispatch limits (`--calibrate-gpu`, above). Without a calibration,
   conservative limits are used and the log reports this.
3. Use the chunked runner for any long simulation. It stops a run that has stalled and keeps
   the log of the failed attempt when the run is repeated.
   ```bash
   python3 scripts/run_chunked_simulation.py --config <config> \
       --total-days <days> --output-dir <dir>
   ```

Large simulations with radiation exchange are handled differently when GPU memory is limited.
When the copy of the view-factor matrix leaves less than 2 GiB of free memory on each GPU, the
radiosity calculation and the GPU thermal solver are placed on separate devices. The
`[Radiosity] ACTIVE RUNG` line of the log reports this, and a 30-minute run takes roughly
8 minutes longer. GPU radiosity is not used at all when less than 0.25 GiB would remain free.
Both limits can be set under `gpu:` and apply only on macOS with AMD hardware.

A hung run can be recognized because the log stops in the middle of a phase while the process
remains alive. After the restart, `pmset -g log | grep -i "shutdown cause"` reports a negative
SMC code (`-108`/`-128`). `pmset -g log` keeps about a week of history and does not need sudo,
unlike `log show`. The recovery procedure and the reason for each precaution, including which
of them apply only on macOS with AMD hardware, are given in
[docs/dev/gpu-wedge-guardrails.md](docs/dev/gpu-wedge-guardrails.md).

### Embree ISA (performance)

The Embree libraries provided by package managers are usually compiled for SSE2 only. On a
CPU with AVX or AVX-512, the same Embree version compiled for several instruction sets makes
the illumination calculation a few percent faster. It can also change the result of a small
number of occlusion tests at grazing incidence, which changes temperatures by hundredths of a
kelvin. Every run prints the Embree instruction set it uses and records it in
`provenance.txt`. [docs/dev/embree-isa.md](docs/dev/embree-isa.md) describes how to check the
installed library and how to build Embree for several instruction sets.

## 4. SPICE kernel data

The model is linked with the CSPICE library (section 2). To run, it also needs SPICE kernel
data, which are ephemerides, orientation and leapseconds. These are large binary files that
are not kept in the repository. They are downloaded from NASA NAIF with a Python script that
needs no additional packages.

```bash
python3 scripts/fetch_kernels.py                  # default (auto): ~1.23 GB — see note
python3 scripts/fetch_kernels.py --bodies moon    # Moon + Mercury core only, ~44 MB
python3 scripts/fetch_kernels.py --bodies europa  # core + jup365.bsp (1.1 GB) for Europa
python3 scripts/fetch_kernels.py --bodies bennu   # core + Bennu's ephemeris (80.7 MB)
```

Without options the script downloads ~1.23 GB. In its default mode, `--bodies auto`, it
downloads the core set and the kernels required by the metakernels in `data/kernels/`. Two of
the metakernels in the repository name a body kernel. `europa.tm` requires the Galilean
ephemeris `jup365.bsp` (1.1 GB), and `bennu.tm` requires `sb-101955-118_long.bsp` (80.7 MB).
The script lists each of these additions before downloading it. For the Moon and Mercury
alone, pass `--bodies moon`, which downloads `naif0012.tls` (5 KB), `pck00011.tpc` (128 KB),
`de440s.bsp` (31.2 MB), `moon_pa_de440_200625.bpc` (12.3 MB) and `moon_de440_250416.tf`
(19 KB), ~44 MB in total.

`sb-101955-118_long.bsp` is the only kernel here that does not come from NAIF. It is JPL
solution 118 for (101955) Bennu, the final OSIRIS-REx reconstruction (Farnocchia et al., 2021),
and covers 1900-01-01 to 2135-10-01. It is downloaded from `https://ssd.jpl.nasa.gov/ftp/xfr/`,
the address that JPL Horizons gives for Bennu in place of generating an SPK on request. This
address is a transfer directory and not a versioned archive, so the SHA-256 checksum recorded
in the script is the only assurance that the file has not changed. If the checksum does not
match, find out whether JPL has replaced the file before assuming that the download is corrupt.

After the download, `data/kernels/` has the following layout.

```
data/kernels/
├── moon.tm  europa.tm  mercury.tm  bennu.tm   # metakernel templates (shipped with the repo)
├── lsk/naif0012.tls                 # leapseconds (required even in Keplerian mode)
├── pck/pck00011.tpc                 # body orientation/shape constants
├── pck/moon_pa_de440_200625.bpc     # lunar principal-axes orientation
├── fk/moon_de440_250416.tf          # MOON_PA / MOON_ME frame definitions
├── spk/de440s.bsp                   # planetary ephemeris
└── spk/sb-101955-118_long.bsp       # (101955) Bennu, 1900-2135 (--bodies bennu/auto)
```

An interrupted download resumes where it stopped. Every file is verified against a SHA-256
checksum, and files that are already complete are skipped when the script is run again.

Two settings determine where the kernels are found.

- The kernel directory is set when the model is built, as the CMake variable
  `SPICE_KERNELS_DIR` (default `<source-dir>/data/kernels`). If you keep the kernels elsewhere,
  configure with `-DSPICE_KERNELS_DIR=/abs/path`.
- The metakernels in the repository use paths relative to the repository root
  (`PATH_VALUES = ('data/kernels')`), which SPICE interprets relative to the directory from
  which the program is started. Run `ptm_thermal_sim` from the repository root, as all
  examples and tutorials do, or change `PATH_VALUES` to an absolute path.

For small bodies with an ephemeris from JPL Horizons (`generic_body` configurations), the
ephemerides are stored under `data/horizons_cache/`, and only `lsk/naif0012.tls` is needed
from the kernels above.

Science datasets (lunar topography, property maps, albedo maps, validation meshes) are
downloaded in the same way with `python3 scripts/fetch_data.py --core` (about 1.4 GB).
`data/README.md` gives the catalog, and `$PTM_DATA_ROOT` in `docs/ENVIRONMENT.md` keeps the
data outside the repository. The quickstart example needs only the SPICE kernels.

## 5. Python environment (QA, notebooks, plots)

The model and its tests need no Python packages. The output checker (`scripts/ptm_qa.py`), the
tutorial notebooks and the plotting tools need the following environment.

```bash
python3 -m venv .venv && source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 -m pip install -e .          # the ptm package, used by every tutorial and notebook
python3 -m pip install scipy         # tutorial 1 and the fitting examples
```

The `requirements.txt` in the repository root contains numpy, matplotlib, PyYAML, JupyterLab,
and SpiceyPy, which tutorial 6 and notebook 04 need. The second line installs the `ptm` package
from the repository, which reads the output and provides the `ptm` command. On a Mac these
three lines take about a minute and about 450 MB of disk. Two other files serve the rest of the
repository. `scripts/requirements.txt` adds astropy and the other packages used by the analysis
and figure scripts, and `gui/requirements.txt` adds PyQt6, PyVista and the other packages of the
desktop GUI. Tutorial 5 also needs the GDAL Python bindings, described under "DEM / DTM tools"
in section 2.

## 6. Verify the build

```bash
# Quick radiation check (analytic bowl-crater view factors and reciprocity), CPU, seconds:
cmake --build build --target test_bowl_crater -j && ./build/bin/test_bowl_crater
#   -> [test_bowl_crater] PASS: 12 of 12 verdicts passed
# Spin-up substep size (not spin-up length):
cmake --build build --target test_equilibration_dt -j && ./build/bin/test_equilibration_dt
#   -> [test_equilibration_dt] PASS
# Unknown configuration keys and flags are reported:
python3 scripts/test_unknown_config_keys.py
```

A GPU build adds four GPU radiosity checks to `test_bowl_crater` (16 of 16). The tests are
standalone programs, and `ctest` is not used.

## 7. First run

Continue with the [tutorials](docs/tutorials/README.md). The first one runs the quickstart
example.
```bash
./build/bin/ptm_thermal_sim --config examples/quickstart_moon_diurnal.yaml
```

Start from the configurations in `examples/`, which are written to be run from a clone of the
repository. One of them, `europa_high_res.yaml`, needs a GPU build. The configurations in
`configs/` are the Europa terrain cases used by the GUI (five Galileo stereo DEMs at eight
latitudes) and one global Europa case. Most of them set `thermal.use_gpu: true`, and a CPU
build stops on them with `FATAL [gpu]` unless the key is set to `false`. The terrain meshes and
view factors they name are downloaded with
`python3 scripts/fetch_data.py --dataset europa-terrain` (`data/external/README.md`).

## 8. Notes

- `BUILD_TESTS=ON` refers to a top-level `tests/` directory that does not exist. Build the
  test programs individually as in section 6.
- `rebuild.sh` is a macOS developer script that always configures the same build (GPU on,
  native instruction set, shared libraries). Use section 3 for any other build.
  `./rebuild.sh --dry-run` prints the commands without running them.
- On Apple Silicon, the portable CPU build (`-DPTM_ARCH_NATIVE=OFF`) is the one built by
  continuous integration on `macos-latest` (arm64), and this release was also tested on an
  Apple M5 machine. Use the `MacM1_OSX_clang_64bit` CSPICE package (section 2).
- The GPU version has been used mainly on macOS with AMD hardware. The precautions for it are
  described in [docs/dev/gpu-wedge-guardrails.md](docs/dev/gpu-wedge-guardrails.md).

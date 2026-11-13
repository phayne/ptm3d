# SPICE kernels

PTM3D loads SPICE data from this directory. Its location is set when the model is built, as
the CMake variable `SPICE_KERNELS_DIR`, which defaults to `<repo>/data/kernels`.

The repository contains only the small metakernel templates (`moon.tm`, `europa.tm`,
`mercury.tm`, `bennu.tm`). The binary kernels they name (~44 MB for the Moon set) are not
stored in the repository. Download them with one of the following commands.

```bash
python3 scripts/fetch_kernels.py                  # default (auto): ~1.23 GB
python3 scripts/fetch_kernels.py --bodies moon    # Moon + Mercury core, ~44 MB
python3 scripts/fetch_kernels.py --bodies europa  # core + jup365.bsp (1.1 GB)
python3 scripts/fetch_kernels.py --bodies bennu   # core + Bennu SPK (80.7 MB)
```

Without options the script uses `--bodies auto`, which also downloads the kernels required
by the metakernels in this directory. `europa.tm` and `bennu.tm` are both here, so the default
also downloads `jup365.bsp` (1.1 GB) and `sb-101955-118_long.bsp` (80.7 MB), for a total of
~1.23 GB. Pass `--bodies moon` if you need only the Moon and Mercury set.

After the download the directory has the following layout.

```
data/kernels/
├── moon.tm  europa.tm  mercury.tm  bennu.tm   # metakernels (shipped)
├── lsk/naif0012.tls                 # leapseconds (required, all modes)
├── pck/pck00011.tpc                 # body orientation/shape constants
├── pck/moon_pa_de440_200625.bpc     # lunar principal-axes orientation
├── fk/moon_de440_250416.tf          # MOON_PA / MOON_ME frame definitions
└── spk/de440s.bsp                   # planetary ephemeris
    (+ spk/jup365.bsp for Europa, spk/sb-101955-118_long.bsp for Bennu)
```

All of these come from NASA NAIF except the Bennu kernel, which comes from JPL Solar System
Dynamics. The header of `bennu.tm` gives its URL, its coverage (1900-01-01 to 2135-10-01), the
orbit solution it contains, and the reason the NAIF copy of the same solution is not used.

Run the model from the repository root. The metakernels use `PATH_VALUES = ('data/kernels')`,
which SPICE interprets relative to the directory from which the program is started. To run
from another directory, change `PATH_VALUES` to the absolute path of this directory.

See `INSTALL.md`, section "SPICE kernel data".

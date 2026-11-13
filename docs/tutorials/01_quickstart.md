# Tutorial 1: a global Moon over one lunar day

**Configuration:** [`examples/quickstart_moon_diurnal.yaml`](../../examples/quickstart_moon_diurnal.yaml)
**Data:** the Moon SPICE kernels (44 MB)
**Time:** about 10 seconds on four cores, under 1 GB of memory

## What the simulation does

The simulation follows one lunar day (29.53 Earth days) over the whole Moon, represented by a
sphere of 1,280 triangular faces, each about 480 km across.

1. **Spin-up.** The regolith column beneath each face is integrated for 80 lunar days under the
   real solar forcing, sampled 48 times per day, before the output window begins.
2. **Illumination.** At each of the 240 output times, about 3 hours apart, the model finds the
   Sun's position with SPICE, traces shadow rays from every face, and adds the sunlight
   reflected by Earth and Earth's thermal emission.
3. **Heat conduction.** Each column has the lunar regolith profile, with density and
   conductivity increasing with depth and conductivity and heat capacity depending on
   temperature.
4. **Statistics.** The minimum, mean and maximum of temperature and of each flux are recorded
   for every face over the output window.

## Check the build (optional)

```bash
./build/bin/test_bowl_crater
```

This compares the view factors of a bowl-shaped crater with the analytic solution and ends with
`[test_bowl_crater] PASS: 12 of 12 verdicts passed`. It is a quick test of the radiation
exchange calculation, which tutorial 2 uses.

## Run it

```bash
python3 scripts/fetch_kernels.py --bodies moon
./build/bin/ptm_thermal_sim --config examples/quickstart_moon_diurnal.yaml
```

Look for these lines in the log.

```
Mesh created with 642 vertices and 1280 faces
Equilibration (integrated): 80.00 day/night cycles = 1 pass(es) x 80.00-cycle forcing block  [default]
    lit   (Qs>0)     n=     1280   tau1 = 9.56 cycles   E0 = 99.088 K   cycles needed: median 66.2, p95 68.8
[Statistics] Exported to output/quickstart_moon_diurnal/statistics
Simulation completed in 10.521 seconds
```

### The spin-up estimate

The `cycles needed` line estimates how many lunar days of spin-up each face needs for its
column to settle within `thermal.equilibration_target_K`, 0.1 K by default. The estimate fits
one relaxation time, `tau1`, to each column and tends to ask for more cycles than a convergence
test shows are needed. This simulation integrates 80 cycles, more than the 69 the estimate asks
for, so no warning is printed.

The surface settles much faster than the deeper column. With a spin-up of one lunar day instead
of 80, the maximum temperatures change by less than 0.25 K, but the nighttime minimum changes by
up to 4.9 K and the mean temperature by up to 2.6 K. Run with `--equilibration-orbits 1` to see
the warning and the difference.

## Check the result

The statistics files hold a 4-byte face count followed by one 4-byte float per face. This reads
them with the Python standard library.

```bash
python3 - <<'PY'
import struct, array
def load(p):
    b = open(p, 'rb').read()
    n = struct.unpack('<i', b[:4])[0]
    a = array.array('f'); a.frombytes(b[4:4 + 4 * n])
    return sorted(a)
for name in ('temperature_max', 'temperature_min', 'temperature_avg'):
    v = load(f'output/quickstart_moon_diurnal/statistics/{name}.bin')
    print(f"{name}: min={v[0]:.1f}  median={v[len(v)//2]:.1f}  max={v[-1]:.1f} K")
PY
```

The expected output is below.

```
temperature_max: min=182.0  median=369.8  max=385.0 K
temperature_min: min=66.8  median=92.6  max=95.0 K
temperature_avg: min=107.8  median=208.2  max=215.7 K
```

Each number has a simple physical explanation. The equatorial noon maximum of 385 K is the
radiative equilibrium temperature near 1 AU for a normal-incidence Bond albedo of 0.12. The
nighttime value near 93 K is the lunar nighttime plateau, set by the heat the regolith stores
during the day. The polar maximum near 182 K reflects sunlight at grazing incidence.

The output directory, about 10 MB, is described in [tutorial 3](03_outputs.md).

## Going further

- Set `mesh.resolution: 4` (5,120 faces) for a sharper map.
- Set `simulation.timesteps: 480` for a 1.5-hour output cadence.
- Open `output/quickstart_moon_diurnal` in the GUI ([tutorial 4](04_gui.md)).
- Convert the output for ParaView with
  `python src/python/vtk_converter.py output/quickstart_moon_diurnal output/qs_paraview`
  (needs the `vtk` package).

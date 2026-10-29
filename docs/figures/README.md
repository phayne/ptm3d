# Figures (`docs/figures/`)

Images used by the README and the documentation pages. Where a figure has a
generator in this repository, it is listed here; regenerate from the repository
root.

| Figure | Used by | Generator |
|---|---|---|
| `ptm3d_workflow.{png,svg}`, `architecture.{png,svg}`, `thermal_profile.{png,svg}` | not on a doc page | `docs/figures/generate_diagrams.py` |
| `v2/bowl_crater_ingersoll.png` | VALIDATION.md, radiation-exchange.md | `src/python/bowl_crater_figure.py` (reads `build/bin/test_bowl_crater` output) |
| `v2/bowl_crater_3d.png` | not on a doc page | `src/python/bowl_crater_3d_figure.py`, `src/python/bowl_crater_ptm_figure.py` |
| `v2/energy_conservation.png` | VALIDATION.md, thermal-model.md | `src/python/energy_conservation_figure.py` |
| `v2/equilibration_convergence.png` | thermal-model.md | `src/python/equilibration_figure.py` |
| `lunar_validation.png`, `v2/lunar_validation.png` | thermal-model.md | `src/python/lunar_validation_figure.py` |
| `thermal_wave_profiles.png`, `v2/thermal_wave_profiles.png` | thermal-model.md | `src/python/plot_thermal_wave.py` |
| `v2/radiation_exchange_schematic.{png,pdf}` | not on a doc page (paper figure) | `src/python/radiation_exchange_schematic.py`, from `examples/crater85N_illumination.yaml` (see the script's docstring) |
| `conductivity/*` | martinez-conductivity-model.md | `scripts/plot_conductivity_models.py` |
| `validation_fig03..09_*.png` | VALIDATION.md, USER_GUIDE.md | Manuscript figures of the Moon3D paper, shown as published |

The lunar-validation, energy-conservation, equilibration and thermal-wave
generators read a global lunar run (`output/lunar_validation*`); its
configuration is not yet shipped in `examples/`. `v2/performance_scaling.png`
(thermal-model.md) is a measurement of particular machines and is shown as
published; the README's three hero images
(`moon3d_testbed_*`, `europa_res7_*`, `moon_global_res9_*`) likewise.

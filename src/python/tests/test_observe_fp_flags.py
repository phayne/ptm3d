#!/usr/bin/env python3
"""Floating-point status checks of the band and disk calculations of ptm.observe.

    python3 src/python/tests/test_observe_fp_flags.py

On Apple silicon, the BLAS library behind NumPy's ``@`` can leave false
divide-by-zero, overflow and invalid flags after a correct product, which
NumPy then reports as "... encountered in matmul" warnings. The band radiance
and its temperature derivative are therefore written with ``np.einsum``, which
does not call BLAS. These checks:

1. evaluate the band radiance, its derivative, the band brightness
   temperature, the face weights and the whole band calculation on a
   realistic temperature field (faces at 20 K on the night side, a day side
   in radiative equilibrium, a few faces at 400 K) with every floating-point
   exception raised as an error and every RuntimeWarning turned into an
   error, and require finite results;
2. require the band brightness temperature of the field to return the
   temperatures it was made from;
3. require that those functions, and the globe and shape plots that show
   their results, contain no ``@`` product outside ``with np.errstate(...)``,
   so the BLAS path does not come back unnoticed. The disk weights keep one
   guarded ``@``, so that the single-frequency results stay bit for bit, and
   check its result for non-finite values instead.

Check 1 passes with either form on processors whose BLAS keeps the flags
clean; check 3 is the one that fails if ``@`` returns. No SPICE kernels or
model runs are needed.
"""

import ast
import inspect
import sys
import textwrap
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src" / "python"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ptm import observe, plots  # noqa: E402
from test_observe import icosphere  # noqa: E402

FAILS = []
R_KM, D_KM = 1560.8, 6.3e8


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def field(mesh):
    """Day side in radiative equilibrium (130 K at the subsolar point), night
    side at 20 K, and the 20 faces nearest the subsolar point at 400 K."""
    c = mesh["vertices"][mesh["faces"]].mean(axis=1)
    cos_i = c[:, 0] / np.linalg.norm(c, axis=1)
    T = np.where(cos_i > 0, 130.0 * np.clip(cos_i, 0.0, None) ** 0.25, 0.0)
    T = np.maximum(T, 20.0)
    T[np.argsort(-cos_i)[:20]] = 400.0
    return T, 50.0 * np.clip(cos_i, 0.0, None)


def strict(fn, *args, **kw):
    """Call fn with every floating-point exception and RuntimeWarning raised."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with np.errstate(all="raise"):
            return fn(*args, **kw)


def test_strict_evaluation():
    print("1-2. band and disk calculations with floating-point exceptions raised")
    mesh = icosphere(5, R_KM)
    T, solar = field(mesh)
    print(f"     {T.size} faces, {T.min():.0f} to {T.max():.0f} K")
    for band in (observe.Band.top_hat(8.0, 14.0), observe.Band.top_hat(3.0, 50.0)):
        L = strict(band.radiance, T)
        check(f"{band.label}: band radiance finite and positive", np.all(np.isfinite(L))
              and np.all(L > 0))
        dL = strict(band._dradiance_dT, T)
        check(f"{band.label}: d(radiance)/dT finite and positive",
              np.all(np.isfinite(dL)) and np.all(dL > 0))
        Tb = strict(band.brightness_temperature, L)
        err = np.max(np.abs(Tb - T) / T)
        check(f"{band.label}: brightness temperature returns T", err < 1e-10,
              f"max rel. error {err:.1e}")
    band = observe.Band.top_hat(8.0, 14.0)
    for lat, lon in ((0.0, 0.0), (10.0, 120.0), (-30.0, 180.0)):
        geom = observe.explicit_geometry(lat, lon, D_KM)
        mu, domega = strict(observe.face_weights, mesh, geom, length_unit="km")
        check(f"face weights from ({lat:g}, {lon:g}) finite",
              np.all(np.isfinite(mu)) and np.all(np.isfinite(domega)) and domega.sum() > 0)
        em = strict(observe.band_emission, mesh, T, geom, band=band, emissivity=0.9,
                    direct_solar=solar, length_unit="km")
        ok = all(np.isfinite(v) for v in (em.flux_w_m2, em.thermal_flux_w_m2,
                                          em.reflected_flux_w_m2, em.tb_k, em.tb_thermal_k))
        check(f"band emission from ({lat:g}, {lon:g}) finite", ok,
              f"T_b {em.tb_k:.3f} K")


def _is_errstate(item):
    f = item.context_expr.func if isinstance(item.context_expr, ast.Call) else None
    return isinstance(f, ast.Attribute) and f.attr == "errstate"


def matmul_lines(fn):
    """Lines of ``@`` products in fn that are not inside ``with np.errstate(...)``."""
    src = textwrap.dedent(inspect.getsource(fn))
    first = inspect.getsourcelines(fn)[1]
    tree = ast.parse(src)
    guarded = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.With) and any(_is_errstate(i) for i in node.items):
            guarded.update(id(n) for n in ast.walk(node))
    return [first + n.lineno - 1 for n in ast.walk(tree)
            if isinstance(n, ast.BinOp) and isinstance(n.op, ast.MatMult)
            and id(n) not in guarded]


def test_no_blas_products():
    print("3. no unguarded @ product in the band, disk and globe functions")
    for fn in (observe.Band.radiance, observe.Band._dradiance_dT, observe._disk_weights,
               plots.globe, plots._draw_arc, plots.shape, plots.mesh_angular_radius):
        lines = matmul_lines(fn)
        check(f"{fn.__qualname__}: no unguarded @", not lines,
              f"@ at line(s) {lines}" if lines else "")


def main():
    print(f"NumPy {np.__version__}")
    test_strict_evaluation()
    test_no_blas_products()
    if FAILS:
        print(f"FAIL: {len(FAILS)} check(s): {', '.join(FAILS)}")
        return 1
    print("PASS: all checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())

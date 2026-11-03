"""When does lateral heat conduction matter for a mesh? The Python half of the advisor.

The engine prints the same advisory before every thermal run
(``src/cpp/thermal/lateral_advisor.h``; silence it with
``thermal.lateral_advisor: false``). This module answers the same question
without running the engine, so a script or the GUI can ask it before launch::

    from ptm.lateral_advisor import advise_config
    advice = advise_config("my_run.yaml")
    print(advice.verdict)            # 'RECOMMENDED', 'BORDERLINE', 'not needed', 'UNMEASURED'
    print("\\n".join(advice.lines))

or from the shell: ``ptm advise my_run.yaml``.

The criterion (Moon3D paper, Hayne et al., in revision, Section 2; after Hayne
et al. 2021): a feature of size D is conductively isolated when its own local
conductance to radiative and conductive loss,

    h = 4 eps sigma T^3 + Gamma sqrt(omega),   Gamma = sqrt(k rho c_p),

outweighs the lateral conductance across it, i.e. when D exceeds

    D* = 8 k / (pi h).

On lunar regolith D* is about 0.9 mm on warm sunlit surfaces, 6.5 mm at 110 K,
2 cm at 70 K and 4 cm at 50 K. Lateral conduction is RECOMMENDED when facets
smaller than D*(50 K) cover more than 1% of the surface area. Facet size is the
mean edge length, the measure the engine's lateral edge-activation rule uses.

Everything here mirrors the engine and is tested against the same numbers
(``src/python/tests/test_lateral_advisor.py``); the verdict rule, thresholds,
and log text are identical, so the GUI and the engine log cannot disagree. What
it does NOT mirror, stated rather than approximated silently:

* ``--thermal-inertia`` on the engine command line (it is not in the config);
* per-face thermal-property maps (a NOTE is added when one is configured);
* a generic body's ``generic_body.config_file`` (only the inline
  ``generic_body:`` block is read; a NOTE is added when a config_file is set).

GUI: not wired in yet. TODO(gui): show ``advise_config(...).verdict`` and its
first line as a "lateral conduction recommended" note on the thermal card. It
does not belong in ``gui_core.preflight``, whose contract is stat-only (it reads
no mesh data, because a cloud placeholder blocks on read); call it off the GUI
thread after the mesh preview has read the mesh.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "ADVISOR_COLD_T", "ADVISOR_ICE_T", "ADVISOR_SUNLIT_T",
    "RECOMMEND_AREA_FRACTION", "BORDERLINE_FLOOR", "BORDERLINE_FACTOR",
    "SurfaceMaterial", "Advice",
    "critical_length", "surface_conductivity", "surface_heat_capacity",
    "facet_sizes", "evaluate", "format_advice",
    "surface_material_from_config", "advise_config",
]

# Same constants as lateral_advisor.h; test_lateral_advisor.py checks they match.
ADVISOR_COLD_T = 50.0      # K, cold trap (the CO2 regime)
ADVISOR_ICE_T = 110.0      # K, water-ice stability threshold
ADVISOR_SUNLIT_T = 300.0   # K, warm sunlit surface
RECOMMEND_AREA_FRACTION = 0.01
BORDERLINE_FLOOR = 0.001
BORDERLINE_FACTOR = 2.0
ROCK_CAVEAT_FACTOR = 100.0

_SIGMA = 5.670374419e-8    # the literal ptm::lateral::conductionLength uses
_MARTINEZ = "martinez2021"


@dataclass
class SurfaceMaterial:
    """The surface node's material, as the engine's batch solver sets it at z = 0."""
    ks: float                       # W/m/K
    rhos: float                     # kg/m^3
    cp0: float                      # J/kg/K
    cp_coeff: tuple = ()            # c_p = cp0 + sum_j cp_coeff[j] T^(j+1)
    emissivity: float = 0.95
    chi: float = 2.7
    day_s: float = 29.530588 * 86400.0
    model: str = "hayne2017"
    label: str = ""

    @property
    def R350(self):
        return self.chi / 350.0 ** 3


def surface_heat_capacity(m, T):
    if not m.cp_coeff:
        return m.cp0
    r = m.cp_coeff[-1]
    for c in reversed(m.cp_coeff[:-1]):
        r = r * T + c
    return m.cp0 + T * r


# Martinez & Siegler (2021) Eqs. 9-10, transcribed from conductivity_model.h
# (test_lateral_advisor.py re-reads the constants from the header and compares).
_MARTINEZ_EQ9 = (-0.203297, -11.472, 22.5793, -14.3084, 3.41742,
                 0.01101, -2.8049e-5, 3.35837e-8, -1.40021e-11)
_MARTINEZ_EQ10 = (5.0821e-6, 5.1e-3, 2.0022e-13, 1.953e-10)
_MARTINEZ_TMIN = 2.7


def _martinez_k(T, rho):
    t = max(T, _MARTINEZ_TMIN)
    ac, bc, cc, dc, ec, fc, gc, hc, ic = _MARTINEZ_EQ9
    it = 1.0 / t
    kam = (ac + bc * it ** 4 + cc * it ** 3 + dc * it ** 2 + ec * it
           + fc * t + gc * t ** 2 + hc * t ** 3 + ic * t ** 4)
    a1, a2, b1, b2 = _MARTINEZ_EQ10
    return (a1 * rho - a2) * kam + (b1 * rho - b2) * t ** 3


def surface_conductivity(m, T):
    if m.model == _MARTINEZ:
        return _martinez_k(T, m.rhos)
    return m.ks * (1.0 + m.R350 * T * T * T)


def critical_length(m, T):
    """D*(T) = 8 k / (pi h), h = 4 eps sigma T^3 + sqrt(k rho c_p) sqrt(omega), in meters."""
    k = surface_conductivity(m, T)
    omega = 2.0 * math.pi / m.day_s
    h = (4.0 * m.emissivity * _SIGMA * T ** 3
         + math.sqrt(k * m.rhos * surface_heat_capacity(m, T)) * math.sqrt(omega))
    return 8.0 * k / (math.pi * h)


# --- geometry ------------------------------------------------------------------

def facet_sizes(vertices, faces, to_m=1.0):
    """(mean edge length [m], area [m^2]) per valid face, as numpy arrays."""
    import numpy as np
    V = np.asarray(vertices, dtype=np.float64)
    F = np.asarray(faces, dtype=np.int64)
    if V.size == 0 or F.size == 0:
        return np.zeros(0), np.zeros(0)
    ok = np.all((F >= 0) & (F < len(V)), axis=1)
    F = F[ok]
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    e = (np.linalg.norm(b - a, axis=1) + np.linalg.norm(c - a, axis=1)
         + np.linalg.norm(c - b, axis=1)) / 3.0
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    good = np.isfinite(e) & np.isfinite(area) & (area > 0.0)
    return e[good] * to_m, area[good] * to_m * to_m


# --- verdict -------------------------------------------------------------------

@dataclass
class Advice:
    verdict: str = "UNMEASURED"
    dstar_cold: float = 0.0
    dstar_ice: float = 0.0
    dstar_sunlit: float = 0.0
    frac_cold: float = 0.0
    frac_ice: float = 0.0
    frac_sunlit: float = 0.0
    frac_cold_x2: float = 0.0
    size_median: float = 0.0
    size_p05: float = 0.0
    size_p95: float = 0.0
    n_faces: int = 0
    lines: list = field(default_factory=list)

    @property
    def recommended(self):
        return self.verdict == "RECOMMENDED"

    def as_dict(self):
        return {k: getattr(self, k) for k in (
            "verdict", "dstar_cold", "dstar_ice", "dstar_sunlit", "frac_cold",
            "frac_ice", "frac_sunlit", "frac_cold_x2", "size_median",
            "size_p05", "size_p95", "n_faces", "lines")}


def _nearest_rank(sorted_values, q):
    if len(sorted_values) == 0:
        return 0.0
    return float(sorted_values[int(q * (len(sorted_values) - 1))])


def evaluate(sizes, areas, material):
    """The advisor's numbers and verdict (no text). Mirrors evaluateLateralAdvice()."""
    import numpy as np
    r = Advice()
    r.dstar_cold = critical_length(material, ADVISOR_COLD_T)
    r.dstar_ice = critical_length(material, ADVISOR_ICE_T)
    r.dstar_sunlit = critical_length(material, ADVISOR_SUNLIT_T)
    sizes = np.asarray(sizes, dtype=np.float64)
    areas = np.asarray(areas, dtype=np.float64)
    r.n_faces = int(sizes.size)
    total = float(areas.sum()) if areas.size else 0.0
    finite = all(math.isfinite(d) for d in (r.dstar_cold, r.dstar_ice, r.dstar_sunlit)) \
        and r.dstar_cold > 0.0
    if r.n_faces == 0 or not total > 0.0 or not finite:
        return r
    r.frac_cold = float(areas[sizes < r.dstar_cold].sum()) / total
    r.frac_ice = float(areas[sizes < r.dstar_ice].sum()) / total
    r.frac_sunlit = float(areas[sizes < r.dstar_sunlit].sum()) / total
    r.frac_cold_x2 = float(areas[sizes < BORDERLINE_FACTOR * r.dstar_cold].sum()) / total
    s = np.sort(sizes)
    r.size_median = _nearest_rank(s, 0.50)
    r.size_p05 = _nearest_rank(s, 0.05)
    r.size_p95 = _nearest_rank(s, 0.95)
    if r.frac_cold > RECOMMEND_AREA_FRACTION:
        r.verdict = "RECOMMENDED"
    elif r.frac_cold > BORDERLINE_FLOOR or r.frac_cold_x2 > RECOMMEND_AREA_FRACTION:
        r.verdict = "BORDERLINE"
    else:
        r.verdict = "not needed"
    return r


# --- text (identical to formatLateralAdvice) -------------------------------------

def _fmt_length(m):
    a = abs(m)
    if not math.isfinite(m):
        return "%g" % m
    if a < 0.01:
        return "%.3g mm" % (m * 1e3)
    if a < 1.0:
        return "%.3g cm" % (m * 1e2)
    if a < 1000.0:
        return "%.3g m" % m
    return "%.3g km" % (m / 1e3)


def _fmt_ratio(x):
    if not (x > 0.0) or not math.isfinite(x):
        return "%g" % x
    if x < 10.0:
        return "%.2g" % x
    scale = 10.0 ** (math.floor(math.log10(x)) - 2.0)
    return "%.0f" % (round(x / scale) * scale)


def _fmt_pct(f):
    p = 100.0 * f
    if p == 0.0:
        return "0%"
    if p < 0.1:
        return "<0.1%"
    if p < 10.0:
        return "%.1f%%" % p
    return "%.0f%%" % p


def format_advice(r, material, lateral_on=False, unit_label="meters", notes=()):
    """The advisory as log lines, each prefixed ``[lateral-advisor] ``."""
    P = "[lateral-advisor] "
    if r.verdict == "UNMEASURED":
        return [P + "WARNING: could not measure the facet sizes of this mesh (or D* is not "
                    "finite for this material), so the lateral-conduction check did NOT run."]
    L = []
    dc = _fmt_length(r.dstar_cold)
    if r.verdict == "RECOMMENDED":
        if lateral_on:
            L.append(P + "thermal.lateral_conduction is ON, and appropriate for this mesh: "
                         "%s of the surface area is in facets smaller than D*(%g K) = %s."
                     % (_fmt_pct(r.frac_cold), ADVISOR_COLD_T, dc))
        else:
            L.append(P + ("RECOMMENDED: %s of the surface area is in facets smaller than D*(%g K) = %s, so this mesh resolves "
                          "features small enough to exchange significant heat laterally. Without thermal.lateral_conduction, "
                          "small cold traps on this mesh are likely too cold, and cold-trap areas can be overestimated, most "
                          "strongly at the coldest thresholds.")
                     % (_fmt_pct(r.frac_cold), ADVISOR_COLD_T, dc))
    elif r.verdict == "BORDERLINE":
        L.append(P + "%s%s of the surface area is in facets smaller than D*(%g K) = %s, and "
                     "%s in facets smaller than %g D*. "
                     "Lateral conduction is likely a small correction for water-ice (110 K) cold traps on this mesh, but "
                     "can still remove the coldest (CO2-class) cold traps, which on lunar regolith sit at 25-45 K where "
                     "D* is 2-3 times D*(50 K); enable it if those matter."
                 % ("thermal.lateral_conduction is ON (borderline for this mesh): "
                    if lateral_on else "BORDERLINE: ",
                    _fmt_pct(r.frac_cold), ADVISOR_COLD_T, dc,
                    _fmt_pct(r.frac_cold_x2), BORDERLINE_FACTOR))
    else:
        ratio = r.size_median / r.dstar_cold
        if lateral_on:
            L.append(P + "thermal.lateral_conduction is ON but unnecessary for this mesh: the "
                         "median facet (%s) is %s times D*(%g K) = %s, so facets are "
                         "conductively isolated. It adds cost for a negligible effect."
                     % (_fmt_length(r.size_median), _fmt_ratio(ratio), ADVISOR_COLD_T, dc))
        else:
            L.append(P + "not needed: the median facet (%s) is %s times D*(%g K) = %s, so "
                         "facets are conductively isolated and exchange heat by radiation."
                     % (_fmt_length(r.size_median), _fmt_ratio(ratio), ADVISOR_COLD_T, dc))

    klaw = ("k(T) Martinez & Siegler 2021" if material.model == _MARTINEZ
            else "k(T) Hayne et al. 2017 with chi = %g" % material.chi)
    label = " (%s)" % material.label if material.label else ""
    L.append(P + "  D* = 8k/(pi h), h = 4 eps sigma T^3 + Gamma sqrt(omega), from this run's "
                 "surface material%s: %s, k_s = %.3g W/m/K, rho_s = %.4g kg/m^3, eps = %.3g, "
                 "solar day %.4g d"
             % (label, klaw, material.ks, material.rhos, material.emissivity,
                material.day_s / 86400.0))
    L.append(P + "  D* = %s at %g K (sunlit), %s at %g K (water-ice threshold), %s at %g K "
                 "(cold trap)"
             % (_fmt_length(r.dstar_sunlit), ADVISOR_SUNLIT_T, _fmt_length(r.dstar_ice),
                ADVISOR_ICE_T, dc, ADVISOR_COLD_T))
    L.append(P + "  facet size (mean edge, %d faces, mesh read as %s): median %s, "
                 "5th-95th percentile %s - %s; area in facets below D*: %s (%g K), "
                 "%s (%g K), %s (%g K)"
             % (r.n_faces, unit_label, _fmt_length(r.size_median), _fmt_length(r.size_p05),
                _fmt_length(r.size_p95), _fmt_pct(r.frac_sunlit), ADVISOR_SUNLIT_T,
                _fmt_pct(r.frac_ice), ADVISOR_ICE_T, _fmt_pct(r.frac_cold), ADVISOR_COLD_T))
    if r.verdict != "not needed":
        L.append(P + "  " +
                 "For reference (Moon3D crater, facets ~2 mm; Hayne et al., in revision): at 110 K, cold traps "
                 "smaller than 1 cm kept about half their area with lateral conduction, those of 1-5 cm about three "
                 "quarters, and the dominant region (equivalent diameter about 91 cm) 99%; facets on its margin "
                 "cooled into the class, leaving the total 110 K area nearly unchanged (+0.5%). At the 54 K CO2 "
                 "threshold every cold trap smaller than 5 cm vanished and the total area fell by 69%.")
        if not lateral_on:
            L.append(P + "  " +
                     "To enable: thermal.lateral_conduction: true, lateral_scheme: monolithic, lateral_picard: 2, "
                     "lateral_activation_c: -1 (the Moon3D configuration), use_gpu: false. Lateral conduction runs on the "
                     "CPU SIMD solver only and solves a coupled system across the mesh at every substep, so expect a "
                     "substantially longer thermal stage.")
    elif r.size_median < ROCK_CAVEAT_FACTOR * r.dstar_cold:
        L.append(P + "  NOTE: this assumes the configured regolith-like material. In rock (k ~ "
                     "0.1-1 W/m/K), D* grows to decimeters or meters, comparable to these facets.")
    for n in notes:
        L.append(P + "  NOTE: " + n)
    L.append(P + "  Silence with thermal.lateral_advisor: false. Criterion and guidance: "
                 "docs/thermal-model.md, \"When to enable lateral conduction\".")
    return L


# --- config resolution -----------------------------------------------------------

def _get(config, dotted, default=None):
    node = config
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return default if node is None else node


def surface_material_from_config(config):
    """(SurfaceMaterial, notes) resolved the way ptm_thermal_sim resolves it.

    Order, as in main(): body table (or the generic-body defaults plus the
    inline ``generic_body.surface`` block) -> ``thermophysical.ks_override`` /
    ``emissivity_override`` -> ``thermophysical.constant_thermal_inertia``
    collapse (cp = TI^2/(ks rhos), chi = 0) -> conductivity law.
    """
    from .body_table import BODY_CAPABILITIES, GENERIC_BODY_SURFACE_DEFAULTS
    notes = []
    body = str(_get(config, "simulation.body", "Moon"))
    gb = _get(config, "generic_body", None)
    if isinstance(gb, dict) and gb:
        surf = dict(GENERIC_BODY_SURFACE_DEFAULTS)
        given = gb.get("surface") or {}
        for key, positive_only in (("emissivity", True), ("ks", True), ("rhos", True),
                                   ("cp0", False), ("chi", False)):
            v = given.get(key)
            if v is not None and (v > 0 if positive_only else v >= 0):
                surf[key] = float(v)
        spin_h = float((gb.get("spin") or {}).get("period_hours", -1) or -1)
        a_au = float((gb.get("orbit") or {}).get("semi_major_axis_au", -1) or -1)
        day_s = 0.0
        if spin_h > 0:
            sid = spin_h * 3600.0
            orb = (a_au ** 1.5) * 365.25 * 86400.0 if a_au > 0 else 0.0
            day_s = (1.0 / abs(1.0 / sid - 1.0 / orb)
                     if orb > 0 and abs(sid - orb) > 1.0 else sid)
        if day_s <= 0:
            notes.append("generic body without generic_body.spin.period_hours; D* uses a 1-day "
                         "solar day (D* depends on it only through the Gamma sqrt(omega) term).")
            day_s = 86400.0
        if gb.get("config_file"):
            notes.append("generic_body.config_file is not read here; D* uses the inline "
                         "generic_body block and the engine defaults.")
        mat = SurfaceMaterial(ks=surf["ks"], rhos=surf["rhos"], cp0=surf["cp0"], cp_coeff=(),
                              emissivity=surf["emissivity"], chi=surf["chi"], day_s=day_s,
                              label="%s surface layer" % (gb.get("name") or body))
    else:
        props = BODY_CAPABILITIES.get(body)
        if props is None:
            raise ValueError("simulation.body %r is not an implemented body and the config "
                             "has no generic_body block" % body)
        mat = SurfaceMaterial(ks=props["ks"], rhos=props["rhos"], cp0=props["cp0"],
                              cp_coeff=tuple(props["cpCoeff"]), emissivity=props["emissivity"],
                              chi=2.7, day_s=props["solar_day_days"] * 86400.0,
                              label="%s surface layer" % body)
    tp = _get(config, "thermophysical", {}) or {}
    if float(tp.get("ks_override", -1) or -1) > 0:
        mat.ks = float(tp["ks_override"])
    if float(tp.get("emissivity_override", -1) or -1) > 0:
        mat.emissivity = float(tp["emissivity_override"])
    ti = float(tp.get("constant_thermal_inertia", -1) or -1)
    if ti > 0 and mat.ks > 0 and mat.rhos > 0:
        mat.cp0 = ti * ti / (mat.ks * mat.rhos)
        mat.cp_coeff = ()
        mat.chi = 0.0
    raw = str(tp.get("conductivity_model", "") or "")
    model = "".join(c for c in raw.lower() if c not in "_- ")   # parseConductivityModel's rule
    if model in ("martinez", "martinez2021", "martinezsiegler2021"):
        mat.model = _MARTINEZ
    elif model not in ("", "hayne", "hayne2017", "default"):
        notes.append("thermophysical.conductivity_model %r is not recognized; D* uses Hayne et "
                     "al. 2017 (the run itself will refuse this value)." % raw)
    if _get(config, "thermal.thermal_inertia_map") or _get(config, "thermal.use_default_property_map"):
        notes.append("a per-face thermal-property map is configured; D* above uses the body-wide "
                     "surface properties, and D* scales roughly with k, so faces with a higher "
                     "thermal inertia have a larger D*.")
    return mat, notes


def _read_mesh(path):
    """(vertices, faces) from PLY or OBJ; None when the format is not readable here."""
    import numpy as np
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".ply":
        from .results import read_ply_mesh_minimal
        m = read_ply_mesh_minimal(path)
        return m["vertices"], m["faces"]
    if suffix == ".obj":
        V, F = [], []
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                if line.startswith("v "):
                    V.append([float(t) for t in line.split()[1:4]])
                elif line.startswith("f "):
                    idx = [int(t.split("/")[0]) for t in line.split()[1:]]
                    idx = [i - 1 if i > 0 else len(V) + i for i in idx]
                    for k in range(1, len(idx) - 1):      # fan-triangulate polygons
                        F.append([idx[0], idx[k], idx[k + 1]])
        return np.asarray(V, dtype=np.float64), np.asarray(F, dtype=np.int64)
    return None


def _icosphere_sizes(resolution, radius_m):
    """Facet statistics of the engine's procedural icosphere, without building it."""
    import numpy as np
    n = 20 * 4 ** int(resolution)
    area = 4.0 * math.pi * radius_m ** 2 / n
    edge = math.sqrt(4.0 * area / math.sqrt(3.0))
    return np.full(1, edge), np.full(1, area * n), n


def advise_config(config_or_path, mesh_units=None, repo_root=None):
    """The advisory for one run config (dict or YAML path). Never runs the engine."""
    import numpy as np
    if not isinstance(config_or_path, dict):
        import yaml
        cfg_path = Path(config_or_path).expanduser().resolve()
        config = yaml.safe_load(cfg_path.read_text()) or {}
        base = cfg_path.parent
    else:
        config, base = config_or_path, Path.cwd()
    from .body_table import BODY_CAPABILITIES
    mat, notes = surface_material_from_config(config)
    lateral_on = bool(_get(config, "thermal.lateral_conduction", False)) and \
        not bool(_get(config, "thermal.lateral_zero_coupling", False))

    mesh_file = _get(config, "mesh.import_file") or _get(config, "generic_body.shape.mesh_file")
    body = str(_get(config, "simulation.body", "Moon"))
    radius_km = (BODY_CAPABILITIES.get(body) or {}).get("radius_km", 0.0)
    if mesh_file:
        p = Path(str(mesh_file)).expanduser()
        candidates = [p] if p.is_absolute() else [base / p, Path(repo_root or Path.cwd()) / p]
        p = next((c for c in candidates if c.exists()), candidates[0])
        if not p.exists():
            raise FileNotFoundError("mesh not found: %s" % p)
        mesh = _read_mesh(p)
        if mesh is None:
            r = evaluate([], [], mat)
            r.lines = format_advice(r, mat, lateral_on) + [
                "[lateral-advisor]   NOTE: %s meshes are not readable here (PLY and OBJ are); "
                "the engine prints this advisory at run time." % p.suffix]
            return r
        V, F = mesh
        V = np.asarray(V, dtype=np.float64)
        F = np.asarray(F, dtype=np.int64)
        centers = V[F].mean(axis=1) if len(F) else np.zeros((0, 3))
        median = float(np.median(np.linalg.norm(centers, axis=1))) if len(centers) else 0.0
        heuristic_m = median > 10000.0
        corroborated = False
        if radius_km > 0 and median > 0:
            says_km = 0.5 <= median / radius_km <= 2.0
            says_m = 0.5 <= (median / 1000.0) / radius_km <= 2.0
            corroborated = (says_km or says_m) and (says_m == heuristic_m)
        if mesh_units in ("m", "km"):
            to_m = 1000.0 if mesh_units == "km" else 1.0
            unit_label = "%s, as requested" % ("kilometers" if mesh_units == "km" else "meters")
            corroborated = True
        else:
            km_confirmed = (not heuristic_m) and corroborated
            to_m = 1000.0 if km_confirmed else 1.0
            unit_label = ("kilometers, confirmed by the body radius" if km_confirmed
                          else "meters, confirmed by the body radius" if corroborated
                          else "meters, not confirmed by a body radius")
        sizes, areas = facet_sizes(V, F, to_m)
        r = evaluate(sizes, areas, mat)
        if not corroborated and r.verdict != "UNMEASURED":
            alt = evaluate(*facet_sizes(V, F, 1000.0), mat)
            if alt.verdict != r.verdict:
                notes.insert(0, "the mesh length unit is not confirmed by a body radius (a locally "
                                "centered patch?). Read as kilometers instead, the verdict would be "
                                "'%s'." % alt.verdict)
        skipped = len(F) - len(sizes)
        if skipped > 0:
            notes.append("%d degenerate or invalid faces were left out of the facet-size "
                         "statistics." % skipped)
    else:
        res = int(_get(config, "mesh.resolution", 6))
        sizes, areas, n = _icosphere_sizes(res, radius_km * 1000.0)
        r = evaluate(sizes, areas, mat)
        r.n_faces = n
        unit_label = "meters, procedural icosphere of resolution %d" % res
    r.lines = format_advice(r, mat, lateral_on, unit_label, notes)
    return r

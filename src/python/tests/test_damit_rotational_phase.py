#!/usr/bin/env python3
"""A body built from a DAMIT model carries the model's rotational phase.

DAMIT gives the orientation of a shape at time t as
``r_ecl = R_z(lambda) R_y(90 - beta) R_z(phi0 + 2 pi (t - t0) / P) r_ast``.
PTM3D builds a generic body's frame from the pole, the J2000 equinox projected
onto the body's equator, and ``W(t) = prime_meridian_deg + 360 (t - spin.epoch)
/ P`` (src/cpp/orbital/body_fixed_frame.cpp). ``damit.rotation`` maps one onto
the other. These checks compare the two frames axis by axis, with the DAMIT
matrix built here from its documented definition, for real spin states with
different poles and epochs, two constructed poles that exercise the engine's
special cases, and every route that writes a body from DAMIT (the client
parser, BodyConfigWriter and the GUI card's fetch function). No network.

Run:  python3 src/python/tests/test_damit_rotational_phase.py
"""
from __future__ import annotations

import logging
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "src" / "python"))
sys.path.insert(0, str(REPO_ROOT / "gui"))

# spin.txt as served by DAMIT (generated_files/open/AsteroidModel/<id>/spin.txt),
# copied 2026-11-12. Only the first two lines matter here.
SPIN_TXT = {
    "433 Eros (asteroid page 1714)": "17 11 5.27025528\n2451545 32.64\n0.43 -0.29 0.22 1 36\n",
    "model 103": "103 27 7.209531\n2434752.77932 0\n0.1\n",
    "model 104": "126 40 16.80061\n2436250 0\n0.1 20.46 0.32 12.18\n",
    "model 1500": "190 -78 3.547454\n2448947 0\n0.1 0.7 0.04 -0.79\n",
    "model 2000": "20 -55 5.59006\n2450989 0\n0.1 0.92 0.21 -0.36\n",
}
# Constructed: a pole within 26 deg of the equinox (the engine then seeds its
# reference axis with ecliptic +y), and an exactly polar one.
CONSTRUCTED = {
    "pole near the equinox": "8 5 9.0\n2455000.5 123.4\n",
    "ecliptic north pole": "0 90 12.0\n2450000.0 300.0\n",
}

# Engine constant (orbital::SpinState::OBLIQUITY_J2000_RAD).
EPS_ENGINE = math.radians(23.4392911)


def rot_z(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rot_y(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def damit_axes(lam, beta, period_h, jd0, phi0, jd):
    """Columns are the shape-frame axes in the ecliptic, from DAMIT's definition."""
    phi = math.radians(phi0) + 2.0 * math.pi * (jd - jd0) * 24.0 / period_h
    return rot_z(math.radians(lam)) @ rot_y(math.radians(90.0 - beta)) @ rot_z(phi)


def engine_axes(spin, jd):
    """Columns are the body-fixed axes in the ecliptic, as BodyFixedFrame builds them
    from the generic_body.spin keys (pole RA/Dec, prime_meridian_deg, epoch, period)."""
    ra, dec = math.radians(spin["pole_ra_deg"]), math.radians(spin["pole_dec_deg"])
    x, y, z = math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)
    p = np.array([x, y * math.cos(EPS_ENGINE) + z * math.sin(EPS_ENGINE),
                  -y * math.sin(EPS_ENGINE) + z * math.cos(EPS_ENGINE)])
    p /= np.linalg.norm(p)
    seed = np.array([1.0, 0.0, 0.0])
    if abs(p @ seed) > 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    e1 = seed - (seed @ p) * p
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(p, e1)
    epoch = spin.get("epoch", "JD 2451545.0")
    assert epoch.startswith("JD "), epoch
    t_sec = (jd - float(epoch[3:])) * 86400.0
    w = (math.radians(spin.get("prime_meridian_deg", 0.0))
         + 2.0 * math.pi / (spin["period_hours"] * 3600.0) * t_sec)
    xb = math.cos(w) * e1 + math.sin(w) * e2
    return np.column_stack([xb, np.cross(p, xb), p])


def worst_axis_angle_deg(a, b):
    """Largest angle between matching columns (atan2 keeps small angles exact)."""
    sin = np.linalg.norm(np.cross(a.T, b.T), axis=1)
    cos = np.sum(a * b, axis=0)
    return float(np.degrees(np.arctan2(sin, cos)).max())


def parse(text, damit_id=0):
    from damit.damit_client import DAMITClient, DAMITShapeModel
    model = DAMITShapeModel(damit_id=damit_id)
    DAMITClient()._parse_spin(text, model)
    return model


# --------------------------------------------------------------------------

def test_spin_txt_second_line_is_read():
    s = parse(SPIN_TXT["433 Eros (asteroid page 1714)"]).spin
    assert (s.lambda_ecl_deg, s.beta_ecl_deg, s.period_hours) == (17.0, 11.0, 5.27025528)
    assert s.has_phase and s.epoch_jd == 2451545.0 and s.phi0_deg == 32.64, s
    s = parse(SPIN_TXT["model 103"]).spin
    assert s.epoch_jd == 2434752.77932 and s.phi0_deg == 0.0, s
    print("    [parse] t0 and phi0 read from line 2 of spin.txt")


def test_frames_agree_with_the_damit_matrix():
    """Every axis within 1e-6 deg, from t0 to decades away (the period is kept exact)."""
    worst = 0.0
    for name, text in {**SPIN_TXT, **CONSTRUCTED}.items():
        s = parse(text).spin
        spin = s.generic_body_spin()
        for days in (0.0, 0.37, -376.3, 3652.5, -18262.0):
            jd = s.epoch_jd + days
            sep = worst_axis_angle_deg(
                damit_axes(s.lambda_ecl_deg, s.beta_ecl_deg, s.period_hours,
                           s.epoch_jd, s.phi0_deg, jd),
                engine_axes(spin, jd))
            assert sep < 1e-6, (name, days, sep, spin)
            worst = max(worst, sep)
    print(f"    [frames] {len(SPIN_TXT) + len(CONSTRUCTED)} spin states x 5 dates: "
          f"largest axis difference {worst:.1e} deg")


def test_the_unconverted_phase_is_detected():
    """The same comparison fails by tens of degrees for what was written before
    (prime meridian 0 at J2000), so the check above can tell the two apart."""
    s = parse(SPIN_TXT["433 Eros (asteroid page 1714)"]).spin
    spin = s.generic_body_spin()
    old = {k: spin[k] for k in ("period_hours", "pole_ra_deg", "pole_dec_deg")}
    jd = 2451171.15447                     # 1998-12-23T15:41:23 UTC, as TDB
    sep = worst_axis_angle_deg(
        damit_axes(17.0, 11.0, 5.27025528, 2451545.0, 32.64, jd), engine_axes(old, jd))
    assert sep > 10.0, sep
    print(f"    [frames] prime meridian 0 at J2000 is {sep:.1f} deg off for Eros")


def test_eros_values():
    spin = parse(SPIN_TXT["433 Eros (asteroid page 1714)"]).spin.generic_body_spin()
    from damit.rotation import phase_offset_deg
    assert abs(phase_offset_deg(17.0, 11.0) - (-86.661)) < 1e-3
    assert abs(spin["prime_meridian_deg"] - (32.64 - 86.661) % 360.0) < 1e-3, spin
    assert spin["epoch"] == "JD 2451545.000000" and spin["period_hours"] == 5.27025528
    print(f"    [eros] W0 = {spin['prime_meridian_deg']:.4f} deg at {spin['epoch']}")


def test_a_model_without_a_phase_says_so_and_keeps_w0_zero():
    records = []

    class _Catch(logging.Handler):
        def emit(self, record):
            records.append(record.getMessage())

    handler = _Catch(level=logging.WARNING)
    for name in ("damit.damit_client", "horizons.body_config_writer"):
        logging.getLogger(name).addHandler(handler)
    try:
        s = parse("17 11 5.27025528\n", damit_id=7).spin
        assert not s.has_phase and "prime_meridian_deg" not in s.generic_body_spin()
        assert any("no epoch and initial rotation angle" in m for m in records), records
        spin = _write(s)["spin"]
        assert spin["prime_meridian_deg"] == 0.0 and "epoch" not in spin, spin
        assert sum("no epoch and initial rotation angle" in m for m in records) == 2, records
    finally:
        for name in ("damit.damit_client", "horizons.body_config_writer"):
            logging.getLogger(name).removeHandler(handler)
    print("    [no phase] warned twice (parser, writer); W0 0 with no epoch, as before")


def _physical(**kw):
    from horizons.horizons_client import PhysicalParams
    phys = PhysicalParams()
    phys.diameter_km = 16.84
    phys.rotation_period_hours = 5.27
    phys.pole_ra_deg = 11.37
    phys.pole_dec_deg = 17.22
    for k, v in kw.items():
        setattr(phys, k, v)
    return phys


def _write(damit_spin=None, spin_overrides=None, physical=None):
    import yaml
    from horizons.horizons_client import OrbitalElements
    from horizons.body_config_writer import BodyConfigWriter
    el = OrbitalElements()
    el.semi_major_axis_au = 1.458
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "body.yaml"
        BodyConfigWriter().write_config(
            str(path), "433 Eros", el, physical or _physical(),
            shape_type="custom_mesh", shape_mesh_file="shape.obj",
            damit_spin=damit_spin, spin_overrides=spin_overrides)
        return yaml.safe_load(path.read_text())["generic_body"]


def test_writer_takes_the_whole_spin_from_damit():
    s = parse(SPIN_TXT["433 Eros (asteroid page 1714)"]).spin
    spin = _write(s)["spin"]
    assert spin == s.generic_body_spin(), spin
    jd = 2451171.15447
    sep = worst_axis_angle_deg(
        damit_axes(17.0, 11.0, 5.27025528, 2451545.0, 32.64, jd), engine_axes(spin, jd))
    assert sep < 1e-6, sep
    # An explicit value wins over DAMIT.
    spin = _write(s, spin_overrides={"prime_meridian_deg": 10.0})["spin"]
    assert spin["prime_meridian_deg"] == 10.0 and spin["epoch"] == "JD 2451545.000000"
    print("    [writer] period, pole, W0 and epoch from DAMIT; spin_overrides win")


def test_writer_unchanged_without_damit():
    spin = _write()["spin"]
    assert spin == {"period_hours": 5.27, "pole_ra_deg": 11.37,
                    "pole_dec_deg": 17.22, "prime_meridian_deg": 0.0}, spin
    spin = _write(physical=_physical(rotation_period_hours=0.0, pole_ra_deg=-1.0,
                                     pole_dec_deg=-1.0))["spin"]
    assert spin == {"period_hours": 6.0, "prime_meridian_deg": 0.0}, spin
    print("    [writer] a body with no DAMIT spin is written as before")


class _Model:
    damit_id = 1714
    vertices = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, -1.0, -1.0]]
    faces = [[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]]
    num_vertices = 4
    num_faces = 4

    def __init__(self, text):
        self.spin = parse(text).spin


class _Client:
    text = SPIN_TXT["433 Eros (asteroid page 1714)"]

    def search_asteroid(self, query):
        return [{"damit_id": 1714, "number": 433, "name": "Eros"}]

    def get_model(self, damit_id):
        return _Model(self.text)


def test_gui_card_fetch_carries_the_phase():
    """The function behind the Asteroid card's "DAMIT shape model..." button,
    and the configuration the card then holds and saves."""
    import yaml
    from gui_core.body_fetch import fetch_damit_shape
    from gui_core.config_manager import ConfigManager
    with tempfile.TemporaryDirectory() as root:
        result = fetch_damit_shape(root, name="433 Eros", naif_id="20000433",
                                   radius_km=8.41, client_factory=_Client)
        assert result.ok, result.message
        expected = parse(_Client.text).spin.generic_body_spin()
        for key in ("pole_ra_deg", "pole_dec_deg", "prime_meridian_deg", "epoch"):
            assert result.values[f"generic_body.spin.{key}"] == expected[key], key
        assert any("rotational phase came with the model" in w for w in result.warnings)
        # From the cache: the sidecar keeps t0 and phi0.
        again = fetch_damit_shape(root, name="433 Eros", naif_id="20000433",
                                  client_factory=lambda: None)
        assert again.from_cache
        assert again.values["generic_body.spin.prime_meridian_deg"] == expected["prime_meridian_deg"]

        cm = ConfigManager()
        cm.set_value("generic_body.name", "433 Eros")
        for key, value in result.values.items():
            cm.set_value(key, value)
        out = Path(root) / "saved.yaml"
        cm.save_to_yaml(str(out))
        spin = yaml.safe_load(out.read_text())["generic_body"]["spin"]
        for key in ("prime_meridian_deg", "epoch"):
            assert spin[key] == expected[key], (key, spin)
    print("    [gui] card values and saved configuration carry W0 and epoch; cache keeps them")


def test_fetch_and_cache_body_takes_the_spin_from_damit():
    """data_manager.fetch_and_cache_body, the route of the GUI's HORIZONS page, with
    stand-ins for the two services. HORIZONS lists a rounded period and its own pole;
    neither may reach a body that uses the DAMIT shape."""
    import yaml
    import damit
    import horizons
    import data_manager
    from horizons.horizons_client import OrbitalElements

    class _Body:
        name = "433 Eros (A898 PA)"
        naif_id = "433"
        body_type = "asteroid"

        def __init__(self):
            self.elements = OrbitalElements()
            self.elements.semi_major_axis_au = 1.458
            self.physical = _physical()

    class _Horizons:
        def get_body(self, query):
            return _Body()

    class _Damit(_Client):
        def __init__(self, **kw):
            pass

    saved = horizons.HorizonsClient, damit.DAMITClient
    horizons.HorizonsClient, damit.DAMITClient = _Horizons, _Damit
    try:
        with tempfile.TemporaryDirectory() as d:
            path = data_manager.fetch_and_cache_body("433", cache_dir=d,
                                                     ephemeris_source="keplerian")
            gb = yaml.safe_load(Path(path).read_text())["generic_body"]
    finally:
        horizons.HorizonsClient, damit.DAMITClient = saved
    assert gb["shape"]["type"] == "custom_mesh", gb["shape"]
    assert gb["spin"] == parse(_Client.text).spin.generic_body_spin(), gb["spin"]
    print("    [data_manager] period, pole and phase from DAMIT, not from HORIZONS")


def test_gui_cached_body_config_keeps_the_epoch():
    """Loading a cached body_config.yaml must bring spin.epoch with the W0 it belongs to."""
    from gui_core import body_lookup
    assert "spin.epoch" in body_lookup._COPIED
    print("    [gui] a cached body_config.yaml copies spin.epoch with prime_meridian_deg")


CHECKS = [v for k, v in sorted(globals().items()) if k.startswith("test_")]


def main():
    failed = 0
    for check in CHECKS:
        name = check.__name__
        try:
            check()
        except Exception as exc:                             # noqa: BLE001
            failed += 1
            import traceback
            print(f"FAIL {name}: {exc!r}")
            traceback.print_exc()
        else:
            print(f"PASS {name}")
    print(f"\n{len(CHECKS) - failed}/{len(CHECKS)} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

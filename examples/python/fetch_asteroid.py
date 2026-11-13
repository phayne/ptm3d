#!/usr/bin/env python3
"""Fetch a numbered asteroid from JPL Horizons and DAMIT and prepare it for PTM3D.

    python3 examples/python/fetch_asteroid.py 433 --damit-id 1714 --start 1993-01-01 \
        --stop 1999-01-01 --epoch 1998-12-23T15:40:00 --outdir data/asteroids/eros

Writes four files to --outdir. ephemeris.bsp is a Horizons SPK covering the dates given.
metakernel.tm names that SPK and the standard kernels. shape.obj is the DAMIT shape, in the
units DAMIT gives, or rescaled with --radius-km. generic_body.yaml is the generic_body block to paste into a configuration,
with the orbit, the spin state carried to --epoch, and the shape. Use --shape sphere for a body
that has no DAMIT model, or to run a sphere. The sphere has the diameter, pole and period that
Horizons lists, and no shape.obj is written.

The surface properties are not written. They are yours to choose.
"""
import argparse
import math
import os
import sys
import urllib.request

import numpy as np
import spiceypy as spice

from damit import DAMITClient
from horizons import HorizonsClient
from horizons.horizons_client import ecliptic_pole_to_equatorial


def prime_meridian(lam, beta, period_h, jd0, phi0, jd_ref):
    """Prime meridian W (degrees) in the PTM3D body-fixed frame at the Julian date jd_ref.

    lam and beta are the DAMIT ecliptic pole, period_h is the sidereal period, and jd0 and phi0
    are the DAMIT reference epoch and phase. The PTM3D frame has its x axis along the J2000
    vernal equinox, projected onto the equator of the body. The DAMIT phase is measured from
    u = (sin b cos l, sin b sin l, -cos b) instead, so W = phi + delta, where delta is the angle
    about the pole from the first direction to the second. Returns W and delta.
    """
    l, b = math.radians(lam), math.radians(beta)
    pole = np.array([math.cos(b) * math.cos(l), math.cos(b) * math.sin(l), math.sin(b)])
    u = np.array([math.sin(b) * math.cos(l), math.sin(b) * math.sin(l), -math.cos(b)])
    seed = np.array([1.0, 0.0, 0.0])
    if abs(pole @ seed) > 0.9:          # the same choice the model makes for a pole near +x
        seed = np.array([0.0, 1.0, 0.0])
    e1 = seed - (seed @ pole) * pole
    e1 /= np.linalg.norm(e1)
    delta = math.degrees(math.atan2(np.cross(e1, u) @ pole, e1 @ u))
    phi = phi0 + 360.0 * (jd_ref - jd0) * 24.0 / period_h
    return (phi + delta) % 360.0, delta


def mesh_size(vertices, faces):
    """Radius of the sphere with the volume of the mesh, and the signed volume of the mesh."""
    v = np.asarray(vertices, dtype=float)
    f = np.asarray(faces, dtype=int)
    vol = np.einsum('ij,ij->i', v[f[:, 0]], np.cross(v[f[:, 1]], v[f[:, 2]])).sum() / 6.0
    return (3.0 * abs(vol) / (4.0 * math.pi)) ** (1.0 / 3.0), vol


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('number', type=int, help='asteroid number, for example 433')
    ap.add_argument('--start', required=True, help='first date of the SPK (ISO, UTC)')
    ap.add_argument('--stop', required=True, help='last date of the SPK (ISO, UTC)')
    ap.add_argument('--epoch', required=True, help='time at which the spin and orbit are set (ISO, UTC)')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--shape', choices=('damit', 'sphere'), default='damit')
    ap.add_argument('--damit-id', type=int, default=None,
                    help='DAMIT model number, the last number in the address of the model page. '
                         'It is searched for by asteroid number if omitted, which can fail when the site is slow')
    ap.add_argument('--radius-km', type=float, default=None,
                    help='rescale the DAMIT shape so that a sphere of equal volume has this radius in km. '
                         'Needed only when the model is not in km')
    ap.add_argument('--lsk', default='data/kernels/lsk/naif0012.tls', help='leap-seconds kernel')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    spice.furnsh(args.lsk)
    jd_ref = 2451545.0 + spice.str2et(args.epoch) / 86400.0      # the epoch as a TDB Julian date

    hz = HorizonsClient()
    body = hz.get_body(str(args.number))
    phys = body.physical
    print(f'Horizons: {body.name}, diameter {phys.diameter_km} km, period {phys.rotation_period_hours} h, '
          f'geometric albedo {phys.albedo}, pole RA {phys.pole_ra_deg} Dec {phys.pole_dec_deg}')
    spk = os.path.join(args.outdir, 'ephemeris.bsp')
    hz.download_spk(str(args.number), args.start, args.stop, spk)
    print(f'SPK: {os.path.getsize(spk) / 1e3:.0f} kB')
    el = hz.get_orbital_elements(f'{args.number};', f'{jd_ref:.6f}')
    ids = spice.cell_int(10)
    spice.spkobj(spk, ids)

    if args.shape == 'damit':
        dm = DAMITClient()
        model = dm.get_model(args.damit_id) if args.damit_id else dm.get_model_by_number(args.number)
        if model is None:
            sys.exit(f'DAMIT has no model for asteroid {args.number}, or the search failed. Give --damit-id.')
        links = dm._get_model_links(model.damit_id)
        rows = [r.split() for r in
                urllib.request.urlopen(links['spin_txt_url'], timeout=30).read().decode().splitlines() if r.strip()]
        lam, beta, period_h = (float(x) for x in rows[0][:3])
        jd0, phi0 = (float(x) for x in rows[1][:2])
        radius, vol = mesh_size(model.vertices, model.faces)
        print(f'DAMIT model {model.damit_id}: {model.num_vertices} vertices, {model.num_faces} faces, '
              f'equal-volume radius {radius:.3f} in the units of the file, signed volume {vol:.1f}')
        print(f'spin: lambda {lam}, beta {beta}, P {period_h} h, JD0 {jd0}, phi0 {phi0}')
        scale = args.radius_km / radius if args.radius_km else 1.0
        radius *= scale
        obj = os.path.join(args.outdir, 'shape.obj')
        with open(obj, 'w') as fh:
            for x, y, z in model.vertices:
                fh.write(f'v {x * scale:.6f} {y * scale:.6f} {z * scale:.6f}\n')
            for a, b, c in model.faces:
                fh.write(f'f {a + 1} {b + 1} {c + 1}\n')
        w0, delta = prime_meridian(lam, beta, period_h, jd0, phi0, jd_ref)
        ra, dec = ecliptic_pole_to_equatorial(lam, beta)
        print(f'pole RA {ra:.4f}, Dec {dec:.4f}; delta {delta:.3f}; W at JD {jd_ref:.6f} = {w0:.3f}')
        shape = (f'    type: custom_mesh\n    mesh_file: {obj}\n    radius_km: {radius:.3f}\n')
    else:
        radius, period_h, ra, dec, w0 = (phys.diameter_km / 2.0, phys.rotation_period_hours,
                                         phys.pole_ra_deg, phys.pole_dec_deg, 0.0)
        shape = f'    type: spherical\n    radius_km: {radius:.3f}\n'

    mk = os.path.join(args.outdir, 'metakernel.tm')
    with open(mk, 'w') as fh:
        fh.write("KPL/MK\n\n\\begindata\n"
                 f"  PATH_VALUES  = ( 'data/kernels', '{args.outdir}' )\n"
                 "  PATH_SYMBOLS = ( 'KERNELS', 'BODY' )\n"
                 "  KERNELS_TO_LOAD = ( '$KERNELS/lsk/naif0012.tls', '$KERNELS/pck/pck00011.tpc',\n"
                 "                      '$KERNELS/spk/de440s.bsp', '$BODY/ephemeris.bsp' )\n"
                 "\\begintext\n")

    block = (f"generic_body:\n"
             f"  name: {body.name}\n"
             f"  ephemeris_source: spice            # or keplerian, to use the elements below\n"
             f"  spice_metakernel: {mk}\n"
             f"  naif_id: '{ids[0]}'\n"
             f"  shape:\n{shape}"
             f"  orbit:\n"
             f"    semi_major_axis_au: {el.semi_major_axis_au:.8f}\n"
             f"    eccentricity: {el.eccentricity:.8f}\n"
             f"    inclination_deg: {el.inclination_deg:.6f}\n"
             f"    lon_ascending_node_deg: {el.lon_ascending_node_deg:.6f}\n"
             f"    arg_perihelion_deg: {el.arg_perihelion_deg:.6f}\n"
             f"    mean_anomaly_deg: {el.mean_anomaly_deg:.6f}\n"
             f"    epoch: 'JD {jd_ref:.6f}'\n"
             f"  spin:\n"
             f"    period_hours: {period_h}\n"
             f"    pole_ra_deg: {ra:.4f}\n"
             f"    pole_dec_deg: {dec:.4f}\n"
             f"    prime_meridian_deg: {w0:.3f}\n"
             f"    epoch: 'JD {jd_ref:.6f}'\n")
    out = os.path.join(args.outdir, 'generic_body.yaml')
    with open(out, 'w') as fh:
        fh.write(block)
    print(f'wrote {out}')


if __name__ == '__main__':
    sys.exit(main())

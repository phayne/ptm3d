#!/usr/bin/env python3
"""Fetch a numbered asteroid from JPL Horizons and DAMIT and prepare it for PTM3D.

    python3 examples/python/fetch_asteroid.py 433 --damit-id 1714 --start 1993-01-01 \
        --stop 1999-01-01 --epoch 1998-12-23T15:40:00 --outdir data/asteroids/eros

Writes four files to --outdir. ephemeris.bsp is a Horizons SPK covering the dates given.
metakernel.tm names that SPK and the standard kernels. shape.obj is the DAMIT shape, in the
units DAMIT gives, or rescaled with --radius-km. generic_body.yaml is the generic_body block to paste into a configuration,
with the orbit at --epoch, the spin state of the DAMIT model, and the shape. Use --shape sphere for a body
that has no DAMIT model, or to run a sphere. The sphere has the diameter, pole and period that
Horizons lists, and no shape.obj is written.

The surface properties are not written. They are yours to choose.
"""
import argparse
import math
import os
import sys

import numpy as np
import spiceypy as spice

from damit import DAMITClient
from horizons import HorizonsClient


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
    ap.add_argument('--epoch', required=True, help='time at which the orbital elements are given (ISO, UTC)')
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
        try:
            model = dm.get_model(args.damit_id) if args.damit_id else dm.get_model_by_number(args.number)
        except RuntimeError as exc:
            sys.exit(f'The DAMIT search did not complete ({exc}). Try again, or give --damit-id.')
        if model is None:
            sys.exit(f'DAMIT has no model for asteroid {args.number}.')
        sp = model.spin
        radius, vol = mesh_size(model.vertices, model.faces)
        print(f'DAMIT model {model.damit_id}: {model.num_vertices} vertices, {model.num_faces} faces, '
              f'equal-volume radius {radius:.3f} in the units of the file, signed volume {vol:.1f}')
        print(f'spin: lambda {sp.lambda_ecl_deg:g}, beta {sp.beta_ecl_deg:g}, P {sp.period_hours} h, '
              + (f'JD0 {sp.epoch_jd}, phi0 {sp.phi0_deg:g}' if sp.has_phase else 'no phase'))
        scale = args.radius_km / radius if args.radius_km else 1.0
        radius *= scale
        obj = os.path.join(args.outdir, 'shape.obj')
        with open(obj, 'w') as fh:
            for x, y, z in model.vertices:
                fh.write(f'v {x * scale:.6f} {y * scale:.6f} {z * scale:.6f}\n')
            for a, b, c in model.faces:
                fh.write(f'f {a + 1} {b + 1} {c + 1}\n')
        spin = sp.generic_body_spin()
        lines = (f"    period_hours: {spin['period_hours']}\n"
                 f"    pole_ra_deg: {spin['pole_ra_deg']:.6f}\n"
                 f"    pole_dec_deg: {spin['pole_dec_deg']:.6f}\n")
        if sp.has_phase:
            lines += (f"    prime_meridian_deg: {spin['prime_meridian_deg']:.6f}\n"
                      f"    epoch: '{spin['epoch']}'\n")
            print(f"pole RA {spin['pole_ra_deg']:.4f}, Dec {spin['pole_dec_deg']:.4f}; "
                  f"W0 {spin['prime_meridian_deg']:.4f} at {spin['epoch']}")
        else:
            lines += "    prime_meridian_deg: 0.000\n"
            print(f"WARNING: DAMIT model {model.damit_id} gives no epoch and initial rotation angle, so "
                  f"prime_meridian_deg is 0 at J2000 and the longitudes of the shape are arbitrary.")
        shape = (f'    type: custom_mesh\n    mesh_file: {obj}\n    radius_km: {radius:.3f}\n')
    else:
        radius = phys.diameter_km / 2.0
        lines = (f"    period_hours: {phys.rotation_period_hours}\n"
                 f"    pole_ra_deg: {phys.pole_ra_deg:.4f}\n"
                 f"    pole_dec_deg: {phys.pole_dec_deg:.4f}\n"
                 f"    prime_meridian_deg: 0.000\n"
                 f"    epoch: 'JD {jd_ref:.6f}'\n")
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
             f"  spin:\n{lines}")
    out = os.path.join(args.outdir, 'generic_body.yaml')
    with open(out, 'w') as fh:
        fh.write(block)
    print(f'wrote {out}')


if __name__ == '__main__':
    sys.exit(main())

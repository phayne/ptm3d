#!/usr/bin/env python3
"""Compare the sub-solar point of a finished run with the DAMIT rotation model.

    python3 examples/python/check_damit_phase.py output/asteroid_eros \
        --metakernel data/asteroids/eros/metakernel.tm --naif 20000433 \
        --spin 17 11 5.27025528 2451545.0 32.64

The five numbers after --spin are the DAMIT ecliptic pole longitude and latitude (degrees), the
sidereal period (hours), the reference epoch (Julian date) and the phase at that epoch (degrees),
as listed on the DAMIT page of the model. The Sun direction comes from SPICE, the body-fixed
direction from the DAMIT rotation matrix, and both are compared with the Sun vector that the run
wrote to body_positions.csv. Run it from the repository root, where the metakernel paths resolve.
"""
import argparse
import math

import numpy as np
import spiceypy as spice

from ptm import Results


def rot_z(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rot_y(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir')
    ap.add_argument('--metakernel', required=True)
    ap.add_argument('--naif', required=True, help='NAIF id of the body in the SPK, for example 20000433')
    ap.add_argument('--spin', nargs=5, type=float, required=True,
                    metavar=('LAMBDA', 'BETA', 'PERIOD_H', 'JD0', 'PHI0'))
    ap.add_argument('--steps', type=int, default=8, help='number of output steps to compare')
    args = ap.parse_args()
    lam, beta, period_h, jd0, phi0 = args.spin

    spice.furnsh(args.metakernel)
    run = Results(args.run_dir, mesh_units='km')
    et = run.times().et
    positions = run.body_positions()
    worst = 0.0
    for i in np.linspace(0, len(et) - 1, args.steps).astype(int):
        sun_ecliptic, _ = spice.spkpos('SUN', et[i], 'ECLIPJ2000', 'LT+S', args.naif)
        phi = math.radians(phi0 + 360.0 * (2451545.0 + et[i] / 86400.0 - jd0) * 24.0 / period_h)
        to_ecliptic = rot_z(math.radians(lam)) @ rot_y(math.radians(90.0 - beta)) @ rot_z(phi)
        predicted = to_ecliptic.T @ sun_ecliptic
        predicted /= np.linalg.norm(predicted)
        row = positions[i]
        modeled = np.array([row['sun_x'], row['sun_y'], row['sun_z']])
        modeled /= np.linalg.norm(modeled)
        sep = math.degrees(math.acos(np.clip(predicted @ modeled, -1.0, 1.0)))
        worst = max(worst, sep)
        lat, lon = (math.degrees(math.asin(predicted[2])), math.degrees(math.atan2(predicted[1], predicted[0])))
        run_lat, run_lon = run.subsolar_latlon(i)
        print(f'step {i:4d}  {run.times().utc[i]}  sub-solar point: DAMIT {lat:7.2f}, {lon:8.2f}  '
              f'run {run_lat:7.2f}, {run_lon:8.2f}  separation {sep:.4f} deg')
    print(f'largest separation {worst:.4f} deg')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Spin-up arm runner (Ada 2026-11-08). One engine run from the n4990 W6h radoff base (or W6h_prod for rad-on).

usage: arm.py NAME --win 2030-09-01 --eq 8 [--base radoff|radon] [--exe build] [--set thermal.seed_seasonal=auto ...]
              [--env K=V ...] [--fine N --state FILE]   (P/I protocol: warm start at W - N*P, integration offset N*P)
Takes /tmp/ptm_cpu_heavy.lock for the run; appends a provenance line to STATUS.txt.
"""
import sys, os, json, time, argparse, subprocess, hashlib, datetime as dt
import yaml

K = '/Users/paha3326/projects/spinup_2026-11-08'
W = '/Users/paha3326/projects/ptm3d_spinup'
GIT = '/Library/Developer/CommandLineTools/usr/bin/git'
BASES = {'radoff': '/Users/paha3326/projects/ptm405_2026-11-07/w6h/D_cpu_radoff.yaml',
         'radon': '/Users/paha3326/projects/ptm846_followup_2026-11-08/gates_m2/tip/W6h_prod.yaml'}
VF = '/Users/paha3326/projects/moon3d_R1/ptm733_tests/csr_up/vf.csr'
SYN = 29.530588 * 86400.0          # engine planet.day for the Moon (planetary_properties.h)
WINDAYS = 29.535040430555554
CLEAN = ['PTM_SEED_WINDOW', 'PTM_SEED_SEASONAL', 'PTM_GPU_RETHERMAL', 'PTM_CPU_RETHERMAL', 'PTM_TIMING_JSON',
         'PTM_GPU_DISPATCH_TARGET_S', 'PTM_MAX_TIMESTEPS_PER_DISPATCH', 'PTM_RAD_POOL_BATCH_BUDGET_MB',
         'PTM_CN_SUBSTEP_DIV', 'PTM_ENERGY_LEDGER', 'PTM_THERMAL_TIME_OUTER', 'PTM_EQ_PERFACE_FLUX',
         'PTM_RAM_ILLUM_MAX_GB', 'PTM_PRIMARY_DISK_STATS', 'PTM_FLUX_REPLAY', 'PTM_FLUX_REPLAY_PERIODIC']
PROD = dict(OMP_NUM_THREADS='8', PTM_ILLUM_TRUE_NORMALS='1', PTM_ILLUM_RISE_EPS='1', PTM_ILLUM_RISE_KEPS='0.2',
            PTM_ILLUM_RISE_HABS='0', PTM_ILLUM_SEAL_HOLES='1', PTM_HMATRIX_VF_REPLACE=VF)


def setkey(c, path, val):
    ks = path.split('.'); d = c
    for k in ks[:-1]: d = d.setdefault(k, {})
    d[ks[-1]] = val if val in ('on', 'off', 'yes', 'no') else yaml.safe_load(val)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name'); ap.add_argument('--win', required=True); ap.add_argument('--eq', type=float, default=8)
    ap.add_argument('--base', default='radoff'); ap.add_argument('--exe', default='build')
    ap.add_argument('--set', action='append', default=[]); ap.add_argument('--env', action='append', default=[])
    ap.add_argument('--fine', type=float, default=None); ap.add_argument('--state', default=None)
    ap.add_argument('--nwin', type=float, default=1.0, help='output span in lunation windows (truth state dumps)')
    ap.add_argument('--final-state', default=None)
    ap.add_argument('--state-at', type=float, default=None, help='truth dump: window ends at W - N*P')
    a = ap.parse_args()
    out = f'{K}/runs/{a.name}'
    c = yaml.safe_load(open(BASES[a.base]))
    c['output']['directory'] = out
    c['thermal']['compute_ice_stability'] = True
    w0 = dt.datetime.fromisoformat(a.win)
    DTW = WINDAYS * 86400.0 / 120.0           # window cadence [s]
    if a.fine is not None:      # warm start at W - nfs*DTW (nfs = N_fine lunations on the window cadence)
        nfs = int(round(a.fine * SYN / DTW))
        t0 = w0 - dt.timedelta(seconds=nfs * DTW)
        c['simulation']['start_time'] = t0.isoformat(timespec='microseconds')
        c['simulation']['duration_days'] = (nfs + 120) * DTW / 86400.0
        c['simulation']['timesteps'] = nfs + 120
        c['thermal']['integration_offset_steps'] = nfs
        c['thermal']['integration_steps'] = 120
        c['thermal']['equilibration_orbits'] = 0.0
        c['thermal']['initial_state_file'] = a.state
    elif a.state_at is not None:  # truth state dump: one window ENDING at W - nfs*DTW, final state written there
        nfs = int(round(a.state_at * SYN / DTW))
        t0 = w0 - dt.timedelta(seconds=(nfs + 119) * DTW)   # last sample (= the final state's epoch) at W - nfs*DTW
        c['simulation']['start_time'] = t0.isoformat(timespec='microseconds')
        c['simulation']['duration_days'] = WINDAYS
        c['simulation']['timesteps'] = 120
        c['thermal']['equilibration_orbits'] = float(a.eq)
    else:
        c['simulation']['start_time'] = w0.isoformat(timespec='seconds')
        c['simulation']['duration_days'] = WINDAYS * a.nwin
        c['simulation']['timesteps'] = int(round(120 * a.nwin))
        c['thermal']['equilibration_orbits'] = float(a.eq)
    if a.final_state: c['thermal']['final_state_file'] = a.final_state
    for s in a.set:
        k, v = s.split('=', 1); setkey(c, k, v)
    if a.base == 'radoff':
        assert not c['radiation']['enabled']
    assert not any((c.get(k) or {}).get('use_gpu') for k in ('thermal', 'radiation', 'illumination'))
    yaml.safe_dump(c, open(f'{out}.yaml', 'w'), sort_keys=False)
    env = {k: v for k, v in os.environ.items() if k not in CLEAN}
    env.update(PROD)
    env['PTM_ENERGY_LEDGER'] = f'{out}.ledger.csv'
    for e in a.env:
        k, v = e.split('=', 1); env[k] = v
    exe = f'{W}/{a.exe}/bin/ptm_thermal_sim'
    lk = '/tmp/ptm_cpu_heavy.lock'
    while True:
        try: os.mkdir(lk); break
        except FileExistsError: time.sleep(10)
    try:
        open(f'{lk}/info', 'w').write(f'pid={os.getpid()} owner=Ada spinup {a.name} (n4990 CPU short run) since '
                                      f'{dt.datetime.utcnow().isoformat()}Z\n')
        subprocess.run(['rm', '-rf', out])
        t = time.time()
        with open(f'{out}.log', 'w') as lf:
            rc = subprocess.call(['nice', '-n', '10', exe, '--config', f'{out}.yaml'], cwd=W, env=env,
                                 stdout=lf, stderr=subprocess.STDOUT)
        wall = time.time() - t
    finally:
        subprocess.run(['rm', '-rf', lk])
    md5 = hashlib.md5(open(exe, 'rb').read()).hexdigest()
    sha = subprocess.run([GIT, '-C', W, 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run([GIT, '-C', W, 'status', '--porcelain', '--untracked-files=no'], capture_output=True,
                           text=True).stdout.strip()
    nstat = len(os.listdir(f'{out}/statistics')) if os.path.isdir(f'{out}/statistics') else 0
    line = (f'{dt.datetime.utcnow().isoformat(timespec="seconds")}Z {a.name} win={a.win} eq={a.eq} fine={a.fine} '
            f'base={a.base} set={a.set} env={a.env} rc={rc} wall={wall:.0f}s nstat={nstat} exe={a.exe}:{md5} '
            f'sha={sha}{"+dirty" if dirty else ""}')
    open(f'{K}/STATUS.txt', 'a').write(line + '\n'); print(line)
    sys.exit(rc)


if __name__ == '__main__':
    main()

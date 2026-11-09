#!/bin/bash
# Default-identity gates for feat/spinup-seed-fixes (Ada 2026-11-08): base = unmodified 961ceb8 build
# (exe/ptm_thermal_sim.961ceb8, md5 043a1c88) vs tip = this worktree's build, DEFAULT env (no diagnostics set),
# plus tip with PTM_DUMP_STATE_ET set (the dump must not perturb products). Ledger off. CPU only.
set -u
K=/Users/paha3326/projects/spinup_2026-11-08; G=$K/gates; M=$G/manifests; mkdir -p $G $M
C=/Users/paha3326/projects/moon3d_R1/runs/ptm727_control/gates; B=/Users/paha3326/projects/moon3d_R1/bench
PY=/Users/paha3326/anaconda3/envs/europa_viz/bin/python3; MF=/Users/paha3326/projects/ptm846_implicit_2026-11-08/tools/manifest.py
VF968=$B/smoke_harness/B_bl_n968_sweep_nr2e4/vf.csr; VF4990=/Users/paha3326/projects/moon3d_R1/ptm733_tests/csr_up/vf.csr
W6H=/Users/paha3326/projects/seed_tier2_2026-11-03/gpucheck/W6h_auto_cpu.yaml
W=/Users/paha3326/projects/ptm3d_spinup
exe(){ case $1 in base) echo $K/exe/ptm_thermal_sim.961ceb8;; *) echo $W/build/bin/ptm_thermal_sim;; esac; }
CLEAN="env -u PTM_SEED_WINDOW -u PTM_SEED_SEASONAL -u PTM_TIMING_JSON -u PTM_CPU_RETHERMAL -u PTM_GPU_RETHERMAL -u PTM_PRIMARY_DISK_STATS -u PTM_HMATRIX_VF_REPLACE -u PTM_ILLUM_TRUE_NORMALS -u PTM_ILLUM_RISE_EPS -u PTM_ILLUM_RISE_KEPS -u PTM_ILLUM_RISE_HABS -u PTM_ILLUM_SEAL_HOLES -u PTM_ILLUM_RISE_EPS2 -u PTM_ILLUM_TOPO_SELFSKIP -u PTM_RAD_POOL_BATCH_BUDGET_MB -u PTM_THERMAL_TIME_OUTER -u PTM_EQ_PERFACE_FLUX -u PTM_RAM_ILLUM_MAX_GB -u PTM_CN_SUBSTEP_DIV -u PTM_ENERGY_LEDGER -u PTM_SEED_STATE_FILE -u PTM_DUMP_STATE_ET -u PTM_DUMP_ABSORBED -u PTM_DUMP_SEASONAL_TABLE -u PTM_SEED_SEASONAL_DC_ONLY"
PROD="PTM_ILLUM_TRUE_NORMALS=1 PTM_ILLUM_RISE_EPS=1 PTM_ILLUM_RISE_KEPS=0.2 PTM_ILLUM_RISE_HABS=0 PTM_ILLUM_SEAL_HOLES=1"
ST=$G/STATUS.txt; log(){ echo "$(date -u +%FT%TZ) $*" >> $ST; }
log "ARM base exe=$(md5 -q $(exe base)) | tip $(cd $W && /Library/Developer/CommandLineTools/usr/bin/git rev-parse --short HEAD)+$(cd $W && /Library/Developer/CommandLineTools/usr/bin/git status --porcelain -uno | wc -l | tr -d ' ')dirty exe=$(md5 -q $(exe tip))"
mk(){ $PY - "$@" <<'PYEOF'
import sys, yaml
src, out, y, m, ice = sys.argv[1:6]
c = yaml.safe_load(open(src)); c['output']['directory'] = out
if m == 'cpu':
    c.setdefault('illumination', {})['use_gpu'] = False; c['thermal']['use_gpu'] = False
    c.setdefault('radiation', {})['use_gpu'] = False; c['output']['skip_time_series'] = False
if ice == 'ice': c['thermal']['compute_ice_stability'] = True
g = [k for k in ('thermal', 'radiation', 'illumination') if (c.get(k) or {}).get('use_gpu')]
if g: sys.exit(f'REFUSE: use_gpu true in {g}')
open(y, 'w').write(yaml.safe_dump(c, sort_keys=False))
PYEOF
}
runarm(){ local n=$1 a=$2 src=$3 mode=$4 ice=$5 env=$6 extra=$7; local O=$G/$a/$n; mkdir -p $G/$a; rm -rf $O
  mk $src $O $G/$a/$n.yaml $mode $ice || { log "  $n $a MK-FAIL"; return; }
  local s=$(date +%s)
  ( cd $W && $CLEAN OMP_NUM_THREADS=8 $env $extra nice -n 10 $(exe $a) --config $G/$a/$n.yaml > $G/$a/$n.log 2>&1 ); local rc=$?
  $PY $MF make $O $M/${n}__$a.json > /dev/null
  log "  $n $a rc=$rc wall=$(( $(date +%s)-s ))s files=$(find $O -type f | wc -l | tr -d ' ')"
}
gate(){ local n=$1 src=$2 mode=$3 ice=$4 env=$5
  log "GATE $n"
  until mkdir /tmp/ptm_cpu_heavy.lock 2>/dev/null; do sleep 15; done
  echo "pid=$$ owner=Ada spinup gate $n since $(date -u +%FT%TZ)" > /tmp/ptm_cpu_heavy.lock/info
  runarm $n base $src $mode $ice "$env" ""
  runarm $n tip $src $mode $ice "$env" ""
  runarm $n dump $src $mode $ice "$env" "PTM_DUMP_STATE_ET=1.0:$G/dump_$n.bin"
  rm -rf /tmp/ptm_cpu_heavy.lock
  for a in tip dump; do log "  $n base-vs-$a $($PY $MF cmp $M/${n}__base.json $M/${n}__$a.json)"; done
}
exe dump >/dev/null
gate G3_def  $C/G3_moon_cpu_simd_new.yaml none none "PTM_HMATRIX_VF_REPLACE=$VF968"
gate G3_prod $C/G3_moon_cpu_simd_new.yaml none none "PTM_HMATRIX_VF_REPLACE=$VF968 $PROD"
gate G4_def  $C/G4_generic_justitia_new.yaml none none ""
gate G5_def  $C/G5_europa_new.yaml none none ""
gate G3_ice  $C/G3_moon_cpu_simd_new.yaml none ice "PTM_HMATRIX_VF_REPLACE=$VF968"
gate G5_ice  $C/G5_europa_new.yaml none ice ""
gate W6h_prod $W6H cpu none "$PROD PTM_HMATRIX_VF_REPLACE=$VF4990"
log "DONE"

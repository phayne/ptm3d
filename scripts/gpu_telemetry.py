#!/usr/bin/env python3
"""Per-die GPU telemetry on macOS without sudo: temperature, power, activity, VRAM.

`powermetrics` needs root. The IOAccelerator registry does not: each GPU die
publishes "PerformanceStatistics" with "Temperature(C)", "Total Power(W)",
"GPU Activity(%)" and "inUseVidMemoryBytes". This samples it at a fixed
interval and appends CSV rows, so a GPU run (or a pause ramp-down trial) can be
judged on measured die temperature and power instead of guesses.

  log mode:        gpu_telemetry.py log --out tel.csv [--interval 2] [--duration S] [--until-file F]
  summarize mode:  gpu_telemetry.py summarize tel.csv [--from ISO] [--to ISO]

Dies are keyed by their IORegistryEntryID (each die appears twice in the
registry; the first entry per id is kept). The registry does not say which
OpenCL device index a die is; in a single-device run the busy die is the one
with activity. Read-only: it never changes a setting.
"""
import argparse
import csv
import os
import plistlib
import subprocess
import sys
import time
from datetime import datetime, timezone

FIELDS = ("utc", "die_id", "temperature_c", "power_w", "activity_pct", "vram_bytes")


def sample():
    raw = subprocess.run(["ioreg", "-a", "-r", "-d", "1", "-c", "IOAccelerator"],
                         capture_output=True, timeout=30).stdout
    rows, seen = [], set()
    for entry in plistlib.loads(raw) if raw else []:
        ps = entry.get("PerformanceStatistics")
        die = entry.get("IORegistryEntryID")
        if not isinstance(ps, dict) or die in seen:
            continue
        seen.add(die)
        rows.append((die, ps.get("Temperature(C)"), ps.get("Total Power(W)"),
                     ps.get("GPU Activity(%)"), ps.get("inUseVidMemoryBytes")))
    return sorted(rows)


def do_log(a):
    new = not os.path.exists(a.out)
    end = time.time() + a.duration if a.duration else None
    with open(a.out, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(FIELDS)
        while True:
            t = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            try:
                for r in sample():
                    w.writerow((t,) + r)
            except Exception as exc:  # noqa: BLE001 -- a failed sample is logged, never silent
                print(f"[gpu_telemetry] WARNING: sample failed at {t}: {exc}", file=sys.stderr)
            f.flush()
            if (end and time.time() >= end) or (a.until_file and os.path.exists(a.until_file)):
                return 0
            time.sleep(a.interval)


def do_summarize(a):
    per = {}
    with open(a.csv) as f:
        for r in csv.DictReader(f):
            if (a.from_ and r["utc"] < a.from_) or (a.to and r["utc"] > a.to):
                continue
            d = per.setdefault(r["die_id"], {"n": 0, "tmax": None, "pmax": None, "amax": None})
            d["n"] += 1
            for key, col in (("tmax", "temperature_c"), ("pmax", "power_w"), ("amax", "activity_pct")):
                if r[col] not in ("", "None"):
                    v = float(r[col])
                    d[key] = v if d[key] is None else max(d[key], v)
    if not per:
        print("no samples in the window")
        return 1
    for die, d in sorted(per.items()):
        print(f"die {die}: {d['n']} samples, peak {d['tmax']} C, peak {d['pmax']} W, "
              f"peak activity {d['amax']}%")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    lg = sub.add_parser("log")
    lg.add_argument("--out", required=True)
    lg.add_argument("--interval", type=float, default=2.0)
    lg.add_argument("--duration", type=float, default=None, help="seconds; default: until stopped")
    lg.add_argument("--until-file", default=None, help="stop when this file exists")
    sm = sub.add_parser("summarize")
    sm.add_argument("csv")
    sm.add_argument("--from", dest="from_", default=None, help="ISO UTC lower bound")
    sm.add_argument("--to", default=None, help="ISO UTC upper bound")
    a = ap.parse_args()
    return do_log(a) if a.mode == "log" else do_summarize(a)


if __name__ == "__main__":
    sys.exit(main())

"""Smoke test of the P23 / P27 / P18 drivers at sizes and lengths far too small to inform
any pass rule (P = 8-32 with 20P measurement rounds; Nowak-May L <= 40 for 60
generations; voter L <= 24). Everything is written to a temporary directory; the real
analysis/outputs/repro_v2 is not touched. Exercises run -> check -> analyze end to end.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

PARTS = {a for a in sys.argv[1:] if a in ("P23", "P27", "P18")} or {"P23", "P27", "P18"}
THEORY_CACHE = next((a for a in sys.argv[1:] if a.endswith(".json")), None)
tmp = Path(tempfile.mkdtemp(prefix="repro_v2_smoke_"))
common.RAW_DIR = tmp / "raw"
common.OUT_DIR = tmp
common.DEVIATIONS_FILE = tmp / "no_deviations.json"

import p18_voter as p18  # noqa: E402
import p23_minority_game as p23  # noqa: E402
import p27_nowak_may as p27  # noqa: E402

for mod in (p23, p27, p18):
    mod.RAW = common.RAW_DIR
p27.INV_FILE = common.RAW_DIR / "P27_breakpoint_invariance.json"
p18.THEORY_FILE = common.RAW_DIR / "P18_theory_voter_exact_crw.json"

# --- P23 at P = 8 / 16 / 32
if "P23" in PARTS:
  p23.GRID = {8: [31, 25, 23, 21, 19, 15, 11, 9], 16: [61, 49, 47, 45, 41, 35, 25, 17]}
  p23.T1_PS, p23.PRIMARY_P = (8, 16), 16
  p23.BURN_MULT, p23.MEAS_MULT = 4, 20
  p23.T2_POINTS = [(32, 31), (32, 15), (16, 15), (16, 7)]
  p23.T2_GATING = [(32, 31), (32, 15)]
  p23.T2_BURN_MULT, p23.T2_MEAS_MULT = 5, 20
  p23.cmd_t1("random", list(range(4)), p23.BURN_MULT, p23._parse_points(None), p23.raw_t1("random", p23.BURN_MULT))
  p23.cmd_stationarity("random", True, list(range(4)))
  p23.cmd_t1("first", list(range(2)), p23.BURN_MULT, p23._parse_points(None), p23.raw_t1("first", p23.BURN_MULT))
  p23.cmd_stationarity("first", True, list(range(2)))
  p23.cmd_t2("linear", list(range(2)), p23.raw_t2("linear"))
  p23.cmd_t2("sign", list(range(2)), p23.raw_t2("sign"))
  p23.cmd_analyze()

# --- P27 at L = 40 / 20 / 10, 60 generations
if "P27" in PARTS:
  p27.L_GATE, p27.L_FS, p27.L_SMALL = 40, 20, 10
  p27.GENS, p27.BURN = 60, 20
  p27.INV_L, p27.INV_GENS = 30, 30
  p27.cmd_invariance()
  for (L, b) in [(40, 1.9), (20, 1.9), (10, 1.9), (40, 1.8)]:
      p27.cmd_run(L, b, list(range(10)))
  p27.cmd_analyze()

# --- P18 at L = 48 / 32 (consensus by t = 3000 in many seeds: uninformative by design)
if "P18" in PARTS:
  if THEORY_CACHE and Path(THEORY_CACHE).exists():          # reuse a cached exact-theory file
      import shutil
      common.RAW_DIR.mkdir(parents=True, exist_ok=True)
      shutil.copy(THEORY_CACHE, p18.THEORY_FILE)
  p18.L_GATE, p18.L_FS = 48, 32
  p18.cmd_run(48, "von_neumann", list(range(20)))
  p18.cmd_precision()
  p18.cmd_run(32, "von_neumann", list(range(10)))
  p18.cmd_run(32, "moore", list(range(10)))
  p18.cmd_analyze()

for name in sorted(PARTS):
    d = json.loads((tmp / f"{name}.json").read_text())
    print(name, "keys:", sorted(d.keys()))
print("smoke outputs in", tmp)

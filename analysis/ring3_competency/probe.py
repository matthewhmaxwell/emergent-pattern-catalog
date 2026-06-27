"""Black-box harness for the agent-observer. The scientist interacts ONLY through this CLI.

  python probe.py info
  python probe.py run <A|B|C> [--start x,y] [--walls "x,y;x,y;..."]

Returns JSON: {"path": [[x,y],...], "reached": bool, "end": [x,y], "steps": n}.
No-walls run = passive observation; with --walls = intervention. Vary --start to probe
convergence from different initial conditions. (Do not open _systems.py — that's the answer key.)
"""
import sys
import json

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from _systems import SYSTEMS, N, TARGET


def _kv():
    d = {}
    a = sys.argv[3:]
    for i in range(0, len(a) - 1, 2):
        d[a[i].lstrip("-")] = a[i + 1]
    return d


if len(sys.argv) < 2:
    print("usage: probe.py info | run <A|B|C> [--start x,y] [--walls 'x,y;...']"); sys.exit(1)

if sys.argv[1] == "info":
    print(json.dumps({"grid": N, "target": list(TARGET), "systems": ["A", "B", "C"],
                      "cells": "(x,y) with 0<=x,y<%d; default start 0,0" % N}))
    sys.exit(0)

sysid = sys.argv[2]
kv = _kv()
start = tuple(int(t) for t in kv.get("start", "0,0").split(","))
walls = []
if kv.get("walls"):
    walls = [tuple(int(t) for t in p.split(",")) for p in kv["walls"].split(";") if p.strip()]
res = SYSTEMS[sysid](walls, start)
res["end"] = res["path"][-1]
res["steps"] = len(res["path"]) - 1
print(json.dumps(res))

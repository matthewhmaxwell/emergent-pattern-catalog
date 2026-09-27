"""Pilot T5 (first pass): recovery benchmark — run every hand-encoded textbook model through the census filter.
Pilot only; nothing here is census data."""
import sys, time, json
from census.bench_models import BENCH
from census import sim, filter as FL

EXPECT = {"voter": "P18", "voter p=0.3": "P18", "voter p=0.1": "P18", "majority coarsen": "P18", "Schelling": "P1", "cyclic CA spirals": "P12",
          "BTW sandpile": "P14", "Nowak-May PD": "P27", "lattice Kuramoto": "P9", "Vicsek": "P5",
          "MIPS (quorum)": "P2", "net Kuramoto": "P9", "co-evolving voter": "P34", "Greenberg-Hastings": "P13",
          "Langton ant(s)": None, "Gray-Scott": "P3", "chemotaxis (KS)": None, "RPS agents": "P12"}
only = sys.argv[1:]
rows = {}
for name, p in BENCH.items():
    if only and name not in only: continue
    t = time.time(); out = sim.run(p, seed0=1); ts = time.time() - t
    t = time.time(); ev = FL.evaluate(out, p, battery_all=True); tf = time.time() - t
    summ = {v: (e["verdict"], e["pattern"], e["emergent_seeds"],
                [(r.get("pattern"), r.get("verdict"), r.get("em_score"), r.get("em_kind"),
                  r.get("fired"), r.get("confirmed")) for r in e["seeds"]]) for v, e in ev.items()}
    rows[name] = {"expect": EXPECT.get(name), "sim_s": round(ts, 1), "filter_s": round(tf, 1), "views": summ}
    print(f"\n### {name}  (expect {EXPECT.get(name)})  sim {ts:.1f}s  filter {tf:.1f}s", flush=True)
    for v, (verd, pat, em, seeds) in summ.items():
        print(f"  view {v:8s} -> {verd:12s} {pat or '':5s} emergent_seeds={em}", flush=True)
        for sd in seeds: print("      seed:", sd, flush=True)
json.dump(rows, open("census/pilot_recovery_v0.json", "w"), indent=1, default=str)

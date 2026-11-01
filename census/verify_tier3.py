"""Tier-3 interventional controls for the undocumented candidate qd:060, as fixed in VERIFICATION_TIER3.md (written and
committed before this script ran). Census settings, the verification seeds of qd:060, 10 seeds, frozen fingerprint.

  python -m census.verify_tier3        # writes census/runs/verify/tier3/qd060_controls.json and prints the summary
"""
import dataclasses, json, os
import numpy as np
from census import verify as V, sim, knockout as KO

PID, VIEW = "qd:060", "C"
OUT = os.path.join(V.OUT, "tier3")


def _fps(out, prog):
    return [s["fp"] for s in KO.per_seed_view_results(out, prog).get(VIEW, [])]


def main():
    os.makedirs(OUT, exist_ok=True)
    it = {x["pid"]: x for x in json.load(open(f"{V.OUT}/inventory.json"))}[PID]
    cl = {x["pid"]: x for x in json.load(open(f"{V.OUT}/verify.json"))}[PID]
    p = V.grammar(it["grammar"]).decode(it["bits"]); seed = V.fresh_seed(it["bits"]); S = V.S_SEEDS
    V._set_size("census")
    full_out = sim.run(p, seed0=seed, S=S)
    full = _fps(full_out, p)
    # agents' distinct occupied cells at the last frame of the full program (sets the pinned-control agent count)
    pos = full_out["pos"][:, -1]
    occ = [len({(int(y) % sim.W, int(x) % sim.W) for y, x in pos[s]}) for s in range(S)]
    n_occ = int(np.median(occ))
    res = {"full": full, "occupied_cells_last_frame": occ, "pinned_agents": n_occ}
    # N: null (structure destroyed after every step)
    res["null"] = _fps(sim.run(p, seed0=seed, S=S, null_shuffle=True), p)
    # S: provably symmetric, field exactly uniform so the gradient is zero everywhere at every step
    orig = sim.State.__init__
    def uniform_init(self, *a, **k):
        orig(self, *a, **k)
        if hasattr(self, "F"): self.F[:] = 0.0; self.F[:, 0] = 1.0
    sim.State.__init__ = uniform_init
    try: res["symmetric"] = _fps(sim.run(p, seed0=seed, S=S), p)
    finally: sim.State.__init__ = orig
    # P: pinned flippers, speed 0, as many agents as distinct occupied cells at the end of the full run
    pp = dataclasses.replace(p, agent=(0.0, p.agent[1])); na = sim.NA; sim.NA = n_occ
    try: res["pinned"] = _fps(sim.run(pp, seed0=seed, S=S), pp)
    finally: sim.NA = na
    # reproduction test per robust feature (plan: |median(full) - median(control)| <= 3 s, s pooled)
    summary = {}
    for f in cl["robust"]:
        x = V._col({VIEW: full}, f, VIEW); row = {"full": float(np.median(x))}
        for c in ("null", "symmetric", "pinned"):
            y = V._col({VIEW: res[c]}, f, VIEW)
            if y is None: row[c] = None; continue
            row[c] = {"median": float(np.median(y)), "reproduces": bool(abs(np.median(x) - np.median(y)) <= 3 * V._scale(x, y))}
        summary[f] = row
    # sub-program 1+2 accounts for the activity decline? (same sign, at least half the full program's size)
    rec = json.load(open(f"{V.OUT}/sims/{PID.replace(':', '__')}.json"))
    sp12 = [sp for sp in rec["subprograms"] if sp["rules"] == [0, 1]][0]
    a12 = float(np.median(V._col(sp12["fps"], "g_activity_trend", VIEW)))
    af = summary["g_activity_trend"]["full"]
    res["summary"] = summary
    res["subprogram_12_activity_trend"] = a12
    res["subprogram_12_accounts_for_decline"] = bool(np.sign(a12) == np.sign(af) and abs(a12) >= 0.5 * abs(af))
    order = ["g_largest_domain", "g_n_domains", "g_corrlen", "g_reversals"]
    res["pinned_reproduces_ordering"] = all(summary[f]["pinned"] and summary[f]["pinned"]["reproduces"] for f in order if f in summary)
    res["H0_supported"] = res["pinned_reproduces_ordering"] and res["subprogram_12_accounts_for_decline"]
    json.dump(res, open(f"{OUT}/qd060_controls.json", "w"), indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("occupied_cells_last_frame", "pinned_agents", "summary", "subprogram_12_activity_trend",
                                          "subprogram_12_accounts_for_decline", "pinned_reproduces_ordering", "H0_supported")}, indent=1))


if __name__ == "__main__":
    main()

"""Ring-2 Stage-0: baseline-validate the model-free novelty tripwire on EXISTING
material (the blind-spot probe corpus: known-emergent controls + nulls), before any
sourcing. Success = the tripwire is QUIET on the known corpus:
  - nulls fall below the complexity floor (not flagged complex), and
  - known-emergent systems are CLASSIFIED by the named lenses (not tripped).
Anything complex-but-unclassified among the KNOWNS is a residual blind spot (a
Tier-3 lens candidate even within our own corpus) and is reported honestly.
"""
from __future__ import annotations

from analysis.blind_spot_probes import PROBES
from epc.phase2a.novelty_tripwire import novelty_tripwire, C_THR, PSI_THR, EM_THR


def run():
    rows = []
    for probe in PROBES:
        try:
            p = probe(0)
            tw = novelty_tripwire(p["history"])
            rows.append({"name": probe.__name__, "truth": p.get("truth"),
                         "fam": p.get("family"), **tw})
        except Exception as e:
            rows.append({"name": probe.__name__, "error": repr(e)[:90]})

    ok = [r for r in rows if "error" not in r]
    nulls = [r for r in ok if r["truth"] == "null"]
    ems = [r for r in ok if r["truth"] == "emergent"]
    maxC = max((r["C"] or 0) for r in nulls) if nulls else 0.0
    maxpsi = max(r["psi"] for r in nulls) if nulls else 0.0
    print(f"NULL calibration: max C={maxC:.3f}  max psi={maxpsi:.3f}  "
          f"-> suggest C_THR={maxC + 0.03:.3f}  PSI_THR={max(0.05, maxpsi + 0.02):.3f}")
    print(f"current constants: C_THR={C_THR}  PSI_THR={PSI_THR}  EM_THR={EM_THR}\n")

    print(f"{'probe':<24}{'truth':<9}{'C':>7}{'psi':>7}{'em':>6}{'cx':>4}{'cls':>4}{'TRIP':>5}  reason")
    for r in rows:
        if "error" in r:
            print(f"{r['name']:<24} ERROR {r['error']}"); continue
        print(f"{r['name']:<24}{r['truth']:<9}{(r['C'] if r['C'] is not None else 0):>7.3f}"
              f"{r['psi']:>7.3f}{r['em_score']:>6.2f}{'Y' if r['is_complex'] else '.':>4}"
              f"{'Y' if r['classified'] else '.':>4}{'Y' if r['tripped'] else '.':>5}  {r['reason']}")

    null_trip = [r['name'] for r in nulls if r["tripped"]]
    known_trip = [r['name'] for r in ems if r["tripped"]]
    print(f"\nnulls tripped (want 0):            {len(null_trip)}  {null_trip}")
    print(f"known-emergent classified:         {sum(1 for r in ems if not r['tripped'])}/{len(ems)}")
    print(f"complex-but-unclassified KNOWNS:   {known_trip}  (residual blind spots, if any)")
    quiet = (len(null_trip) == 0) and (len(known_trip) == 0)
    print(f"\nBRIDGE BASELINE: {'QUIET on the known corpus — calibrated.' if quiet else 'NOT yet quiet — see trips above.'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())

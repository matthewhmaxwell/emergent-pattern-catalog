"""Validate the live Ring-2 descriptor: run it across the probe corpus + the 3 new
substrate families and confirm it produces clean, substrate-AWARE fingerprints (the right
lenses fire per substrate; no crashes; named em + tripwire consistent with their own
validations). This is the integration check that the admitted lenses are wired live."""
import numpy as np
from analysis.blind_spot_probes import PROBES
from epc.phase2a.ring2_descriptor import ring2_descriptor
from epc.models.reaction_diffusion import reaction_diffusion
from epc.models.kuramoto_lattice import kuramoto_lattice
from epc.models.lenia import lenia

rows = []
for p in PROBES:
    pr = p(0)
    rows.append((p.__name__, pr.get("truth"), ring2_descriptor(pr["history"], {})))
rows.append(("RD_stripes", "emergent", ring2_descriptor(reaction_diffusion(0, F=0.030, k=0.057, N=96, steps=8000, record=24), {})))
rows.append(("kuramoto_sync", "emergent", ring2_descriptor(kuramoto_lattice(0, N=40, K=3.0, steps=3000, record=120), {})))
rows.append(("lenia_glider", "emergent", ring2_descriptor(lenia(0, N=64, steps=400, record=80, n_creatures=1), {})))

ABBR = {"structure_factor": "SF", "persistent_homology": "PH", "graph_structure": "GS",
        "directed_info_flow": "DTE"}
print(f"{'system':<24}{'truth':<10}{'em':>5} {'kind':<14}{'cx':>3}{'trp':>4}  lenses_fired")
for nm, t, d in rows:
    fired = "+".join(ABBR.get(x, x) for x in d["lenses_fired"]) or "(none)"
    cx = "Y" if d["mf_complex"] else "-"; trp = "Y" if d["tripped"] else "-"
    print(f"{nm:<24}{str(t):<10}{(d['em_score'] or 0):>5.2f} {str(d['em_kind']):<14}{cx:>3}{trp:>4}  {fired}")

# coverage: which lenses fire on how many systems; any system with NO Tier-2 lens?
from collections import Counter
cnt = Counter(x for nm, t, d in rows for x in d["lenses_fired"])
none_t2 = [nm for nm, t, d in rows if not d["lenses_fired"]]
print(f"\nlens fire counts (of {len(rows)} systems): {dict(cnt)}")
print(f"systems with NO Tier-2 lens firing: {none_t2 or 'none'}")
print(f"nulls tripped (want 0): {[nm for nm,t,d in rows if t=='null' and d['tripped']]}")
print(f"known-emergent classified: {sum(1 for nm,t,d in rows if t=='emergent' and d['classified'])}"
      f"/{sum(1 for nm,t,d in rows if t=='emergent')}")

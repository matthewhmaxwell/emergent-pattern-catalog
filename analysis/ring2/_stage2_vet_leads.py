"""Vet the 4 Stage-2 Lenia tripwire leads. Params are reproducible from the seed (gen(i)
used rng(1000+i) in the exact order nc/mu/sigma). For each: re-run, inspect the mass
trajectory + final field (is it a sustained structure, or a dead/dying transient?), and the
descriptor. A genuine lead = sustained non-trivial structure that no named lens classifies;
a mundane lead = a death/explosion transient that is temporally complex but spatially trivial."""
import numpy as np
from epc.models.lenia import lenia
from epc.phase2a.ring2_descriptor import ring2_descriptor

for i in [0, 19, 29, 34]:
    rng = np.random.default_rng(1000 + i)
    nc = int(rng.integers(1, 6)); mu = float(rng.uniform(0.10, 0.20)); sig = float(rng.uniform(0.012, 0.020))
    h = lenia(i, N=64, mu=mu, sigma=sig, steps=400, record=60, n_creatures=nc)
    mass = np.array([float(f["field"].sum()) for f in h])
    last = h[-1]["field"]
    # centroid path (is anything moving / is it localized?)
    coms = []
    for f in h:
        A = f["field"]; m = A.sum()
        if m > 1e-6:
            y, x = np.indices(A.shape); coms.append((float((A * y).sum() / m), float((A * x).sum() / m)))
    pl = sum(np.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(coms, coms[1:])
             if np.hypot(b[0] - a[0], b[1] - a[1]) < 10) if len(coms) > 1 else 0.0
    d = ring2_descriptor(h, {})
    # fraction of grid that is "on" (active) at the end
    active = float((last > 0.05).mean())
    print(f"#{i}  nc={nc} mu={mu:.3f} sig={sig:.4f}")
    print(f"    mass: start {mass[0]:.0f} -> mid {mass[len(mass)//2]:.0f} -> end {mass[-1]:.0f}  "
          f"(end/start {mass[-1]/(mass[0]+1e-9):.2f})")
    print(f"    final field: mean {last.mean():.3f} std {last.std():.3f} max {last.max():.3f} active_frac {active:.3f}")
    print(f"    centroid path {pl:.0f} | em={d['em_score']} C={d['mf_C']} struct={d['mf_struct']} "
          f"cx={d['mf_complex']} TRIPPED={d['tripped']}")
    if last.std() < 0.02 or active < 0.005:
        verdict = "DEAD/uniform final state -> MUNDANE (death transient is temporally complex, spatially trivial)"
    elif active > 0.9:
        verdict = "near-FULL grid (saturated/exploded) -> MUNDANE"
    else:
        verdict = "SUSTAINED non-trivial structure -> inspect further (possible genuine lead)"
    print(f"    VERDICT: {verdict}\n")

"""Substrate run 3: exercise the lens battery on Lenia (autonomous moving creatures) —
the qualitative regime RD/Kuramoto lacked, and ASAL's own search space. Sharp question:
does the emergence indicator's NAMED channels recognise a lone gliding creature, or does
it slip through (a locomotion blind spot)? If complex-but-unclassified, the tripwire should
fire — a (known-phenomenon) novelty lead illustrating the bridge working as intended.
structure_factor on the field; persistent_homology on creature-occupancy positions."""
import numpy as np
from epc.models.lenia import lenia
from epc.metrics.structure_factor import structure_factor_peak
from epc.metrics.persistent_homology import persistent_homology
from epc.phase2a.novelty_tripwire import novelty_tripwire

REGIMES = [
    ("single_glider", dict(n_creatures=1)),
    ("multi_glider",  dict(n_creatures=5)),
    ("soup",          dict(soup=True)),
]


def stats(h):
    masses, coms = [], []
    for f in h:
        A = f["field"]; m = A.sum(); masses.append(m)
        if m > 1e-6:
            y, x = np.indices(A.shape); coms.append((float((A * y).sum() / m), float((A * x).sum() / m)))
        else:
            coms.append(None)
    pl = sum(np.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(coms, coms[1:])
             if a and b and np.hypot(b[0] - a[0], b[1] - a[1]) < 10)
    return masses[-1] / (masses[0] + 1e-9), pl


print(f"{'regime':<14}{'alive':>6}{'move':>7}{'sk_peak':>9}{'h1_max':>8}{'em':>6} {'em_kind':<12}{'cplx':>5}{'trip':>5}")
for name, kw in REGIMES:
    h = lenia(0, N=64, steps=400, record=80, **kw)
    alive, pl = stats(h)
    sf = structure_factor_peak(h)                       # field view
    pos_hist = [{"positions": np.argwhere(f["field"] > 0.3).astype(float)} for f in h]
    ph = persistent_homology(pos_hist)                  # creature-occupancy positions view
    tw = novelty_tripwire(h)
    sk = sf.get("sk_peak") if sf else None
    h1 = ph.get("h1_max") if ph else None
    print(f"{name:<14}{alive:>6.2f}{pl:>7.0f}"
          f"{(f'{sk:9.1f}' if sk is not None else '      n/a')}"
          f"{(f'{h1:8.3f}' if h1 is not None else '     n/a')}"
          f"{tw['em_score']:>6.2f} {str(tw['em_kind']):<12}"
          f"{('Y' if tw['is_complex'] else '-'):>5}{('Y' if tw['tripped'] else '-'):>5}")
print("\nalive = final/initial mass; move = centroid path length. If a gliding creature is"
      " complex but em<0.5 (unclassified) -> tripwire fires = the bridge flagging structure the"
      " named lenses have no word for (locomotion gap). ASAL substrate: our interpretable axes.")

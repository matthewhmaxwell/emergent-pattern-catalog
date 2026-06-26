"""Substrate run 2: exercise the lens battery on a Kuramoto 2D lattice, OOD from the
catalog corpus and especially stressing directed_info_flow (causal direction) — the
newest admitted lens, which RD could not reach. Expectation: travelling/spiral WAVES have
directional phase coupling (higher directionality) vs symmetric global SYNC vs uncoupled
INCOHERENCE. Also report the Kuramoto order parameter r (sanity) + emergence/tripwire."""
import numpy as np
from epc.models.kuramoto_lattice import kuramoto_lattice, order_parameter
from epc.metrics.directed_info_flow import directed_transfer_entropy
from epc.phase2a.novelty_tripwire import novelty_tripwire

REGIMES = [
    ("incoherent", dict(K=0.2, omega_spread=1.5, init="random")),
    ("sync",       dict(K=3.0, omega_spread=0.0, init="random")),
    ("spiral",     dict(K=1.5, omega_spread=0.0, init="spiral")),
    ("plane_wave", dict(K=1.5, omega_spread=0.0, init="plane")),
]

print(f"{'regime':<12}{'r':>6}{'mte_grid':>9}{'dir_grid':>9}{'mte_row':>9}{'dir_row':>9}{'em':>6} {'em_kind':<12}")
res = []
for name, kw in REGIMES:
    h = kuramoto_lattice(0, N=40, steps=3000, record=120, **kw)
    r = order_parameter(h)
    N = h[0]["phases"].shape[0]
    grid = directed_transfer_entropy(np.array([np.sin(f["phases"]).ravel() for f in h]))      # spatially-blind subsample
    row = directed_transfer_entropy(np.array([np.sin(f["phases"])[N // 2, :] for f in h]))     # ordered along x (flow axis)
    tw = novelty_tripwire(h)
    res.append((name, row["directionality"] if row else None))
    g = lambda d, k: (f'{d[k]:9.3f}' if d else '      n/a')
    print(f"{name:<12}{r:>6.2f}{g(grid,'mean_te')}{g(grid,'directionality')}"
          f"{g(row,'mean_te')}{g(row,'directionality')}{tw['em_score']:>6.2f} {str(tw['em_kind']):<12}")

print("\nKEY: dir_grid = spatially-blind ravel subsample; dir_row = components ordered along x"
      " (the plane wave's propagation axis). If dir_row separates the travelling wave from"
      " sync/incoherent while dir_grid does not, the lens DOES see propagation -- but only when"
      " components are ordered along the flow axis (scope refinement, not a failure).")
pw = [d for n, d in res if n == "plane_wave" and d is not None]
sx = [d for n, d in res if n in ("sync", "incoherent") and d is not None]
if pw and sx:
    print(f"dir_row  plane_wave {pw[0]:.3f}  vs  sync+incoherent max {max(sx):.3f}  -> "
          f"{'SEES propagation when ordered' if pw[0] > max(sx) else 'still mixed'}")

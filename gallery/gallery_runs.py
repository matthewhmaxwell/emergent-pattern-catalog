"""Per-pattern gallery run overrides: finer/longer recording so the emergence
is captured at enough resolution to play smoothly (frame-adequacy fixes)."""
from __future__ import annotations
import numpy as np

def _p9():   # Kuramoto: shorter+finer so incoherent->synced fills the clip
    from epc.models.kuramoto import KuramotoModel, KuramotoParams
    m = KuramotoModel(KuramotoParams(N=300, K=8.0, gamma=0.5, dt=0.05, seed=0))
    return [m.run(n_steps=2500, record_every=5)], m.get_metadata()

def _p29():  # Ant trail: finer snapshots so the network forms gradually
    from epc.models.trail_network import AntTrailModel, AntTrailParams
    m = AntTrailModel(AntTrailParams(n_nodes=7, n_agents=40, grid_size=100,
        alpha=2.0, beta=3.0, deposition_rate=10.0, evaporation_rate=0.05,
        n_steps=600, snapshot_interval=4, seed=42))
    return [m.simulate()], m.get_metadata()

def _p3():   # Gray-Scott: longer run, fine recording so the labyrinth spreads
    from epc.models.gray_scott import GrayScott
    m = GrayScott(rows=64, cols=64, feed_rate=0.037, kill_rate=0.060, seed=0); m.setup()
    h = []
    for t in range(3000):
        st = m.step()
        if t % 30 == 0:
            h.append({"field": np.asarray(st["field"], dtype=np.float32), "step": t})
    return [h], {}

GALLERY_RUNS = {"p3": _p3, "p9": _p9, "p29": _p29}

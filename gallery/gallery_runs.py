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

def _p11():  # Lotka-Volterra: smaller lattice -> larger, more visible anti-correlated population cycles
    from epc.models.lotka_volterra_lattice import LotkaVolterraLattice
    m = LotkaVolterraLattice(rows=40, cols=40, predation_rate=4.0,
        prey_reproduction_rate=1.0, predator_death_rate=1.0, seed=1)
    return [m.run(n_steps=1200)], m.get_metadata()

# NOTE: P24 keeps the validated PROPORTIONAL homeostat (build_p24). The IntegralHomeostat was
# tried for zero-offset return but is unstable at every gain (oscillatory divergence); proportional
# rejects ~80% of the perturbation with a small residual offset — legitimate regulation, shown honestly.

def _p16():  # Hopfield: store RECOGNISABLE 10x10 glyphs so recall is a VISIBLE clean-up of a corrupted letter
    rng = np.random.default_rng(3)
    def g(rows): return np.array([[1.0 if c == "#" else -1.0 for c in r] for r in rows]).ravel()
    glyphs = [
        g(["..######..", ".##....##.", "##......##", "##......##", "##......##",
           "##......##", "##......##", "##......##", ".##....##.", "..######.."]),   # O
        g(["##########", "##########", "....##....", "....##....", "....##....",
           "....##....", "....##....", "....##....", "....##....", "....##...."]),   # T
        g(["#........#", "##......##", ".##....##.", "..##..##..", "...####...",
           "...####...", "..##..##..", ".##....##.", "##......##", "#........#"]),   # X
        g(["....##....", "....##....", "....##....", "##########", "##########",
           "##########", "....##....", "....##....", "....##....", "....##...."]),   # +
    ]
    pats = np.array(glyphs); N = pats.shape[1]
    W = sum(np.outer(p, p) for p in pats) / N
    np.fill_diagonal(W, 0.0)
    target = pats[0].copy(); s = target.copy()
    s[rng.choice(N, size=int(0.28 * N), replace=False)] *= -1     # 28%-corrupted cue
    sp = pats.astype(int)
    def rec(t, conv):
        return {"state": s.copy().astype(int), "overlap": float((s * target).mean()), "step": t,
                "trial": 0, "cue_pattern_idx": 0, "converged": conv, "stored_patterns": sp}
    hist = [rec(0, False)]
    order = rng.permutation(N); k = 0
    for t in range(1, 60):
        for _ in range(6):                                       # 6 async neuron updates per frame -> gradual clean-up
            i = order[k % N]; k += 1
            s[i] = 1.0 if (W[i] @ s) >= 0 else -1.0
        ov = float((s * target).mean()); hist.append(rec(t, ov > 0.999))
        if ov > 0.9999 and t > 8: break
    for t in range(len(hist), len(hist) + 4): hist.append(rec(t, True))   # hold the recalled glyph
    return [hist], {}

def _p12():  # RPS: larger lattice + higher mobility (toward M_c~4.5e-4) -> active swirling cyclic domains, not frozen coarsening
    from epc.models.rps_spatial import RPSSpatialModel
    return [RPSSpatialModel(rows=90, cols=90, mobility=3e-4, seed=2).run(n_steps=400)], {}

def _p2():   # MIPS: higher packing (phi=0.7) so the active gas clearly phase-separates into dense+dilute domains
    from epc.models.active_brownian_particles import ActiveBrownianParticles
    N = 1000; phi = 0.7; box = float(np.sqrt(N * np.pi / 4.0 / phi))
    m = ActiveBrownianParticles(n_particles=N, box_size=box, v0=1.0, D_r=0.01,
        rho_star=4.0, r_cg=1.0, dt=0.05, init_mode="uniform", seed=1)
    return [m.run(n_steps=2000)], m.get_metadata()

GALLERY_RUNS = {"p3": _p3, "p9": _p9, "p29": _p29, "p11": _p11, "p16": _p16, "p12": _p12, "p2": _p2}

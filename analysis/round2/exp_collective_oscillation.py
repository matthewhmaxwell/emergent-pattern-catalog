"""Round-2 experiment: collective-oscillation channel co-design.

Pairs a genuinely SYNCHRONIZING probe (mean-field Kuramoto above critical coupling)
with the coordination-gated complexity detector (mpr_emergence). Confirms the
detector fires on emergent collective oscillation (K>Kc) and NOT on the trivial
independent-oscillator null (K=0) — the distinction that broke the naive
"is the macro periodic?" approach.
"""
from __future__ import annotations

import numpy as np

from epc.phase2a.emergence import generic_emergence
from epc.phase2a.info_channels import mpr_emergence


def kuramoto(seed: int, N: int = 60, K: float = 4.0, steps: int = 260,
            dt: float = 0.05, spread: float = 0.4):
    """Mean-field Kuramoto. K>~2*spread synchronizes (collective oscillation);
    K=0 is independent rotators (no emergence). Frames emit cos(phase) as a
    'state' vector so the generic micro_macro extractor sees per-oscillator series."""
    rng = np.random.default_rng(seed)
    omega = rng.normal(1.0, spread, N)
    phi = rng.uniform(-np.pi, np.pi, N)
    H, micro = [], []
    for _ in range(steps):
        z = np.mean(np.exp(1j * phi))
        phi = phi + dt * (omega + K * np.abs(z) * np.sin(np.angle(z) - phi))
        c = np.cos(phi)
        H.append({"state": c.copy()})
        micro.append(c.copy())
    return H, np.array(micro)


def run():
    print(f"{'condition':<22}{'coord_z':>9}{'C_obs':>8}{'generic_em':>12}")
    print("-" * 51)
    for K, label, truth in [(4.0, "synced(emergent)", "emergent"),
                            (0.0, "uncoupled(null)", "null")]:
        for s in (0, 1, 2):
            H, micro = kuramoto(s, K=K)
            z, c = mpr_emergence(micro)
            em = generic_emergence(H)["score"]
            # sanity: order parameter at the end
            print(f"K={K:<4} {label:<14} s{s}  z={z:>6.2f}  C={c:>6.3f}  em={em:>6.3f}")
    print("\nWANT: synced z>3 (fires), uncoupled z<3 (does not). Then the "
          "coordination-gate temporal channel is safe to wire.")


if __name__ == "__main__":
    run()

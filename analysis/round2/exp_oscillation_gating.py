"""Round-2: does a spectral-peak oscillation channel, GATED to exclude phase-kind
observations, cleanly catch scalar/field limit cycles without false-firing on the
T2c nulls? Diagnoses the per-source frame 'kind' + oscillation_score so we can set
a safe gate + threshold."""
from __future__ import annotations

import numpy as np

from epc.phase2a.info_channels import micro_macro, oscillation_score
from analysis import t2c_systems as T
from analysis import blind_spot_probes as B


def frame_kind(history):
    f = history[-1]
    for k in ("grid", "field", "velocities", "positions", "phases", "theta"):
        if k in f:
            return k
    for k in ("state", "opinions", "wealth", "attendance", "x", "fraction_on",
              "task_assignments", "cell_types"):
        if k in f:
            return f"vec:{k}"
    return "?"


def osc_of(history):
    M, _ = micro_macro(history)
    if M is None:
        return None
    return oscillation_score(M.mean(1))


print(f"{'source':<26}{'truth':<9}{'kind':<14}{'osc(mean-micro)':>16}")
print("-" * 65)
# T2c nulls (must NOT fire)
for name, fn in [("null_spatial_noise", T.null_spatial_noise),
                 ("null_random_walk", T.null_random_walk),
                 ("null_uncoupled_phases", T.null_uncoupled_phases),
                 ("null_frozen_noise", T.null_frozen_noise)]:
    h, _ = fn(0)
    print(f"{name:<26}{'null':<9}{frame_kind(h):<14}{(osc_of(h) or float('nan')):>16.2f}")
# emergent temporal probes (should fire if non-phase)
for name, fn in [("limit_cycle", B.limit_cycle), ("spatiotemporal_chaos", B.spatiotemporal_chaos),
                 ("traveling_wave", B.traveling_wave)]:
    p = fn(0)
    print(f"{name:<26}{'emergent':<9}{frame_kind(p['history']):<14}{(osc_of(p['history']) or float('nan')):>16.2f}")
print("\nGATE PLAN: run oscillation channel only when kind is NOT phases/theta "
      "(phase sync already handled by the order-parameter channel).")

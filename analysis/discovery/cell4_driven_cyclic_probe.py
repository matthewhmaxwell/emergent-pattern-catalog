"""Cell-4 probe: ext_forced x nontransitive (driven cyclic dominance).

Compare autonomous spatial RPS (P12) vs the SAME system under external periodic
cyclic forcing. If the driven version reads the same as autonomous (no distinct
new emergent structure the instrument flags), the cell is an evidenced non-prediction.
"""
import numpy as np
from epc.models.rps_spatial import RPSSpatialModel
from analysis.battery_profile import profile_observation
from epc.phase2a.emergence import generic_emergence


def run_rps(seed, driven, L=80, gens=300, rec=10, tau=40, strength=0.06, mobility=5e-4):
    m = RPSSpatialModel(rows=L, cols=L, mobility=mobility, init_mode="random", seed=seed)
    hist = [m.setup()]
    rng = np.random.default_rng(seed * 991 + 1)
    for g in range(1, gens + 1):
        if driven and g % tau == 0:
            fav = (g // tau) % 3 + 1                 # rotating favored species (external drive)
            mask = rng.random((L, L)) < strength
            m._grid[mask] = fav
        st = m.step()
        if g % rec == 0:
            hist.append(st)
    return hist


def show(name, hist):
    em = generic_emergence(hist)
    r = profile_observation(hist, {}, match_min_tier="confirmation")
    v = r["verdict"]
    fired = [(row["pattern_id"], row["detector_tier"]) for row in r["profile"]
             if row.get("fired_detected")]
    print(f"{name:<16} verdict={v['verdict']:<24} pid={v.get('pattern_id','-'):<5} "
          f"em={em['score']:.2f} kind={em.get('kind')} fired={fired[:5]}", flush=True)


if __name__ == "__main__":
    print("compare: does external cyclic forcing create a DISTINCT emergent pattern?", flush=True)
    for s in range(2):
        show(f"autonomous s{s}", run_rps(s, driven=False))
        show(f"driven s{s}", run_rps(s, driven=True))

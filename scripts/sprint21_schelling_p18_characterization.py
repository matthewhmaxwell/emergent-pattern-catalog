"""Sprint 21 — Schelling × P18 content-level characterization.

Closes Sprint 20 carry-forward #20 ("Schelling × P18 content-level
negative test"). Empirically tests whether the P18 detector rejects
Schelling segregation on metric grounds at canonical Schelling
parameters (L = 64, density = 0.9, threshold = 0.375), without
relying on the metadata-flag-based P1 exclusion.

Findings (5-seed audit):
  - All 5 seeds remain at SCREENING tier (4 detected = False, 1 detected
    = True at SCREENING because wall_final_qtr ≈ 0.36 > 0.30
    confirmation ceiling). The metric-only rejection holds.
  - The Sprint 20 detector docstring claim "Schelling P1 saturates to
    wall ~0.02" is empirically wrong: Schelling's three-state grid
    yields wall_final_qtr ≈ 0.36 across all seeds because empty cells
    contribute boundary segments under the Moore neighborhood.
  - The Sprint 20 §6.10 claim "Schelling moran_growth below 0.20 floor"
    is also wrong: moran_growth values were 0.23-0.30, all above 0.20.
  - The actual rejection mechanism is two-pronged: moran_final_qtr ≤
    0.30 fails screening (4 seeds), or wall_final_qtr ≥ 0.30 fails the
    confirmation ceiling (the 5th seed).

A separate exploratory check at Schelling threshold = 0.5
(strong-segregation, sometimes cited in textbook expositions) found
that all 5 seeds reach P18 DEFINITIVE with P1 marked "inconclusive".
This false-positive is recorded as Sprint 21 carry-forward #20b in
REPLICATION_NOTES.md.

Output: schelling_p18_characterization.json with full per-seed records.

Run: python scripts/sprint21_schelling_p18_characterization.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from epc.detectors.p18_consensus import P18ConsensusDetector
from epc.models.schelling import run_schelling


def characterize_seed(seed: int, L: int = 64, n_steps: int = 200) -> dict:
    """Run Schelling at one seed, push through P18, capture all metrics.

    We also compute the wall trajectory metrics directly (independent of
    detector), so they're available even when the detector early-exits at
    screening with an empty secondary_metrics dict.
    """
    from epc.detectors.p18_consensus import _moran_i_moore, _wall_density_moore
    import numpy as np

    history = run_schelling(grid_size=L, density=0.9, threshold=0.375,
                            n_steps=n_steps, seed=seed)

    detector = P18ConsensusDetector(n_permutations=199, seed=42)
    # Schelling-style metadata: no 'update' key with copy/imitation/voter.
    metadata = {'threshold': 0.375, 'density': 0.9}
    result = detector.detect(history, model_metadata=metadata)

    primary = result.primary_metric or {}
    secondary = result.secondary_metrics or {}

    # Direct trajectory diagnostic — always computed even when detector
    # early-exits at screening. Mirrors the detector's internal calc.
    grids = [np.asarray(s['grid']).astype(np.int8) for s in history]
    moran_traj = np.array([_moran_i_moore(g) for g in grids])
    wall_traj = np.array([_wall_density_moore(g) for g in grids])
    n_qtr = max(1, len(grids) // 4)

    direct = {
        'direct_wall_initial': float(wall_traj[0]),
        'direct_wall_final': float(wall_traj[-1]),
        'direct_wall_final_qtr_mean': float(wall_traj[-n_qtr:].mean()),
        'direct_wall_decay': float(wall_traj[0] - wall_traj[-1]),
        'direct_moran_initial': float(moran_traj[0]),
        'direct_moran_final': float(moran_traj[-1]),
        'direct_moran_final_qtr_mean': float(moran_traj[-n_qtr:].mean()),
    }

    return {
        'seed': seed,
        'L': L,
        'n_steps': n_steps,
        'detected': result.detected,
        # Primary metrics
        'moran_spearman_early': primary.get('moran_spearman_early'),
        'moran_initial': primary.get('moran_initial'),
        'moran_final': primary.get('moran_final'),
        'moran_growth': primary.get('moran_growth'),
        'moran_final_qtr_mean': primary.get('moran_final_qtr_mean'),
        # Secondary metrics (empty if screening failed)
        'wall_initial': secondary.get('wall_initial'),
        'wall_final': secondary.get('wall_final'),
        'wall_decay': secondary.get('wall_decay'),
        'wall_spearman_early': secondary.get('wall_spearman_early'),
        'wall_final_qtr_mean': secondary.get('wall_final_qtr_mean'),
        'minority_fraction_final': secondary.get('minority_fraction_final'),
        # Null
        'null_p_value': result.null_p_value,
        # Tier
        'tier': result.tier.name if result.tier else None,
        'confidence': result.confidence,
        # Exclusions
        'exclusion_results': dict(result.exclusion_results) if result.exclusion_results else {},
        # Direct trajectory diagnostic (always populated)
        **direct,
    }


def evaluate_gates(rec: dict) -> dict:
    """Manually check which gates Schelling fails, for diagnostic purposes.

    Falls back to direct trajectory metrics when detector secondaries are
    empty (the screening-fail case).
    """
    det = P18ConsensusDetector
    msp = rec['moran_spearman_early'] or 0.0
    mfm = rec['moran_final_qtr_mean'] or 0.0
    mgr = rec['moran_growth'] or 0.0
    # Wall metrics may be missing if screening failed; fall back to direct.
    wsp = rec['wall_spearman_early']  # only computed by detector
    wfm = rec['wall_final_qtr_mean'] if rec['wall_final_qtr_mean'] is not None \
        else rec['direct_wall_final_qtr_mean']
    wdc = rec['wall_decay'] if rec['wall_decay'] is not None \
        else rec['direct_wall_decay']
    mfq = mfm
    mff = rec['minority_fraction_final']  # may be None
    p = rec['null_p_value']

    gates = {}
    # Screening
    gates['screen_moran_spearman_>0.70'] = msp > det.SCREENING_MORAN_SPEARMAN_MIN
    gates['screen_moran_final_>0.30'] = mfm > det.SCREENING_MORAN_FINAL_MIN
    gates['screen_moran_growth_>0.20'] = mgr > det.SCREENING_MORAN_GROWTH_MIN
    # Confirmation (assumes screening pass + null pass)
    gates['conf_null_p_<0.01'] = p is not None and p < 0.01
    gates['conf_wall_spearman_<-0.40'] = (
        wsp is not None and wsp < det.CONFIRMATION_WALL_SPEARMAN_MAX
    )
    gates['conf_wall_final_<0.30'] = wfm < det.CONFIRMATION_WALL_FINAL_MAX
    gates['conf_wall_decay_>0.15'] = wdc > det.CONFIRMATION_WALL_DECAY_MIN
    # Definitive (assumes confirmation pass)
    gates['def_moran_final_in_[0.30,0.75]'] = (
        det.DEFINITIVE_MORAN_FINAL_MIN <= mfq <= det.DEFINITIVE_MORAN_FINAL_MAX
    )
    gates['def_wall_final_>0.05'] = wfm > det.DEFINITIVE_WALL_FINAL_MIN
    gates['def_minority_>0.05'] = (
        mff is not None and mff > det.DEFINITIVE_MINORITY_FINAL_MIN
    )

    return gates


def fmt(v) -> str:
    if v is None:
        return 'None'
    if isinstance(v, float):
        return f'{v:.4f}'
    return str(v)


def main() -> None:
    print("=" * 80)
    print("Schelling × P18 content-level characterization (Sprint 21)")
    print("L=64, density=0.9, threshold=0.375, n_steps=200, n_perms=199")
    print("=" * 80)

    records = []
    for seed in [0, 1, 2, 3, 4]:
        print(f"\n--- Seed {seed} ---")
        rec = characterize_seed(seed, L=64, n_steps=200)
        records.append(rec)
        print(f"  detected:                    {rec['detected']}")
        print(f"  tier:                        {rec['tier']}")
        print(f"  null_p_value:                {fmt(rec['null_p_value'])}")
        print("  -- detector primary --")
        print(f"  moran_spearman_early:        {fmt(rec['moran_spearman_early'])}")
        print(f"  moran_initial:               {fmt(rec['moran_initial'])}")
        print(f"  moran_final:                 {fmt(rec['moran_final'])}")
        print(f"  moran_final_qtr_mean:        {fmt(rec['moran_final_qtr_mean'])}")
        print(f"  moran_growth:                {fmt(rec['moran_growth'])}")
        print("  -- detector secondary (empty if screening failed) --")
        print(f"  wall_initial:                {fmt(rec['wall_initial'])}")
        print(f"  wall_final:                  {fmt(rec['wall_final'])}")
        print(f"  wall_final_qtr_mean:         {fmt(rec['wall_final_qtr_mean'])}")
        print(f"  wall_decay:                  {fmt(rec['wall_decay'])}")
        print(f"  wall_spearman_early:         {fmt(rec['wall_spearman_early'])}")
        print(f"  minority_fraction_final:     {fmt(rec['minority_fraction_final'])}")
        print("  -- direct trajectory diagnostic --")
        print(f"  direct_wall_initial:         {fmt(rec['direct_wall_initial'])}")
        print(f"  direct_wall_final:           {fmt(rec['direct_wall_final'])}")
        print(f"  direct_wall_final_qtr_mean:  {fmt(rec['direct_wall_final_qtr_mean'])}")
        print(f"  direct_wall_decay:           {fmt(rec['direct_wall_decay'])}")
        print(f"  direct_moran_initial:        {fmt(rec['direct_moran_initial'])}")
        print(f"  direct_moran_final:          {fmt(rec['direct_moran_final'])}")
        print(f"  direct_moran_final_qtr_mean: {fmt(rec['direct_moran_final_qtr_mean'])}")

        gates = evaluate_gates(rec)
        print("  -- gate-by-gate (using detector's reported metrics) --")
        for g, ok in gates.items():
            mark = "PASS" if ok else "FAIL"
            print(f"    {mark}  {g}")

    # Summary
    print("\n" + "=" * 80)
    print("Summary across 5 seeds")
    print("=" * 80)
    tiers = [r['tier'] for r in records]
    detecteds = [r['detected'] for r in records]
    print(f"  tiers:        {tiers}")
    print(f"  detecteds:    {detecteds}")
    print(f"  unique tiers: {set(tiers)}")
    print(f"  any P18 detected: {any(detecteds)}")

    print("\n  First failing gate per seed (in tier order):")
    for r in records:
        gates = evaluate_gates(r)
        first_fail = next((g for g, ok in gates.items() if not ok), None)
        print(f"    seed {r['seed']}: detected={r['detected']}, tier={r['tier']}, "
              f"first_fail={first_fail}")

    # Save full records
    out = Path(__file__).resolve().parent / 'schelling_p18_characterization.json'
    out.write_text(json.dumps(records, indent=2, default=str))
    print(f"\n  Full records → {out}")


if __name__ == '__main__':
    main()

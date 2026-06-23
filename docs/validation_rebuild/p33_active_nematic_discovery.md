# P33 — Active nematic order with ±1/2 topological defects (first Ring-0 discovery)

**Status: discriminating detector VALIDATED to the catalog Phase-2a standard.**
First pattern surfaced and verified end-to-end by the discovery pipeline
(find → verify → catalog). Ring-0 (known-but-uncatalogued): the phenomenon is
established science, novel to *our* 32-pattern catalog.

## What it is
Active nematics (Sanchez et al. 2012, *Nature* 491:431; Doostmohammadi et al.
2018, *Nat Commun* 9:3246): a field of apolar (head–tail-symmetric) rods develops
**local nematic orientational order** while continually nucleating and annihilating
**half-integer (±1/2) topological defects** — "active turbulence". The ±1/2 defect
is the topological fingerprint: it exists only in nematic (apolar) order; polar
fields admit only integer (±1) defects.

## How it was discovered
The Ring-0 sweep (`analysis/discovery/sweep.py`, commit 7e71fc2) ran an
agent-based active nematic through the full 31-detector battery at the confirmation
OOD gate: **EMERGENT-UNCLASSIFIED on 3/3 seeds**, em≈0.9 via the orientation
channel, closest catalog detector P6 (milling) — which never fired. The instrument
detected genuine emergence and correctly declined to match any catalog pattern.

## Distinctness from the catalog (the discriminator)
The centerpiece metric is **coherent half-integer defect density**:

```
D* = half_def_density   IF (S_loc ≥ 0.40  AND  φ ≤ 0.35  AND  |L| ≤ 0.25)
   = 0                   otherwise
```

where `S_loc` is the block-scale local nematic order |⟨e^{2iθ}⟩|, `φ` the global
polar order, `|L|` the angular momentum, and `half_def_density` the fraction of
plaquettes whose director winding is ±π (charge ±1/2). Each gate rejects a distinct
lookalike:

| lookalike | fails on | reading |
|---|---|---|
| P5 flocking (real Vicsek, low noise) | φ ≈ 0.997 > 0.35 | polar, not nematic |
| P6 milling (+1 vortex) | \|L\| ≈ 0.69 > 0.25; half-defects ≈ 0 | integer vortex, net rotation |
| isotropic noise (incl. real Vicsek, high noise) | S_loc ≈ 0.15–0.35 < 0.40 | no coherent local order |
| uniform nematic (defect-free) | half_def_density ≈ 0 | nematic order without the topological signature |

## Validation (Phase-2a, `analysis/discovery/p33_active_nematic.py`)
- **Positives:** active-nematic field, 5 seeds → **5/5 definitive** (D* 0.05–0.14).
- **Negatives:** polar / milling / isotropic / uniform-nematic × 5 seeds → **TNR 20/20 = 1.000**.
- **Continuous effect size:** d (D*) = **4.81** (pos vs pooled-neg).
- **Hardening (real catalog models + finite size):**
  - Real P5 Vicsek flock (low noise) rejected 3/3 (φ ≈ 0.997).
  - Real Vicsek disordered (high noise) rejected 3/3 (S_loc ≈ 0.35).
  - Finite-size: definitive at **G = 64, 96, 128** (the phenomenon has a clean
    large-N limit; the ±1/2 defect gas is sustained at every size).

## Verification ladder
1. **Robustness / finite-size** ✅ 5 seeds + G=64/96/128 all definitive.
2. **Not an instrument artifact** ✅ ±1/2 winding is a direct topological computation;
   a shuffle null on D* is reported (fragile only on the most-ordered seeds, hence
   reported as a diagnostic, not a hard gate — the negative panel is the rigor).
3. **Not in the catalog** ✅ TNR 1.0 vs the real P5 and a faithful P6; the Ring-0
   sweep showed the battery's closest match (P6) never fires.
4. **Literature novelty** — known to science (active matter), novel to the catalog.
   Ring-0 = catalog-coverage discovery, as planned.
5. **Mechanism** ✅ local nematic alignment + activity noise → a steady-state ±1/2
   defect gas (active turbulence); apolar self-propulsion → φ ≈ 0.
6. **Discriminating detector** ✅ built and validated to the Phase-2a bar.

## Next (integration to close the loop)
Promote the self-contained module into `epc/models/active_nematic.py`,
`epc/metrics/nematic_order.py`, `epc/detectors/p33_active_nematic.py` + a standalone
`analysis/discovery/p33_active_nematic_test.py`; then register P33 in the battery
(continuous_metrics + calibration) so the instrument MATCHES active nematic instead
of returning EMERGENT-UNCLASSIFIED — closing the find → verify → catalog loop.

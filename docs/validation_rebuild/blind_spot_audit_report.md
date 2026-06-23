# Blind-spot / detection-recall audit of the emergence instrument

**Date:** 2026-06-23. **Question:** can the catalog instrument detect all the
emergent behavior it might encounter — and which additional channels most improve
its recall? Method: ~13 unambiguously-emergent/null probe systems spanning the six
emergence families (including deliberately-suspected blind spots) run through three
detection channels; recall measured against ground truth, thresholds for the new
channels calibrated against the null probes. (3-seed robustness; `analysis/blind_spot_probes.py`,
`analysis/blind_spot_audit.py`, `analysis/emergence_channels.py`.)

## Channels tested
- **CURRENT** — the instrument's `generic_emergence` indicator (4 morphology-specific
  channels: clustering, spatial autocorrelation, phase/orientation order, distributional entropy).
- **Ψ_CE** — Rosas/Mediano synergy / causal-emergence criterion
  `Ψ_CE = I(Vₜ;Vₜ′) − Σᵢ I(Xᵢ,ₜ;Vₜ′) > 0` (best over a generic macro-feature bank).
- **MPR-C** — Martín-Plastino-Rosso statistical complexity on Bandt-Pompe ordinal
  patterns (structure orthogonal to entropy).

Both new estimators were validated on synthetic signals before use (noise→Ψ_CE≈0;
redundant→Ψ_CE≪0; XOR-parity→Ψ_CE≈+1; noise/sine/chaos→correct complexity-entropy
behavior). That validation caught a real estimator bug — quantile-binning silently
collapsed binary data and zeroed Ψ_CE — which would have broken every discrete probe.

## Result

| | recall on emergent probes | null false-positives |
|---|---|---|
| CURRENT (4 channels) | ~6/11 robust (0.64 single-seed) | 0.00 |
| + Ψ_CE + MPR-C | ~9–10/11 | Ψ_CE 0.00; **MPR-C non-zero (see caveat)** |

**Families the current instrument MISSES:** connectivity, temporal-oscillation,
transient-wave, chaos.

**Recovered by the new channels:**
- **Ψ_CE → connectivity (percolation, 2/3) and synergy (XOR, 3/3).** Clean: exactly
  0 on every non-synergistic system; zero null false-positives. The synergy class
  ("greater-than-the-parts", XOR-type) is invisible to all correlation/clustering/
  order-parameter channels by construction — Ψ_CE is the only thing that sees it.
- **MPR-C → temporal-oscillation (limit cycle, 3/3) and chaos** — with two conditions
  the audit surfaced (below).

**Still blind to ALL three channels: the transient wave** (an SIR-style front that
crosses and leaves a uniform final state). Its signature is in the full trajectory,
not the final frame; needs a dedicated moving-front / full-run detector.

## Honest caveats the audit exposed about the NEW channels
1. **MPR-C macro choice is decisive.** Chaos was missed on the *mean-field* summary
   (averaging 64 chaotic sites → near-constant, washes the chaos out) but clearly
   caught on a *single-site* series (C = 0.33 vs 0.14). MPR-C must be applied to
   per-component / multiple candidate macros (max), not a global average.
2. **MPR-C false-fires on smooth drift.** A bounded random walk's running mean is a
   Brownian path — temporally autocorrelated, so it scores as "structured" (null FP
   1–2/3). MPR-C therefore needs the same **shuffle/phase-randomized surrogate-null**
   discipline the rest of the indicator uses, not a fixed threshold.
3. **Ψ_CE needs the right candidate macro** (it is sufficient-not-necessary, so it can
   false-negative under redundancy); a generic bank (mean/std/parity/majority) worked
   here but is not guaranteed complete.
4. **Some "catches" are via a secondary channel.** The vortex mill is flagged via its
   *clustering* (it is also a cohesive disk), lanes via *nematic* (antiparallel =
   nematic axis) — the system is correctly flagged emergent, but not because rotation
   or banding per se is detected. A discovery user would get the right verdict, the
   wrong reason.

## Bottom line / recommendation
- **The instrument cannot detect all emergent behavior** (current recall ~0.64), and
  no finite channel set makes it complete — emergence is open-ended.
- **Add Ψ_CE now** (highest-confidence: closes synergy + connectivity, zero FP).
- **Add MPR-C with a surrogate-null wrapper and per-component macros** (closes
  temporal-oscillation + chaos).
- **Build a dedicated transient/front (full-trajectory) detector** for the remaining
  blind class.
- For the DISCOVERY phase: treat any "emergent-but-unclassified" hit as a candidate
  *only after* it clears a surrogate null on the firing channel, and report the
  detection channel so the blind classes above are explicitly flagged as low-recall.

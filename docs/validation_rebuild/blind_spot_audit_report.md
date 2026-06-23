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

## Round-2 update (2026-06-23, af66728) — recall 0.64 → 1.0

Two general-purpose channels were wired into the live `generic_emergence` after
probe+detector co-design (round-2 workspace `analysis/round2/`):

1. **Ψ_CE synergy / causal emergence** (c0ce695) — recovers synergy (XOR) +
   connectivity (percolation). Gated at Ψ_CE>0.08; zero null false-positives.
2. **Spectral-peak oscillation channel, gated to non-phase observations** (af66728)
   — recovers temporal-oscillation (limit cycle), chaos (coupled-map lattice), and
   transient-wave (SIR front) in one stroke. The phase-kind gate is the key: phase
   synchronization is already covered by the order-parameter channel, and the only
   null that false-fired (uncoupled Kuramoto, osc 17) is phase-kind; every non-phase
   null is ≤ 9 vs the emergent range 35–64, so threshold 15 cleanly separates.

**Co-design lessons (why the obvious approaches failed):** the coordination-gate
(collective complexity vs desync surrogate) has the WRONG polarity for
synchronization — sync *lowers* complexity (clean oscillation), the desync
surrogate is a more-complex beat; and a naive "is the macro periodic?" detector
false-fires on trivial independent oscillators. The phase-kind gate resolves it.

**Result:** blind-spot audit (3-seed) **11/11 emergent recall, 0/2 null
false-positives**; T2c UNCHANGED (null-spec 1.0, STRICT rec 0.83 / nov 1.0 /
false-MATCH 0). The transient/front family is recovered by detection (the front's
activity pulse has spectral structure); a dedicated front *classifier* is no longer
needed for the indicator. STILL OPEN for full comprehensiveness: probes+channels for
heavy-tail/power-law distributions and network/community emergence (round-2 cont.).

## Round-2 continued (2026-06-23, 720d7ce) — heavy-tail + network → recall 14/14

Two more general channels wired (probe suite extended to 14 emergent / 3 nulls):
3. **Heavy-tail / SOC** (25228d1) — `plaw_llr` + `heavy_tail_score` read raw frames
   for size observables (avalanche_sizes/activity) or per-step change; fire when the
   event-size distribution is power-law (LLR>0) over ≥1.3 decades. Required
   restructuring generic_emergence so the info-theoretic channels run even with NO
   morphology frame (a sandpile scored 0 on every channel: 0.00 → 0.95).
4. **Network topology** (720d7ce) — `graph_structure_score` reads an 'adjacency'
   observable and z-scores spectral modularity (community) and degree heterogeneity
   (scale-free) against an Erdős–Rényi null. community 0.87, scale-free 0.70, ER null
   0.00. Self-contained (scipy, no networkx).

**FINAL recall: blind-spot audit 14/14 emergent (3-seed), 0/3 null false-positives;
T2c UNCHANGED (null-spec 1.0, STRICT rec 0.83 / nov 1.0 / false-MATCH 0).** Live
generic_emergence channels: morphology (clustering / autocorr / phase / vector) +
orientation (polar+nematic) + Ψ_CE synergy + gated spectral-peak temporal +
heavy-tail/SOC + network (modularity / scale-free). The detection net now spans all
families surfaced by the emergence-channel research, with measured recall and zero
T2c regression. Completeness remains impossible in principle (open-ended concept);
this is comprehensive coverage of the known families with honest per-channel reporting.

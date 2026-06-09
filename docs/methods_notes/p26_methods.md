# P26 Methods Note — Stochastic Resonance

**Pattern:** P26 — Stochastic resonance
**Canonical model:** Bistable double-well with subthreshold periodic signal + noise
**Detector:** `epc/detectors/p26_stochastic_resonance.py`
**Model:** `epc/models/stochastic_resonance.py`
**Reproduction sprint:** Sprint 74 (`analysis/reproductions/p26_collins.py`)

---

## 1. Pattern and canonical reference

P26 is the phenomenon where noise improves system performance at an
intermediate optimal level — too little noise yields no signal detection,
too much noise destroys the signal, but a sweet spot amplifies it. The
defining signature is the **inverted-U performance-vs-noise curve**:
performance is low at zero noise, peaks at an intermediate noise level,
and declines at high noise.

**Distinctness from P24 (homeostasis):** P24 is about maintaining a
variable near a set-point despite perturbation. P26 is about noise
*improving* signal detection — the noise is the mechanism, not the
disturbance. P24 has active feedback; P26 has no feedback, just
constructive noise-signal interaction.

**Primary reference:** Gammaitoni, L., Hänggi, P., Jung, P., & Marchesoni, F.
(1998). Stochastic resonance. *Rev. Mod. Phys.* 70(1), 223–287.

**Secondary references:**
- Collins, J. J., Chow, C. C., & Imhoff, T. T. (1995). Stochastic resonance
  without tuning. *Nature*, 376(6537), 236–238.
- Benzi, R., Sutera, A., & Vulpiani, A. (1981). The mechanism of stochastic
  resonance. *J. Phys. A*, 14, L453–L457.

---

## 2. Model implementation

### 2.1 BistableDoubleWell (canonical)

Overdamped particle in a symmetric double-well potential:

```
V(x) = -a·x²/2 + b·x⁴/4
dx/dt = a·x - b·x³ + A·sin(2πft) + √(2D)·ξ(t)
```

where D is noise intensity and ξ(t) is unit white noise. Wells at
x = ±√(a/b), barrier height = a²/(4b).

At zero noise, x stays in one well with tiny intra-well oscillation
(amplitude ~ A/(2a)). At optimal noise, the Kramers escape rate
matches the signal frequency, producing inter-well hopping synchronized
with the signal. At high noise, hopping is random and unsynchronized.

### 2.2 ThresholdUnit (T1b cross-model)

```
output = 1 if (signal(t) + noise(t)) > threshold else 0
```

Signal is subthreshold (amplitude < threshold). At intermediate noise,
the signal+noise sum crosses the threshold preferentially when the
signal is near its peak, producing correlated detections. This is an
independent SR implementation with analytically tractable behavior.

### 2.3 Multi-trial design

Both models run multiple independent trials (n_trials, default 20) per
noise level with alternating starting wells. The detector computes
the coherent response averaged across all trials, yielding a clean
performance estimate even for single realizations.

### 2.4 Parameters

| Parameter | Canonical value | Role |
|-----------|----------------|------|
| a         | 4.0            | Well depth coefficient |
| b         | 1.0            | Quartic confinement |
| barrier   | 4.0            | Barrier height a²/(4b) |
| A         | 1.0            | Subthreshold signal amplitude |
| f         | 0.005          | Signal frequency |
| dt        | 0.01           | Integration timestep |
| n_steps   | 20000          | Steps per trial (1 signal period) |
| n_trials  | 20             | Independent trials per noise level |
| n_levels  | 15             | Noise levels in sweep |

---

## 3. Detector design

### 3.1 Observation contract (T1a)

The detector reads a **noise-sweep observation bundle** via
`extract_observation_bundle()`. The bundle contains five aligned arrays:
`time`, `x` (system output), `signal` (driving signal), `noise_level`,
and `noise_level_idx`. Any system producing these fields can be tested
by the P26 detector without modification.

### 3.2 Performance metric

The coherent response: |⟨x · signal⟩| — the absolute value of the
time-averaged product of x and the driving signal. This metric is
sensitive to both the **amplitude** and **phase** of x's response:

- Zero noise: x stays in one well → tiny oscillation → small ⟨x·s⟩
- Optimal noise: inter-well hopping synchronized with signal → large ⟨x·s⟩
- High noise: random hopping → ⟨x·s⟩ averages to ~0

Unlike Pearson correlation, the coherent response is sensitive to
response amplitude, not just phase. At zero noise, x tracks the signal
phase perfectly but with negligible amplitude → small metric. This is
exactly the SR signature.

### 3.3 Null model

Time-shuffle at the peak noise level: randomly permute the temporal
order of x while preserving its marginal distribution. Under H0
(no synchronization), the shuffled ⟨x·signal⟩ ≈ 0. The observed
peak should far exceed this baseline.

This directly tests whether the coherent response requires genuine
temporal synchronization between x and the driving signal.

### 3.4 Tier criteria

- **Screening:** Performance at some nonzero noise > performance at
  zero noise (gain > 0.02).
- **Confirmation:** Inverted-U shape — interior peak with both rise
  (gain > 0.02) and fall (decline > 0.02), null p < 0.05.
- **Definitive:** Gain > 0.05, decline > 0.02, Cohen's d > 1.0,
  p ≤ 0.005, metadata confirms subthreshold signal.

### 3.5 Key design decisions

1. **Coherent response over FFT-based SNR:** The FFT SNR metric is
   degenerate at zero noise (infinite SNR, zero noise floor) and
   requires many signal periods for convergence. The coherent response
   gives a clean, finite value at all noise levels and works for both
   continuous (double-well) and binary (threshold unit) outputs.

2. **Multi-trial averaging:** A single realization of the double-well
   at any noise level has high variance in the coherent response. Running
   20 independent trials per level and averaging reduces this variance
   by √20 ≈ 4.5×, producing a clean inverted-U curve.

3. **Time-shuffle null over phase-shuffle:** Circular shifting preserves
   the periodicity of x, leaving ⟨x·s⟩ nonzero. Time-shuffling fully
   destroys temporal structure, giving a proper null baseline near zero.

4. **Metadata gate at definitive:** If metadata indicates
   `has_subthreshold_signal=False`, definitive tier is denied. This
   prevents false positives from suprathreshold systems where signal
   detection is already maximal without noise.

---

## 4. Reproduction results

### 4.1 dim1 (Sprint 74)

Bistable double-well (a=4, b=1, barrier=4, A=1.0, f=0.005):
- Peak noise: D = 1.5
- Peak coherent response: 0.918
- Zero-noise coherent response: 0.063
- **Gain over zero: 0.855** (tolerance: > 0.05) — **PASS**
- **Decline after peak: 0.811** (tolerance: > 0.02) — **PASS**
- **Interior peak:** idx=7 of 15 — **PASS**
- Detector tier: **DEFINITIVE** (confidence=0.90) — **PASS**

### 4.2 dim2 (Sprint 74)

20-seed campaign (n_trials=20, n_steps=20000 per seed):
- Peak coherent response: 0.897 ± 0.049 (CV=5.4%)
- Gain over zero: 0.833 ± 0.049 (CV=5.8%)
- Peak noise location: D = 1.27 ± 0.25
- All 20 seeds: **DEFINITIVE** (fraction=1.00)

See `analysis/outputs/p26_multiseed.json`.

### 4.3 T1b cross-model (Sprint 74)

Threshold unit (threshold=1.0, amplitude=0.7) tested with same detector:
- DEFINITIVE detection (p ≤ 0.005)
- Confirms detector recognizes the *phenomenon* (stochastic resonance),
  not the specific implementation (bistable double-well vs threshold unit).

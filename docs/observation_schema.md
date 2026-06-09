# Observation Schema — T1a Canonical Observation Contracts

**Purpose:** Document the stable observation-bundle schemas that detectors
consume. Each bundle defines a set of aligned arrays that any system — native
model, independent ABM, or real data — must provide for a detector to operate.
Detectors read bundles via adapter functions, NOT model objects directly.

**Instrument roadmap:** This file is the Track 1a deliverable. See
`docs/instrument_roadmap.md` for context.

---

## Scalar-regulated-variable bundle (P24)

**Detector:** `epc/detectors/p24_homeostasis.py`
**Adapter:** `extract_observation_bundle(history) → dict`

| Key           | Type              | Description |
|---------------|-------------------|-------------|
| `time`        | `ndarray(T,)`     | Timestamps (monotonically increasing) |
| `x`           | `ndarray(T,)`     | Regulated variable value at each time |
| `setpoint`    | `ndarray(T,)`     | Target set-point value at each time |
| `perturbation`| `ndarray(T,)`     | External perturbation magnitude at each time |

**Input format:** list of dicts, each containing keys `'time'`, `'x'`,
`'setpoint'`, `'perturbation'` (all numeric scalars). The adapter extracts
and stacks these into aligned numpy arrays.

**Usage example:**
```python
from epc.detectors.p24_homeostasis import P24HomeostasisDetector

# Any system producing dicts with these four keys works:
history = [{'time': t, 'x': x_val, 'setpoint': sp, 'perturbation': p}
           for t, x_val, sp, p in zip(times, xs, setpoints, perts)]

detector = P24HomeostasisDetector(n_permutations=199, seed=42)
result = detector.detect(history)
```

**Native models:** `ProportionalHomeostat`, `IntegralHomeostat`
(`epc/models/homeostasis.py`) — both produce history dicts matching this schema.

---

## Noise-sweep-timeseries bundle (P26)

**Detector:** `epc/detectors/p26_stochastic_resonance.py`
**Adapter:** `extract_observation_bundle(history) → dict`

| Key              | Type              | Description |
|------------------|-------------------|-------------|
| `time`           | `ndarray(T,)`     | Timestamps (monotonically increasing within each noise level) |
| `x`              | `ndarray(T,)`     | System output at each timestep |
| `signal`         | `ndarray(T,)`     | Driving signal value at each timestep |
| `noise_level`    | `ndarray(T,)`     | Noise intensity D for this timestep |
| `noise_level_idx`| `ndarray(T,)`     | Integer index of the noise level (for grouping) |

**Input format:** list of dicts, each containing keys `'time'`, `'x'`,
`'signal'`, `'noise_level'`, `'noise_level_idx'` (all numeric scalars).
The adapter extracts and stacks these into aligned numpy arrays.

**Usage example:**
```python
from epc.detectors.p26_stochastic_resonance import P26StochasticResonanceDetector

# Any system producing dicts with these five keys works:
history = [{'time': t, 'x': x_val, 'signal': s, 'noise_level': D, 'noise_level_idx': idx}
           for t, x_val, s, D, idx in zip(times, xs, signals, noise_levels, indices)]

detector = P26StochasticResonanceDetector(n_permutations=199, seed=42)
result = detector.detect(history)
```

**Native models:** `BistableDoubleWell`, `ThresholdUnit`
(`epc/models/stochastic_resonance.py`) — both produce history dicts matching
this schema. The bundle concatenates all noise levels and trials into a single
flat sequence; the detector groups by `noise_level_idx` to compute the
per-level coherent response.

---

## Attendance/choice time series bundle (P23)

**Detector:** `epc/detectors/p23_anticoordination.py`
**Adapter:** `extract_observation_bundle(history) → dict`

| Key           | Type              | Description |
|---------------|-------------------|-------------|
| `round`       | `ndarray(T,)`     | Round number (or time index) |
| `attendance`  | `ndarray(T,)`     | Number of agents choosing side 1 (or attending) at each round |
| `n_agents`    | `int`             | Total number of agents |
| `capacity`    | `float`           | Capacity threshold (N/2 for symmetric MG, stated capacity for El Farol) |

**Input format:** list of dicts, each containing keys `'attendance'` and
`'n_agents'` (both integer). Optional: `'round'` (int), `'capacity'` (int).
The adapter extracts and stacks these into aligned numpy arrays.

**Usage example:**
```python
from epc.detectors.p23_anticoordination import P23AnticoordinationDetector

# Any system producing dicts with attendance and n_agents works:
history = [{'round': t, 'attendance': a, 'n_agents': N}
           for t, a in enumerate(attendance_series)]

detector = P23AnticoordinationDetector(n_permutations=199, seed=42)
result = detector.detect(history)
```

**Native models:** `MinorityGame`, `ElFarolBar`
(`epc/models/minority_game.py`) — both produce history dicts matching this
schema. The Minority Game reports binary-choice attendance (side 1 count);
El Farol reports bar attendance count.

---

## Attractor-network state-vector bundle (P16)

**Detector:** `epc/detectors/p16_associative_memory.py`
**Adapter:** `extract_observation_bundle(history) → dict`

| Key               | Type              | Description |
|-------------------|-------------------|-------------|
| `state`           | `ndarray(N,)`     | Binary neuron states (±1) at each timestep |
| `step`            | `int`             | Update sweep number within the current trial |
| `trial`           | `int`             | Trial index (which cue presentation) |
| `cue_pattern_idx` | `int`             | Index of the stored pattern used as cue |
| `overlap`         | `float`           | Overlap m = (1/N) Σ ξ_i s_i with the cued pattern |
| `stored_patterns` | `ndarray(P, N)`   | The P stored template patterns |
| `converged`       | `bool`            | Whether the network has reached a fixed point |

**Input format:** list of dicts, each containing the keys above. The adapter
groups records by `trial` and extracts per-trial state trajectories, overlaps,
and convergence status. The stored patterns array is read from the first record.

**Usage example:**
```python
from epc.detectors.p16_associative_memory import P16AssociativeMemoryDetector

# Any system producing dicts with these keys works:
history = [{'state': s, 'step': t, 'trial': tr, 'cue_pattern_idx': idx,
            'overlap': m, 'stored_patterns': patterns, 'converged': conv}
           for ...]

detector = P16AssociativeMemoryDetector(n_permutations=199, seed=42)
result = detector.detect(history)
```

**Native models:** `HopfieldNetwork`, `BooleanGRN`
(`epc/models/hopfield.py`) — both produce history dicts matching this schema.

---

## Future bundles

As new detectors are added with T1a contracts, their observation bundles
will be documented here. Planned:
- Agent-position-velocity bundle (continuous_2d detectors: P5, P6, P7, P17, P19)
- Grid-state bundle (lattice_2d detectors: P1, P13, P15, P22, etc.)
- Oscillator-phase bundle (P9, P10)

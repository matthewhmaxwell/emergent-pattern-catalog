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

## Future bundles

As new detectors are added with T1a contracts, their observation bundles
will be documented here. Planned:
- Agent-position-velocity bundle (continuous_2d detectors: P5, P6, P7, P17, P19)
- Grid-state bundle (lattice_2d detectors: P1, P13, P15, P22, etc.)
- Oscillator-phase bundle (P9, P10)

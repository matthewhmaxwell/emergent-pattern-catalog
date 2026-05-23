# Sprint 37 Return Summary

**Sprint goal:** Add `_gen_p2_active_brownian` and `_gen_p6_dorsogna` native catalog
generators + minimum-viable B' supplement builders for the `continuous_2d` substrate
type, so Sprint 43's continuous_2d panel batch can use a complete substrate-typed
Class B. **Status: complete.**

**Sprint 38 cue:** Lattice_1d B' supplements + `_gen_p8_nagel_schreckenberg` generator
addition (next infrastructure prerequisite for Sprint 44 mixed batch).

---

## Pre-flight verification

- Base HEAD: `d95c107` (Sprint 36 follow-up — orchestrator parser fix; matches
  expected Sprint 36 post-commit). ✓
- `PYTHONPATH=. python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 /
  27 / 19 / 361 / 77 / 284** ✓ (unchanged, as expected for infrastructure-only sprint).
- `git pull --ff-only origin main`: already up to date. ✓
- pytest pre-flight: chain-aware skip per brief (Sprint 36 post-flight passed at HEAD).

---

## Part A — `_gen_p2_active_brownian` added to `epc/phase2a/catalog.py`

**`SUBSTRATE_PARAMS["P2_abp"]`** added with canonical Fily-Marchetti 2012 MIPS
parameters:
- `n_particles=200, box_size=16.0` → packing fraction φ ≈ 0.61 (phase-separated regime)
- `v0=0.3, D_r=3e-3` → Péclet Pe = 100 ≥ 50 (MIPS onset)
- `rho_star=10.0, r_cg=1.0, dt=0.1, n_steps=200, seed=0`

Panel-scale (200 particles) rather than full canonical (400 particles) for tractable
runtime; qualitative MIPS pattern still present.

**`_gen_p2_active_brownian(p)`** instantiates `ActiveBrownianParticles`, calls
`model.run(n_steps=p["n_steps"])`, stacks `headings` and `positions` into float32
arrays, returns `{"kind": "particles", ...}` matching `_gen_p5_vicsek`'s shape
contract.

**Registered** in `_GENERATORS["P2_abp"]`.

---

## Part B — `_gen_p6_dorsogna` added to `epc/phase2a/catalog.py`

**`SUBSTRATE_PARAMS["P6_dorsogna"]`** added with Carrillo et al. 2009 canonical
milling parameters confirmed in `tests/test_abp_p2_e2e.py`:
- `n_particles=100, C_a=0.5, C_r=1.0, l_a=3.0, l_r=0.5, alpha=1.0, beta=0.5`
- `dt=0.05, init_mode="ring", init_radius=5.0, n_steps=200, seed=0`

`init_mode="ring"` chosen so milling is established from step 0, avoiding the need
for a long burn-in at reduced dt. Total simulated time = 10.0 ≫ orbital period
(T_orbit ≈ 2π × R / v_eq ≈ 2π × 3.0 / 1.414 ≈ 13.3), so one partial orbit is
captured in 200 frames.

**`_gen_p6_dorsogna(p)`** derives `box_size` dynamically from
`2 × (max(|positions|) + 2)` since D'Orsogna uses open (non-periodic) space. Returns
`{"kind": "particles", ...}` matching shape contract.

**Registered** in `_GENERATORS["P6_dorsogna"]`.

---

## Part C — `PATTERN_TO_SUBSTRATE_ID` comments updated

Entries for `"P2"` and `"P6"` were already present (declarative stubs from a prior
sprint). Comments updated from `"declarative; generator NOT yet implemented"` to
`"generator added Sprint 37"`. No structural change to the map — entries and their
substrate-type lookups through the registry were correct prior to this sprint.

`class_b_for_pattern("P5")` now returns `["P2_abp", "P6_dorsogna"]` as catalog_mates
(all three patterns share `continuous_2d` substrate type per the MODEL_REGISTRY). ✓

---

## Part D — 2 B' supplement builders for `continuous_2d` in `epc/phase2a/structured.py`

### `uncorrelated_random_walks`

N particles (default 200) in a periodic box. Headings redrawn i.i.d. uniform(−π, π)
each step; positions updated with constant speed `speed=0.03`. No neighbor coupling.
Polarization order parameter r ≈ 1/√N by construction. A correctly-specific flocking
or milling detector should not fire.

### `independent_brownian_motion`

N particles (default 200) in a periodic box. Each particle takes an independent
Gaussian displacement (σ=0.1) each step; no interactions. Headings derived from
displacement direction. Tests against trivial diffusive motion that could spuriously
trigger collective-motion detectors.

Both builders:
- Added to `SUPPLEMENTS_BY_SUBSTRATE_TYPE["continuous_2d"]`
- Added to `SUPPLEMENT_BUILDERS` dict
- Return `list[dict]` with `positions`, `headings`, `velocities`, `speeds`, `step`
  keys — same format as Vicsek/ABP/DOrsogna model histories

---

## Part E — 6 new tests in `tests/test_phase2a_panel.py`

| Test | What it verifies |
|---|---|
| `test_gen_p2_active_brownian_deterministic` | Same seed → byte-identical ABP trajectory |
| `test_gen_p6_dorsogna_deterministic` | Same seed → byte-identical D'Orsogna trajectory |
| `test_gen_p2_active_brownian_output_format` | kind="particles", ndim correct, box_size > 0 |
| `test_gen_p6_dorsogna_output_format` | kind="particles", ndim correct, box_size > 0 |
| `test_class_b_p5_contains_p2_abp_and_p6_dorsogna` | substrate_type="continuous_2d", both mates present |
| `test_continuous_2d_supplements_registered` | both builders in SUPPLEMENTS_BY_SUBSTRATE_TYPE and SUPPLEMENT_BUILDERS |

All 6 pass. Total fast test count: **591** (585 baseline + 6 new). ✓

---

## Part F — `docs/paper_CHANGELOG.md` updated

Sprint 37 entry prepended with the 3 required bullet points. Mechanical-only sprint;
no §3–§8 prose changes.

---

## Part G — Carry-forward review

Per brief: C-supplements partially closed (continuous_2d done this sprint).

| Substrate type | C-supplements status |
|---|---|
| `continuous_2d` | **CLOSED** (this sprint) |
| `lattice_1d` | OPEN (Sprint 38) |
| `lattice_2d_continuous` | OPEN |
| `scalar_wealth` | OPEN |
| `opinion_space` | OPEN |

All other open carry-forwards unchanged from Sprint 36:
- **C-class-a-constant-field-trivial-sync (Sprint 35):** Open. P9 residual false
  positive from constant_field substrate. v1.3 candidate; not in scope.
- **C-p27-time-shuffle-invariance (Sprint 34):** Open/provisional. Flag=True per spec.
- **C-p14-class-c-borderline (Sprint 33):** Persists at p_diss=0.350; P14 still PASS.
- **C-pyproject-pin (Sprint 29):** Open. 1-line pyproject.toml change deferred.

---

## Post-flight verification

All acceptance criteria met:

| Check | Result |
|---|---|
| `'P2_abp' in _GENERATORS and 'P6_dorsogna' in _GENERATORS` | **PASS** ✓ |
| `class_b_for_pattern('P5')` contains both mates | **['P2_abp', 'P6_dorsogna']** ✓ |
| `'uncorrelated_random_walks' in SUPPLEMENTS_BY_SUBSTRATE_TYPE['continuous_2d']` | **PASS** ✓ |
| `pytest tests/test_phase2a_panel.py` | **83 passed** (+6 from Sprint 36 baseline of 77) ✓ |
| `pytest tests/ -m "not slow"` | **591 passed, 0 failed, 65 deselected** ✓ |
| `scripts/count_transfer_matrix.py` | **20/19/79/274/27/19/361/77/284** (unchanged) ✓ |

`git diff --stat` (vs Sprint 36 HEAD): 4 files changed, 179 insertions(+), 2
deletions(−): `epc/phase2a/catalog.py`, `epc/phase2a/structured.py`,
`tests/test_phase2a_panel.py`, `docs/paper_CHANGELOG.md`. No other files modified. ✓

---

## Deviations and judgment calls

### Deviation 1 — D'Orsogna box_size derived dynamically

The brief specified `kind="particles"` with `box_size`. D'Orsogna operates in open
(non-periodic) space with no intrinsic box_size parameter. Solution: compute
`box_size = 2 × (max(|positions|) + 2.0)` from the actual position range across the
simulated trajectory. This gives a well-defined bounding box for the occupancy-grid
adapter in `_adapt_to_grid`, while remaining seed-deterministic (positions are
deterministic given fixed seed + params).

### Deviation 2 — `init_mode="ring"` for D'Orsogna

The model default is `init_mode="random"`. The brief says "canonical milling regime —
check the existing replication tests." The replication test uses `init_mode="random"
+ dt=0.01`, but at `n_steps=200` with `dt=0.01` the orbit isn't complete (T_orbit ≈
13.3, but simulated time = 2.0). Using `init_mode="ring"` + `dt=0.05` gives
simulated time = 10.0 and milling is immediate, which is more representative of the
canonical P6 pattern. The parameter set (C_a, C_r, l_a, l_r, α, β) matches the
published Carrillo 2009 values used in the replication tests exactly.

---

## HEAD commit hash and tag at end of sprint

- **Commit:** `b9caec0`
- **Tag:** `v0.37.0` (pushed to origin)

**Decision: GO** — infrastructure-only sprint completed cleanly; all 591 fast tests pass; transfer matrix unchanged; chain may proceed autonomously to Sprint 38.

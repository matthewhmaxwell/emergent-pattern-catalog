# P7 Methods Note — Lane Formation in Counterflow

**Pattern:** P7 — Lane formation in counterflow
**Canonical model:** Counterflow social-force (Helbing & Molnár 1995)
**Detector:** `epc/detectors/p7_lane_formation.py`
**Model:** `epc/models/lane_formation.py`
**Reproduction sprint:** Sprint 65 (`analysis/reproductions/p7_helbing1995.py`)

---

## 1. Pattern and canonical reference

P7 is the spontaneous segregation of two oppositely-moving populations into
same-direction lanes within a shared corridor. The phenomenon arises from
short-range lateral repulsion between agents: head-on encounters push agents
sideways, and once local directional purity is established it self-reinforces
(same-direction neighbors don't push laterally). The result is a dynamic flow
structure (not a built path like P29 trails) that differs from P5 flocking
(single heading consensus) by maintaining TWO persistent opposing streams.

**Primary reference:** Helbing, D. & Molnár, P. (1995). Social force model for
pedestrian dynamics. *Phys. Rev. E* 51(5), 4282–4286.

**Secondary references:**
- Helbing, D., Farkas, I. J., & Vicsek, T. (2000). Simulating dynamical
  features of escape panic. *Nature* 407, 487–490.
- Nowak, S. & Schadschneider, A. (2012). Quantitative analysis of pedestrian
  counterflow in a cellular automaton model. *Phys. Rev. E* 85(6), 066128.

---

## 2. Model equations and update rule

The implementation uses a first-order (overdamped) social-force model:

```
v_i(t+dt) = v_desired_i + (F_social_i + F_wall_i) * dt / tau
x_i(t+dt) = x_i(t) + v_i(t+dt) * dt
```

**Desired velocity:** `v_desired_i = v0 * (+1, 0)` for population 0 (rightward),
`v_desired_i = v0 * (-1, 0)` for population 1 (leftward).

**Social force** (pairwise repulsion):
```
F_ij = A * exp(-d_ij / B) * n_ij
F_social_i = sum_{j != i, d_ij < r_cut} F_ij
```
where `d_ij` is the inter-agent distance (periodic in x), `n_ij` is the unit
vector from j to i, A is repulsion amplitude, B is decay length, and r_cut is
the interaction cutoff.

**Wall force** (y-boundaries):
```
F_wall_bottom = A_w * exp(-y_i / B_w) * (0, +1)
F_wall_top = A_w * exp(-(H - y_i) / B_w) * (0, -1)
```

**Boundary conditions:** Periodic in x (flow axis), reflecting walls at y=0
and y=H (lateral axis).

---

## 3. Parameters

| Parameter | Symbol | Default | Role |
|-----------|--------|---------|------|
| n_agents | N | 200 | Total agents (split evenly between populations) |
| corridor_width | W | 20.0 | Flow-axis length (periodic) |
| corridor_height | H | 4.0 | Lateral width (walls) |
| desired_speed | v0 | 1.0 | Target speed magnitude |
| repulsion_amplitude | A | 5.0 | Social force strength |
| repulsion_range | B | 0.3 | Social force decay length |
| tau | τ | 0.5 | Velocity relaxation time |
| dt | Δt | 0.05 | Integration timestep |
| interaction_radius | r_cut | 2.0 | Force computation cutoff |
| wall_repulsion_amplitude | A_w | 10.0 | Wall force strength |
| wall_repulsion_range | B_w | 0.2 | Wall force decay length |

The density `ρ = N/(W*H)` is the primary control parameter. At the canonical
regime ρ = 2.5 agents/m², lanes form reliably within ~50 time units.

---

## 4. Observable extraction

**Lane order parameter** (primary): Following Nowak & Schadschneider (2012),
the corridor is divided into `n_bins` lateral strips (perpendicular to flow).
For each strip containing agents:
```
phi_strip = |n_right - n_left| / (n_right + n_left)
```
The lane order parameter is the mean over occupied strips:
```
phi_lane = mean(phi_strip)
```
Random mixing → φ ≈ 0 (for large N/n_bins); perfect lanes → φ = 1.

**Head-on encounter rate:** Count pairs from opposing populations within
encounter_radius = 0.5, normalized per agent.

**Throughput:** Mean speed in desired direction:
```
throughput = mean(v_x_i * sign_i)
```
where sign_i = +1 for population 0, −1 for population 1.

---

## 5. Deviations from canonical

1. **Overdamped (first-order) dynamics:** Helbing-Molnár 1995 uses second-order
   (inertial) dynamics with mass. Our implementation omits inertia (`m → 0`
   limit) for simplicity. This accelerates lane formation but does not qualitatively
   change the steady-state lane structure.

2. **Circular (isotropic) social force:** The original model uses elliptical
   (anisotropic) force fields accounting for pedestrians' visual awareness.
   We use isotropic repulsion; lanes still form via the same lateral-exclusion
   mechanism.

3. **No body radius / hard-core exclusion:** Agents are point particles with
   soft exponential repulsion. This simplification means agents can occasionally
   overlap at very close range rather than experiencing hard-core contact forces.

4. **Reflecting vs. periodic lateral walls:** Helbing 1995 uses a corridor
   with explicit walls; our implementation matches this with reflecting boundaries.

---

## 6. Reproduction cross-links and limitations

**Reproduction:** `analysis/reproductions/p7_helbing1995.py` →
`analysis/outputs/p7_helbing1995_reproduction.json`
- 5 seeds at canonical regime (ρ=2.5)
- φ rises from ~0.50 (mixed) to ~0.92 (steady-state lanes)
- Throughput fraction of v0: ~0.998 in steady state

**Multi-seed:** `analysis/reproductions/p7_multiseed.py` →
`analysis/outputs/p7_multiseed.json`
- 20 seeds: φ = 0.897 ± 0.091, CV = 10.2%

**Limitations:**
- The overdamped approximation means formation timescales differ from the
  published inertial model (lanes form faster without inertia).
- At very high density (ρ > 6), the soft-core repulsion may permit unphysical
  overlaps that an excluded-volume model would prevent.
- The lane order metric has finite-size bias: at low N/n_bins, random fluctuations
  inflate φ. The canonical regime (N=200, n_bins=10) gives ~20 agents/bin,
  sufficient for statistical reliability.
- dim4 (Phase-2a panel) pending Sprint 66.

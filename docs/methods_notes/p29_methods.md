# P29 Methods Note — Trail / Network Formation

## 1. Pattern definition
P29 detects emergent efficient transport networks: agents collectively build
connective pathways linking source/food nodes that approach the efficiency
of the minimum spanning tree (MST) while retaining fault tolerance through
redundant connections. The emergent network balances total path length
(efficiency) against robustness (multiple routes between nodes).

**Canonical reference:** Tero, A. et al. (2010). Rules for biologically
inspired adaptive network design. *Science*, 327(5964), 439-442.

**Distinctness from catalog neighbours:**
- P4 (territoriality): P4 makes EXCLUSIVE domains via scent avoidance;
  P29 makes CONNECTIVE transport networks via trail reinforcement.
- P7 (lane formation): P7 lanes are dynamic flow structures that vanish
  when agents stop; P29 trails are persistent built paths reinforced by
  environmental modification (pheromone/conductance).
- P1 (aggregation): P1 is spatial clustering; P29 is spatial connectivity
  (linking distributed nodes, not grouping agents).

## 2. Model implementations

### 2.1 AntTrailModel (canonical)
Ant colony optimization on a complete graph of N food/source nodes.
Each edge carries a pheromone level. At each step:

1. **Evaporation**: all edge pheromones decay by factor (1 - evaporation_rate).
2. **Agent movement**: each ant at node i chooses next node j with probability
   proportional to pheromone(i,j)^alpha / distance(i,j)^beta. The alpha
   exponent controls pheromone-following; beta controls distance preference.
3. **Deposition**: traversed edge (i,j) receives pheromone += deposition_rate / dist(i,j).

Over many iterations, frequently-used short edges accumulate pheromone
(positive feedback) while infrequent edges decay → a sparse, efficient
network emerges near the MST.

Default parameters: n_nodes=7, n_agents=40, alpha=1.0, beta=2.0,
deposition_rate=10.0, evaporation_rate=0.02, n_steps=500.

### 2.2 PhysarumModel (T1b cross-model)
Tero et al. (2010) flux-reinforcement on a complete graph. Each edge
has a conductance D_ij evolving by:

    dD_ij/dt = |Q_ij|^gamma - decay_rate * D_ij

where Q_ij is the electrical-analogue flux (solved via Kirchhoff's
current law for each source-sink pair). gamma > 1 gives positive
feedback: high-flow edges strengthen, low-flow edges decay.

The model uses a qualitatively different mechanism (global flow
optimization vs local pheromone following) to achieve the same
emergent outcome: an efficient transport network.

Default parameters: n_nodes=7, gamma=1.8, decay_rate=0.01, n_steps=2000.

### 2.3 NoReinforcementModel (negative control)
Agents on the same complete graph choose next nodes uniformly at random
(no pheromone bias, no distance bias). Uniform deposition (constant per
traversal) produces edge weights uncorrelated with distance → no efficient
network structure emerges. Correctly rejected by the P29 detector.

## 3. Detector design

### 3.1 T1a observation contract
The detector reads a **trail-network observation bundle** via
`extract_observation_bundle()`. The bundle contains:

| Key | Type | Description |
|-----|------|-------------|
| `node_positions` | list[ndarray(N,2)] | Node coordinates per snapshot |
| `edge_weights` | list[ndarray(N,N)] | Edge weight matrix per snapshot |
| `pheromone_fields` | list[ndarray(G,G)] | 2D visualization field |
| `steps` | ndarray(T,) | Step numbers |
| `n_nodes` | int | Number of food/source nodes |
| `grid_size` | int | Domain size |

### 3.2 Primary metric
**Weight-distance correlation (Spearman)**: rank correlation between edge
weight and 1/distance across all edges. For reinforced networks, short edges
accumulate more weight → strong positive correlation (>0.5). For random
traffic, no systematic correlation (~0).

This metric captures the fundamental signature of trail formation: the
emergent preference for efficient (short) routes over random paths.

### 3.3 Null model
**Edge-weight shuffle**: permute edge weights across the complete graph's
edges, preserving the marginal weight distribution but destroying the
weight-distance correlation. The p-value is the fraction of null shuffles
with correlation ≥ observed.

### 3.4 Tier criteria

| Tier | Criteria |
|------|----------|
| Screening | corr > 0.1, connectivity ≥ 0.6, p < 0.10 |
| Confirmation | corr > 0.3, length/MST < 2.0, p < 0.05, d > 1.0 |
| Definitive | corr > 0.5, length/MST < 2.0, ft > 0, p ≤ 0.005, metadata |

### 3.5 Design decisions
1. **Spearman over Pearson**: rank correlation is robust to the heavy-tailed
   pheromone distributions typical of ACO convergence.
2. **Median threshold for network extraction**: edges above median positive
   weight define the "strong" network for length/MST and connectivity.
3. **Definitive requires metadata**: `has_pheromone_reinforcement=True`
   confirms the mechanism is stigmergic, not coincidental distance bias.

## 4. Reproduction results

### dim1 (Sprint 88)
Physarum (grid layout, seed=42):
- length/MST = 1.3536 (tolerance [1.0, 1.5]: PASS)
- fault_tolerance = 1.000 (> 0: PASS)
- weight_dist_corr = 0.846 (> 0.5: PASS)
- Detector tier: **DEFINITIVE** (confidence=0.90, p=0.005, d=2.60)

### dim2 (Sprint 88)
20-seed campaign (Physarum, random layout):
- length/MST: 1.548 ± 0.112 (CV=7.2%)
- weight-dist corr: 0.647 ± 0.096
- Detected: 19/20 (95%), Definitive: 1/20
See `analysis/outputs/p29_multiseed.json`.

### T1b cross-model (Sprint 88)
AntTrailModel (grid layout): CONFIRMATION+ detection (corr=0.818, p=0.005).
PhysarumModel (grid layout): DEFINITIVE detection (corr=0.846, p=0.005).
Confirms detector recognizes the *phenomenon* (efficient transport network),
not the specific implementation (ACO vs flux-reinforcement).

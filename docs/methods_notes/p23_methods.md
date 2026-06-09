# P23 — Anti-coordination / Emergent Load Balancing: Methods Note

## Pattern summary

Agents repeatedly self-distribute across options to avoid overcrowding,
producing collective utilization near capacity without central assignment.
The group neither converges to one shared choice nor fragments into stable
camps — it dynamically balances.

## Canonical model

**Minority Game** (Challet & Zhang, 1997). N agents (odd N), each holding
S strategies mapping the last m winning sides (a binary history of length
m) to a binary choice {0, 1}. The minority side wins each round; agents
accumulate virtual points for strategies that would have predicted the
winning side, and play their highest-scoring strategy. The control
parameter is α = 2^m / N.

**El Farol Bar** (Arthur, 1994). N agents decide whether to attend a bar
with capacity C. Each agent holds K linear predictors of next attendance
from recent history. Agents attend if their best predictor forecasts
attendance below capacity. Predictor scores decay with prediction error.

## Detection approach

**Primary metrics:**
- Scaled variance σ²/N of attendance after burn-in, compared to the
  random-choice baseline p̂(1−p̂) where p̂ = mean(attendance)/N.
- Lag-1 autocorrelation of attendance (negative = anti-persistence).

**Null model:** Random-choice surrogate. For each permutation, generate
i.i.d. Binomial(N, p̂) attendance series of the same length. Compare
observed σ²/N and lag-1 autocorrelation against this null distribution.
This tests whether reduced variance and temporal anti-persistence require
strategic adaptation or could arise from independent random choice.

**Tier criteria:**
- Screening: Mean attendance within 20% of capacity.
- Confirmation: σ²/N < random baseline (p < 0.01 vs surrogate) OR
  negative lag-1 autocorrelation (p < 0.01 vs surrogate).
- Definitive: BOTH conditions satisfied simultaneously.

**Burn-in:** 20% of rounds discarded as transient (agents' initial
strategy scores are zero; the system needs time to develop adaptive
behavior).

## Key design decisions

1. **Surrogate null over temporal shuffle.** A temporal shuffle preserves
   the marginal distribution (and hence the variance) of the attendance
   series, so it cannot test whether variance reduction is real. The
   random-choice surrogate tests the right null hypothesis: whether the
   observed statistics could arise without strategic adaptation.

2. **Separate p-values for variance and autocorrelation.** Some regimes
   (e.g., the efficient phase near α ≈ 0.3) show strong variance
   reduction but weak or absent negative autocorrelation. Others may
   show the reverse. Using OR for confirmation (either metric alone
   suffices) and AND for definitive (both required) captures this.

3. **p̂-adaptive baseline.** Using p̂(1−p̂) rather than the fixed 0.25
   baseline handles asymmetric games (El Farol with capacity ≠ N/2)
   correctly.

## Savit curve (dim1 anchor)

The canonical Minority Game result is the Savit curve: σ²/N as a
function of α = 2^m / N. Our reproduction (N=101, m=1..11, 10 seeds
per point, 3000 rounds, burn-in 600) shows:

- Interior minimum at α ≈ 0.32 (m=5) with σ²/N ≈ 0.077.
- Random baseline σ²/N = 0.25 crossed between m=8 and m=9 (α ≈ 2.5–5).
- Symmetric (maladaptive) phase at small α: σ²/N ≈ 1.45 at m=1.

This matches the qualitative shape reported by Savit, Manuca & Riolo
(1999, Fig. 1).

## Multi-seed stability (dim2)

25 seeds at α ≈ 0.63 (m=6, N=101): σ²/N mean = 0.075 ± 0.006
(CV = 0.087). All 25 seeds detected (23 confirmation, 2 definitive).

## Limitations

- The lag-1 autocorrelation is not consistently negative across all
  efficient-phase seeds; some seeds show weak positive ac1. This limits
  the rate of definitive-tier detection at this α. The variance reduction
  path (confirmation tier) is robust across all seeds.
- The El Farol model uses linear predictors, which are simpler than the
  original Arthur (1994) formulation with diverse predictor families.
  This is sufficient for T1b cross-model validation.

## References

- Arthur, W. B. (1994). Inductive reasoning and bounded rationality.
  AER, 84(2), 406–411.
- Challet, D. & Zhang, Y.-C. (1997). Emergence of cooperation and
  organization in an evolutionary game. Physica A, 246, 407–418.
- Savit, R., Manuca, R. & Riolo, R. (1999). Adaptive competition,
  market efficiency, and phase transitions. PRL, 82(10), 2203–2206.

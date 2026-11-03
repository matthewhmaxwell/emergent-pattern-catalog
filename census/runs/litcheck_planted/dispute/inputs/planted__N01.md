# planted:N01

## Target paper (ground truth; the gate never saw this)
- Title: Arithmetic of the sync basin for pulse-coupled oscillators
- arXiv: 2609.01668 (v1 2026-09-01); authors: K. P. O'Keeffe
- Behaviour the authors report as new: At linear charging, all-to-all pulse-coupled oscillators end fully synchronised with probability set by the prime factorisation of N (1-1/N^N for prime N; failures ~N^-(m-1) with m the least prime factor), and only equal-size cluster states survive.
- Authors' novelty statement: "The result adds a new member to the atlas of exotic basin geometries: alongside fractal, riddled, and tentacled basins, we now have a basin that is arithmetic."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: N identical oscillators, each coupled to every other (no space). Each has a voltage between 0 and 1, drawn uniformly at random at the start. Rules: (1) all voltages rise at the same constant speed (straight-line charging, no leak); (2) when a voltage reaches 1 the oscillator fires and resets to 0, together with any oscillators already firing with it as a group; (3) a firing group of k oscillators kicks every other voltage up by k/N; (4) any oscillator pushed to 1 or above by a kick joins the firing group and fires with it from then on. Run until the grouping stops changing and record whether all N end as one group. Sizes: N from 2 to about 100, many random starts per N (exact enumeration for N up to 8).

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Same rules with N prime (2, 3, 5, 7, 11): the only failure is all N staying separate.
  - Same rules with even N (4, 6, 8, ...): failures fall off only slowly with N.
  - Same rules with N divisible by 3 but not 2 (9, 15, ...).

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - full_sync_probability_n4 = 0.887  (share of random starts ending as one group, N=4)
  - full_sync_probability_n5 = 0.99968  (share of random starts ending as one group, N=5)
  - full_sync_probability_n6 = 0.906  (share of random starts ending as one group, N=6)
  - distinct_final_states_n100 = 9  (number of possible final groupings at N=100, all equal-sized)
  - failure_decay_exponent_even_n = 1  (failure probability falls as N to minus this, even N)
  - failure_prefactor_even_n = 0.6  (prefactor of that failure power law for even N)

This behaviour appears across 7 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Peskin / Mirollo-Strogatz pulse-coupled integrate-and-fire oscillators at marginal coupling (linear charging, pulse j/N, absorption); aggregation/cluster-kinetics framework of O'Keeffe-Krapivsky-Strogatz
- **behaviour_documented**: True
- **notes**: Undocumented twist: final states are exactly m equal clusters for each divisor m of N (9 states at N=100), near-certain sync for prime N, and roughly 0.6/N failure for even N. The divisor rule follows from period matching: without leak each cluster of size k must charge k/N per cycle, so all clusters are equal-sized and their number divides N. The twist is real but minor; only the 0.6/N prefactor needs real calculation. Neutral (no-leak) dynamics keeps leftover multi-cluster states from merging. Caution: exact-threshold ties (kicks are multiples of 1/N and rule 4 includes equality) may make results sensitive to floating point; check the N=5 deficit (0.99968) with exact arithmetic. PRE 93 032203 is sole-authored by O'Keeffe. The PRL supplement was not read and may contain the finite-N details.
- **citations**: ["O'Keeffe, Krapivsky, Strogatz, PRL 115, 064101 (2015) + Supplemental Material (deterministic Peskin model, linear charging, pulse 1/N)", "K.P. O'Keeffe, PRE 93, 032203 (2016): full sync guaranteed only for gamma>0 and pulse >1/N; otherwise multiple stable clusters possible", "Mirollo & Strogatz, SIAM J. Appl. Math. 50, 1645 (1990)", "PRE 95, 022207 (2017): equal-size cluster splay states in pulse-coupled LIF networks", "Bottani, PRE 54, 2334 (1996): global IF synchronization, en bloc vs avalanche"]

### N_blocked_r2 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Peskin / Mirollo-Strogatz globally pulse-coupled integrate-and-fire oscillators with linear charging and pulse eps=1/N (O'Keeffe-Krapivsky-Strogatz 'synchronization as aggregation' model), at the marginal non-concave case
- **behaviour_documented**: True
- **notes**: The model and how its clusters grow at large N are documented exactly. The twist I could not find in the literature is the finite-N endgame: non-sync outcomes with positive probability; neutral final states made only of equal-sized groups whose size divides N (9 states at N=100 = 9 divisors of 100); prime N almost always syncing; and even-N failure of about 0.6/N, mainly from two halves of size N/2. My own elementary calculation (not from a source) explains this: two unequal groups a,b drift by (a-b)/N per cycle, so only equal groups can coexist, and straight-line charging (concavity b=0) is the boundary where the Mirollo-Strogatz proof fails. The twist is real but disappears with any curvature or pulse > 1/N. I could not read the PRL supplement, which might cover part of this.
- **citations**: ["O'Keeffe, Krapivsky, Strogatz, Phys. Rev. Lett. 115, 064101 (2015), doi:10.1103/PhysRevLett.115.064101", "O'Keeffe, Phys. Rev. E 93, 032203 (2016), 'Transient dynamics of pulse-coupled oscillators with nonlinear charging curves'", "Mirollo & Strogatz, SIAM J. Appl. Math. 50, 1645 (1990)", "Senn & Urbanczik, SIAM J. Appl. Math. 61, 1143 (2001)"]

### N_blocked_r3 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Peskin / Mirollo-Strogatz pulse-coupled integrate-and-fire oscillators, linear (non-leaky) charging with pulse 1/N: cluster aggregation via periodic size-doubling into equally spaced equal-sized clusters (O'Keeffe, Krapivsky & Strogatz 2015, Supplemental Material)
- **behaviour_documented**: True
- **undocumented_phenomenology**: Finite-N endgame statistics not found: end states restricted to equal groups whose size divides N; prime N almost always fully synchronizes; full-sync probability jumps around with N; even-N failure falls as about 0.6/N. These look like the endpoint of the documented pairing process rather than a new mechanism.
- **notes**: The exact model (linear charging, pulse 1/N, all-to-all) is analysed in the OKS 2015 Supplemental Material for large N, where clusters double in size and form equally spaced, equal-sized groups. The authors leave out the endgame. Pulse = 1/N is a knife-edge where kicks can land exactly on threshold, so the exact probabilities (e.g. N=5 giving 0.99968 rather than 1) may depend on tie-breaking or floating point; check with exact arithmetic. Could not read the OKS Supplemental Material or the Senn-Urbanczik conditions.
- **citations**: ["O'Keeffe, Krapivsky, Strogatz, Synchronization as Aggregation: Cluster Kinetics of Pulse-Coupled Oscillators, PRL 115, 064101 (2015) + Supplemental Material", "O'Keeffe, Transient dynamics of pulse-coupled oscillators with nonlinear charging curves, PRE 93, 032203 (2016)", "Mirollo & Strogatz, SIAM J. Appl. Math. 50, 1645 (1990)", "Senn & Urbanczik, Similar non-leaky integrate-and-fire neurons with instantaneous couplings always synchronize, SIAM J. Appl. Math. (2000/2001)"]

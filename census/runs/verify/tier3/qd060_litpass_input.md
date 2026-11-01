RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 3.14159 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent advances the type of the cell under it by one (cyclically); (3) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_activity_trend, g_corrlen, g_reversals, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 4.4659 with interaction, 6.3483 with interaction knocked out (difference -1.8824, robust spread 0.2163); difference at 4x area -1.9168, at 3x run length -1.8939
  - g_activity_trend (activity rising): median -2.666 with interaction, 0.0 with interaction knocked out (difference -2.666, robust spread 0.0847); difference at 4x area -2.9408, at 3x run length -2.7764
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.864 with interaction, 0.0569 with interaction knocked out (difference 0.8071, robust spread 0.0184); difference at 4x area 0.8517, at 3x run length 0.8192

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 3 rules):
  - rules 1 only: g_n_domains 6.3483; g_activity_trend 0.0; g_corrlen 1.0; g_reversals 0.0; g_largest_domain 0.0569
  - rules 2 only: g_n_domains 6.3953; g_activity_trend 0.0001; g_corrlen 1.0; g_reversals 1.0; g_largest_domain 0.0386
  - rules 3 only: g_n_domains 4.7664; g_activity_trend -0.5022; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.8235
  - rules 1, 2 only: g_n_domains 6.3578; g_activity_trend -3.5294; g_corrlen 1.0; g_reversals 1.0; g_largest_domain 0.0505
  - rules 1, 3 only: g_n_domains 4.7664; g_activity_trend -0.5022; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.8235
  - rules 2, 3 only: g_n_domains 4.8402; g_activity_trend -0.0038; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.3419
No sub-program came within 3 robust spreads of the full program on every feature above.

CONTROL RESULTS (census settings, same 10 seeds; median of each feature):
- Distinct cells occupied by the 400 agents at the last frame, per seed: [10, 20, 6, 34, 6, 34, 11, 13, 6, 4] (median 10).
- g_n_domains: full program 4.466; null (spatial structure shuffled every step) 6.322; uniform field (gradient exactly zero) 4.84; 10 agents pinned in place (speed 0) 4.783
- g_activity_trend: full program -2.666; null (spatial structure shuffled every step) 0.005687; uniform field (gradient exactly zero) -0.003804; 10 agents pinned in place (speed 0) -0.5245
- g_corrlen: full program 2; null (spatial structure shuffled every step) 1; uniform field (gradient exactly zero) 2; 10 agents pinned in place (speed 0) 2
- g_reversals: full program 1; null (spatial structure shuffled every step) 1; uniform field (gradient exactly zero) 1; 10 agents pinned in place (speed 0) 1
- g_largest_domain: full program 0.864; null (spatial structure shuffled every step) 0.7196; uniform field (gradient exactly zero) 0.3419; 10 agents pinned in place (speed 0) 0.8169
Notes from the code: the field starts at 1.0 with 20 random 3x3 patches at 0.5 plus 1% noise and has no source or decay, so it only smooths out; an agent's heading is set to the local gradient direction when the gradient is non-zero; 'activity' is the share of cells that changed between recorded frames two steps apart; g_activity_trend is log(late activity / early activity); g_reversals counts flip-backs and is saturated at 1.0 whenever cells change with two types.

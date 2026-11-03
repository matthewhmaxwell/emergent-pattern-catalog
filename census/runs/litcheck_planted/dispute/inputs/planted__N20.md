# planted:N20

## Target paper (ground truth; the gate never saw this)
- Title: Nonreciprocity Propels Nematic Defects
- arXiv: 2610.39037 (v1 2026-10-28); authors: Yunpeng Cao, Yu Du, Yao Li
- Behaviour the authors report as new: Vision-cone-type nonreciprocal coupling on a nematic lattice makes +1/2 defects self-propel along their polar axis without active stress, and makes +1 defects split into three +1/2 and one -1/2 defect at strong nonreciprocity.
- Authors' novelty statement: "Nonreciprocal interactions are a generic source of nonequilibrium dynamics, yet their consequences for nematic topological defects remain unexplored."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: a triangular lattice (6 neighbours per site) of linear size 100 with open edges; each site holds a director angle defined modulo 180 degrees. Rules, deterministic (zero noise): each site's angle changes at a rate equal to the sum over its 6 neighbours of w*sin(2*(neighbour angle - own angle)), where w = exp(s*cos(2*phi)) and phi is the angle between the bond to that neighbour and the site's own director, s from 0 to 0.5. A site thus listens more to neighbours lying along its own axis, so influence is not mutual. Start from imposed defect textures (single half-integer defects, or one integer +1 defect with texture phase mu); locate defects by angle winding around triangles.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Continuum version with the corresponding extra torque term: isolated +1/2 defect speed rises with s.

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - plus_half_speed_trend = 1  (sign of change in +1/2 defect self-propulsion speed as s rises)
  - minus_half_defect_drift = 0  (intrinsic translation speed of isolated -1/2 defects)
  - four_defect_split_share_s05 = 0.64  (share of initial phases where +1 splits into four defects, s=0.5)
  - four_defect_split_share_s010 = 0.0  (same share at s=0.10)
  - four_defect_onset_s = 0.15  (asymmetry strength s above which four-defect splitting appears)

This behaviour appears across 8 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Nonreciprocal anisotropic-kernel XY model (Rouzaire-Pearce-Pagonabarraga-Levis 2025), nematic (2-theta, fore-aft-symmetric kernel) variant; flow-free +1/2 defect self-propulsion as in dry/active nematics
- **behaviour_documented**: True
- **notes**: Core behaviour (+1/2 drift that grows with s, no -1/2 drift, shape/phase-dependent defect dynamics) is documented in the polar anisotropic-kernel NRXY and active-nematic literature. The kernel exp(s cos 2phi) is fore-aft symmetric, so this is an axial kernel, not a vision cone; +1/2 motion comes from the comet shape plus the non-mutual dynamics. Twist not found: a threshold (s~0.15) splitting of +1 into 3(+1/2)+1(-1/2), i.e. nonreciprocity creating an extra +-1/2 pair at the core, with a phase-dependent share (0.64 at s=0.5). Moderately significant if robust. Rule out boundary nucleation (open edges), lattice-scale core effects, and finite-size effects (N=100); the Rouzaire et al. full text and supplement were paywalled and not checked.
- **citations**: ["Rouzaire, Pearce, Pagonabarraga, Levis, PRL 134, 167101 (2025), arXiv:2401.12637", "Loos, Klapp, Martynec, PRL 130, 198301 (2023)", "Dopierala, Chate, Shi, Solon, PRL 135, 088302 (2025)", "Brauns et al., PRL 136, 058302 (2026)", "Shankar/Bowick et al.-type active nematic defect self-propulsion; -1/2 non-motile by threefold symmetry (R. Soc. Open Sci. 10, 221229 (2023))", "New J. Phys. abe8a8 (2021): higher-charge defects dissociate into elementary ones"]

### N_blocked_r2 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Non-reciprocal (anisotropic-kernel / vision-cone) XY model, nematic analogue: shape-dependent defect dynamics with self-propelled +1/2 and passive -1/2 defects
- **behaviour_documented**: True
- **notes**: Self-propelled +1/2 and non-drifting -1/2 defects, plus phase-dependent dynamics, are documented (polar NRXY with the same exp-kernel; active nematics generally). -1/2 zero drift is symmetry-forced and uninformative. Twist not found: phase-dependent splitting of a +1 defect into four defects (likely 3x(+1/2)+1x(-1/2)), with onset near s~0.15 and 64% of phases at s=0.5, in a nematic (cos2phi) non-reciprocal kernel model; no published nematic version of this model located. The twist is plausibly real (an anisotropy-driven pair nucleation at the core), but it may be a triangular-lattice or core-discretization or open-boundary artifact; check with a square lattice, larger L, small noise, and defect charge bookkeeping.
- **citations**: ["Rouzaire, Pearce, Pagonabarraga, Levis, PRL 134, 167101 (2025), Nonreciprocal Interactions Reshape Topological Defect Annihilation", "Inescapable Anisotropy of Nonreciprocal XY Models, PRL 135, 088302 (2025)", "Popli, Maitra, Ramaswamy, Ordering and Defect Cloaking in Nonreciprocal Lattice XY Models, PRL 135 (2025)", "Loos, Klapp, Martynec, PRL 130, 198301 (2023)", "The role of fluid flow in the dynamics of active nematic defects, New J. Phys. (2021)", "Active Solids: Topological Defect Self-Propulsion Without Flow, PRL 136, 058302 (2026)"]

### N_blocked_r3 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Nonreciprocal (vision-cone/anisotropic-kernel) XY model, nematic version; motile +1/2 and passive -1/2 defects as in dry active nematics
- **behaviour_documented**: True
- **notes**: +1/2 speed rising with s and zero -1/2 drift follow from symmetry and match documented active-nematic behaviour; texture-dependent defect dynamics under an anisotropic non-mutual kernel is documented by Rouzaire et al. for the polar case. The twist not found: above s~0.15 a +1 defect splits into four defects (three +1/2 plus one -1/2, i.e. a ±1/2 pair is created) for a phase-dependent share of initial textures (0.64 at s=0.5), instead of the standard split into two +1/2. Plausibly meaningful (non-mutual coupling creating defect pairs) but needs checks for edge seeding (open boundaries, L=100), brief pairs that annihilate, and lattice/detection artifacts. Only abstracts and snippets were searched; a nematic follow-up to Rouzaire et al. may already contain this.
- **citations**: ["Rouzaire, Pearce, Pagonabarraga, Levis, PRL 134, 167101 (2025), arXiv:2401.12637", "Loos, Klapp, Martynec, PRL 130, 198301 (2023)", "Inescapable Anisotropy of Nonreciprocal XY Models, PRL (2025), doi:10.1103/r3dx-7lrd", "Popli, Maitra, Ramaswamy, PRL 135, 088303 (2025)", "Brauns et al., PRL 136, 058302 (2026)", "Shankar et al., PRL 121, 108002 (2018)", "Santhosh et al., R. Soc. Open Sci. 10, 221229 (2023), doi:10.1098/rsos.221229", "Thijssen et al., Phys. Rev. Research 2, 042008 (2020)"]

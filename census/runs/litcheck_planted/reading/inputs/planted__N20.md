# planted:N20

**Reported behavior (what the authors report as new):** Vision-cone-type nonreciprocal coupling on a nematic lattice makes +1/2 defects self-propel along their polar axis without active stress, and makes +1 defects split into three +1/2 and one -1/2 defect at strong nonreciprocity.

**Rule on the card:** World: a triangular lattice (6 neighbours per site) of linear size 100 with open edges; each site holds a director angle defined modulo 180 degrees. Rules, deterministic (zero noise): each site's angle changes at a rate equal to the sum over its 6 neighbours of w*sin(2*(neighbour angle - own angle)), where w = exp(s*cos(2*phi)) and phi is the angle between the bond to that neighbour and the site's own director, s from 0 to 0.5. A site thus listens more to neighbours lying along its own axis, so influence is not mutual. Start from imposed defect textures (single half-integer defects, or one integer +1 defect with texture phase mu); locate defects by angle winding around triangles.

**Observed features on the card:**
  - plus_half_speed_trend = 1 (sign of change in +1/2 defect self-propulsion speed as s rises)
  - minus_half_defect_drift = 0 (intrinsic translation speed of isolated -1/2 defects)
  - four_defect_split_share_s05 = 0.64 (share of initial phases where +1 splits into four defects, s=0.5)
  - four_defect_split_share_s010 = 0.0 (same share at s=0.10)
  - four_defect_onset_s = 0.15 (asymmetry strength s above which four-defect splitting appears)

## Reply R02

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Nonreciprocal (vision-cone / anisotropic-kernel) XY model, nematic variant; nonreciprocity-induced defect self-propulsion", "notes": "No source found for the nematic (2-theta) nonreciprocal variant or for the +1 -> four-defect splitting threshold (s~0.15) and phase-dependent share. Those parts may be new, or may sit in papers I could not read. Most key papers were paywalled, so only abstracts were checked. -1/2 zero drift follows from threefold symmetry. Four defects from one +1 imply charges 3x(+1/2) + 1x(-1/2)."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The model family and the general mechanism are published. I found no source that describes this exact nematic (2θ) system or the four-defect splitting of a +1 defect.

**1. Model family.** This is a nonreciprocal XY-type lattice spin model, in the dry active matter / nonequilibrium lattice spin class. Each coupling depends on the bond direction measured in the receiving spin's own frame. Here the director is nematic (mod 180°), so the coupling is sin 2Δθ and the weight is exp(s cos 2φ). The published relatives use polar spins and a vision-cone or kernel weight.

**2. What the literature describes**
- **Vision-cone nonreciprocal XY (earliest clear source).** Loos, Klapp & Martynec (PRL 130, 198301, 2023) study a 2D XY model where each spin interacts only with neighbours in a "vision cone" around its current orientation. They report that defects propagate in a directional manner, breaking parity and time-reversal symmetry. I could only read the abstract and snippets.
- **Closest match to your rule.** Rouzaire, Pearce, Pagonabarraga & Levis (PRL 134, 167101, 2025) weight the XY coupling by an anisotropic kernel. They find that defect shape matters in addition to topological charge, and that nonreciprocal coupling twists the spin field and selects specific defect shapes. This resembles the role of your texture phase μ for +1 defects. I could only read the abstract.
- **Lattice analysis.** Popli, Maitra & Ramaswamy (PRL 135, 088303, 2025) treat square and triangular lattices. They find defect interactions are highly anisotropic with respect to the mean ordering direction. A companion paper is "Inescapable Anisotropy of Nonreciprocal XY Models" (PRL, 2025, DOI 10.1103/r3dx-7lrd). It shows nonreciprocity is akin to self-propulsion and that the late state is a dynamic foam of asters.
- **Review-level source.** "Hamiltonian description of non-reciprocal interactions" (Nature Physics, 2026) covers vision-cone XY models. I did not read it in detail.
- **Active nematic benchmarks.** In fluid active nematics, a +1/2 defect is self-propelled and a −1/2 defect is not self-propelled under uniform activity. This matches your finding of zero drift for −1/2 defects. By symmetry, a threefold-symmetric −1/2 core cannot have a preferred direction.
- **Self-propulsion without flow.** Brauns et al. (PRL 136, 058302, 2026) propose +1/2 defects that move by local remodeling of the nematic texture without advection. This mechanism can unbind defect pairs and stabilize +1 defects. It is a conceptual relative of your result, but it is an elastic active solid, not nonreciprocal alignment.

**3. What differs, and whether it matters**
- I found no paper on a nonreciprocal nematic (head-tail symmetric) lattice model with exp(s cos 2φ) weights. I also found none on the matching continuum extra-torque term where +1/2 defects speed up with s. A nematic version is the natural extension of the polar NRXY work. +1/2 self-propulsion that grows with s, and zero drift for −1/2, are therefore expected rather than surprising.
- I found no report that a +1 defect splits into four defects above a threshold in s (about 0.15 here), with a phase-dependent share (0.64 at s=0.5, 0 at s=0.1). I only searched for it briefly. The ordinary passive splitting of integer defects into half-integer pairs, for example in biaxial or thin-cell nematics, is a different mechanism. Four defects from one +1 must carry total charge +1, so three +1/2 and one −1/2. This is a quantitative detail I did not find in any source.
- The polar-model papers show that the nonreciprocal twist reshapes defects. The nematic four-way splitting is plausibly a new case of that, but I could not confirm it.

**4. Confidence: medium-low.** The core mechanism and the symmetry results are well supported. The specific threshold and splitting statistics are unconfirmed. Several key papers were paywalled, so I read only abstracts. I also ran only 8 searches, which is a limited sweep for the splitting claim.

**

## Reply R12

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Non-reciprocal (anisotropic-kernel / vision-cone) XY model, nematic analogue: shape-dependent defect dynamics with self-propelled +1/2 and passive -1/2 defects", "notes": "Self-propelled +1/2 and non-drifting -1/2 defects, plus phase-dependent dynamics, are documented (polar NRXY with the same exp-kernel; active nematics generally). -1/2 zero drift is symmetry-forced and uninformative. Twist not found: phase-dependent splitting of a +1 defect into four defects (likely 3x(+1/2)+1x(-1/2)), with onset near s~0.15 and 64% of phases at s=0.5, in a nematic (cos2phi) non-reciprocal kernel model; no published nematic version of this model located. The twist is plausibly real (an anisotropy-driven pair nucleation at the core), but it may be a triangular-lattice or core-discretization or open-boundary artifact; check with a square lattice, larger L, small noise, and defect charge bookkeeping."}

**Reply text:**

**Verdict: KNOWN-VARIANT (medium confidence).** I agree with the first pass's family, but I'd tighten its reasoning. Two of the five features are documented or expected, two are only sign or symmetry statements, and one is a real twist I couldn't find anywhere.

**What is documented**

- **The polar version of this exact model exists.** Rouzaire, Pearce, Pagonabarraga and Levis (PRL 2025) weight the 2D XY coupling by an anisotropic kernel. They find that in addition to the topological charge q, the actual shape of the defects becomes crucial to faithfully describe their dynamics, and that nonreciprocal coupling twists the spin field, selecting specific defect shapes. A dependence on the texture phase μ is therefore expected, not new.
- **Non-reciprocity acts like activity.** For lattice non-reciprocal XY models, nonreciprocity is akin to a self-propulsion. Related vision-cone models show defects that propagate in a directional manner, thereby breaking the parity and time-reversal symmetry.
- **Moving +1/2 and stationary −1/2 is the standard active-nematic picture.** +1/2 defects self-propel, while self-propulsion is absent for the −1/2 defects due to their threefold symmetry. A recent model even gets +1/2 self-propulsion without any flow, by local remodeling of the nematic texture without advection. So `minus_half_defect_drift = 0` is forced by symmetry and isn't evidence of anything. `plus_half_speed_trend = +1` (faster as s rises) is the expected direction.
- **+1 defects breaking up is generic.** Higher-charge defects are energetically unstable to disassociation into elementary defects. Splitting into two +1/2 is the default expectation.

**The twist I couldn't find**

1. A **nematic (apolar, cos 2φ) version** of the anisotropic-kernel non-reciprocal model. I found only polar non-reciprocal XY / vision-cone papers.
2. **A +1 defect splitting into four defects** (by charge conservation, presumably three +1/2 and one −1/2). This happens for 64% of initial phases at s = 0.5, never at s = 0.10, and starts around s ≈ 0.15. The nearest documented relatives are pair creation in active nematics and the "minimal four +1/2" states on spheres and shells. Neither is this mechanism.

**Does the twist matter?** Plausibly, yet it's the feature most at risk of being an artifact:
- **Lattice effects:** lattice non-reciprocal XY models are known to show strong lattice anisotropy (the lattice always induces anisotropy on large scales). The six-fold triangular lattice could favour particular split geometries.
- **Unresolved defect cores:** the defect cores are about one lattice spacing wide, so the extra ±1/2 pair could be a discretization effect.
- **Edges and noise:** the open edges at L = 100 and the zero-noise deterministic dynamics could bias which initial phases split.

The obvious checks are a square lattice, larger L, a weak noise floor, and confirming the charges of the four defects.

**Trivial or artifact?** Overall, no. Defect motility tied to the anisotropy strength is a real emergent effect, not a frozen, diverged or inert state. The four-defect onset is the only feature with a meaningful artifact risk.

**

## Reply R19

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Nonreciprocal (vision-cone/anisotropic-kernel) XY model, nematic version; motile +1/2 and passive -1/2 defects as in dry active nematics", "notes": "+1/2 speed rising with s and zero -1/2 drift follow from symmetry and match documented active-nematic behaviour; texture-dependent defect dynamics under an anisotropic non-mutual kernel is documented by Rouzaire et al. for the polar case. The twist not found: above s~0.15 a +1 defect splits into four defects (three +1/2 plus one -1/2, i.e. a ±1/2 pair is created) for a phase-dependent share of initial textures (0.64 at s=0.5), instead of the standard split into two +1/2. Plausibly meaningful (non-mutual coupling creating defect pairs) but needs checks for edge seeding (open boundaries, L=100), brief pairs that annihilate, and lattice/detection artifacts. Only abstracts and snippets were searched; a nematic follow-up to Rouzaire et al. may already contain this."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass. Most of the behaviour is well documented. The one new piece is that a +1 defect splits into four defects above s ≈ 0.15, depending on the texture phase. I couldn't find that described anywhere, but it isn't enough to call the whole result undocumented.

**What's already documented**
- **Defects in non-mutual XY models.** The closest polar relative is Rouzaire, Pearce, Pagonabarraga and Levis (PRL 2025). They weight the XY coupling by an anisotropic kernel, so that besides the topological charge, the actual shape of the defects becomes crucial to their dynamics. They also report that nonreciprocal coupling twists the spin field, selecting specific defect shapes, dramatically altering the pair annihilation process. Your rule looks like the nematic version of their kernel: a 2φ dependence, so a site listens equally to both ends of its axis. Defects whose behaviour depends on the texture phase is exactly their main point.
- **Defects moving in one direction.** In the vision-cone XY model, defects propagate in a directional manner, thereby breaking the parity and time-reversal symmetry of the spin dynamics (Loos, Klapp & Martynec, PRL 2023). Other work shows that non-mutual alignment is akin to a self-propulsion.
- **+1/2 moves, −1/2 doesn't.** This is the standard active-nematic picture. The +1/2 is comet-shaped with a head–tail polarity, while the −1/2 has threefold symmetry, and the −1/2 remains a non-motile particle. Self-propelled +1/2 defects also appear without any fluid flow, in nematic active solids (Brauns et al., PRL 2026). Your results (+1/2 speed rising with s, zero drift for −1/2) follow from symmetry and are expected.
- **+1 splitting in two.** In nematics, two half-integer defects have less elastic energy than one full-integer defect, so any full-integer defect should be inherently unstable. Splitting into two +1/2 defects is the expected baseline.

**The twist I couldn't find**
Above a threshold s ≈ 0.15, a +1 defect splits into four defects for about 64% of starting phases at s = 0.5, and never at s = 0.10. Because total charge is conserved, four defects must be three +1/2 plus one −1/2, so an extra ±1/2 pair is being created. I found no paper describing this or its dependence on the starting phase.

Does it matter? Plausibly yes, a little: it would mean non-mutual coupling can create defect pairs from a single integer defect, not just change how existing defects move. But it fits easily within the Rouzaire-type "shape matters" picture. Note that my searches only reached abstracts and snippets, and the nematic follow-ups to these papers could easily contain it.

**Real pattern or artifact?**
The motility results look like real physics with no sign of artifacts. The four-defect splitting is plausible but not yet verified. Before treating it as real, I'd check that:
- the extra ±1/2 pair survives over time rather than appearing briefly and annihilating;
- the result doesn't change with lattice size (100 is small, and the edges are open, so the edges could be seeding the extra pair);
- it isn't caused by the winding-detection method or by pinning to the triangular lattice;
- the onset at 0.15 holds up with a finer scan of s.

**

## Reply R35

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Nonreciprocal (anisotropic-kernel / vision-cone) XY model, nematic analogue; nonreciprocity-induced +1/2 defect self-propulsion with immotile -1/2 (dry active-nematic defect phenomenology)", "notes": "Core behaviour is documented: +1/2 self-propulsion growing with s, -1/2 immotile by symmetry, +1 instability. Not found: the phase-dependent, above-threshold (s~0.15) splitting of +1 into 3(+1/2)+1(-1/2). Most likely the standard two-way split plus one extra nucleated +-1/2 pair (consistent with twisting/nucleation effects in the polar nonreciprocal XY papers). Should be checked against boundary effects, system size, time step and noise before calling it new. No nematic (2-theta, cos2phi kernel) nonreciprocal lattice paper located; only abstracts were checked."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium confidence.** I agree with the first pass's family assignment, and I'm narrowing its claim about what is new. The core behaviour is documented. The new part is how a +1 defect breaks up, and that may be a lattice or boundary effect rather than new physics. I found no paper on this exact nematic (2θ, fore-aft-symmetric kernel) model.

**What is documented**
- **The model family.** A 2D XY model where each spin only listens to neighbours within a cone around its own direction is published. In it, "defects propagate in a directional manner, thereby breaking the parity and time-reversal symmetry of the spin dynamics" (Loos, Klapp & Martynec, PRL 130, 198301, 2023).
- **Your kernel type.** The polar version of your rule is also published: "Nonreciprocity is introduced by weighting the coupling strength of the two-dimensional XY model by an anisotropic kernel." There, defect shape matters as well as charge, and nonreciprocal coupling twists the field and selects specific defect shapes (Rouzaire, Pearce, Pagonabarraga & Levis, PRL 134, 167101, 2025). This is the closest match I found, but it is polar, not nematic.
- **Nucleation near defects.** In a related nonreciprocal lattice XY study, ±1 pairs are shielded from each other, delaying annihilation and "allowing time for the nucleation of further defects" (Popli et al., PRL 2025). That is a plausible route to extra defects appearing near a split +1.
- **Moving +1/2, stationary −1/2.** This is the textbook active-nematic result: the active backflow of a −1/2 defect vanishes at its core because of its threefold symmetry, while +1/2 defects self-propel. It has also been shown without any fluid flow: "self-propelled +1/2 defects that move relative to the elastic medium by local remodeling of the nematic texture without advection" (Brauns et al., PRL 2026).
- **+1 splitting into two +1/2.** This is standard. Two half-integer defects cost less energy than one +1 (any full-integer defect should be inherently unstable), and in active nematics a +1 defect splits into two +1/2 defects.

So your +1/2 speed rising with s and the −1/2 drift of zero are expected, not new.

**The twist I could not find**
Above s ≈ 0.15, a +1 sometimes splits into four defects (three +1/2 and one −1/2) instead of two, and whether it does depends on the starting texture phase μ (64% of phases at s = 0.5, none at s = 0.10). I found no description of this four-defect outcome, its threshold, or its phase dependence. It is probably just the usual two-way split plus one extra ±1/2 pair created near the core. That would fit the field-twisting and nucleation effects described in the polar papers.

Does it matter? Moderately at most. It would only count as new if it holds up when you:
- check whether the extra pair is created at the core or at the open boundary;
- change lattice size and integration time step;
- check how long it lasts (an extra pair that briefly appears and then annihilates is not a distinct outcome);
- add a small amount of noise.

**Genuine or artifact?** The defect motion is genuine: it comes from the non-symmetric influence term, and it is zero at s = 0. The four-defect share comes from a single deterministic run per starting phase on a 100-site open lattice. So it could be sensitive to the lattice geometry or the edges, but I see no sign of a frozen, inert or numerically diverged state.

**

## Reply R38

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Nonreciprocal XY (NRXY) model, nematic analogue; motile +1/2 defects as in dry active nematics", "notes": "No paper found on a nematic (cos 2phi) nonreciprocal lattice model, nor on nonreciprocity-driven splitting of a +1 defect into four +1/2 defects with a texture-phase dependence. Zero -1/2 drift is expected by symmetry. Lattice pinning could make the s ~ 0.15 threshold and the 0.64 share lattice-specific. Only abstracts and snippets were read; Dadhichi 2020 was not checked directly."}

**Reply text:**

**Verdict: KNOWN-VARIANT, with low-to-medium confidence.** The model family and the self-propulsion mechanism are described in the literature. The nematic (apolar) version of this lattice rule is not. I also found no source for the four-defect splitting of a +1 defect.

**1. Model family.** This is a nonreciprocal planar-spin (XY-type) lattice model, a dry active-matter model. Each spin's coupling to a neighbour depends on where that neighbour sits in the spin's own frame. Here the weight exp(s·cos 2φ) has period 180°, so it is the nematic analogue of the polar "nonreciprocal XY" (NRXY) model. The zero-noise version is deterministic gradient-like (Kuramoto-type) dynamics on a triangular lattice.

**2. What the literature says.**
- **Polar NRXY model, same lattices and the same mechanism.**
  - Popli, Maitra and Ramaswamy (Phys. Rev. Lett. 135, 088303, 2025) examine the nonreciprocal planar spin model on square and triangular lattices, building on the original model of Dadhichi et al. (Phys. Rev. E 101, 052601, 2020).
  - Their companion paper, "Inescapable Anisotropy of Nonreciprocal XY Models" (Phys. Rev. Lett. 135, 088302, 2025; I inferred the 088302 number from the ADS bibcode 2025PhRvL.135h8302D), shows that nonreciprocity is akin to a self-propulsion.
  - The same paper reports that, as in flocks, the asymptotic state is a dynamic foam of asters.
  - Popli et al. also find that defect interactions are highly anisotropic with respect to the mean ordering direction.
- **Directional defect motion in nonreciprocal XY models.** Loos, Klapp and Martynec (Phys. Rev. Lett. 130, 198301, 2023; author names taken from the Nature Physics reference list) study a vision-cone version. They find that defects propagate in a directional manner. A separate paper, "Nonreciprocal Interactions Reshape Topological Defect Annihilation" (Phys. Rev. Lett. 134, 167101, 2025), reports that defect shape, not only topological charge, becomes crucial to their dynamics.
- **Motile +1/2 and passive −1/2 defects in active nematics.** This is well established. In active nematics, +1/2 defects are motile, with a self-propelling speed proportional to activity (Shankar et al., Phys. Rev. Lett. 121, 108002, 2018). The −1/2 defect creates no net flow at its position and is not self-propelled in systems with uniform activity. This matches your finding of zero drift for −1/2 defects, which is expected from its threefold symmetry.
- **Flow-free defect self-propulsion.** Brauns et al. (Phys. Rev. Lett. 136, 058302, 2026) show self-propelled +1/2 defects that move by local remodeling of the nematic texture without advection. They are in an elastic active solid, not a nonreciprocal-coupling model. I found no review of these ideas for nonreciprocal nematic lattices.

**3. What differs from your system, and whether it matters.**
- I found no paper on a nematic (cos 2φ) nonreciprocal lattice model. Your +1/2 speed rising with s is therefore a plausible and symmetry-consistent result, but I have no direct source for it.
- The +1 defect splitting into four +1/2 defects above a threshold s ≈ 0.15, with the splitting share depending on the texture phase μ, was not found. Two things are known. One is that integer defects are generically unstable to splitting into half-integer pairs, because the elastic energy scales as m². The other is that active nematics form virtual full-integer defects such as asters and vortices from bound +1/2 defects. Splitting into four, driven by nonreciprocity alone and depending on μ, is the part that looks new. I have not shown that it is.
- Because the system has no noise, a lattice-pinning effect cannot be excluded. Popli et al. report that lattice anisotropy persists at large scales in the polar model, so the 0.15 threshold and the 0.64 share may be lattice-specific.

**4. Confidence and gaps.** My confidence is low-to-medium. The searches were limited. I could not read full texts, only abstracts and snippets. I did not find the Dadhichi 2020 paper's own title or authors in the results, so I give it only as cited by Popli et al.

**

## Reply R57

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Nonreciprocal anisotropic-kernel XY model (Rouzaire-Pearce-Pagonabarraga-Levis 2025), nematic (2-theta, fore-aft-symmetric kernel) variant; flow-free +1/2 defect self-propulsion as in dry/active nematics", "notes": "Core behaviour (+1/2 drift that grows with s, no -1/2 drift, shape/phase-dependent defect dynamics) is documented in the polar anisotropic-kernel NRXY and active-nematic literature. The kernel exp(s cos 2phi) is fore-aft symmetric, so this is an axial kernel, not a vision cone; +1/2 motion comes from the comet shape plus the non-mutual dynamics. Twist not found: a threshold (s~0.15) splitting of +1 into 3(+1/2)+1(-1/2), i.e. nonreciprocity creating an extra +-1/2 pair at the core, with a phase-dependent share (0.64 at s=0.5). Moderately significant if robust. Rule out boundary nucleation (open edges), lattice-scale core effects, and finite-size effects (N=100); the Rouzaire et al. full text and supplement were paywalled and not checked."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass's verdict, though one point in its framing needs fixing. The core behaviour is documented in nearby work. The one twist I couldn't find described is the four-defect splitting above a threshold. I found nothing showing it's already known, but my search was limited, so the gap may be a search gap.

**What the rule is.** It is the anisotropic-kernel nonreciprocal XY model of Rouzaire, Pearce, Pagonabarraga and Levis (PRL 134, 167101, 2025), with the angle doubled to make it nematic. In that paper, nonreciprocity is introduced by weighting the coupling strength of the two-dimensional XY model by an anisotropic kernel. They also show that in addition to the topological charge, the actual shape of the defects becomes crucial to faithfully describe their dynamics. Nonreciprocal coupling twists the spin field, selecting specific defect shapes. Your texture phase mu is exactly this shape variable. The wider model family (spins that weight neighbours by position relative to their own direction) is reviewed in Dopierala et al., who show that nonreciprocity is akin to a self-propulsion. Loos et al. found that vision-cone interactions cause defects to propagate in a directional, parity-broken manner.

**Correction to the first pass.** Calling this a "vision cone" is misleading. A weight of exp(s·cos 2φ) treats front and back equally, so this is an axial "listen along my own axis" kernel, not a cone. That also means the +1/2 self-propulsion can't come from a forward bias in the coupling. It has to come from the comet shape of the +1/2 defect combined with the non-mutual (and so non-energy-minimising) part of the dynamics. That fits a known phenomenology: active-nematic-like +1/2 motion with no fluid flow. Brauns et al. (PRL 2026) report self-propelled +1/2 defects that move relative to the elastic medium by local remodeling of the nematic texture without advection.

**Feature by feature:**
- **+1/2 self-propulsion that gets faster as s rises: KNOWN in kind.** It is the standard active-nematic signature, here coming from nonreciprocity, consistent with the "nonreciprocity is akin to self-propulsion" result.
- **−1/2 defects don't drift: KNOWN and expected.** The negative defect has three equivalent polarizations, so by symmetry no drift direction is picked out. In active nematics the −1/2 is not self-propelled in systems with uniform activity.
- **A +1 breaking up into half-integer defects: KNOWN.** Higher-charge defects are energetically unstable to disassociation into elementary defects. In passive nematics the usual result is two +1/2 defects.
- **The twist I couldn't find described:** above s ≈ 0.15, a +1 splits into four defects (three +1/2 and one −1/2, so the total charge stays +1) in 64% of initial phases at s = 0.5 and in none at s = 0.10. In other words, the nonreciprocity creates an extra ±1/2 pair at the defect core, and whether that happens depends on the texture phase. Rouzaire et al. discuss shape-dependent dynamics, but the summaries I could read (abstract and a conference abstract) don't describe splitting or a threshold. I couldn't read the full paper or its supplement, both paywalled, and I found no nematic version of this model in any paper.

**Does the twist matter?** Moderately. It would mean nonreciprocity creates new defects rather than just moving existing ones, which bears on whether defects multiply and spread. That would echo the metastable ordered phase that is metastable to the nucleation of topological defects.

**Genuine pattern or artifact?** Probably genuine, not trivial. The drift trends are smooth, show the right symmetry, and switch on at a threshold. Three caveats:
- The four-defect events could be a lattice-scale effect at the defect core, since Dopierala et al. show lattice anisotropy is "inescapable" in these models.
- With open edges and N = 100, defects may be created near the boundary.
- A 0/1 share sampled over a handful of phases is coarse.

Before claiming the twist as new, check three things: that the extra pair forms at the core and not at an edge, that the share holds up at larger lattice sizes, and that it also happens on a square lattice.

**

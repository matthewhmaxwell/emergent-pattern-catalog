# planted:N08

**Reported behavior (what the authors report as new):** Purely repulsive inertial dumbbells driven by constant torques phase-separate into a dilute liquid and a dense hexatic phase (binodal, spinodal, t^1/3 coarsening), with rigidly rotating clusters; inertia is essential.

**Rule on the card:** World: about 2.1 million dumbbells (two identical beads of diameter 1 joined by a stiff spring) in a square periodic 2D box at area fraction 0.65. Rules, Langevin dynamics with mass: (1) every pair of beads repels on contact (purely repulsive, no attraction); (2) each dumbbell feels a constant torque, equal and opposite forces perpendicular to its axis on its two beads, all spinning the same way, strength set by Peclet number 5 (swept); (3) beads have mass and weak friction (momentum persists; Prandtl number 0.005) plus thermal noise at temperature 0.01. Start homogeneous and watch the local density field.

**Observed features on the card:**
  - onset_torque_peclet = 1.0 (torque strength above which the uniform state splits apart)
  - domain_growth_exponent = 0.333 (exponent of domain size growth with time after splitting)
  - coexisting_states = 2 (coexisting states: dilute liquid and dense six-fold ordered state)
  - example_area_fraction = 0.65 (global area fraction of the illustrated split state)

## Reply R07

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Rotation-induced phase separation (RIPS) of inertial chiral spinners, here with dumbbells and Lifshitz-Slyozov-like coarsening", "notes": "I found no paper on torque-driven (pure-torque) inertial dumbbells with a hexatic dense phase and a growth exponent of 1/3, so this is a variant of the disk RIPS result, not an exact match. I could only read the PRL abstract and fragments, so I could not check whether it reports a hexatic dense phase, a growth exponent or an onset threshold. The Pe about 1 threshold was not found in any readable source. The 1/3 exponent is the generic Lifshitz-Slyozov value and is unsurprising. Check the full text of the PRL and arXiv 2504.08533 for dumbbell or growth-law content before treating anything as new."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** Phase separation of repulsive, inertial spinners with a constant torque is published for disks. I found no paper that reports the exact dumbbell version with a hexatic dense phase and a growth exponent of 1/3.

**1. Model family.** This is chiral active matter: inertial self-spinning particles under underdamped Langevin dynamics with purely repulsive interactions. The relevant sub-families are "rotation-induced phase separation" (RIPS) of chiral spinners and active or passive dumbbell fluids that show liquid–hexatic coexistence.

**2. What the literature says.**
- **Closest source (disks).** Digregorio, Pagonabarraga and Vega Reyes, "Phase Separation in a Chiral Active Fluid of Inertial Self-Spinning Disks", Phys. Rev. Lett. 136, 218301 (published 27 May 2026; arXiv 2504.08533).
  - Systematic particle rotations in a fluid of disk-shaped spinners can spontaneously lead to phase separation.
  - A pressure feedback mechanism drives it when the active rotation is not counterbalanced by translational friction. They call it rotation-induced phase separation (RIPS).
  - The setup matches yours: Langevin dynamics with a thermal bath and a homogeneous, constant torque.
  - The mechanism needs large translational inertia. Their appendix says phase separation is observed only in the regime of large particles' translational inertia, which is consistent with your weak friction.
  - The abstract describes a low-density phase coexisting with a dense chiral liquid. The excerpts I could read do not give a six-fold ordered dense phase, a growth exponent, or dumbbells.
  - The same paper mentions arrested phase separation and an unsteady regime at large spin rates.
- **Earlier related spinner work.** Ding et al., "Odd Response-Induced Phase Separation of Active Spinners", Research 7:0356 (2024).
  - This is repulsive spinners under external shear, so it is a different mechanism. The dense phase there is closer to a solid stripe with a hexagonal lattice structure than to a liquid.
  - Massana-Cid, Levis, Hernández, Pagonabarraga and Tierno, "Arrested phase separation in chiral fluids of colloidal spinners", Phys. Rev. Research 3, L042021 (2021). I saw only its reference-list entry, not the paper, so its content is not verified.
- **Dumbbells, liquid–hexatic coexistence.**
  - Passive dumbbells show a first-order liquid–hexatic transition.
  - Active dumbbells with a transverse (perpendicular) force also phase-separate into a dilute phase and a hexatic dense phase. This is Entropy 27, 692 (2025), "Transverse Self-Propulsion Enhances the Aggregation of Active Dumbbells". Across all activity, it finds phase separation into dilute and dense phases. The dense phase is hexatic.
  - That model uses a constant transverse force per dumbbell. This is related to, but not the same as, equal and opposite perpendicular forces on the two beads, which give a pure torque with no net force. It is also not an inertial system.
- **Coarsening kinetics in 2D active dumbbells.** Caporusso, Cugliandolo, Digregorio, Gonnella and Suma, Soft Matter 20, 4208 (2024).
  - They report growth faster than in active disk models, with a non-universal exponent, in self-propelled (axial-force) dumbbells.
  - The exponent of about 0.6 is a different class from your 1/3.
- **Review.** Cates and Tailleur, "Motility-induced phase separation and coarsening in active matter" (C. R. Physique 16, 2015). It covers active dumbbells, with rotating clusters and hexatic and nematic order. It does not treat torque-driven spinners.

**3. The 1/3 exponent.** A growth exponent of 1/3 is the standard Lifshitz–Slyozov value for coarsening with a conserved density. A 2026 review of two-temperature phase separation (via a Pith summary) reports a high-density growth exponent of 0.33, identified with the Lifshitz–Slyozov class. That is a different system, so I would treat 1/3 here as plausible and unsurprising, not as a new result. I did not find it reported for spinning dumbbells.

**4. What differs from your system, and whether it likely matters.**
- *Dumbbells instead of disks.* This probably does not change the instability. It plausibly explains the six-fold ordered dense phase, since dumbbell fluids freeze through a hexatic phase at area fraction near 0.65 and above. This is untested.
- *Dense phase character.* The PRL abstract says "dense chiral liquid". A hexatic or six-fold ordered dense phase is a difference to check.
- *Onset at Pe about 1.* I did not find this threshold in any source I could read.

**5. Confidence: medium.** The mechanism and the gas–dense coexistence are firmly documented for inertial spinners. The dumbbell and hexatic specifics rest on neighbouring papers, and I could only read abstracts and fragments of the PRL.

**

## Reply R25

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Rotation-induced phase separation (RIPS) in inertial chiral spinner fluids; with dumbbells, liquid-hexatic coexistence", "notes": "Torque-only inertial repulsive spinners phase-separating out of a uniform state is documented as RIPS (disks). Twist not found in the literature: the dense phase here is hexatic/six-fold ordered (RIPS reports a dense chiral liquid with packing near 1 and spin near 0), and coarsening keeps going with exponent 1/3 instead of the arrested phase separation RIPS reports in some regime. The hexatic order probably comes from the passive dumbbell liquid-hexatic transition, so it is plausible but likely minor. 1/3 is the generic conserved-coarsening value. The Pe~1 onset was not checked against a readable source. Artifact checks: low temperature 0.01 plus weak friction means strong heating from the spinning, so rule out a jammed dense state and fit the exponent over a long time range."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass, with one correction: the dense phase is the part that isn't in the source paper. I could only read fragments of the main paper's full text, not the whole thing.

**The behaviour itself is documented.** A 2026 Physical Review Letters paper (Digregorio, Pagonabarraga, Vega Reyes) describes it. It uses inertial disks that repel only on contact and each feel the same constant torque. It shows that systematic particle rotations in a fluid composed of disk-shaped spinners can spontaneously lead to phase separation, starting from a uniform state. The authors call this rotation induced phase separation (RIPS), when the active rotation is not properly counterbalanced by translational friction. That fits your weak friction with momentum carried over. The ingredients match yours: a purely repulsive Weeks-Chandler-Anderson (WCA) potential, translational and rotational inertia, and a constant torque. The dense phase they report has packing fraction close to 1 and spin near zero, alongside a loose phase. They argue the effect should be general: RIPS is expected to appear as a general nonequilibrium phenomenon in chiral fluids. Your dumbbells pass rotation into translation through their shape rather than through a sliding force at contact. That is a different rule, not a different behaviour, and it is the kind of generalisation the paper predicts.

**Each feature separately is also on record:**
- **Six-fold ordered dense phase of dumbbells.** Passive dumbbells already have a first-order liquid–hexatic transition. Dumbbells pushed sideways, which also puts torque on clusters, split into a dilute phase and a dense phase with hexatic ordering at all Péclet numbers.
- **Six-fold order and chirality.** In a spinning-disk fluid with attraction, chirality can also break the solid phase, giving rise to a dense fluid made of rotating hexatic patches.
- **Growth exponent of 1/3.** This is the textbook coarsening value for systems that conserve particle number, so it adds nothing new.
- **Onset near Pe ≈ 1.** A threshold is expected because RIPS depends on spinning beating friction. I could not check the exact number against the paper.

**The twist in what you see:**
1. The RIPS paper calls its dense phase a dense chiral liquid. I found no report of a hexatic or crystalline dense phase in a system driven only by torque. With dumbbells it may simply come from the passive liquid–hexatic tendency. That would make it plausible but probably not important.
2. In at least one regime (described in an appendix), the RIPS paper says RIPS is characterized by arrested phase separation, forming an array of domains. Your domains keep coarsening with exponent 1/3. If that holds up, it is a real difference, though it may just reflect where you are in the parameter space.

**Genuine or artifact?** It looks genuine, not a trivial outcome. There is a nonzero threshold, the passive system stays uniform at 0.65, the domains coarsen, and the two coexisting phases are physically distinct. Three cautions:
- The temperature of 0.01 is very low, and with weak friction the spinning will heat the system far above it. Check that the dense phase isn't just a jammed, frozen region.
- With weak friction, momentum-driven coarsening (exponent 1/2 or 2/3) would also be possible. The value 1/3 may only describe a friction-dominated window, so fit it over a long time range.
- Dumbbell clusters driven by self-propulsion grow faster than disk clusters (Caporusso et al. 2024), so check whether 1/3 is the true long-time value.

## Reply R34

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Rotation-induced phase separation (RIPS) in inertial chiral spinner fluids; here with dumbbells giving liquid-hexatic coexistence", "notes": "Torque-only inertial phase separation into dilute and dense phases is documented (RIPS, disks with contact friction). Twists not found in the literature: (1) the rotation-translation coupling comes from dumbbell shape with purely repulsive contacts, not contact friction; (2) the dense phase is hexatic rather than a chiral liquid (expected for dumbbells, whose passive and transversely driven forms show liquid-hexatic coexistence); (3) the disk RIPS SI reports arrested phase separation (cluster arrays) in some regime, while this system shows continued t^(1/3) coarsening. If the coarsening is truly unarrested over long runs, that is the most notable difference. The Pe ~ 1 onset was not found in any readable source and may just reflect torque-driven energy exceeding the low thermal noise. Suggested checks: kinetic temperature in each phase (cold-dense/hot-dilute splits are known from inertial dumbbell MIPS), particle mobility inside the hexatic phase to rule out jamming, and longer runs to rule out late-time arrest."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass, with medium confidence. The behaviour is documented: a chiral fluid driven only by torque splits from a uniform state into a dilute phase and a dense phase. The twists are in the dense phase's structure, in how rotation couples to translation, and in whether the coarsening continues or stops.

**What is documented**
- **The behaviour itself (2026 PRL; arXiv 2504.08533, by Digregorio, Pagonabarraga and Vega Reyes).** Spinning disks with inertia can spontaneously lead to phase separation, starting out of a homogeneous and hydrostatic stationary state, due to a pressure feedback mechanism. It happens when the active rotation is not properly counterbalanced by translational friction. The authors expect RIPS to appear as a general nonequilibrium phenomenon in chiral fluids. Your system has weak friction and inertia, which is exactly that regime.
- **Inertia is required.** The paper has an appendix titled "Suppression of phase separation for small inertia" and says that otherwise phase separation disappears. That fits your momentum-dominated setup.
- **Dumbbells give a dense hexatic phase in related systems.** Passive dumbbells have a liquid-hexatic phase transition of first-order. Dumbbells pushed sideways (perpendicular to their axis) separate between a dilute phase and a dense phase with hexatic ordering at all Péclet numbers. That also covers Pe = 1, so a hexatic dense phase near your onset is expected.
- **Spinning dimers are an established model.** Torque-driven spinning dimers are known to show spatiotemporal order and active melting with no equilibrium counterparts (van Zuiden et al., PNAS 2016).

**The twists I could not find described**
1. **The rotation–translation coupling.** In the disk RIPS paper it comes from rubbing contact friction between particles: tangential friction and active rotation in disks' collisions. In your system it comes only from the dumbbell shape, with purely repulsive contacts. This supports the paper's claim that RIPS is generic, rather than being a new effect.
2. **The dense phase.** The disk paper's dense phase is a dense chiral liquid; yours is hexatic. That is expected for dumbbells at these densities.
3. **Coarsening.** The disk SI says RIPS is characterized by arrested phase separation, forming an array of clusters, at least in some regimes. Your domains keep growing as t^(1/3). This difference could matter. If your growth truly continues with no plateau, it is the most interesting part of the result.
4. **The Pe ≈ 1 onset.** I could not find this threshold reported anywhere.

**Genuine or artifact?** It looks genuine: two clearly different phases and self-similar growth. Some caveats:
- With weak friction and constant torque the system heats up through collisions. A "cold dense / hot dilute" split is already known from inertial MIPS in dumbbells (motility-induced phase separation, where self-propelled particles cluster). It is worth checking the kinetic temperature of each phase.
- An onset at Pe ≈ 1 may just mark where the torque-driven energy overtakes the very low thermal noise (T = 0.01).
- The hexatic dense phase could be a slowly arrested, locally jammed state. Measure the hexatic correlation length and check that particles still move in the dense phase.
- Before saying "unarrested", run long enough to rule out the cluster-array arrest that the disk paper reports.

**Searches I ran:** RIPS with dumbbells; arXiv 2504.08533; torque-driven dumbbells with inertia; Gonnella/Suma-group dumbbell papers; van Zuiden spinner dimers; self-spinning dumbbells with a hexatic dense phase. I could not read the full PRL text, so its growth-law content is still unchecked.

## Reply R39

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Rotation-induced phase separation (RIPS) in inertial chiral spinner fluids, with liquid-hexatic coexistence as in active/passive dumbbell systems", "notes": "Dilute-liquid / dense-hexatic splitting driven only by spin, with inertia, is documented (RIPS for frictional spinner disks; hexatic dense phase for dumbbells). Twist: rotation is turned into motion by the dumbbell shape instead of contact friction, and coarsening continues as t^(1/3), while RIPS disks and colloidal spinners show arrested coarsening. The continued coarsening matters if it holds at late times. Onset Pe~1 and the 1/3 exponent were not verified against any source; 1/3 is the generic exponent for conserved-density coarsening. Checks to run: Pe=0 should be uniform at phi=0.65 (passive dumbbells have liquid-hexatic coexistence), and measure kinetic temperatures per phase because T=0.01 is small next to torque-driven heating."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass, but on narrower grounds. Splitting into a dilute liquid and a dense six-fold (hexatic) state is documented for closely related systems. What I couldn't find is the combination of torque-only drive, dumbbell shape, no contact friction and coarsening that keeps going as t^(1/3).

**What is documented**
- **Spinning disks with inertia (RIPS).** A 2026 PRL (PRL 136, 218301; arXiv 2504.08533) shows that particles which only spin can split into dense and dilute phases. They call it rotation-induced phase separation (RIPS). It happens when the active rotation is not properly counterbalanced by translational friction. The authors argue that RIPS is expected to appear as a general nonequilibrium phenomenon in chiral fluids. In that model, though, rotation turns into motion through friction between particles: between active rotation and interparticle friction.
- **Dumbbells with a dense hexatic phase.** For self-propelled dumbbells, dumbbells tend to form a single global hexatic and dense phase coexisting with the liquid, and this regime continuously extends to the passive dumbbell case. Dumbbells pushed sideways to their axis (the closest geometry to a torque couple) separate between a dilute phase and a dense phase with hexatic ordering at all Péclet numbers. A side note on scales: in this group's convention, Pe=2σdFa/kBT.
- **Coarsening.** Domain growth in 2D active dumbbells has been measured. The growth observed is faster than in active disk models, and this effect is further enhanced under stronger activity. An exponent of 1/3 is the generic value for coarsening when the number of particles is conserved, so it is unsurprising.
- **Inertial dumbbells.** Phase separation of underdamped dumbbells has been studied, with four distinct kinetic temperatures across the coexisting phases. That work uses self-propulsion, not torque.

**The twist**
1. **How rotation becomes motion.** Here pure repulsion plus the elongated shape converts spin into translation, with no contact friction. This is a sensible stand-in for the friction in the RIPS paper, which covers the mechanism in principle. The specific geometric route is not documented.
2. **Coarsening keeps going.** The supplement of the RIPS paper shows that RIPS is characterized by arrested phase separation, forming an array of separate domains. Separate experiments with magnetic colloidal spinners also find stopped coarsening, there due to fluid-mediated interactions (hydrodynamics dominate over attractive magnetic interactions and impede coarsening). Steady t^(1/3) growth would therefore be a real difference, if it is not a short-time or system-size effect. It matters to whether RIPS is generic. A good test is to run longer times and check for a plateau.

**Genuine pattern or artifact?** It looks genuine: two clearly different coexisting states, a threshold, and a growth law. Some caveats remain:
- **Temperature.** At T=0.01 with weak friction, much of the effective temperature is set by the torque. Pe=1 may simply be where torque-driven heating outweighs the thermal noise. Report kinetic temperatures in each phase.
- **Passive baseline.** Passive dumbbells already show liquid–hexatic coexistence. Check that the Pe=0 run is uniform at area fraction 0.65.
- **Untested numbers.** I found no source giving the onset at Pe≈1 or the 1/3 exponent for this exact system.

I could not open the full PRL text (only abstract and snippets), so the "arrested" claim rests on supplement fragments. I also could not identify any torque-driven dumbbell paper.

**

## Reply R54

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Rotation-induced phase separation (RIPS) in an inertial chiral spinner fluid; dense phase hexatic as in dumbbell liquid-hexatic coexistence", "notes": "Same verdict as the first pass, with stronger evidence. van Zuiden et al. 2016 simulated essentially the same particles (inertial repulsive dimers with force-dipole torque) but reported crystal, liquid and jammed phases, not torque-driven dilute/dense coexistence with coarsening. Disk RIPS (2026) gives the mechanism and predicts it should occur in chiral fluids generally. Twists not found in the literature: (a) torque turns into translation through dumbbell shape, with no friction between particles; (b) the dense phase is hexatic, not a dense chiral liquid with spin near zero; (c) uninterrupted t^(1/3) coarsening, while disk RIPS shows arrest in some regimes. Things to check: whether 0.333 is assumed or fitted over a short window, whether the Pe=1 onset is just a sweep grid point, and whether torque-shifted effective temperature produces equilibrium-like liquid-hexatic coexistence instead of RIPS. Only abstracts and snippets were read."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass, but for somewhat different reasons. The behaviour looks like a real emergent effect, not an artifact. Confidence is medium.

**What the literature covers**

1. **The mechanism is documented for inertial spinning disks.** Digregorio, Pagonabarraga & Vega Reyes (PRL 136, 218301, 2026) show that particles spinning the same way can spontaneously lead to phase separation. They name it rotation induced phase separation (RIPS), when the active rotation is not properly counterbalanced by translational friction. That fits your setup: weak friction, momentum that persists, uniform torque, and a split out of a uniform state. They also say RIPS is expected to appear as a general nonequilibrium phenomenon in chiral fluids. So a dumbbell version is something the theory predicts, not a new effect.

2. **The same particle model has been simulated before, without this finding.** van Zuiden, Paulose, Irvine, Bartolo & Vitelli (PNAS 113, 12919, 2016) used dimers, each consisting of two point particles of mass m connected by a stiff link, with a torque implemented as a force dipole and purely repulsive interactions. That is essentially your rule set (theirs used soft Yukawa repulsion). They reported spinner crystals, a liquid with edge currents, and jamming. In the snippets I saw, they did not report a dilute liquid coexisting with a dense ordered phase, or any coarsening. This is the strongest prior work, and the first pass missed it.

3. **A dilute liquid next to a dense hexatic phase is standard for dumbbells.** Dumbbells pushed sideways (perpendicular to their axis) split into a dilute gas and a hexatically ordered dense phase at all Péclet numbers from Pe=1 to Pe=100. Even without activity, dumbbells show a liquid-hexatic phase transition of first-order, characterized by a coexistence region. Growth of domain size as t^(1/3) is the generic result for density separation where particle number is conserved.

**The twist in the behaviour (not just the rules)**
- **Ordered dense phase.** In the disk version, the dense phase is a dense chiral liquid, with packing near 1 and spin near zero. Here the dense phase is six-fold ordered, probably because the dumbbells lock together.
- **Uninterrupted coarsening.** The disk paper notes that RIPS is characterized by arrested phase separation in at least part of its phase diagram. A clean t^(1/3) growth with no arrest would be a difference.
- **No friction needed between particles.** Disk RIPS relies on active rotation and interparticle friction. In dumbbells, the shape alone transfers spin into translation. This is the twist that plausibly matters most, but it is a difference in mechanism, and the overall behaviour is the same family.

I could not find any paper on torque-driven, frictionless inertial dumbbells that separate into a liquid and a hexatic phase.

**Possible artifacts to check**
- **The exponent is suspiciously exact.** 0.333 looks like it was assumed or fitted over a short window rather than measured. With weak friction, momentum-driven coarsening (exponents of 1/2 to 2/3) might show up at intermediate times.
- **The onset value may just be a grid point.** Pe = 1 is round enough that it could be a value from the sweep rather than a measured threshold.
- **It could be ordinary coexistence in disguise.** Without strong friction, the torque may change the effective temperature. Then the split could just be equilibrium-like liquid–hexatic coexistence at a shifted temperature, rather than RIPS. Passive dumbbells at 0.65 should be in the liquid, which argues against this. Comparing the kinetic temperature of each phase with the passive phase diagram would settle it.
- **Related instabilities.** The RIPS paper itself compares its mechanism to clustering in heated granular fluids.

None of these points to a diverged simulation, a frozen absorbing state or a random colouring. Two coexisting phases with domains that keep growing looks like genuine emergent behaviour.

**

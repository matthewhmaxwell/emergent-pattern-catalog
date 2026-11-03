# planted:N16

**Reported behavior (what the authors report as new):** A self-avoiding ring of active Brownian beads with rotational inertia collapses into a dense globule (over tenfold smaller) at intermediate activity with outward-polarised propulsion, then re-swells at high activity; translational inertia opposes collapse.

**Rule on the card:** World: one closed loop of 2048 beads (a chain joined end to end) in open 3D space. Rules, Langevin dynamics with inertia: (1) consecutive beads are linked by springs and all beads repel at contact (no attraction), so the chain cannot cross itself; (2) each bead pushes itself with constant force along its own heading (activity: Peclet number Pe, swept); (3) headings diffuse randomly but with rotational inertia (reduced moment of inertia J swept 0.001-100, e.g. J=5); (4) beads have translational inertia (reduced mass M=0.3), friction and thermal noise. Measure loop size (radius of gyration) versus Pe, J and bead number.

**Observed features on the card:**
  - size_reduction_factor = 10 (squared loop size shrinks by over this factor versus passive)
  - activity_of_max_compaction = 17 (Peclet number giving the most compact loop)
  - reswelling_activity = 40 (Peclet number above which the loop swells back open)
  - size_vs_length_exponent_j10 = 0 (exponent of loop size versus bead number at J=10 (near zero))

## Reply R09

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "activity-induced collapse (and re-swelling) of active ring polymers; inertial active Brownian polymer", "notes": "Ring compaction by activity with steric repulsion, non-monotonic size vs activity, and inertial ABP physics are each documented. No source found that sweeps rotational inertia J for an ABP ring, so the J-window and the N-independent size at J=10 (exponent ~0, odd for a compact globule) are not confirmed in the literature. Check whether the N-independence is real or a numerical artifact such as chain crossing or bead overlap at high inertia. Search tool returned mostly abstracts, so a paper may have been missed."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium confidence.** The general phenomenon is well described. I found no paper that covers the specific J-controlled version.

**1. Model family**
This is an active polymer, or active Brownian polymer (ABPO), in the active-matter / polymer-physics family. It is a closed ring of self-propelled beads with random-direction (ABP-type) propulsion, excluded volume only, and underdamped dynamics. Both translational and rotational inertia are included, which makes it an "active Langevin" or inertial-ABP system.

**2. What the literature covers**
- **Ring collapse from activity.** Locatelli, Bianco and Malgaretti (2021, PRL 126, 097801) is the earliest clear source I verified. Simulations of isolated active rings with steric repulsion found short rings (N≲100) swell, whereas longer rings (N≳200) collapse, at sufficiently high activity. They stress the role of steric interactions and topology conservation. Their collapsed rings are arrested, and some spin about their axis. Their activity is tangential (polar), not ABP-type.
- **Non-monotonic ring size with activity.** Shrinkage followed by swelling as activity rises is reported in several ring models. These include non-monotonic behavior of the radius of gyration of an active Brownian ring as a function of active force (the chiral active ring, 2D) and partially active rings in 3D. The 3D partially active rings show shrinkage followed by swelling as the number of active monomers increases. In the same vein, Kumar and Thakur (Macromolecules 2023) report polar activity swelling short rings and crumpling larger ones.
- **Collapse with no attraction.** Buglakov, Lelecova and Chertovich (PCCP 2025, vol. 27, p. 1691) studied linear chains of ABP monomers. They found that self-propelled monomer units cause chain collapse where there is no explicit attraction between particles. They link this to motility-induced phase separation (MIPS) and used underdamped dynamics. Linear ABP-polymer collapse and re-swelling are also covered in the review "The physics of active polymers and filaments" (J. Chem. Phys. 153, 040901, 2020).
- **Inertia in active polymers.** Tejedor, Ramírez and Ripoll (Phys. Rev. Research 6, L032002, 2024) study linear polar chains. They find that breakdown of self-similarity occurs both in the absence and presence of inertia, although it is stronger in the latter case. A separate paper, "Effects of inertia on conformation and dynamics of tangentially driven active filaments" (PRE 108, 024606, 2023), finds that inertial collisions increase the persistence length of active polymers.
- **Rotational inertia in single-particle and bulk systems.** Reviews by Löwen (J. Chem. Phys. 152, 040901, 2020) and Sprenger, Caprini, Löwen and Wittmann (J. Phys.: Condens. Matter 2023) cover inertial ABPs. In the Löwen review, long-time diffusion depends strongly on the moment of inertia J and increases with J, showing extra persistence at larger J. Rotational inertia also changes MIPS for self-propelled disks (PMC9800414). This gives a plausible mechanism for a J-dependent window: J changes the effective orientational persistence, and so changes the effective Pe.

**3. What I did not find**
I found no source that sweeps rotational inertia J for a closed ring of ABP beads. That means I found nothing on the following:
- J-dependent compaction.
- A narrow Pe window of re-opening at large J (J=30).
- A ring size nearly independent of N at J=10.

Rotational-inertia work on polymers is silent in my searches. The ring-collapse work I found uses tangential or overdamped activity.

**Differences that likely matter**
- Propulsion here is bead-heading ABP-type, with inertial heading dynamics. In the ring-collapse papers it is tangential and mostly overdamped or low-inertia.
- J is a new control axis. It probably acts through persistence time, but that is my inference, not something a source states.
- Sources report collapse at large N (Locatelli: N≳200 for tangential rings) and swelling of short rings. A ring size independent of N at J=10 (exponent ≈ 0) is not what I would expect from an ordinary compact globule, which should scale roughly as N^(1/3). This is not reported anywhere I found. Check whether it is a real feature, such as a size set by an active length scale, or an artifact. One possible artifact is bead overlap or chain crossing when inertia is high and the repulsion is soft. Checking the time step and the topology (no strand crossings) would settle this.

**Confidence and reasoning**
Medium. The ring and ABP-polymer compaction, the non-monotonic size, and the re-swelling at high Pe are all well documented. The specific J-dependence and the N-independence are not in anything I could read. The search tool returned few abstracts for the inertial-polymer papers, so a paper I missed could cover them.

**

## Reply R11

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Motility-induced collapse (and high-activity re-swelling) of an underdamped active Brownian polymer; non-monotonic Rg of active Brownian rings", "notes": "Collapse, re-swelling and inertia dependence are documented for the linear underdamped ABP chain; non-monotonic Rg is documented for overdamped ABP rings. The ring + inertia combination appears unstudied, but that is a rule-combination gap, not new phenomenology. Not found in the literature: a size-vs-N exponent of ~0 at J=10. It is physically suspicious (a dense globule would give ~1/3) and should be checked for incomplete equilibration, a narrow N range, or a measurement artifact before it is treated as real. The sharp J-dependent re-opening is a quantitative refinement of documented inertia effects. Full texts were not read."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass. The core behaviour is documented: a chain with no attraction collapses when its beads self-propel, inertia changes where that happens, and the chain swells back at high activity. Closing the chain into a loop is not new for this model family either. The one unexplained feature is the size that doesn't grow with bead number (exponent ≈ 0), and that feature looks more like a measurement problem than a discovery. I didn't read any full text, so everything below rests on abstracts and preprint excerpts.

**Where the main behaviour is documented**
- **Collapse without attraction (linear chain).** Buglakov, Lelecova & Chertovich, PCCP 27, 1691 (2025), preprinted as "Monomer activity solely induces collapse…". They demonstrate the ability of self-propelled monomer units to cause chain collapse under conditions where there is no explicit attraction between particles. The collapsed state has a dense condensed core surrounded by elongated loops and stabilized by an effective temperature gradient. They explain it as being similar to motility-induced phase separation in active media.
- **Same inertia model.** They use a chain of active brownian particles in an underdamped regime with m - mass of the particles and J - moment of inertia.
- **Re-swelling.** Their preprint states that activity-induced collapse of the chain of self-propelled monomer units is followed by polymer re-swelling at high-enough activity.
- **Inertia shifts the collapse.** Rotational inertia induces large decorrelation time of velocity direction, reducing effective rotational diffusion. For large J = 16 and m = 16 collapse occurs at lower Pe.

**The loop (ring) version**
- Loops whose beads push in random directions, as here, have been simulated without inertia: R_g of unknotted polymer rings exhibits a nonmonotonic change with Pe, at N = 200 (Macromolecules 2025). A 2D loop shows the same non-monotonic size versus activity (non-monotonic behavior of the radius of gyration of an active Brownian ring).
- Loops pushed along their own backbone (a different activity rule) also collapse, and the outcome depends on length: Short rings swell, whereas longer rings collapse, at sufficiently high activity (Locatelli, Bianco & Malgaretti, PRL 126, 097801, 2021).
- So "loop + inertia" combines known pieces. I found no paper on exactly this combination, but that is a gap in the rule set, not in the behaviour.

**The twist I couldn't find documented**
- **Exponent ≈ 0 at J=10.** A loop whose size doesn't grow with bead number isn't reported anywhere I looked. It also doesn't fit physics: a dense collapsed blob should scale as about N^(1/3), and a core with loops hanging off it should be larger still. Before treating it as real, I would check:
  - whether the runs had reached steady state at N = 2048;
  - how narrow the range of bead numbers was;
  - whether only part of the loop is being measured.
- **Sharp re-opening at J=30.** This is a quantitative detail of the inertia dependence the PCCP paper already describes, so I don't count it as new behaviour.
- **Specific numbers** (squared size about 10× smaller than passive, tightest at Pe≈17, re-swelling above Pe≈40) are plausible for this model family, but I couldn't check them against the paper's figures.

**Genuine or artifact?** The collapse and re-swelling look genuine and have a known mechanism, so this is not a frozen or diverged state. The exponent-zero result is the only part that looks likely to be an artifact.

## Reply R29

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "activity-induced (motility-induced) collapse and re-swelling of an active Brownian polymer with translational and rotational inertia", "notes": "Linear-chain version with identical ABP+inertia (m,J) model is published; closed-ring version with J sweep, sharp J=30 window and N-independent size at J=10 not found. Full text of key paper not read, only preprint fragments."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The same model, with the same inertial ABP monomers, has been reported for a linear chain. I did not find it for a closed ring. Confidence is medium.

**1. Model family**
- This is active matter: a polymer made of underdamped active Brownian particles (ABPs, or "active Langevin" particles). Each bead has translational mass and rotational inertia.
- Beads repel only, so any compaction is motility-induced. It is the polymer analogue of motility-induced phase separation (MIPS).

**2. Closest literature**
- **Earliest and closest source: Buglakov, Lelecova and Chertovich, *Phys. Chem. Chem. Phys.* 27, 1691–1700 (2025), "Motility-induced collapse of an active Brownian particle polymer chain."** The preprint is on ChemRxiv as "Monomer activity solely induces collapse of polymer chain in a good solvent."
  - It shows self-propelled monomers collapsing a chain with no explicit attraction between particles. It also reports that the effect resembles MIPS.
  - The preprint uses the same model as yours. It considers a chain of active Brownian particles in an underdamped regime, with mass m and moment of inertia J.
  - It also reports that collapse is followed by polymer re-swelling at high enough activity. It gives an Rg minimum at Pe = 15, close to your Pe of about 17.
  - On rotational inertia, it says rotational inertia induces a large decorrelation time of the velocity direction, reducing effective rotational diffusion.
  - I could read only the preprint's fragments, not the full paper. I did not verify the J-dependence details or the 10× size reduction.
- **Review: Physics of active polymers, *Rep. Prog. Phys.* (2024), "Non-equilibrium structural and dynamic behaviors of active polymers in complex and crowded environments" (IOPscience, doi 10.1088/1361-6633/ad3e11).**
  - It says the structural and dynamical properties of active Brownian ring polymers are in general similar to those of linear chains.
  - It does not cover inertial rings.
- **Ring-specific sources:**
  - Kumar and Thakur, *Macromolecules* 56, 5229 (2023). They find polar activity tends to swell the shorter rings and crumple the larger ones. This is overdamped, and I did not check whether the activity is tangential.
  - Locally relevant ring work includes the knotted active-ring papers (*Macromolecules* 2025). They show swelling-then-collapse and a shrinkage dip for rings, but these are overdamped, tangentially driven, Brownian-dynamics models.
- **Background on inertial ABPs:** Löwen, *J. Chem. Phys.* 152, 040901 (2020), "Inertial effects of self-propelled particles: From active Brownian to active Langevin motion." It gives the result that MIPS is shifted to higher Péclet numbers and finally destroyed as particle mass increases. This is context for why M and J matter.
- **Reswelling at high activity** is also a general feature. Examples are a polymer in an ABP bath, which collapses at moderate activity and reswells at high activity, and ABP filaments with hydrodynamic interactions.

**3. What differs from your system, and whether it matters**
- **Topology (ring versus linear).** The ring closure probably matters little for the existence of the collapse and reswelling. The review says rings are similar to linear chains for ABP drive. Ring topology could shift the numbers, such as the collapse threshold or the size scaling, because it forbids free ends.
- **Open points I could not confirm:**
  - the sharp J = 30 re-opening window;
  - the size independent of bead number at J = 10 (a compact, bead-number-independent loop is not described for linear chains in the fragments I read);
  - the 10× reduction;
  - whether the condensed state is a single-loop globule or has dense and stretched regions. For linear chains the preprint reports separation within the collapsed state into dense and stretched loop regions.
- I did not find a paper that sweeps J on a closed ring.

**4. Confidence: medium.** The model class and the collapse-then-reswelling behaviour are clearly documented. I could not read the full text of the key paper, and the ring-specific, J-sweep and N-scaling results were not found.

**

## Reply R30

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Motility-induced collapse of an underdamped (inertial) active Brownian particle polymer; non-monotonic Rg(Pe) of active rings", "notes": "Agrees with the first-pass KNOWN-VARIANT label but corrects its note: varying rotational inertia J in an inertial self-propelled-bead chain IS documented (Buglakov et al. 2025, linear chains; collapse at lower Pe for large J and m, MIPS-like core-plus-loops state). Remaining differences: (1) closed-loop topology (minor); (2) loop size independent of bead number (exponent ~0) at J=10, which is not documented and is physically implausible for a non-crossing chain at finite density (a globule gives ~1/3). This most likely comes from beads overlapping or the chain passing through itself at high inertia/Pe, unconverged runs, or a narrow fitting range. Verify with minimum non-neighbour bead distances, a check that the loop's knot type is unchanged, time-step sensitivity and convergence before treating it as real. The collapse and reswelling themselves look genuine, not trivial."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass's label, but one of its findings is wrong. A published study already varies rotational inertia J in this kind of model, though for an open chain rather than a closed loop. The size-vs-length exponent of about 0 at J=10 is the only feature I couldn't find in the literature, and it looks more like a numerical artifact than new physics.

**What is documented**
- **The closest match.** Buglakov, Lelecova and Chertovich (PCCP 2025; preprint on ChemRxiv 2024) simulate an open chain of self-propelled beads with both translational and rotational inertia. They consider a chain of active brownian particles in an underdamped regime with m - mass of the particles and J - moment of inertia. They show the ability of self-propelled monomer units to cause chain collapse under conditions where there is no explicit attraction between particles. They also cover how inertia shifts the collapse: rotational inertia induces large decorrelation time of velocity direction, reducing effective rotational diffusion, and for large J and m, collapse happens at lower Pe. The collapsed state is heterogeneous, characterized by a dense condensed core surrounded by elongated loops. They link it to motility-induced phase separation. So the first pass was wrong that no source varies J; that has been done for linear chains.
- **Why more rotational inertia helps collapse.** A separate study found that the increase of the moment of inertia of each particle favors non-equilibrium phase coexistence, known as motility induced phase separation.
- **Shrink-then-reswell with activity in rings.** For a 2D ring of active Brownian beads, We observe a non-monotonic behavior of the radius of gyration of an active Brownian ring as a function of active force. In 3D, polar activity is reported to swell the shorter rings and crumple the larger ones. A 2025 study of knotted active rings builds on an earlier swelling-collapse transition, though there the activity pushes along the chain's backbone rather than along each bead's own heading. Non-monotonic size with activity also appears in partially active rings.

**What I could not find**
1. The same collapse model with inertia shown specifically for a closed loop. This only swaps the topology, and I'd count it as a minor difference.
2. Loop size that doesn't grow with bead number (exponent ≈ 0 at J=10). I'd treat this as suspect, not as a discovery:
   - A compact globule should give an exponent of about 1/3. The documented core-plus-loops state would not give 0 either.
   - An exponent of 0 out to 2048 beads would mean the loop's size doesn't depend on its length at all. That can't hold for a chain that truly can't cross itself and has finite density.
   - Likely causes: beads overlapping or the chain passing through itself at high speed (large J gives long straight runs at high Pe, with soft contact repulsion and a finite time step); runs that haven't settled for the longest chains; or a fit over too narrow a range of chain lengths.
   - Before taking it seriously, check minimum distances between non-neighbouring beads, whether the loop's knot type stays the same, sensitivity to the time step, and convergence over time.

**Is it a real pattern or an artifact?** The collapse by roughly 10× in squared size, the most compact point around Pe≈17 and the reswelling above Pe≈40 all look like a genuine, documented effect (collapse of a self-propelled chain driven by motility-induced phase separation). The exponent of about 0 is probably an artifact or an under-sampled fit.

## Reply R48

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Motility-induced collapse and re-swelling of an active Brownian particle polymer with translational and rotational inertia (here with ring topology)", "notes": "Collapse then re-swelling with activity in an inertial ABP chain is documented for linear chains with parameters close to these (Pe_min ~15 vs observed ~17). Twists: ring topology (likely minor); sharp J=30 re-opening (inertia shifting the collapse is documented qualitatively); size-vs-N exponent ~0 at J=10, which is physically implausible for a chain whose beads cannot overlap (a dense globule should give at least ~1/3). The exponent most likely reflects a narrow N range, runs not at steady state, or Rg dominated by loops; treat it as a possible artifact rather than a discovery. PCCP full text was not opened; Buglakov numbers come from preprint snippets."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass on the behaviour, but I'd tighten it. The main behaviour (shrinking, then swelling back open as activity rises) is already documented in a model that differs from yours only in using a linear chain instead of a ring. The size-vs-bead-number exponent of about 0 is not documented, and I think it's more likely a measurement problem than a real effect.

**Where the main behaviour is documented**
- **Buglakov, Lelecova & Chertovich, PCCP 27, 1691 (2025).** They show that self-propelled monomer units cause chain collapse under conditions where there is no explicit attraction between particles. The collapsed state is heterogeneous, characterized by a dense condensed core surrounded by elongated loops, and the authors show that this effect is similar to motility-induced phase separation.
- **It is the same kind of model as yours.** The preprint models a chain of active brownian particles in an underdamped regime with m - mass of the particles and J - moment of inertia, keeping inertia in both translational and rotational equations.
- **Your numbers fit theirs.** In their work, collapse is followed by polymer re-swelling at high-enough activity. The most compact state sits at Pe = 15 for N = 1024, and they plot Rg dependence of active chain with N = 1024 on Pe at different J and m. They also report that inertia changes where collapse happens (for large J = 16 and m = 16 collapse occurs at lower Pe). Your maximum compaction at Pe ≈ 17 and re-opening near Pe ≈ 40 for N = 2048 fall in the same range.

**Rings specifically**
- A review notes that the structural and dynamical properties of active Brownian ring polymers are in general similar to that of linear chains. That makes "ring instead of linear chain" a minor twist.
- Activity-driven collapse of long rings is also documented, though with a different kind of push: each bead driven along the chain direction rather than its own heading. There, short rings (N≲100) swell, whereas longer rings (N≳200) collapse, at sufficiently high activity (Locatelli et al., PRL 2021).
- I found no paper on a ring of inertial active Brownian beads with a sweep over rotational inertia J.

**The twists, and whether they matter**
1. **Ring instead of linear chain:** probably unimportant at N = 2048. The collapsed state (dense core plus loops) doesn't need free ends.
2. **Sharp re-opening at J = 30:** not documented. But the effect of inertia on where collapse happens is documented qualitatively, so this is most likely a parameter-dependent shift.
3. **Size not depending on bead number (exponent ≈ 0 at J = 10):** not documented, and physically suspect. Beads that can't overlap must fill space, so even a fully dense globule should grow roughly as N^(1/3). An exponent near zero suggests one of these:
   - a fit over too narrow a range of N;
   - runs that haven't reached steady state (collapsed active rings are known to get stuck);
   - the measured size being set by a few dangling loops rather than the whole ring.
   I would treat it as a likely artifact until it's checked over a wide range of N with careful equilibration.

**Real pattern or artifact?** The collapse itself (a 10× drop in squared size, then swelling back) looks real. It is not a frozen absorbing state or a numerical blow-up, and it matches the published motility-induced-collapse picture. Only the exponent-0 claim looks like a possible artifact.

**Searches I ran:** ring polymers of active Brownian particles with rotational inertia, collapse and re-swelling; Buglakov's inertial active-polymer collapse work; the PCCP abstract and the re-swelling at high Pe; underdamped active Brownian rings; Buglakov/Chertovich work on active rings; recent (2025–26) inertial active-ring collapse papers. I couldn't open the PCCP full text, so the Buglakov numbers come from the preprint snippets. I also didn't verify the Breoni et al. paper on knotted active rings (Macromolecules 2025).

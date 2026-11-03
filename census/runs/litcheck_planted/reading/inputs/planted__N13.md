# planted:N13

**Reported behavior (what the authors report as new):** A noisy voter model hit by random population-wide opinion shocks shows trimodal and asymmetric bimodal stationary opinion distributions in addition to the usual unimodal and bimodal ones.

**Rule on the card:** World: 500 agents, each holding opinion +1 or -1, every agent able to interact with every other. Rules each step: pick an agent; with probability eta (around 0.002, the inverse population; swept) it flips its opinion at random; otherwise it copies the opinion of a randomly chosen other agent. Separately, at rate r (swept), a population-wide shock occurs: opinions are instantly reassigned so that the average opinion jumps to 0.3 (65% holding +1). Measure the long-run distribution of the average opinion.

**Observed features on the card:**
  - distribution_peak_count = 3 (peaks in long-run distribution of the average opinion)
  - middle_peak_position = 0.3 (average opinion at the extra, shock-induced peak)
  - distinct_regimes = 6 (distinct long-run regimes across the noise-versus-shock-rate plane)
  - noise_threshold_without_shocks = 0.002 (noise level separating one-peak and two-peak behaviour without shocks)

## Reply R04

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Mean-field noisy voter (Kirman) model finite-size bimodal-unimodal transition (a_c ~ h/N) under stochastic resetting of the magnetization (resetting-induced NESS)", "notes": "The noise threshold of 0.002 = 1/N simply reproduces the known Kirman/noisy-voter transition. The middle peak at 0.3 is the generic resetting-induced peak at the reset point, put in directly by the shock rule rather than emergent from the interactions; it shows up when the shock rate is comparable to the voter relaxation rate (~1/N per sweep). The closest documented equivalent is the Ising model with magnetization resetting (Magoni et al. 2020), which has a rich (T,r) phase diagram. No paper found on the exact noisy-voter-plus-resetting combination. The only twist is that the voter model's bistability is a finite-size effect with no deterministic pull, which may change regime boundaries but probably not the qualitative picture. The count of 6 regimes depends on how regimes are defined. Sources were checked through search snippets only, not full texts. Not a numerical artifact, but close to trivial."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass. I didn't find any paper on this exact model (a noisy voter model with its average opinion reset to 0.3), but each part of the behaviour is documented. The three-peak result follows directly from combining those parts and is largely expected.

**1. The no-shock baseline is a textbook result.** This is the mean-field noisy voter model, also called the Kirman herding model. On a fully connected population it is known to have a finite-size transition between a bimodal and a unimodal behavior, depending on the relative magnitude of the noise and the interaction parameters. Above the threshold (noise a > h/N), the distribution of n becomes unimodal with a peak at n = N/2. With N = 500 that threshold is 1/N = 0.002, which is exactly your reported value. That confirms the simulation is correct, not that something new was found. It is also a known finite-size effect: it disappears in the limit of an infinite number of agents.

**2. Resetting the average opinion of a two-state spin system has been studied.** Magoni, Majumdar and Schehr (Phys. Rev. Research 2, 033182, 2020) take an Ising system and use a protocol that resets the magnetization to its fixed initial value m0 at a constant rate r. In the resulting NESS, the magnetization acquires a nontrivial distribution, leading to a rich phase diagram in the (T,r) plane. That is structurally the same setup as yours, with temperature playing the role of your noise η. A general feature of resetting is a non-equilibrium steady state (NESS) with weight concentrated near the reset point (Evans, Majumdar and Schehr's 2020 review). So a peak at m = 0.3 is the expected signature of resetting, not something emerging from the opinion dynamics.

**3. Multi-peak distributions from outside interventions are also documented.** The paper "Noisy voter models in switching environments" (Phys. Rev. E 108, 044301, 2023) shows that a population of noisy voters subject to groups of external influencers develops multipeak stationary distributions.

**The twist, and why it matters little.** I found no published phase diagram for noisy voter dynamics with resetting. The difference from the Ising work is that the voter model has no deterministic pull towards either consensus; its two peaks exist only because the population is finite. That changes the exponents and the regime boundaries, but probably not what you see qualitatively. Peaks at ±1 plus a peak at the reset point is what simple reasoning predicts, as follows:
- Voter dynamics take roughly N sweeps to move the average opinion by order one.
- So when the shock rate is comparable to that relaxation rate, the system spends a visible amount of time near 0.3.
- When shocks are much rarer, the distribution goes back to the baseline; when they are much more frequent, it collapses to a single peak at 0.3.

The count of six regimes depends on how regimes are defined. It is probably just the combinations of {one or two baseline peaks} with {no, intermediate or dominant reset peak}.

**Emergent or artifact?** Not a numerical artifact, but close to trivial. The middle peak is imposed directly by the shock rule, not produced by the agents' interactions. The only part that comes from the interactions is the baseline one-peak versus two-peak transition, which is well known.

**What I searched:** noisy voter model with resetting and three peaks; voter model with stochastic resetting; Kirman model with resetting; Ising model with stochastic resetting. Caveat: I relied on search-result abstracts and snippets rather than reading the full papers. I did not check what the Ising-resetting paper's phase diagram looks like in detail, or how it compares to your six regimes.

## Reply R18

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Noisy voter model / Kirman model finite-size bimodal-unimodal transition (threshold ~1/N), plus stochastic resetting of the magnetization", "notes": "The no-shock part (two peaks at consensus below eta~1/N, one central peak above) is textbook. The shock-induced third peak at the reset value 0.3 and the 6-regime map were not found in the literature; the third peak is the expected result of resetting competing with slow consensus formation. A dedicated resetting-plus-voter-model paper may exist that I did not find. The 6-regime count was not checked."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The base model and its transition are well documented. I found no source that adds a population-wide reset to a fixed bias of 0.3. The extra peak at 0.3 looks like a direct consequence of combining two known ingredients.

**1. Model family.** This is the noisy voter model (opinion dynamics and statistical physics), on a complete graph with N = 500. It is the same as the Kirman ant/herding model from economics. The noisy voter model is the network-embedded equivalent of the Kirman model in its original formulation. The shocks are a stochastic-resetting protocol: at rate r the whole state is reset to a fixed magnetization of 0.3.

**2. The behaviour, in two parts.**
- **Part A: two peaks without shocks, and the noise threshold.**
  - This is known. In the fully-connected case there is a finite-size transition between bimodal and unimodal behaviour, depending on the relative size of the noise and the interaction.
  - Bimodal means peaks near ±1 (consensus). Unimodal means a single peak near 0.
  - The threshold scales as 1/N. A related paper gives a bimodal phase for p<1/(N+1), a plain phase at p=1/(N+1), and a unimodal phase for p>1/(N+1).
  - Your threshold of η ≈ 0.002 = 1/500 matches this. The exact prefactor depends on conventions. In your rule a "noise" event picks a random sign, so the effective flip rate is η/2.
- **Part B: a third peak at 0.3 from resets.**
  - I found no paper on the voter or Kirman model with resets to a fixed bias.
  - In the Kirman literature, shocks appear only as the alternative, exogenous explanation for regime switching, not as a modelled ingredient. Kirman said regime switches in other models are the result of exogenous shocks, while in his model they are endogenous.
  - Stochastic resetting is a large field, started by Evans and Majumdar: Diffusion with Stochastic Resetting, Phys. Rev. Lett. 106, 160601 (2011). It is known that resetting breaks detailed balance and drives nonequilibrium steady states.
  - A reset-induced peak at the reset point is the generic outcome when resetting is faster than relaxation. Below η ≈ 1/N, the relaxation toward consensus takes about N steps. I expect the mix of peaks at −1, +0.3 and +1 to change with r. I expect the six regimes in the noise-versus-r plane to be crossovers among the consensus peaks, the reset peak and the central peak. I did not verify this.

**3. Citations checked against search results.**
- Kirman 1993, "Ants, Rationality, and Recruitment", Quarterly Journal of Economics 108: 137–156. This is the earliest clear source for the model. Its bibliographic details are confirmed (result 26-1).
- Peralta, Carro, San Miguel and Toral 2016, "The noisy voter model on complex networks", Scientific Reports (DOI 10.1038/srep24775). It serves as a review of the bimodal–unimodal transition. I did not confirm the author list or year in the search results; the author names and 2016 come from my memory.
- Evans and Majumdar 2011, Phys. Rev. Lett. 106, 160601. This is the founding resetting paper, confirmed in a reference list (result 34-7).
- Liggett-style duality results for the noisy voter model (Granovsky and Madras 1995, Stochastic Processes and their Applications) also appear, but only as background.

**4.

## Reply R23

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Mean-field noisy voter (Kirman herding) model with finite-size bimodal-unimodal transition at eta~1/N, under Poissonian stochastic resetting of the magnetization", "notes": "No-shock part is textbook (threshold eta=1/N=0.002 matches). The third peak at the reset value 0.3 follows directly from the renewal formula (resetting cusp/peak at the reset point); the closest documented analogue is the Ising model with resetting, which has a rich (T,r) phase diagram. I found no voter/Kirman-specific resetting paper, and the 6-regime count is unverified and probably depends on how peaks were detected (combinations of consensus, centre and reset peaks). The pattern is real, but the middle peak is imposed by the shocks rather than produced by the agents' interactions."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass, though I'd describe the twist more narrowly. The behaviour is a documented model plus a documented general mechanism. The specific combination, a noisy voter model with the average opinion reset to a value other than zero, has no paper I could find. Its main outcome, an extra peak at the reset value, follows directly from standard resetting theory.

**Part 1: without shocks (textbook).** This is the noisy voter model, which is the same thing as Kirman's herding model. Economics and chemistry papers do not refer to the voter model and use terms such as "catalytic reaction", "herding" or "Kirman model" instead. For the plain two-opinion version, the only transition is trivial: with no noise the system ends in consensus, and with any noise it keeps moving between states. With 500 agents this appears as a size-dependent switch between one peak and two peaks. In the usual scaling, two peaks appear when noise×N/copy rate is about 0.2, the curve is flat at 1, and there is a single central peak at 4. Your threshold of 0.002 = 1/N is exactly this critical point.

**Part 2: the shocks (generic resetting).** In standard resetting theory, the long-run distribution under resetting is the Laplace transform, with the reset rate as the variable, of the distribution without resetting. A cusp at the reset point is a robust signature of resetting. So when shocks come faster than consensus forms (consensus takes about N steps per agent), a peak at 0.3 is expected. The closest documented analogue is the Ising model with resetting. There, the magnetization is reset to a fixed value m0 at constant rate r, giving a nontrivial distribution and a rich phase diagram in the (temperature, reset rate) plane. Several regimes in the (noise, reset-rate) plane are therefore the expected kind of result.

**What I couldn't find:**
- No paper on the voter or Kirman model with resetting the average opinion to an off-centre value. An exact search for "stochastic resetting" plus "voter model" returned nothing relevant.
- No published count of six regimes. I suspect it is partly an artefact of how regimes were counted. Mixing three peak locations (the two consensus states, the centre, and the reset point) gives several combinations, and the count will depend on the peak-detection settings.

**Does the twist matter?** Only a little. The result follows directly from the renewal formula applied to the known diffusion limit. It would be a modest analytic exercise, not a new phenomenon.

**Real or artefact?** The pattern is real: the system keeps moving rather than freezing, and nothing diverged. However, the middle peak is imposed by the shocks rather than produced by the agents' interactions, so "emergent" applies only to the consensus and centre peaks.

**Searches I ran:** noisy voter model + stochastic resetting; Kirman model + resetting; voter model resetting phase diagram; "voter model" "resetting" opinion dynamics; noisy voter trimodal; Wright-Fisher + resetting; opinion dynamics + resetting; "stochastic resetting" "voter model"; Curie-Weiss + resetting. One caution: the trimodal distributions in the nonlinear, aging and partisan noisy voter papers come from different mechanisms. For example, in the noisy partisan model the long-run distribution goes from two peaks to one through a three-peak stage. Those papers should not be cited as describing this behaviour.

## Reply R33

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Mean-field noisy voter (Kirman) model finite-size bimodal-unimodal transition (noise ~ 1/N) under global stochastic resetting of magnetization (resetting-induced peak at reset point)", "notes": "Noise threshold ~1/N reproduces the known Kirman criterion (1/N vs 2/N depends on whether 'random flip' is forced or a random reassignment). The third peak at the reset value 0.3 is the generic resetting signature (renewal/Laplace-transform form of the NESS). The specific 6-regime count across the noise-shock plane is unverified and likely reflects combinations of uni/bimodal background, presence of the reset peak, and which peak dominates (shock rate vs relaxation rate ~ max(eta, 1/N)); it is not a qualitatively new phenomenon. Not an artifact: the system stays ergodic and the shock peak is real, but it is expected by construction."}

**Reply text:**

**Verdict: KNOWN-VARIANT, and I agree with the first pass.** Both parts of this behaviour are documented. What is put together here is the textbook result for stochastic resetting applied to a well-known model. I found no paper on this exact model. But the behaviour you'd predict from the published pieces matches what was observed, so this isn't a discovery.

**1. The bimodal/unimodal switch with no shocks is a known result.** This is the noisy voter model, also called the Kirman model. Papers describe it as having a finite-size phase transition, upon increasing the relative importance of the noise to the herding, form a bimodal phase where most of the agents share the same opinion to a unimodal phase where almost the same fraction of agent are in opposite states. In the same setup as yours, where one agent per step either flips randomly or copies another agent, the switch happens when noise × N ÷ copying rate is about 1 (< 1 U-shaped (bimodal) ... > 1 peaked at 1/2). With N = 500 that gives a noise level of about 0.002, which is exactly what was observed. Whether the threshold is 1/N or 2/N depends on whether "flip at random" means a forced flip or picking a random opinion (which only changes it half the time). This model is the network-embedded equivalent of the Kirman model in its original, extensive formulation.

**2. The extra peak at 0.3 is the standard effect of resetting.** The long-run distribution under resetting is the no-reset distribution averaged over the random times since the last reset. That always puts probability mass around the reset point. For simple diffusion, the distribution exhibits exponential decay away from the resetting position. In a model that is already bimodal, the reset point therefore adds a third peak at 0.3. That is the expected result, not something new.

**3. Resetting a whole population's average opinion has been studied in a close equivalent.** Magoni, Majumdar and Schehr (Phys. Rev. Research 2020) studied the Ising model with a protocol that resets the magnetization to its fixed initial value m0 at a constant rate r. In the resulting NESS, the magnetization acquires a nontrivial distribution, leading to a rich phase diagram in the (T,r) plane. That is the same structure as your two-parameter sweep, with temperature in place of noise. Gupta and Jayannavar's review of resetting in interacting systems (J. Phys. A 2023) covers this kind of whole-population resetting.

**4. A trimodal outcome isn't new for this model family either.** A noisy voter variant with agent preferences undergoes a discontinuous transition where the probability distribution of the number of agents in one state changes from unimodal to trimodal.

**The twist I couldn't find:** the detailed phase diagram for this specific case, i.e. the "6 regimes" across the noise-versus-shock-rate plane and where the boundaries fall. I searched for: noisy voter model with stochastic resetting; voter model resetting; Kirman model with resetting or shocks; Moran/Wright-Fisher models with resetting; resetting in interacting particle systems; noisy voter resetting 2025. None of these turned up a direct study. It probably doesn't matter much. The regimes are likely just combinations of three things: one or two noise peaks, whether a reset peak shows up, and which peak is tallest. That depends on the shock rate compared with how fast the system relaxes. Nobody has checked the count of six against a source, and it may depend on how peaks are detected and counted.

**Genuine or artifact?** It's a real pattern, not a broken or frozen simulation. The noise stops the population from locking into permanent consensus, and the reset peak exists because of the shocks. But it is what the setup guarantees by construction, a simple superposition of known effects, so it isn't surprising emergent behaviour.

## Reply R50

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Mean-field noisy voter (Kirman) model finite-size bimodal-unimodal transition (a_c=2/(N+2)) + stochastic resetting of the magnetization (renewal-formula NESS, cf. Ising model with stochastic resetting)", "notes": "No-shock threshold matches a_c=2/(N+2) (~0.004 per flip-to-opposite, ~0.002 for random-flip convention at N=500). Third peak at the reset value 0.3 follows directly from the renewal formula when reset rate exceeds the O(1/N) voter drift rate, with consensus peaks at +/-1 surviving while noise < ~1/N; analogous rich (T,r) phase diagram documented for the Ising model with resetting. No paper found specifically on noisy voter/Kirman + global resetting; the 6-regime count is unverified and likely depends on how regimes were classified. Middle peak is an externally imposed imprint of the reset target, not self-organised: genuine but shallow."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass, but I'm adding a closer precedent it missed. Both parts of the behaviour are documented. What I couldn't find is a paper that applies resetting to the noisy voter (Kirman) model itself.

**1. No-shock part: textbook.** The mean-field noisy voter model is the Kirman model. A network version is described as the network-embedded equivalent of the Kirman model in its original, extensive formulation. It undergoes a finite-size phase transition, upon increasing the relative importance of the noise to the herding, form a bimodal phase where most of the agents share the same opinion to a unimodal phase. The exact threshold is known: at [ac=2/(N+2)], the stationary distribution is flat. For N = 500 that is about 0.004 in rate units. Your ~0.002 is the same point once you convert between "flip to a random opinion" and "flip to the opposite opinion" (a random flip changes the opinion only half the time). This part is confirmed exactly.

**2. Shock part: same idea is documented for an equivalent system.** The closest match is Magoni, Majumdar and Schehr (2020), on a mean-field-type magnetization variable that is periodically reset. Their protocol resets the magnetization to its fixed initial value m0 at a constant rate r. In the resulting NESS, the magnetization acquires a nontrivial distribution, leading to a rich phase diagram in the (T,r) plane. The general theory gives the long-run distribution with resetting directly from a standard formula (the "renewal formula"): The stationary distribution for the resetting problem is again related to the Laplace transform of the propagator in the absence of resetting. In plain terms, the long-run distribution is the no-shock evolution starting from 0.3, averaged over how long it has been since the last shock. Your setup is that construction with a 2-state noisy voter model in place of the Ising model.

So the three peaks are expected. Shocks keep refilling the population at 0.3. If shocks come faster than the population drifts towards consensus (that drift takes on the order of N steps per agent), mass builds up at 0.3. When the noise is below ~1/N, the population still reaches consensus between shocks, so the ±1 peaks remain. Similarly, a noise level × shock rate plane splitting into several regions is the same kind of result as Magoni et al.'s rich (T, r) phase diagram.

**Twist I couldn't find described, and whether it matters:**
- I found no paper on the noisy voter or Kirman model with resetting. My searches for voter-model resetting returned nothing, or only unrelated "resetting" of an internal confidence or age variable.
- I also didn't find the specific three-peak picture: a peak at the 0.3 reset value alongside the two consensus peaks.
- This likely matters little. It's a direct, predictable consequence of the renewal formula, not a new mechanism.
- I haven't verified the count of 6 regimes. It probably comes from crossing about 2 no-shock shapes (one peak / two peaks) with about 3 shock levels (none / intermediate / dominant). That makes it partly a classification choice, not a robust fact.

**Genuine or artifact:** It isn't a numerical artifact, frozen state or random colouring. But the middle peak is imposed from outside: it's an imprint of the reset target, not something the agents organise themselves into. I'd call it a genuine but shallow pattern.

**

## Reply R53

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Noisy voter model / Kirman ant model, finite-size bimodal-unimodal transition at noise ~ 1/N; shocks = stochastic resetting", "notes": "No-shock threshold ~1/N is well known (factor-of-2 convention differences). No paper found on noisy voter model with population-wide resets to fixed magnetization; the shock-induced middle peak is mechanistically expected, but the 6-regime noise-vs-shock-rate phase diagram could not be checked against any source. Absence in a limited search is weak evidence of novelty. The Evans-Majumdar-Schehr resetting review was not verified by search."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The noise-driven part is well established. The population-wide shock that resets the average opinion to 0.3 is the part I could not find in the literature.

**1. Model family.** This is the noisy voter model on a complete graph. In economics it is the Kirman "ant" or herding model. It belongs to opinion dynamics and statistical physics (mean-field, finite-size transitions). The shocks are a stochastic-resetting protocol added on top: at rate r the state is reset to a fixed magnetisation of 0.3.

**2. Evidence from the searches.**
- **No shocks (standard result).** Carro, Toral and San Miguel state that, in the fully-connected case, the system is characterized by the existence of a finite-size transition between a bimodal and a unimodal behavior. They add that the noisy voter model becomes the network-embedded equivalent of the Kirman model.
- **Where the threshold sits.** Artime et al. (2019) give the flat-distribution point as ac = 2/(N+2) at the precise value where the exponent is 0. That is about 1/N up to a convention-dependent factor of order 1. Your threshold of 0.002 = 1/500 is therefore consistent with this literature. In Banisch's contrarian version the threshold is p* = 1/(N+1), with a bimodal distribution of opinions below it and a unimodal one above.
- **Earliest sources.**
  - Kirman, "Ants, rationality, and recruitment", Quarterly Journal of Economics 108, 137–156 (1993). I verified the reference only through citing papers, not the paper itself.
  - Granovsky and Madras, "The noisy voter model", Stochastic Processes and their Applications 55 (1995), pp. 23–43. The ideas.repec.org page confirms the title, journal and pages.
- **Review-style source.** Carro, Toral and San Miguel, "The noisy voter model on complex networks", Scientific Reports 6, 24775 (2016), arXiv:1602.06935. Artime et al., "Herding and idiosyncratic choices: Nonlinearity and aging-induced transitions in the noisy voter model", Comptes Rendus Physique 20, 262–274 (2019), also covers it.
- **Three-peak distributions.** Trimodal and "W"/"M" distributions are known in noisy-voter generalisations. One example is "Polarization-induced stress in the noisy voter model" (Physica A, 2024, ScienceDirect S0378437124003492). In that model for α>5, an additional (trimodal) phase emerges. These arise from nonlinear internal mechanisms, not from external resets.
- **With shocks.** I found no paper on a noisy or plain voter model with stochastic resetting or population-wide shocks to a fixed magnetisation. The closest hits were external-field or driven variants: the oscillating-propaganda contrarian voter model (Entropy 24, 1140, 2022) and Kirman with external information (Carro et al., PLoS ONE 2015, "Markets, Herding and Response to External Information"). Neither is a reset protocol.

**3. What differs and whether it matters.** The shock is a Poisson reset to m = 0.3, not a smooth field or a nonlinear internal term. Mechanistically, a third peak at the reset value is the expected outcome. Shocks that are frequent compared with the slow consensus-drift time (~N steps) should leave a peak at 0.3. Rare shocks should leave the two noise-driven peaks at the extremes. Intermediate rates give mixed shapes, which is plausibly where the six regimes come from. That reading is my own inference, not a sourced result. I could not verify the "6 regimes" count or the phase diagram against any source. Stochastic resetting is a large field (Evans, Majumdar and Schehr, 2020), but I did not verify that citation here, so treat it as unconfirmed. The qualitative result is probably not new, but the specific noise-versus-shock-rate phase diagram might not have been published.

**4. Confidence: medium.** The no-shock threshold is solidly sourced and matches your 1/N value. The shock part rests on about eight searches that found nothing, and absence in a search engine is weak evidence.

**

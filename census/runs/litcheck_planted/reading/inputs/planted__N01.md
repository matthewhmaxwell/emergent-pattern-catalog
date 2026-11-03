# planted:N01

**Reported behavior (what the authors report as new):** At linear charging, all-to-all pulse-coupled oscillators end fully synchronised with probability set by the prime factorisation of N (1-1/N^N for prime N; failures ~N^-(m-1) with m the least prime factor), and only equal-size cluster states survive.

**Rule on the card:** World: N identical oscillators, each coupled to every other (no space). Each has a voltage between 0 and 1, drawn uniformly at random at the start. Rules: (1) all voltages rise at the same constant speed (straight-line charging, no leak); (2) when a voltage reaches 1 the oscillator fires and resets to 0, together with any oscillators already firing with it as a group; (3) a firing group of k oscillators kicks every other voltage up by k/N; (4) any oscillator pushed to 1 or above by a kick joins the firing group and fires with it from then on. Run until the grouping stops changing and record whether all N end as one group. Sizes: N from 2 to about 100, many random starts per N (exact enumeration for N up to 8).

**Observed features on the card:**
  - full_sync_probability_n4 = 0.887 (share of random starts ending as one group, N=4)
  - full_sync_probability_n5 = 0.99968 (share of random starts ending as one group, N=5)
  - full_sync_probability_n6 = 0.906 (share of random starts ending as one group, N=6)
  - distinct_final_states_n100 = 9 (number of possible final groupings at N=100, all equal-sized)
  - failure_decay_exponent_even_n = 1 (failure probability falls as N to minus this, even N)
  - failure_prefactor_even_n = 0.6 (prefactor of that failure power law for even N)

## Reply R14

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "Finite-N endgame statistics not found: end states restricted to equal groups whose size divides N; prime N almost always fully synchronizes; full-sync probability jumps around with N; even-N failure falls as about 0.6/N. These look like the endpoint of the documented pairing process rather than a new mechanism.", "standard_name": "Peskin / Mirollo-Strogatz pulse-coupled integrate-and-fire oscillators, linear (non-leaky) charging with pulse 1/N: cluster aggregation via periodic size-doubling into equally spaced equal-sized clusters (O'Keeffe, Krapivsky & Strogatz 2015, Supplemental Material)", "notes": "The exact model (linear charging, pulse 1/N, all-to-all) is analysed in the OKS 2015 Supplemental Material for large N, where clusters double in size and form equally spaced, equal-sized groups. The authors leave out the endgame. Pulse = 1/N is a knife-edge where kicks can land exactly on threshold, so the exact probabilities (e.g. N=5 giving 0.99968 rather than 1) may depend on tie-breaking or floating point; check with exact arithmetic. Could not read the OKS Supplemental Material or the Senn-Urbanczik conditions."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass, but for a more specific reason: this exact model, including its main clustering behaviour, is already in the literature.

**This exact system is documented.** In the paper by O'Keeffe, Krapivsky and Strogatz, the authors say in the Supplemental Material [26], we apply aggregation theory to the deterministic Peskin model [16], assuming the oscillators rise linearly to threshold and fire pulses of size 1/N, where N ≫ 1 is the number of oscillators. That is your rule set: straight-line charging, no leak, and a group of k firing together giving a k/N kick. They state the result in the main text: the oscillators condense into clusters whose size doubles periodically. At the end of the first period, all oscillators synchronize into pairs spaced equally apart. Groups of 4 follow, and so on.

So the core behaviour is documented: equal-sized groups, evenly spaced, built up by pairing. Your odd/even differences and equal-sized end states fit naturally with that pairing picture.

**Non-sync is not a contradiction.** A follow-up paper says Mirollo and Strogatz guaranteed sync only for a leak γ > 0 and pulse size > 1/N. For other values of γ and other pulse sizes, full sync is possible, but not certain to occur. Your case (no leak, pulse exactly 1/N) is outside the theorem. One possible conflict: Senn & Urbanczik proved a result titled "Similar non-leaky integrate-and-fire neurons with instantaneous couplings always synchronize". That only holds under specific conditions, and I could not read their paper to check whether the exact 1/N case is excluded.

**What I could not find (the twist):** I found nothing on the small-N endgame:
- that the end states are equal groups whose size divides N (9 end states at N=100, which has 9 divisors);
- that prime N almost always ends in full sync;
- that full-sync probability jumps around with N (0.887, 0.99968, 0.906 for N=4, 5, 6);
- that failure falls as about 0.6/N for even N.

The O'Keeffe group says outright that it models only the early and middle stages, where N ≫ 1, as opposed to the end game. So the finite-N statistics are a real gap. They are a counting result about where the known pairing process stops, not a new kind of behaviour. I could not read the Supplemental Material, which may say more.

**Genuine pattern or artifact?** The grouping is genuine and matches the published mechanism. But pulse = 1/N is a knife-edge. When m equal groups each kick by 1/m, voltages can land exactly on the threshold of 1, so the "≥ 1 joins" rule and floating-point ties may decide some outcomes. Prime N=5 giving 0.99968 rather than exactly 1 points that way. The exact numbers (the 0.6 prefactor, the sync probabilities) need checking with exact arithmetic before anyone relies on them.

**

## Reply R17

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Peskin model with linear (non-leaky) charging and pulse 1/N; Mirollo-Strogatz absorption dynamics; cluster aggregation / size-doubling", "notes": "Could not read the PRL supplement or Senn-Urbanczik. Finite-N divisor structure of final groupings, N-factorisation-dependent failure probabilities, and the ~0.6/N even-N failure law were not found in any source; they may be new or only in unread supplements. Large-N doubling clusters (O'Keeffe 2015) is consistent with equal-sized final groups."}

**Reply text:**

**Verdict: KNOWN-VARIANT, confidence low-to-medium.** The model and its large-N cluster behaviour are described in the literature. I did not find the finite-N result you observed: final groupings that are equal-sized and tied to the divisors of N, with full-sync probability depending on N's factorisation.

**1. Model family**
- This is a pulse-coupled, all-to-all, identical integrate-and-fire oscillator system.
- It is the Peskin cardiac-pacemaker model, as analysed by Mirollo & Strogatz.
- The specific parameters are a linear charging curve (no leak), pulse size k/N, and absorption.
- Your rules are the degenerate boundary case of Mirollo–Strogatz. Their theorem assumes a monotonic, concave-down charging curve. They assume the curve rises monotonically and concave down, and they prove sync for all N and almost all initial conditions. A linear curve with pulse exactly 1/N sits outside that proof.
- The O'Keeffe 2016 paper gives the same boundary. It states that full sync is guaranteed for γ > 0 and pulse size > 1/N. The case γ = 0 with pulse exactly 1/N is not covered.
- Your failures at N=4, 5 and 6 are consistent with this boundary status. The 1990 theorem should not be read as contradicting them.

**2. Closest sources**
- **Mirollo & Strogatz 1990**, "Synchronization of pulse-coupled biological oscillators", SIAM J. Appl. Math. 50(6):1645–1662. This is the foundational paper. It describes how absorptions create groups that fire in unison with a combined pulse strength proportional to group size. I confirmed the title, authors, year and venue from the search results.
- **O'Keeffe, Krapivsky & Strogatz 2015**, "Synchronization as Aggregation: Cluster Kinetics of Pulse-Coupled Oscillators", Phys. Rev. Lett. 115, 064101. This is the closest match.
  - It treats the deterministic Peskin model with a linear charging curve and pulse ε = 1/N. It reports that the oscillators condense into clusters whose size doubles periodically.
  - Pairs form first, then merge into clusters of size 4, and so on.
  - That result is in the Supplemental Material, which I could not read. It is a mean-field result for N ≫ 1.
  - This fits your "equal-sized final groups" and your slow even-N failure pattern. It does not give finite-N divisor counts, failure probabilities or prefactors.
- **O'Keeffe 2016**, "Transient dynamics of pulse-coupled oscillators with nonlinear charging curves", Phys. Rev. E 93, 032203. It describes synchrony developing through aggregation and generalises the earlier linear-charging results to nonlinear charging curves. It is a follow-up and the nearest thing to a review.
- **Related, not read:** Senn & Urbanczik 2001, "Similar non-leaky integrate-and-fire neurons with instantaneous couplings always synchronize", SIAM J. Appl. Math. 61(4):1143–1155. The title suggests a non-leaky case. I could only see citations of it, not the paper. Its pulse and coupling assumptions may differ from yours, so I did not rely on it.

**3. What is different or unverified**
- The literature gives the large-N, mean-field picture of cluster-size doubling.
- I found no source reporting the finite-N, exact-enumeration outcome. That outcome includes the "number of possible final groupings = number of divisors of N", which gives 9 for N=100.
- It also includes the prime-N result that the only failure is all N staying separate, and the even-N failure law of roughly 0.6/N.
- The number-theoretic structure is plausibly new, or at least buried in the unread supplement. The mean-field doubling picture likely reflects the same mechanism, but I did not verify that.

**Confidence**
- It is not high because I could not read the PRL supplement or Senn–Urbanczik.
- Search coverage of finite-N, number-theoretic statements is weak.

**

## Reply R22

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Peskin / Mirollo-Strogatz pulse-coupled integrate-and-fire oscillators at marginal coupling (linear charging, pulse 1/N); aggregation/scrambler framework", "notes": "Known: the model, that full sync can fail at pulse 1/N with linear charging, and equal-sized cluster states. Not found: final states as m equal clusters for each divisor m of N (9 at N=100), prime-N near-certain sync with only the all-separate state as alternative, and the ~0.6/N even-N failure law. Equal sizes are likely forced by the neutral dynamics of linear charging, so this is an arithmetic refinement, not a new mechanism. The paywalled PRL supplement was not read and may contain it. Check tie handling at exactly 1 for large N."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass on the model family. I don't overturn it, but I'd describe the new part more narrowly and more carefully.

**What is documented**
- **The model itself.** O'Keeffe, Krapivsky and Strogatz (PRL 115, 064101, 2015) say their supplement applies aggregation theory to the deterministic Peskin model, assuming the oscillators rise linearly to threshold and fire pulses of size 1/N, where N ≫ 1. In their model, as in this one, clusters never desynchronize, and a firing cluster absorbs oscillators within a distance proportional to its size. Your rules match this exactly.
- **Full sync can fail at exactly this pulse size with linear charging.** The PRE 93, 032203 manuscript (I still couldn't confirm its author list) restates the earlier result: full sync is guaranteed for γ > 0 and pulse size > 1/N, and for other values of γ and other pulse sizes, full sync is possible, but not certain to occur. Your case (no leak, pulse exactly 1/N) is in that "not certain" region. It also says that in such cases full synchrony is not guaranteed to occur, and so multiple, stable clusters are possible.
- **Final states made of equal-sized clusters are a known type of attractor** in identical pulse-coupled oscillators. A related leaky integrate-and-fire paper reports splay states of clusters with equal sizes greater than 1 and uses the neutrality of the dynamics at zero coupling (PRE 95, 022207). Bottani (PRE 54, 2334, 1996, sole author Samuele Bottani) and Mauroy–Sepulchre (Chaos 18, 037122) also describe convergence to a phase-locked clustered behavior in excitatory integrate-and-fire populations.

**The twist I couldn't find described**
- Small-N, number-theoretic results: the possible final states are exactly m equal clusters of size N/m for each divisor m of N. That gives 9 states at N=100, which is the number of divisors of 100. For prime N, the only alternative to full sync is N separate oscillators, which almost never happens. This explains the jump from 0.887 (N=4) to 0.99968 (N=5) and back down at N=6. For even N, failure falls off as about 0.6/N.
- The PRL papers focus on large-N cluster kinetics, and I couldn't read the paywalled supplement. So this twist may be in there.
- **Does it matter?** Mildly. It's a finite-size, arithmetic refinement of a known failure mode, not a new mechanism. The divisor structure follows from the setup: with straight-line charging, nothing drifts back apart, so the only lasting non-sync states are ones where every cluster gets the same total kick per cycle, which forces equal cluster sizes. The exact even-N law may be new as a number.

**Genuine or artifact?** I think it's genuine. Groupings with unequal clusters keep merging, and the leftover equal-cluster states are stable periodic states, not frozen or diverged ones. One caveat: at pulse exactly 1/N, how ties at exactly 1 are handled ("1 or above") could affect the rarest outcomes. The exact enumeration up to N=8 makes this less of a worry, but floating-point ties at N=100 should be checked.

**Searches I ran:** the O'Keeffe–Krapivsky–Strogatz aggregation work and its supplement; Peskin model with linear charging, pulse 1/N and divisor/equal-cluster final states; non-leaky integrate-and-fire sync probability vs prime/even N; Bottani 1996; Mauroy–Sepulchre clustering; neutral equal-cluster splay states.

## Reply R26

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Peskin / Mirollo-Strogatz pulse-coupled integrate-and-fire oscillators at marginal coupling (linear charging, pulse j/N, absorption); aggregation/cluster-kinetics framework of O'Keeffe-Krapivsky-Strogatz", "notes": "Undocumented twist: final states are exactly m equal clusters for each divisor m of N (9 states at N=100), near-certain sync for prime N, and roughly 0.6/N failure for even N. The divisor rule follows from period matching: without leak each cluster of size k must charge k/N per cycle, so all clusters are equal-sized and their number divides N. The twist is real but minor; only the 0.6/N prefactor needs real calculation. Neutral (no-leak) dynamics keeps leftover multi-cluster states from merging. Caution: exact-threshold ties (kicks are multiples of 1/N and rule 4 includes equality) may make results sensitive to floating point; check the N=5 deficit (0.99968) with exact arithmetic. PRE 93 032203 is sole-authored by O'Keeffe. The PRL supplement was not read and may contain the finite-N details."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass, but I'd narrow its claim. The model and the main behaviour (marginal coupling, full sync not guaranteed, several stable clusters can survive) are documented. What I couldn't find in any readable source is the exact counting of final states for finite N.

**What is documented**
- **The model itself.** The O'Keeffe–Krapivsky–Strogatz paper (PRL 115, 064101, 2015) says that its Supplemental Material applies aggregation theory to the deterministic Peskin model, assuming the oscillators rise linearly to threshold and fire pulses of size 1/N, where N ≫1 is the number of oscillators. That is your system. I couldn't open the supplement, so it may contain more.
- **Cluster growth and firing rules.** In O'Keeffe's follow-up (PRE 93, 032203, 2016 — the search results list him as sole author, which settles one of the first pass's open questions) clusters can only ever increase in size. They can never break apart because (a) the oscillators are identical. An oscillator pushed to threshold by a pulse doesn't fire its own pulse; it joins the firing cluster, which then fires a pulse of total size j/N (so that the pulse has total size *j/N*). That matches your rules 3 and 4.
- **Sync can fail at marginal coupling.** The same paper says Mirollo and Strogatz showed that for γ > 0 and pulse size > 1/N, then full sync is guaranteed for all IC except for a set of measure ... For other values of γ and other pulse sizes, full sync is possible, but not certain to occur. Your case (no leak, pulse exactly 1/N) is one of those "other" cases, and it adds that multiple, stable clusters are possible.
- **Equal-sized cluster states.** These are known in a related integrate-and-fire model (PRE 95, 022207), which reports splay states of clusters with equal sizes greater than 1.

**The twist I couldn't find**
- Exactly one possible final state for each divisor m of N: m equal groups of size N/m. N=100 has 9 divisors, which matches your 9.
- For prime N, sync is almost certain; for even N, failures fall off as about 0.6/N.

This twist matters less than it looks, because the divisor rule follows from a short argument. With no leak, every group has the same period T. In each cycle, a group of size k gets kicks totalling (N−k)/N from the others, so it must charge the remaining k/N on its own. So T = k/N for every group, all groups are the same size, and the number of groups must divide N. Prime N therefore leaves only full sync or N separate oscillators, and even N always allows the two-group state. For two groups there is a continuous range of non-merging layouts (the second group's voltage anywhere in (0, ½) when the first fires), which explains why failure is common for even N. The ~0.6/N prefactor and its exact value are the only parts that would need real calculation, and they look minor.

**Genuine or artifact:** genuine. The non-synced outcomes cover a positive share of random starts and follow from the model's structure. Without leak, the spacing between groups stays fixed between merges, so leftover multi-group states persist rather than being pulled back towards sync. Two cautions:
1. Rule 4 triggers at exactly "1 or above", and kicks are exact multiples of 1/N, so exact ties can occur. Floating-point rounding could turn some borderline merges into non-merges or the reverse. The N=5 value of 0.99968, slightly short of 1, deserves a check with exact arithmetic.
2. The first pass's report of an exact-enumeration value isn't directly comparable to sampled ones without knowing the method.

**

## Reply R32

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Peskin / Mirollo-Strogatz model at marginal pulse strength 1/N with linear charging (scrambler/aggregation framework of O'Keeffe-Krapivsky-Strogatz)", "notes": "Marginal-case failure of full sync is documented (pulse size 1/N, linear charging). The specific divisor-structured final states, the prime-N pattern and the ~0.6/N even-N failure law were not found in any readable source. The PRL supplement was not read and could contain them. PRE 93, 032203 author list and PRE 54, 2334 authors not confirmed from search results."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The model is a well-studied one. I found no source that describes the specific outcome, which is the divisor-structured final groupings and their N-dependent failure rates.

**1. Model family.** This is a pulse-coupled integrate-and-fire oscillator model with all-to-all coupling, which belongs to synchronization and oscillator dynamics. It is the Peskin / Mirollo–Strogatz model with a linear charging curve and pulse size 1/N per firing oscillator. Absorption makes groups irreversible. It is also the "deterministic Peskin model" that O'Keeffe, Krapivsky and Strogatz treat in their supplement.

**2. What the literature says.**
- Earliest source: Mirollo & Strogatz (1990, SIAM J. Appl. Math. 50, 1645–1662).
  - When an oscillator fires, it pulls the others up by a fixed amount or to threshold. For almost all initial conditions, the population ends up firing synchronously.
  - Their proof assumes the charging curve is monotonic and concave down. The oscillators "rise toward threshold with a time-course which is monotonic and concave down".
  - Your rule system has a linear curve, so it is outside the strict hypotheses. It is the boundary case.
- Review-level source: the O'Keeffe, Krapivsky & Strogatz paper (PRL 115, 064101, 2015), plus the follow-up on nonlinear charging (PRE 93, 032203, 2016).
  - The follow-up paper states that full sync is guaranteed for all initial conditions except a measure-zero set when γ > 0 and pulse size > 1/N. For other values of γ and other pulse sizes, full sync is possible but not certain.
  - Your rule system sits exactly at pulse size 1/N with γ = 0. That is the case where full synchrony is not guaranteed, which fits your observed failures.
  - Those papers focus on transients (the cluster-size distribution during aggregation), not on the final state at finite N.
- Related result: a 1996 PRE paper (Phys. Rev. E 54, 2334) on globally coupled integrate-and-fire oscillators. It finds synchronization in two forms: oscillators evolve "en bloc" in phase, or they stay out of phase but relax in stable avalanches. I did not check whether this covers the linear 1/N case.
- Fragments of the PRL text, which I could not read in full, mention that oscillators condense into clusters, and that oscillators synchronize into pairs spaced equally apart. This may be a remark about the splay/cluster endpoint. I could not confirm its context or whether it covers ε = 1/N.

**3. What differs and whether it matters.** I found no source that states any of the following:
- the final states at ε = 1/N are equal-sized groups, so group sizes divide N;
- prime N fails only by all N staying separate;
- failure falls roughly as 0.6/N for even N.

These look like arithmetic consequences of the marginal pulse size 1/N. A group of size k advances every other voltage by exactly k/N, which is why divisor structure is plausible. I think the behaviour is a corollary of the known marginal case, but I have not verified that. The 1996 PRE paper's description of avalanche states may be relevant, and I did not read it.

**4. Citation check.** I verified these against search results:
- Mirollo & Strogatz 1990, SIAM J. Appl. Math. 50(6), 1645–1662 (JSTOR).
- O'Keeffe, Krapivsky & Strogatz 2015, Phys. Rev. Lett. 115, 064101 (PubMed and APS).
- The PRE 93, 032203 paper on nonlinear charging curves (APS and PubMed). I could not confirm its author list from the results, so I have not attributed authors. The PRL page does not mention it.
- Gwon et al. 2023, Chaos 33(3), 033102 (PubMed and AIP), a percolation treatment of the same scrambler-model aggregation.
- The 1996 PRE 54, 2334 paper (APS). I could not confirm its authors from the results.

**5. Confidence: medium-low.**
- The model identification is solid, and the sync-versus-not boundary at pulse size 1/N is documented.
- I could not read the PRL supplement, which is the most likely place for the exact ε = 1/N linear analysis. I also could not check for a finite-N exact-enumeration or number-theoretic treatment.
- I ran 9 searches on the model, the marginal case, cluster and divisor structure, and equal-spaced clusters. None returned a source describing the divisor or prime pattern.
- If that pattern turns out to be absent from the supplement, the verdict could move toward LITERATURE-SILENT for the specific finite-N arithmetic.

## Reply R56

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Peskin / Mirollo-Strogatz globally pulse-coupled integrate-and-fire oscillators with linear charging and pulse eps=1/N (O'Keeffe-Krapivsky-Strogatz 'synchronization as aggregation' model), at the marginal non-concave case", "notes": "The model and how its clusters grow at large N are documented exactly. The twist I could not find in the literature is the finite-N endgame: non-sync outcomes with positive probability; neutral final states made only of equal-sized groups whose size divides N (9 states at N=100 = 9 divisors of 100); prime N almost always syncing; and even-N failure of about 0.6/N, mainly from two halves of size N/2. My own elementary calculation (not from a source) explains this: two unequal groups a,b drift by (a-b)/N per cycle, so only equal groups can coexist, and straight-line charging (concavity b=0) is the boundary where the Mirollo-Strogatz proof fails. The twist is real but disappears with any curvature or pulse > 1/N. I could not read the PRL supplement, which might cover part of this."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass. The model is a documented one. What I couldn't find described is how finite systems end up: the specific final cluster states and the failure probabilities. Those look real rather than an artifact, but they come almost directly from the fact that straight-line charging is a borderline case.

**The model is documented exactly.** O'Keeffe, Krapivsky and Strogatz (PRL 115, 064101, 2015) apply aggregation theory to the deterministic Peskin model, assuming the oscillators rise linearly to threshold and fire pulses of size 1/N, where N ≫ 1. In that model, when a cluster of size j fires it adds a voltage pulse jε to every other oscillator, which is your rules 3 and 4. Their companion paper (PRE 93, 032203, 2016) extends this to curved charging. It notes that clusters can never break apart because the oscillators are identical and all receive the same number of pulses, which explains why your groupings only ever merge. These papers cover how clusters grow in very large systems. They do not cover how small systems finish.

**Why full sync isn't guaranteed here.** The Mirollo–Strogatz proof needs a charging curve that is smooth, monotonic increasing, and concave down, and they introduce a parameter b that measures the extent to which f is concave down. Straight-line charging (b = 0) sits on the boundary where the proof stops working. The PRE paper's own text says full sync is shown for γ > 0 and pulse size > 1/N, which is exactly not your case. A separate result proves synchronization for similar but not identical straight-line oscillators (Urbanczik and Senn (2001)). So it is known that identical straight-line oscillators with pulse 1/N are a delicate case.

**The twist I couldn't find described:**
- Runs fail to fully synchronize with a positive probability at finite N.
- The only possible final states are equal-sized groups whose size divides N. Your 9 states at N=100 match the 9 divisors of 100 exactly.
- Prime N almost always synchronizes (0.99968 at N=5), while even N fails at about 0.6/N.

A quick calculation of my own (not from a source) explains all three. Take two groups of sizes a and b. Each full cycle, their relative offset shifts by (a−b)/N. Unequal groups therefore keep drifting until one absorbs the other. Equal groups don't drift at all, so each neutral final state is a set of equal groups, one for each divisor of N. For even N, the main way to fail is two halves of size N/2 that never merge. For prime N, the only alternative to full sync is all N oscillators firing separately, which needs very wide spacing and is rare. The two-oscillator version of this neutral behaviour is implied by the Mirollo–Strogatz concavity condition. I found no paper that works out the general N-group picture or the failure numbers.

**Does the twist matter? Only somewhat.** It is a correct, number-theoretic description of a borderline case. However, it relies on the charging being exactly straight and the pulse being exactly 1/N. Any curvature (b > 0) or any pulse larger than 1/N brings back guaranteed sync, so the effect disappears under small changes. Two caveats: I couldn't read the PRL supplement, which might mention some of this, and the 0.6 prefactor and the prime/even difference are not confirmed by any source.

**Genuine or artifact?** Genuine. These are real, neutrally stable repeating states of the exact rules, not numerical blow-up or random colouring. One thing to check: the rule that 1 *or above* joins the firing group matters for perfectly evenly spaced starts, but those occur with probability zero.

**

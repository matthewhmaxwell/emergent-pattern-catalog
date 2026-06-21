# -*- coding: utf-8 -*-
"""Educational content for the EPC gallery: a general intro to emergence, and a
per-pattern 'how it emerges' mechanism (the causal story from local rule to
global pattern). Kept separate from registry so it's easy to edit."""

OVERVIEW_HTML = """
<h2>Emergence — complex behavior from simple rules</h2>
<p>A flock of starlings wheels as one. A market settles on a price. An ant colony
finds the shortest path to food. No bird leads, no trader sees the whole market,
no ant holds a map — each follows a few simple <b>local</b> rules, yet coherent,
large-scale order appears that none of them designed or even represents. That is
<b>emergence</b>: global structure and behavior arising from many agents
interacting locally, with no central controller and no blueprint.</p>

<p>The recurring surprise is that <b>the pattern is not "in" any individual</b> — it
is a collective property. You cannot find the flock in a single bird, segregation
in one household's mild preference, or a memory in one neuron. The behavior lives
in the <i>interactions</i>.</p>

<p>This catalog collects <b>32 of these phenomena</b> from across physics, biology,
and society — flocking, phase separation, Turing patterns, synchronization,
epidemics, opinion dynamics, wealth inequality, division of labor,
self-organized criticality, and more. Each is a <b>minimal model</b>: the
simplest rule that still produces the effect, so you can see exactly how little is
needed.</p>

<p><b>How to read each entry:</b> play the animation, note <i>what to look for</i>,
read <i>how it emerges</i> (the mechanism), see the few lines of code that drive
it, and watch the validated detector recognize the pattern and report its
calibrated confidence.</p>

<p>The through-line: <b>complexity does not require a complex cause.</b> A handful of
local rules — copy your neighbor, avoid the crowd, react to a threshold — is
enough to generate the rich collective behavior that fills the living and
physical world.</p>
"""

MECHANISM = {
"P1": "Each agent moves only if too few of its immediate neighbors share its type — a mild local preference, not a wish to segregate. But every move changes other agents' neighborhoods, triggering more moves; the feedback cascades until the grid sorts into large single-type blocks. Strong global segregation emerges from weak individual preference, with no one intending it.",
"P2": "Particles slow down where they are crowded. Slower particles linger in dense regions, making those regions denser still — a positive feedback. With no attraction between particles at all, this density-dependent speed alone splits the system into dense and dilute phases, like a gas condensing purely from how motion responds to crowding.",
"P3": "Two species diffuse at different rates: a slow, self-amplifying 'activator' and a fast-spreading 'inhibitor'. Activation builds a local peak while the faster inhibitor suppresses the surroundings, stopping peaks from merging. This short-range activation with long-range inhibition selects one characteristic spacing, turning a near-uniform medium into a stationary striped/spotted pattern — Turing's insight that diffusion can create structure, not just smear it away.",
"P4": "Each animal deposits scent as it moves and steps away from cells marked by others. Because everyone retreats from foreign scent, overlapping ranges are pushed apart until each settles into its own patch with sharp borders. Exclusive territories emerge from a simple 'avoid the neighbors' smell' rule — no fighting, no negotiation.",
"P5": "Each agent steers toward the average heading of its nearby neighbors, plus a little noise. Alignment is contagious: a locally agreed direction spreads neighbor to neighbor until the whole group moves together. Below a noise threshold this local copying tips into global polar order — a flock — with no leader and no shared goal.",
"P6": "Self-propelled particles attract at a distance and repel up close. Balanced, these forces cannot settle into a still clump or a straight march, so the group resolves the tension by circulating. A rotating mill (vortex) emerges purely from the interplay of propulsion with attraction–repulsion.",
"P7": "In two opposing streams, each person just sidesteps to avoid collisions. Following someone already heading your way is easier than crossing the opposing flow, so people fall in behind like-movers, and this self-reinforces. Orderly counter-flowing lanes emerge from nothing but local collision-avoidance.",
"P8": "Each driver accelerates in open road, brakes for the car ahead, and occasionally dawdles. One random slow-down forces the follower to brake, then the next, and the disturbance travels backward as a growing jam. Stop-and-go waves emerge above a critical density — congestion with no crash or bottleneck causing it.",
"P9": "Each oscillator nudges its rhythm toward its neighbors'. When coupling beats the spread of natural frequencies, a few that happen to align pull others in, recruiting the rest in a cascade until the whole population ticks together. Spontaneous synchrony — fireflies, applause, pacemaker cells — emerges from each unit simply adjusting to those around it.",
"P10": "Identical oscillators, coupled more strongly to nearby ones, surprisingly split: part of the population locks into synchrony while the rest stays incoherent — and both states persist side by side. This coexistence of order and disorder among identical, symmetrically-coupled units was long thought impossible; it emerges from the non-local coupling alone.",
"P11": "On a lattice, prey reproduce into empty cells and predators eat adjacent prey, then starve without food. Locally, predators boom where prey are dense, crash the prey, then crash themselves — letting prey recover. These local cycles couple across space into traveling waves and sustained population oscillations from simple eat / reproduce / die rules.",
"P12": "Three strains play rock–paper–scissors: each beats one neighbor and loses to another. No strain can win outright, so locally each chases the type it dominates and flees the type that dominates it. This intransitive loop organizes space into rotating spiral waves that let all three coexist indefinitely — diversity maintained by a cycle, not by balance.",
"P13": "Each cell rests, fires when enough neighbors fire, then briefly cannot fire again (a refractory period). A firing cell excites its neighbors, but the refractory trail behind it blocks back-propagation, so activity only travels forward — into self-sustaining spiral and target waves. The same rule drives heart tissue and the Belousov–Zhabotinsky reaction.",
"P14": "Grains drop on a pile; when a site exceeds a threshold slope it topples, spilling onto neighbors and sometimes triggering a chain reaction. The pile self-organizes to a critical state where one more grain can set off an avalanche of ANY size — the avalanche sizes follow a scale-free power law with no characteristic scale. Criticality emerges spontaneously, with no parameter tuned to reach it (self-organized criticality).",
"P15": "Cells live or die by counting live neighbors — three rules, nothing more. Yet stable structures, oscillators, and 'gliders' that travel and collide appear, and their interactions can carry and transform information. Universal computation emerges from a trivially simple update rule: complexity with no complex cause.",
"P16": "Neurons flip to match the weighted majority of their inputs, where the weights encode several stored patterns. From a corrupted cue the network rolls 'downhill' in an energy landscape whose valleys are the memories, settling into the nearest one. Content-addressable recall — reconstructing the whole from a fragment — emerges from simple threshold updates.",
"P17": "Each agent senses the environmental gradient poorly and noisily, but also tends to move with its neighbors. Averaging over many noisy individuals cancels the error, so the group climbs the gradient far more accurately than any member could alone. 'Many wrongs' become a collective right — accurate navigation from inaccurate individuals.",
"P18": "Each agent simply copies the opinion of a random neighbor. Copying makes local agreement self-reinforcing: single-opinion patches grow and merge, their boundaries wandering, until one opinion fills the whole population. Consensus emerges from imitation alone — no persuasion, no leaders, no goal.",
"P19": "A small minority knows the preferred direction; everyone else just stays with the group. Because the informed few consistently pull the same way while the uninformed cancel out, the whole group is steered accurately — and the bigger the group, the smaller the informed fraction needed. Effective leadership emerges with no signaling or designated leaders.",
"P20": "Each cell releases a shared signal and switches behavior only once the surrounding signal crosses a threshold. At low density the signal disperses harmlessly; once enough cells accumulate it, they all flip on together, sharply and with hysteresis. A population-scale decision emerges from each cell reading nothing but a local concentration.",
"P21": "Each person averages opinions only with others already close to their own (bounded confidence). Moderates between two leanings get pulled both ways and thin out, while like-minded clusters tighten. A single continuous spread fragments into separated camps — polarization emerging from being open-minded, but only toward the similar.",
"P22": "Each infected cell infects susceptible neighbors with some probability, then recovers and becomes immune. Infection can only advance into still-susceptible territory, so it sweeps outward as a wavefront leaving immunity behind. The epidemic curve and its spatial wave emerge from a purely local infect-then-recover rule.",
"P23": "Agents repeatedly choose a side and win only if they picked the minority. Since everyone wants to be where others are not, they learn to spread out, and attendance self-balances around capacity with smaller swings than random choosing would give. Efficient load-balancing emerges with no communication and no coordinator.",
"P24": "A controller measures how far a variable has drifted from its set-point and pushes back in proportion. Negative feedback automatically cancels disturbances, holding the variable steady against perturbations. Stable self-regulation — the principle behind body temperature and the thermostat — emerges from 'correct the error' alone.",
"P25": "Many different starting states flow 'downhill' on a canalized landscape toward the same attractor. Because trajectories are funnelled into shared valleys, the destination is set by the landscape, not the starting point. Robust convergence to one outcome from diverse beginnings — Waddington's developmental canalization — emerges from the basins of the dynamics.",
"P26": "A signal too weak to cross a threshold can be detected once noise is added: the right amount of noise occasionally bumps the system over in time with the signal. Too little noise and nothing crosses; too much and the signal is swamped — so detection peaks at an intermediate noise level. Noise helping a signal, counter-intuitively, emerges from threshold dynamics.",
"P27": "Cooperators and defectors play their neighbors and then copy whichever nearby strategy scored best. A lone cooperator is exploited, but cooperators that cluster shield one another and out-earn the defectors around them, so clusters persist and spread. Cooperation survives in a world that rewards defection — purely because interaction is local (spatial structure).",
"P28": "Two agents trade by staking a random fraction of the poorer one's wealth on a fair coin flip. Even though every trade is fair, the multiplicative nature of gambling means losses hurt the poor proportionally more, so wealth drifts toward whoever is already ahead. Extreme inequality condenses from repeated fair exchanges — no cheating required.",
"P29": "Agents lay attractant on the paths they use, and the marking evaporates over time. Heavily-used short routes are reinforced faster than they fade, while redundant ones disappear, so the network self-prunes toward a few efficient connections. Near-optimal transport networks (ant trails, slime mould) emerge from reinforce-and-evaporate, with no surveyor.",
"P30": "A catalyst stitches substrate particles into links that form a boundary, while those links also decay and must be continually replaced. When production keeps pace with decay, a thin closed membrane self-assembles and maintains itself around the catalyst. A self-producing, self-bounded unit — a minimal 'cell' — emerges from local make-and-decay reactions.",
"P31": "Agents sometimes accept a locally worse move now because it unlocks a better outcome later. Foregoing immediate reward, applied locally and repeatedly, lets the collective reach configurations that purely greedy steps cannot. Goal-directed-looking patience emerges from simple per-agent trade-offs rather than any global plan.",
"P32": "Each agent takes whichever task it senses most strongly, and its threshold for that task falls with practice, so it specialises. Because specialists handle their task before others' thresholds are tripped, stable, distinct roles differentiate across the group. Division of labor — like castes in social insects — emerges from response thresholds plus learning, with no assignment.",
}

METHODS_HTML = """
<h2>How these models are validated</h2>
<p>The 32 phenomena in this catalog are <b>established, literature-validated
science</b> — Schelling segregation, Turing patterns, Kuramoto synchronisation,
and so on. The science is not what is new here. What is new is a set of
<b>minimal executable models</b>, each paired with a <b>detector</b> — an
algorithm that decides, from the raw simulation output alone, whether the named
pattern is present — and a test of whether that detector actually works.</p>
<p>"Validated", on every pattern page, means the detector passed the tests
below. It does <b>not</b> mean the science is novel, and it does <b>not</b> mean
the detector is infallible. Where it has limits, we show them.</p>

<h3>1 · Negative-control discrimination — the main test</h3>
<p>A detector that just says "yes" is worthless. So each detector faces a
<b>panel of negatives it must reject</b>, in three families:</p>
<ul>
<li><b>Synthetic nulls</b> — structureless or random inputs with no pattern at all.</li>
<li><b>Catalog look-alikes</b> — the <i>other</i> patterns in this catalog,
several of which superficially resemble the target.</li>
<li><b>Failed regimes</b> — the target's own model run in parameter settings
where the phenomenon does not form.</li>
</ul>
<p>Across the 31 patterns with a panel, detectors rejected <b>636 of 636</b>
negative controls — a true-negative rate of <b>1.0</b>, every panel at 100%.
Each pattern page shows its own count ("rejected N / N negative controls"). This
is the evidence that a detector recognizes <i>its</i> pattern, not merely "something happened".</p>

<h3>2 · Cross-model generalization — an instrument, not a lookup table</h3>
<p>A detector could cheat by memorising quirks of the one model it was built on.
To check, for seven phenomena we built a <b>second, independent implementation</b>
(for example, associative memory via a Boolean gene-regulatory network instead of
a Hopfield net) and tested the detector on the model it was <i>not</i> tuned on.
All seven still ranked their own pattern first; six fired a firm match. The
detectors track the <b>phenomenon</b>, not the implementation.</p>

<h3>3 · Self-recognition — necessary but weak</h3>
<p>Running all 32 detectors over each positive, the correct detector wins for
<b>30 of 31</b> patterns (one cross-fire between two grid-coarsening patterns).
We report it for completeness, but it is <b>in-sample</b> evidence — the weakest
of the three tests, and not a validation on its own.</p>

<h3>Honest caveats</h3>
<ul>
<li><b>A withdrawn effect size.</b> An earlier draft reported large Cohen's
<i>d</i> values, some infinite. Those were an <b>artefact</b> of computing
<i>d</i> over a discrete confidence score: a cleanly discriminating detector has
near-zero within-group variance, which inflates <i>d</i> without meaning. We have
withdrawn those numbers. On the continuous metric, clean separation often makes
the effect size <b>undefined</b> rather than infinite, and we say so.</li>
<li><b>"Generic emergence" is not the verdict.</b> Each page also shows a coarse,
cross-pattern "generic emergence index". It is a rough screen, <b>not</b> the
pattern-specific result — a detector can fire a definitive match while this index
is low (associative memory is one). Trust the verdict and the negative-control
panel, not this number.</li>
<li><b>Tiers differ.</b> A <b>definitive</b> match is stronger than a
<b>confirmation</b> or <b>screening</b> one; a few patterns reach only screening
on a single seed.</li>
<li><b>Known limits, surfaced not hidden.</b> One pattern (activity-induced phase
separation) rejects its negatives but does not yet reach a firing tier on its own
positive — shown as "unclassified". The delayed-gratification / sorting pattern is
validated by a separate multi-run test, not a panel. An integral-controller
homeostat is recognized but stays below threshold.</li>
<li>All of this validates the faithfulness of our <b>detectors and test
harness</b> — not any new claim about nature.</li>
</ul>

<h3>Reproduce</h3>
<p>Everything is open. The models, detectors, the Phase-2a negative-control
panels and the full validation note (section&nbsp;4A) are in the source
repository:
<a href="https://github.com/matthewhmaxwell/emergent-pattern-catalog" target="_blank" rel="noopener">github.com/matthewhmaxwell/emergent-pattern-catalog</a>
— specifically <a href="https://github.com/matthewhmaxwell/emergent-pattern-catalog/blob/validation-rebuild/docs/paper_section4A_validation_rebuild.md" target="_blank" rel="noopener">the validation note</a>
and the <a href="https://github.com/matthewhmaxwell/emergent-pattern-catalog/tree/validation-rebuild/analysis/outputs" target="_blank" rel="noopener">per-pattern panel outputs</a>.</p>
"""

# Where each phenomenon shows up in the real world (natural / technical / both) — grounds the
# "simple rules -> complex behavior" point in concrete instances.
WHERE = {
"P1": "Residential and school segregation, and the clustering of like-with-like in social networks — sharp divides emerging from only mild individual preferences (society); echoed in self-sorting of some shaken granular and colloidal mixtures (physics).",
"P2": "Swimming bacteria and cells in dense tissue piling into clusters (biology); synthetic active matter — self-propelled Janus colloids and microswimmers that jam into dense patches with no glue (lab/technology).",
"P3": "Animal coat and skin markings — leopard spots, zebra stripes, fish pigmentation — plus sand ripples and banded desert vegetation (nature); the Belousov–Zhabotinsky reaction and patterning in chemical reactors (lab).",
"P4": "Animal home ranges marked by scent or song — wolves, songbirds, reef fish hold non-overlapping territories (nature); spatial coverage allocation such as cell-tower or sensor footprints (technology).",
"P5": "Starling murmurations, fish schools and insect swarms (nature); drone swarms and swarm robotics that align direction from local rules alone (technology).",
"P6": "Fish bait-balls and torus mills, circling army-ant 'death spirals', and bacterial vortices (nature); rotating/orbiting robot swarms (technology).",
"P7": "Pedestrians spontaneously forming lanes in dense counterflow, and two-way ant trails (society/nature); self-organized lanes in traffic and particle sorting in microfluidic channels (technology).",
"P8": "Phantom highway jams that form with no accident and crawl backward through the traffic (society); congestion waves in data networks and clogging in granular and pipe flow (technology/physics).",
"P9": "Synchronously flashing fireflies, chirping crickets, and pacemaker cells in heart and brain (nature); coupled lasers, power-grid generators locking to one frequency, and audiences clapping in unison (technology/society).",
"P10": "Unihemispheric sleep — dolphins and many birds keep one brain hemisphere awake while the other sleeps (nature); arrays of coupled oscillators, lasers and power grids where synced and desynced groups coexist (technology).",
"P11": "Boom-and-bust predator–prey cycles such as the lynx–hare record, plus host–parasite and disease cycles (nature); chemical oscillators and ecosystem/epidemic models (technical).",
"P12": "Rock-paper-scissors strategy cycling in side-blotched lizards and in toxin-producing / resistant / sensitive microbial strains, and reef space competition (nature); a workhorse of evolutionary game theory (technical).",
"P13": "Spiral electrical waves in heart tissue (the basis of dangerous arrhythmias) and in neural tissue, plus chemical BZ waves and slime-mould signaling (nature/lab); reaction–diffusion 'unconventional' computing (technology).",
"P14": "Earthquakes, landslides, forest fires and neuronal avalanches in the brain — power-law distributed with no characteristic size (nature); cascading failures in power grids and financial markets (technology/society).",
"P15": "Chiefly computational: a few local rules give rise to gliders, logic gates and universal computation — the cleanest proof that simple rules can compute anything (technology/theory); echoed in some excitable chemical and physical media (nature).",
"P16": "Memory recall and pattern completion in the brain's attractor networks — recognizing a face from a glimpse (nature); content-addressable memory, denoising, and the modern 'Hopfield' attention layers inside transformers (technology).",
"P17": "Fish schools and cell collectives that track a faint light, temperature or chemical gradient far better than any member alone (nature); distributed sensor networks and source-seeking robot swarms (technology).",
"P18": "Adoption of a shared convention, language or opinion across a population, and neutral gene fixation by drift (society/biology); distributed consensus protocols and spin-system relaxation (technology/physics).",
"P19": "A few informed individuals steering a whole fish school or bee swarm to a resource or nest site (nature); leader–follower coordination and information-limited navigation in robot swarms (technology).",
"P20": "Bacterial quorum sensing — light, biofilms and virulence switched on only above a density threshold (biology); threshold-triggered distributed systems and social tipping points / bandwagons (technology/society).",
"P21": "Political and social polarisation into hardened camps with a hollowed-out middle (society); recommendation-driven echo chambers and filter bubbles online (technology).",
"P22": "Epidemics of infectious disease and the spread of rumours, memes and fads as traveling waves (nature/society); computer-virus and malware propagation and viral marketing cascades (technology).",
"P23": "The El Farol bar problem and market speculation (the minority choice wins) and drivers picking less-crowded routes (society); distributed load balancing across servers and links (technology).",
"P24": "Body-temperature, blood-sugar and pH regulation, and ecological balance (biology); thermostats, cruise control, PID controllers and datacenter autoscaling (technology).",
"P25": "Embryos reaching the same body plan from very different starts, and regeneration after injury — developmental robustness (biology); robust attractor-based control and error-correcting computation (technology).",
"P26": "Sensory neurons using noise to detect otherwise-imperceptible signals — prey detection in paddlefish and crayfish, human hearing and balance — and noise-paced ice-age cycles (nature); deliberately added noise (dithering) to sharpen detectors and analog-to-digital converters (technology).",
"P27": "Evolution of cooperation in structured populations — microbial cooperators, cells in tissue, social networks — where cooperators survive by clustering (biology/society); sustaining cooperation in spatial and networked multi-agent systems (technology).",
"P28": "Extreme wealth and income inequality (Pareto, 'rich get richer') and city-size distributions (society); preferential attachment that concentrates links on a few hub nodes in networks (technology).",
"P29": "Slime-mould and ant transport networks, fungal mycelium, leaf veins and blood vessels — efficient networks built with no central plan (nature); network design and transport-route optimization, including Physarum-inspired algorithms (technology).",
"P30": "The living cell's self-producing membrane and metabolism, and candidate origin-of-life protocells (biology); artificial-life and self-repairing / self-replicating systems (technology).",
"P31": "Cells reaching target arrangements through purely local moves in morphogenesis and tissue sorting — a window on basal cognition (biology); sorting algorithms that show unexpected 'competencies' like delayed gratification (Levin's algorithmic 'Platonic space') (technology/theory).",
"P32": "Division of labor and castes in ant and bee colonies via response thresholds, and cell differentiation (biology); task allocation and work scheduling in swarm robotics and distributed systems (technology).",
}

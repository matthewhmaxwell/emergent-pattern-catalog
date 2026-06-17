# -*- coding: utf-8 -*-
"""Educational content for the EPC gallery: a general intro to emergence, and a
per-pattern 'how it emerges' mechanism (the causal story from local rule to
global pattern). Kept separate from registry so it's easy to edit."""

OVERVIEW_HTML = """
<h2>Emergence — complex behaviour from simple rules</h2>
<p>A flock of starlings wheels as one. A market settles on a price. An ant colony
finds the shortest path to food. No bird leads, no trader sees the whole market,
no ant holds a map — each follows a few simple <b>local</b> rules, yet coherent,
large-scale order appears that none of them designed or even represents. That is
<b>emergence</b>: global structure and behaviour arising from many agents
interacting locally, with no central controller and no blueprint.</p>

<p>The recurring surprise is that <b>the pattern is not "in" any individual</b> — it
is a collective property. You cannot find the flock in a single bird, segregation
in one household's mild preference, or a memory in one neuron. The behaviour lives
in the <i>interactions</i>.</p>

<p>This catalog collects <b>32 of these phenomena</b> from across physics, biology,
and society — flocking, phase separation, Turing patterns, synchronization,
epidemics, opinion dynamics, wealth inequality, division of labour,
self-organized criticality, and more. Each is a <b>minimal model</b>: the
simplest rule that still produces the effect, so you can see exactly how little is
needed.</p>

<p><b>How to read each entry:</b> play the animation, note <i>what to look for</i>,
read <i>how it emerges</i> (the mechanism), see the few lines of code that drive
it, and watch the validated detector recognise the pattern and report its
calibrated confidence.</p>

<p>The through-line: <b>complexity does not require a complex cause.</b> A handful of
local rules — copy your neighbour, avoid the crowd, react to a threshold — is
enough to generate the rich collective behaviour that fills the living and
physical world.</p>
"""

MECHANISM = {
"P1": "Each agent moves only if too few of its immediate neighbours share its type — a mild local preference, not a wish to segregate. But every move changes other agents' neighbourhoods, triggering more moves; the feedback cascades until the grid sorts into large single-type blocks. Strong global segregation emerges from weak individual preference, with no one intending it.",
"P2": "Particles slow down where they are crowded. Slower particles linger in dense regions, making those regions denser still — a positive feedback. With no attraction between particles at all, this density-dependent speed alone splits the system into dense and dilute phases, like a gas condensing purely from how motion responds to crowding.",
"P3": "Two species diffuse at different rates: a slow, self-amplifying 'activator' and a fast-spreading 'inhibitor'. Activation builds a local peak while the faster inhibitor suppresses the surroundings, stopping peaks from merging. This short-range activation with long-range inhibition selects one characteristic spacing, turning a near-uniform medium into a stationary striped/spotted pattern — Turing's insight that diffusion can create structure, not just smear it away.",
"P4": "Each animal deposits scent as it moves and steps away from cells marked by others. Because everyone retreats from foreign scent, overlapping ranges are pushed apart until each settles into its own patch with sharp borders. Exclusive territories emerge from a simple 'avoid the neighbours' smell' rule — no fighting, no negotiation.",
"P5": "Each agent steers toward the average heading of its nearby neighbours, plus a little noise. Alignment is contagious: a locally agreed direction spreads neighbour to neighbour until the whole group moves together. Below a noise threshold this local copying tips into global polar order — a flock — with no leader and no shared goal.",
"P6": "Self-propelled particles attract at a distance and repel up close. Balanced, these forces cannot settle into a still clump or a straight march, so the group resolves the tension by circulating. A rotating mill (vortex) emerges purely from the interplay of propulsion with attraction–repulsion.",
"P7": "In two opposing streams, each person just sidesteps to avoid collisions. Following someone already heading your way is easier than crossing the opposing flow, so people fall in behind like-movers, and this self-reinforces. Orderly counter-flowing lanes emerge from nothing but local collision-avoidance.",
"P8": "Each driver accelerates in open road, brakes for the car ahead, and occasionally dawdles. One random slow-down forces the follower to brake, then the next, and the disturbance travels backward as a growing jam. Stop-and-go waves emerge above a critical density — congestion with no crash or bottleneck causing it.",
"P9": "Each oscillator nudges its rhythm toward its neighbours'. When coupling beats the spread of natural frequencies, a few that happen to align pull others in, recruiting the rest in a cascade until the whole population ticks together. Spontaneous synchrony — fireflies, applause, pacemaker cells — emerges from each unit simply adjusting to those around it.",
"P10": "Identical oscillators, coupled more strongly to nearby ones, surprisingly split: part of the population locks into synchrony while the rest stays incoherent — and both states persist side by side. This coexistence of order and disorder among identical, symmetrically-coupled units was long thought impossible; it emerges from the non-local coupling alone.",
"P11": "On a lattice, prey reproduce into empty cells and predators eat adjacent prey, then starve without food. Locally, predators boom where prey are dense, crash the prey, then crash themselves — letting prey recover. These local cycles couple across space into travelling waves and sustained population oscillations from simple eat / reproduce / die rules.",
"P12": "Three strains play rock–paper–scissors: each beats one neighbour and loses to another. No strain can win outright, so locally each chases the type it dominates and flees the type that dominates it. This intransitive loop organises space into rotating spiral waves that let all three coexist indefinitely — diversity maintained by a cycle, not by balance.",
"P13": "Each cell rests, fires when enough neighbours fire, then briefly cannot fire again (a refractory period). A firing cell excites its neighbours, but the refractory trail behind it blocks back-propagation, so activity only travels forward — into self-sustaining spiral and target waves. The same rule drives heart tissue and the Belousov–Zhabotinsky reaction.",
"P14": "Grains drop on a pile; when a site exceeds a threshold slope it topples, spilling onto neighbours and sometimes triggering a chain reaction. The pile self-organises to a critical state where one more grain can set off an avalanche of ANY size — the avalanche sizes follow a scale-free power law with no characteristic scale. Criticality emerges spontaneously, with no parameter tuned to reach it (self-organized criticality).",
"P15": "Cells live or die by counting live neighbours — three rules, nothing more. Yet stable structures, oscillators, and 'gliders' that travel and collide appear, and their interactions can carry and transform information. Universal computation emerges from a trivially simple update rule: complexity with no complex cause.",
"P16": "Neurons flip to match the weighted majority of their inputs, where the weights encode several stored patterns. From a corrupted cue the network rolls 'downhill' in an energy landscape whose valleys are the memories, settling into the nearest one. Content-addressable recall — reconstructing the whole from a fragment — emerges from simple threshold updates.",
"P17": "Each agent senses the environmental gradient poorly and noisily, but also tends to move with its neighbours. Averaging over many noisy individuals cancels the error, so the group climbs the gradient far more accurately than any member could alone. 'Many wrongs' become a collective right — accurate navigation from inaccurate individuals.",
"P18": "Each agent simply copies the opinion of a random neighbour. Copying makes local agreement self-reinforcing: single-opinion patches grow and merge, their boundaries wandering, until one opinion fills the whole population. Consensus emerges from imitation alone — no persuasion, no leaders, no goal.",
"P19": "A small minority knows the preferred direction; everyone else just stays with the group. Because the informed few consistently pull the same way while the uninformed cancel out, the whole group is steered accurately — and the bigger the group, the smaller the informed fraction needed. Effective leadership emerges with no signalling or designated leaders.",
"P20": "Each cell releases a shared signal and switches behaviour only once the surrounding signal crosses a threshold. At low density the signal disperses harmlessly; once enough cells accumulate it, they all flip on together, sharply and with hysteresis. A population-scale decision emerges from each cell reading nothing but a local concentration.",
"P21": "Each person averages opinions only with others already close to their own (bounded confidence). Moderates between two leanings get pulled both ways and thin out, while like-minded clusters tighten. A single continuous spread fragments into separated camps — polarization emerging from being open-minded, but only toward the similar.",
"P22": "Each infected cell infects susceptible neighbours with some probability, then recovers and becomes immune. Infection can only advance into still-susceptible territory, so it sweeps outward as a wavefront leaving immunity behind. The epidemic curve and its spatial wave emerge from a purely local infect-then-recover rule.",
"P23": "Agents repeatedly choose a side and win only if they picked the minority. Since everyone wants to be where others are not, they learn to spread out, and attendance self-balances around capacity with smaller swings than random choosing would give. Efficient load-balancing emerges with no communication and no coordinator.",
"P24": "A controller measures how far a variable has drifted from its set-point and pushes back in proportion. Negative feedback automatically cancels disturbances, holding the variable steady against perturbations. Stable self-regulation — the principle behind body temperature and the thermostat — emerges from 'correct the error' alone.",
"P25": "Many different starting states flow 'downhill' on a canalized landscape toward the same attractor. Because trajectories are funnelled into shared valleys, the destination is set by the landscape, not the starting point. Robust convergence to one outcome from diverse beginnings — Waddington's developmental canalization — emerges from the basins of the dynamics.",
"P26": "A signal too weak to cross a threshold can be detected once noise is added: the right amount of noise occasionally bumps the system over in time with the signal. Too little noise and nothing crosses; too much and the signal is swamped — so detection peaks at an intermediate noise level. Noise helping a signal, counter-intuitively, emerges from threshold dynamics.",
"P27": "Cooperators and defectors play their neighbours and then copy whichever nearby strategy scored best. A lone cooperator is exploited, but cooperators that cluster shield one another and out-earn the defectors around them, so clusters persist and spread. Cooperation survives in a world that rewards defection — purely because interaction is local (spatial structure).",
"P28": "Two agents trade by staking a random fraction of the poorer one's wealth on a fair coin flip. Even though every trade is fair, the multiplicative nature of gambling means losses hurt the poor proportionally more, so wealth drifts toward whoever is already ahead. Extreme inequality condenses from repeated fair exchanges — no cheating required.",
"P29": "Agents lay attractant on the paths they use, and the marking evaporates over time. Heavily-used short routes are reinforced faster than they fade, while redundant ones disappear, so the network self-prunes toward a few efficient connections. Near-optimal transport networks (ant trails, slime mould) emerge from reinforce-and-evaporate, with no surveyor.",
"P30": "A catalyst stitches substrate particles into links that form a boundary, while those links also decay and must be continually replaced. When production keeps pace with decay, a thin closed membrane self-assembles and maintains itself around the catalyst. A self-producing, self-bounded unit — a minimal 'cell' — emerges from local make-and-decay reactions.",
"P31": "Agents sometimes accept a locally worse move now because it unlocks a better outcome later. Foregoing immediate reward, applied locally and repeatedly, lets the collective reach configurations that purely greedy steps cannot. Goal-directed-looking patience emerges from simple per-agent trade-offs rather than any global plan.",
"P32": "Each agent takes whichever task it senses most strongly, and its threshold for that task falls with practice, so it specialises. Because specialists handle their task before others' thresholds are tripped, stable, distinct roles differentiate across the group. Division of labour — like castes in social insects — emerges from response thresholds plus learning, with no assignment.",
}

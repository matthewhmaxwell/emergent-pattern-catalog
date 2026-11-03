# planted:N15

**Reported behavior (what the authors report as new):** Self-propelled particles that only steer toward neighbours inside a narrow finite-range vision cone form a system-spanning, continuously restructuring fractal-like network of particle streams (dimension ~1.5) that never coarsens.

**Rule on the card:** World: 10,000 point particles in a 100x100 periodic box (density 1), with no collisions. Rules: (1) each moves forward at constant speed 0.1 along its heading (small inertia); (2) it sees only particles within distance R (5 to 10) inside a vision cone of half-angle pi/8 around its heading; (3) it turns toward the positions of the visible particles, turning rate 3 times the average sine of the angle between its heading and the bearing to each (no copying of their headings); (4) its heading gets rotational noise of strength 0.08. Seeing is one-sided: being seen does not imply seeing back. Run 10 million steps of size 0.005.

**Observed features on the card:**
  - pattern_fractal_dimension = 1.5 (fractal dimension of particle pattern (1 = lines, 2 = uniform))
  - void_area_scaling_exponent = 3 (areas of empty voids scale as vision range to this power)
  - commonest_junction_degree = 3 (number of particle streams meeting at the most common junction)
  - max_common_junction_degree = 4 (largest junction degree seen; five or more essentially absent)

## Reply R01

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Filament/nematic-lane network state of the minimal cognitive flocking (Barberis-Peruani vision-cone, position-based attraction) model", "notes": "Interconnecting stream/nematic-band networks are documented for vision-cone position-attraction particles; honeycomb-like junction networks appear in iABP mixtures; similar networks appear in Physarum-type forward-sensing agents. The twist is quantitative only: D_f ~1.5, void area ~R^3, mostly 3-way junctions. Mostly 3-way junctions is generic for 2D networks. D_f comes from barely one decade of scale. The R^3 exponent comes from only a factor-of-2 range of R and is dimensionally suggestive of a second length (turning radius or persistence length); worth checking over a wider R range before treating it as meaningful. Barberis-Peruani full text not verified."}

**Reply text:**

**Verdict: KNOWN-VARIANT, which agrees with the first pass, but for somewhat different reasons.** The pattern itself, a network of particle streams built by vision-cone attraction, is documented. The specific numbers (fractal dimension, void scaling, junction counts) are not, and I don't think they are strong enough to count as a new phenomenon.

**The model is already in the literature.** Barberis & Peruani (PRL 2016) studied exactly this kind of system: active particles, with no memory, that interact by a short-ranged, position-based, attractive force, which acts inside a vision cone (VC), and lack velocity-velocity alignment. Depending on parameters they report aggregates or millinglike patterns, the formation of moving—locally polar—files with particles at the front of these structures acting as effective leaders, and the self-organization of particles into macroscopic nematic structures. A follow-up with topological rather than distance-based vision describes how various small nematic bands grow in size and interconnect. That is the same "streams joining into a web" picture.

**Junction-dominated networks are also documented in close relatives:**
- In binary mixtures of vision-steered particles with half-angles π/8 and π/4, Negi et al. (NJP 2025) report locally diffuse honeycomb-like structures and the prevalence of the honeycomb lattice structures. Honeycombs are made of 3-way junctions.
- In Jones's agent-based slime-mould (Physarum) model, agents sense ahead and steer toward what they sense. That population spontaneously forms complex and dynamic transport networks. The main difference is that those agents follow a chemical trail rather than looking at each other directly.

**What's new is only the quantification, and it is weak:**
- **Mostly 3-way junctions, rarely more than 4:** this is what almost any coarsening 2D network does, because 3-way junctions are the stable ones. It is not a distinctive signature.
- **Fractal dimension 1.5:** in a 100×100 box with streams spaced about 5–10 apart, box counting covers barely one decade of scale. The number most likely describes "a mesh of lines with some width", not a true fractal.
- **Void area ∝ R³:** the vision range was only varied from 5 to 10, a factor of 2, so the exponent is poorly constrained. If a single length set the void size, area would grow as R², so a cube suggests a second length is involved, such as the turning-limited radius or the noise persistence length. That would be interesting if it held up over a wider range of R, but as it stands it is fragile.

**Artifact or real pattern?** It looks like a genuine emergent pattern, not a trivial outcome. Streams with junctions and voids whose size tracks the vision range is a real structure, consistent with the known stream and nematic-band states. It is not a frozen state, a random scatter, or a numerical blow-up.

**Gaps in this check:**
- I could not read the full Barberis & Peruani text or its supplement, so I could not confirm the exact sine-of-bearing turning term. I also could not see whether their figures already show this web-like state.
- I did not find a separate Peruani 2017 paper to check, as the first pass suggested.
- Searches for fractal dimension, void scaling and junction-degree statistics in vision-cone systems found nothing.

## Reply R21

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Position-based vision-cone ('minimal cognitive flocking') active particles forming files/nematic streams; network morphology analogous to Jones-type Physarum agent transport networks", "notes": "The behaviour (vision-cone, non-reciprocal, position-based attraction producing files/streams that form a network) is documented. The twist is the quantified network shape: fractal dimension ~1.5, void area ~R^3, mostly degree-3 junctions with degree 5+ absent. Mostly degree-3 junctions is generic for coarsening 2D networks. R^3 void scaling is the only possibly meaningful item, but it rests on a factor-2 range of R in a small box and may reflect unfinished coarsening. The fractal dimension covers only about one decade of scale. I could not read the Barberis-Peruani full text to check whether their nematic regime looks like a network."}

**Reply text:**

**Verdict: KNOWN-VARIANT. I agree with the first pass, but only partly.** The model class and its tendency to form streams and chains are well documented. The specific network-shape measurements are not. I don't think that gap makes this a discovery.

**What is documented**
- **Barberis & Peruani (PRL 117, 248001, 2016)** is essentially this model. Their particles interact by a short-ranged, position-based, attractive force, which acts inside a vision cone (VC), and lack velocity-velocity alignment. They report aggregates or millinglike patterns, the formation of moving—locally polar—files with particles at the front of these structures acting as effective leaders, and the self-organization of particles into macroscopic nematic structures. Streams of particles moving both ways, which is what your network edges look like, fit their "nematic structures" regime. I couldn't read the full text to check whether any of their snapshots look like a network.
- **Negi, Winkler & Gompper (Soft Matter 2022)** use the same kind of steering: a limited turning rate in response to an instantaneous visual input of the positions of neighbors within a vision cone and a cut-off radius. The patterns they get depend on selected parameter values, and, in particular, the presence of excluded-volume interactions. So dropping collisions, as your system does, is a known way to change the outcome.
- **PRE 110, 064602 (2024)** shows that particles attracted to what they see in a vision cone form chains: active particles with nonreciprocal attraction form particle chains, which eventually become closed loops. That paper needed obstacles to keep the loops alive.
- **The closest match for the network shape itself is Jones's Physarum agent model** (Artificial Life, 2010). It is a different mechanism: agents sense ahead and follow a chemical trail that they lay down themselves. It does produce self-organised transport networks, with the circular contraction of network lacunae and bifurcation of network paths to maintain network connectivity. Your particles steer toward particles ahead of them instead of toward a trail, so the "sense ahead, steer toward density" loop is the same idea.

**What I couldn't find (the twist)**
I found no paper reporting a fractal dimension of about 1.5, void areas scaling as R³, or junction-degree statistics for a vision-cone model like this. These look like new measurements of a known pattern rather than a new behaviour. Specifically:
- **Mostly 3-way junctions with almost no 5-way ones** is what you'd expect from almost any 2D network that coarsens. Y-shaped junctions are stable, and 4-way ones are plausibly two streams crossing.
- **Void area ∝ R³** is the one result that might matter. It would mean the void size grows as R^1.5 rather than in proportion to R. But it rests on R between 5 and 10 only, a factor of 2, in a box just 10–20 R wide. That is too narrow to establish an exponent. The voids may also still be coarsening rather than at a steady size.
- **A fractal dimension of 1.5** measured over about one decade of scale in a 100×100 box is weak evidence of real fractal structure. It may just be a mix of 1D streams and empty space.

**Real pattern or artifact?** It looks real. At density 1 with R between 5 and 10, each particle sees roughly 10–40 others in its cone, so interactions are strong. A network of streams is the expected kind of outcome, not random noise, a frozen state or a numerical blow-up. The exponents, though, are not yet trustworthy.

**Searches I ran:** vision-cone particles plus networks, junctions and fractal dimension; the Barberis–Peruani model and its patterns; vision-perception active particles forming network-like lanes; Negi et al. without excluded volume; narrow-cone nematic patterns; vision-cone networks, cells and voids; Jones's Physarum network scale and junctions; Physarum-like networks without a chemical trail; vision-cone chains forming branching networks.

**To strengthen it:** check whether the void size settles over time, use a wider range of R or a bigger box, and compare how the pattern scale depends on R against Jones's results for sensor distance.

## Reply R37

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Minimal cognitive flocking model (vision-cone, position-based attraction, no alignment); iABP with visual perception", "notes": "The model class is clearly documented. A fractal (D~1.5) stream network with R^3 void-area scaling and mostly 3-way junctions was not found in any abstract or snippet returned. Full texts were not read, so a description in a figure, supplement or thesis, or under other terminology, cannot be excluded. Reported small-angle regimes in the Negi-Gompper papers (homogeneous or dilute phases, small clusters) differ in parameters: they include excluded volume, overdamped dynamics or alignment. Barberis-Peruani's polar files and nematic structures are the nearest reported pattern."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The model family is well documented. I found no source that reports the specific fractal network pattern.

**1. Model family**
- This is active matter: self-propelled particles with non-reciprocal, vision-cone, position-based steering. Other names for it are "minimal cognitive flocking model" and "intelligent active Brownian particles (iABPs) with visual perception".
- Your rule 3 turns each particle toward the positions of visible neighbours, with no heading alignment. That is the defining feature of this family.
- The 2016 Barberis–Peruani paper describes active particles with no memory, a short-ranged, position-based attractive force acting inside a vision cone, and no velocity-velocity alignment.
- Negi, Winkler and Gompper's iABP models use an orientational response with limited maneuverability to an instantaneous visual input of the positions of neighbors within a vision cone and a cut-off radius.
- A 2025 roadmap calls this steering non-additive and non-reciprocal.

**2. Is the specific behaviour reported?**
- The behaviour is a pattern of particle streams with a fractal dimension of about 1.5. It has void areas scaling as R³ and mostly 3-way (Y) junctions, with 4-way junctions as the maximum.
- I found no source describing this.
- Closest earlier source: Barberis & Peruani 2016 (PRL 117, 248001). Their vision cone breaks Newton's third law and produces aggregates or milling-like patterns, moving locally polar files with effective leaders, and macroscopic nematic structures with long-ranged nematic order. The "files" and "nematic structures" are stream-like. The abstract gives no fractal dimension, void scaling or junction statistics. I read only the abstract, not the full paper.
- Review: the Negi–Winkler–Gompper line. Soft Matter 18, 6167 (2022) is the main paper. PRResearch 6, 013118 (2024) adds alignment. The 2025 roadmap (J. Phys.: Condens. Matter) summarises the field.
- What those papers report at small vision angles differs from your result:
  - When the vision angle is low, very few particles are detected within the VC, no clusters can form, and particles are distributed homogeneously.
  - Small vision angles imply the formation of small clusters, while large vision angles lead to more complex clusters.
  - Their small-angle regimes are not network-like. Their parameters also differ from yours: they include excluded volume, strong overdamping, or alignment, and they use a different torque normalisation.
- Other vision-cone models are not equivalent:
  - Lavergne et al. (Science 2019) switch propulsion on and off depending on how many particles are in the cone.
  - The Vicsek-type vision-cone models align headings.
  - The XY and lattice-gas vision-cone models (Loos et al. PRL 2023; Du et al. 2025) are on a lattice and fixed in space. Their long-range order and aster clusters are not stream networks.
  - The topological-vision paper (PRResearch 7, 023234) uses Voronoi neighbours and finds gas, ordered, nematic-band and polar-band phases.

**3. What differs, and does it matter?**
- Your system has point particles with no collisions, a speed of 0.1 with small inertia, a narrow cone (π/8 and π/12), a range of 5 to 10 at density 1, and a mean-sine turning rule. That puts it in the Barberis–Peruani / Negi–Gompper family but in a regime they have not described.
- The likely difference is that your long runs (10⁷ steps) and narrow cone reach a coarsened, network-like steady state that earlier work did not measure. Earlier work used shorter runs, different parameters, or only looked for clusters and bands.
- The measured exponents may well be new even though the mechanism is known. The void-area scaling with R³ is the least likely to have been reported. I could not find the fractal dimension of 1.5, the junction-degree statistics or the R³ scaling anywhere.

**4. Verified citations**
- Barberis & Peruani, PRL 117, 248001 (2016): confirmed via the APS page and arXiv 1912.07929.
- Negi, Winkler & Gompper, Soft Matter 18, 6167 (2022): confirmed via the RSC page.
- Negi, Winkler & Gompper, Phys. Rev. Research 6, 013118 (2024): confirmed via the APS page.
- Goh, Winkler & Gompper, New J. Phys. 24, 093039 (2022), "Noisy pursuit and pattern formation of self-steering active particles": the title and journal came from a reference list, not from the paper itself. It is context only.
- 2025 motile active matter roadmap, J. Phys.: Condens. Matter: confirmed via the IOP page. I did not read the roadmap section on vision-based steering in full.

**5. Confidence: low to medium.**
- I ran only about 8 searches, and the search tool mostly returned abstracts and reference lists. I could not read full texts for narrow-cone, long-run, point-particle simulations.
- The pattern could be buried in a figure, supplement or thesis, or described in other terms such as "cobweb", "filamentary", "foam" or "nematic". Nothing I found describes it, but I have not ruled that out.

**

## Reply R40

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "non-reciprocal vision-cone (cognitive) flocking model; vision-based self-steering active particles", "notes": "Model class is firmly established; no source found reporting a fractal (D~1.5) branching stream network, void-area ~ R^3 scaling or degree-3-dominated junctions for finite-range, collision-free, weak-noise vision-cone attraction. Full texts of Barberis & Peruani 2016 and Brigatti & Peruani 2026 were not read and could contain a network regime; the search budget ran out before I could check."}

**Reply text:**

**Verdict: KNOWN-VARIANT, with low-to-medium confidence.** The model family is documented. I found no source that reports the specific network statistics seen in the simulation.

**1. Model family.** This is active matter, specifically self-propelled particles with non-reciprocal "vision-cone" interactions, also called minimal cognitive flocking or self-steering "intelligent active Brownian particles" (iABPs).
- Particles steer toward the positions of neighbours inside a forward cone.
- There is no velocity alignment.
- Being seen does not imply seeing back, so Newton's third law is broken.

**2. Closest literature.**
- **Earliest clear source (same rules):** Barberis & Peruani, PRL 117, 248001 (2016), "Large-Scale Patterns in a Minimal Cognitive Flocking Model: Incidental Leaders, Nematic Patterns, and Aggregates."
  - The model uses active particles with no memory. They interact through a short-ranged, position-based attractive force acting inside a vision cone, and they lack velocity-velocity alignment.
  - Depending on parameters, they report aggregates or milling-like patterns, moving polar files with effective leaders, and macroscopic nematic structures with long-ranged nematic order.
  - None of the snippets I saw mention fractal dimension, voids or junctions.
- **Same rule class, with excluded volume:** Negi, Winkler & Gompper, Soft Matter (2022), "Emergent collective behavior of active Brownian particles with visual perception."
  - Particles respond with limited maneuverability to the positions of neighbours inside a vision cone and a cut-off radius.
  - The structures they report are worms, worm-aggregate coexistence and hexagonally close-packed structures.
  - The same group's roadmap has a section on this model class: "The 2025 motile active matter roadmap," J. Phys.: Condens. Matter, section 15. It notes that the vision-induced steering is non-additive and non-reciprocal. I am treating it as the review.
- **Closest in spirit:** Brigatti & Peruani, Soft Matter 22(13), 2474–2481 (2026), "Zoology of collective patterns modulated by non-reciprocal, long-range interactions."
  - The particles have distance-independent interactions restricted only by a field of view of angle β.
  - In 2D they find a closed nematic ring filament, closed and open polar filaments, and a structure with two singular topological points.
  - They also find an absence of gas phases, and a close link between topology and transport.
  - I could only read the abstract. I could not check whether it reports network-like states, fractal dimension, void scaling or junction statistics.

**3. What differs, and whether it matters.**
- **Interaction range:** The simulated system has a finite range R of 5 to 10. The 2026 paper has unlimited range, and Negi et al. use a cutoff but with different weighting.
- **Collisions:** The simulation has none. Negi et al. use excluded volume, which probably drives their dense clusters.
- **Noise and speed:** The simulation has weak rotational noise, a slow speed and a "mean sine" steering law. Other papers use different torque normalisations.
- **Reported states:** The prior work reports clusters, worms, files, nematic structures and filaments. I found no report of a space-spanning, fractal, branching stream network with a void-area scaling of R³ and a degree-3 junction majority.
- **Likely effect:** The network state is plausibly a distinct regime of the same model class, reached at finite R, low density, no collisions and weak noise. It is not obviously a different phenomenon.

**4. Searches and verification.**
- I verified the titles, authors, years and venues of Barberis & Peruani (2016), Negi et al. (2022) and Brigatti & Peruani (2026) against the search results. Negi et al. is confirmed by title and abstract. The year comes from the search result's age metadata and the publisher URL, not from a visible date. The roadmap's title and section appear in the results, but I did not confirm its authors or exact year.
- I ran about 10 searches. Their themes were:
  - vision cone, non-reciprocal and fractal network terms;
  - bearing-only steering, filament networks and junctions;
  - Barberis & Peruani details;
  - pursuit steering and web-like networks;
  - cosmic-web analogies;
  - the Brigatti & Peruani paper.
- I found no source describing fractal dimension, void-area scaling or junction-degree statistics for this system.
- I could not read full texts for the Barberis & Peruani or Brigatti & Peruani papers, and I could not search further. A network phase may be shown in a figure or supplement without being named in the abstract.

**5. Confidence: low to medium.** The model class and its non-reciprocal vision-cone mechanism are well established. The specific quantitative network statistics are unconfirmed in either direction. A close reading of Brigatti & Peruani (2026) and Barberis & Peruani (2016) could upgrade this to KNOWN. If those papers show no network state, the verdict would move toward LITERATURE-SILENT for the quantitative features.

## Reply R44

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Minimal cognitive flocking model (Barberis-Peruani): vision-cone, position-based non-reciprocal attraction producing polar files / filaments / nematic band networks", "notes": "Rules essentially match the Barberis-Peruani 2016 model; filament/file/band patterns are documented there. The undocumented part is only the quantitative network measurements: D~1.5, degree-3 junctions dominating, voids scaling as R^3. Degree-3 dominance and D between 1 and 2 are generic for 2-D stream networks. The R^3 exponent is fitted over only a factor-2 range in R in a small box and may be a finite-size effect. Full texts of the 2016 and 2026 papers were not read; if they show the network explicitly, the verdict becomes KNOWN."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The first pass's verdict stands, but I'd put less weight on the twist it flagged. This rule set is essentially the Barberis & Peruani 2016 model, and a network of streams is a documented feature of that model family. What I couldn't find is the specific measurements: D≈1.5, voids scaling as R³, and junctions of degree 3 to 4. Those look like a description of known behaviour, not a new phenomenon.

**The model is already published.** The 2016 paper studies particles with no memory that interact by a short-ranged, position-based, attractive force, which acts inside a vision cone (VC), and lack velocity-velocity alignment. That matches all four rules here: constant speed, a vision cone with a finite range, turning toward the positions of seen particles, and noise. Neither version has collisions.

**What the paper reports seeing.** Depending on parameters, the emergence of aggregates or millinglike patterns, the formation of moving—locally polar—files with particles at the front of these structures acting as effective leaders, and the self-organization of particles into macroscopic nematic structures. A later paper sums up this class: such particles aggregate into clusters, polar filaments, and nematic bands. A web of thin streams in a very narrow cone (π/8), with weak noise and point particles that never collide, fits these filament and band regimes.

**Related work by the same group:**
- Saavedra & Peruani (2024): vision-cone attraction makes particles form particle chains, which eventually become closed loops. That study includes obstacles, which this run doesn't.
- Brigatti & Peruani (2026) found a nematic closed filament in the form of a ring, polar filaments with topological defects, and open polar filaments. Their interactions have no distance limit, though, so they don't explain how the network spacing depends on R.

**Why the twist probably doesn't matter much:**
- **Degree-3 junctions:** Y-shaped junctions dominate in almost any 2-D network formed by streams merging. Merging three streams is far more likely than four, and five or more almost never happens. This is expected, not distinctive.
- **D≈1.5:** box-counting a line network gives values between 1 and 2 when the mesh size is close to the measurement window. That makes it a crossover number, not real fractal scaling.
- **Voids scaling as R³:** this is the only odd number, since simple geometry predicts R². But it was fitted over only R = 5 to 10, a factor of 2, in a 100×100 box. At R=10 the voids are probably limited by the box size. I would treat it as unreliable unless it holds across a wider range of R and a larger box.

**Real pattern or artifact?** It looks like a real emergent pattern of filaments and junctions, not random scatter, a numerical blow-up or a frozen state. One caveat: because nothing pushes the particles apart, streams can collapse into lines of essentially zero width, which inflates line-like statistics.

**Limits of this check:** I didn't read the full texts of the 2016 PRL or the Brigatti paper. If their figures show a branching network explicitly, the verdict would move to KNOWN.

## Reply R58

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "No source reports fractal (D~1.5) branching stream networks, void-area ~ R^3 scaling, or degree-3-dominated junction statistics for vision-cone attraction models; the R^3 exponent implies a second length scale (turning radius or interparticle spacing) or unfinished coarsening.", "standard_name": "Minimal cognitive flocking / vision-cone position-based (non-reciprocal) attraction model (Barberis & Peruani 2016); worm/polar-file and interconnected nematic-band regimes; iABP family (Negi, Winkler, Gompper)", "notes": "Same model class as Barberis-Peruani. Files led by incidental leaders, and small nematic bands that interconnect during a transient, are documented; 3-way merging junctions are expected from follower convergence. The twist is the quantitative network statistics; the R^3 void scaling plausibly matters but needs finite-size and run-length checks to exclude a slow coarsening transient. Full texts of Barberis-Peruani and Brigatti-Peruani 2026 were not read."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass, though for a sharper reason. Merging streams and interconnected bands are documented for this exact class of model. What I couldn't find anywhere is the network statistics you measured: a fractal dimension of about 1.5, void area scaling as R³, and junctions that are mostly 3-way.

**Your system is the documented model.** Barberis & Peruani (PRL 2016) studied particles with no memory, that interact by a short-ranged, position-based, attractive force, which acts inside a vision cone (VC), and lack velocity-velocity alignment. That matches your rules (2)–(4): no heading-copying, one-sided vision, steering toward bearings. Their phase diagram plots vision angle β vs angular noise intensity D_θ, with i) gas phase, ii) aggregate phase, iii) worm phase, and iv) nematic phase. They report the formation of moving—locally polar—files with particles at the front of these structures acting as effective leaders. A narrow cone (half-angle π/8) and weak noise (0.08) should put your run in the files/"worm" region rather than the bands region, which they place at larger values of D_θ and β. I'm inferring that placement because I didn't read their numbers.

**Network-like structure is partly documented.** In the nematic phase, after a complex transient where various small nematic bands grow in size and interconnect, the system reaches a steady state, with one or several bands. Related work also describes particle chains: active particles with nonreciprocal attraction form particle chains, which eventually become closed loops (here only around obstacles). Point-particle vision-cone models also show a quasi-periodic, honeycomb-lattice-like structure in binary mixtures. And the single-species point-particle model gives clusters, polar filaments, and nematic bands.

**What I couldn't find (the twist):**
- **Fractal dimension ~1.5 of the whole pattern.** Nobody reports a fractal measure for these structures.
- **Void area ∝ R³.** This is not reported, and it is odd. If R were the only length scale, void area should go as R². An exponent of 3 means another length matters, probably the turning radius (speed/turn rate), the spacing between particles, or a pattern that hasn't finished coarsening. This plausibly matters, but only if it survives a test across system sizes and run lengths.
- **Mostly 3-way junctions, 5+ absent.** This is unreported but expected: followers converging on leaders' files naturally give tree-like Y-merges. I see it as a descriptor of known behaviour, not a new phase.

**Genuine or artifact?** Probably a genuine emergent pattern, not trivial. Without collisions there is no excluded volume, so streams can collapse onto lines without limit, which makes D < 2 easy to get. The thing to check is whether the network is a steady state or a slow transient on the way to the documented bands or aggregates. Barberis & Peruani saw an interconnecting network only as a transient in the bands regime. The run (50,000 time units, roughly 5,000 length units travelled per particle) is long, but I can't rule out that the network is still coarsening.

**Gaps in my search:** I ran 10 searches covering networks, fractals, junctions, voids, branched filaments, lanes and the Brigatti & Peruani paper. I couldn't locate that 2026 paper, and I only saw snippets of Barberis & Peruani, not the full text or supplement. Either could still contain a network regime. Negi et al. (Soft Matter 2022, NJP 2025) and the PRResearch 2024 alignment-plus-vision paper report worms, clusters and honeycomb structures, but no fractal or junction statistics.

# Feature check input

## Feature definitions (census fingerprint, computed from one simulation run)

- f_std: standard deviation of the field over the grid (last frame)
- f_mean: mean of the field (last frame)
- f_peak_k: wavenumber of the peak of the radial power spectrum
- f_peak_sharp: share of spectral power in that peak (sharp pattern)
- f_n_blobs: log(1 + number of connected regions above mean + 1 sd)
- f_change: mean absolute change of the field per recorded frame, late run
- f_moran: Moran's I spatial autocorrelation of the field (spatial structure)
- a_local_align: local heading alignment within distance 1.5, minus the random expectation
- a_polar: global polar order |mean heading vector|, late run
- a_polar_std: fluctuation of global polar order
- a_nn_ratio: mean nearest-neighbour distance / random expectation (below 1 = clustered)
- a_density_cv: coefficient of variation of counts on a 16 x 16 grid (density inhomogeneity)
- a_ang_mom: crude global angular momentum about the centroid
- a_step: mean displacement per recorded frame
- p_r: global phase order |mean exp(i theta)|, late run
- p_r_std: fluctuation of global phase order
- p_local_r: local phase order over 4 lattice neighbours
- p_vortex_density: share of plaquettes with a phase winding (vortices)
- p_freq_spread: spread of mean phase velocities

## Card K01

**Rule system:** World: two continuous concentration fields u and v on a 64x64 square grid with wrap-around edges, started at their uniform steady state (u=b, v=b/(a+b^2)) plus small random noise. Rules, integrated forward in time: u changes at rate -u + a*v + u^2*v + D1*(Laplacian of u) + c1*div(u grad v); v changes at rate b - a*v - u^2*v + D2*(Laplacian of v) + c2*div(v grad u). The last terms make each species drift along the other's gradient, with strength and sign set by c1 and c2. Parameters: a=0.015, b=10, D1=0.01, D2=1000*D1. Setting: c1=c and c2=-c (opposite signs, so one species is drawn toward the other while the other is pushed away), c=0.17.

**Features on the card (value):** f_std = 13.029, f_mean = 10.0, f_n_blobs = 2.614, f_moran = 0.967, f_peak_sharp = 0.768, f_peak_k = 0.047

**What the source paper reports as new about this system:** In a two-species reaction-diffusion (glycolysis-type) model with mutual chemotaxis, steady Turing pattern amplitude varies non-monotonically with chemotactic strength (re-entrant uniform/patterned transitions) and spots turn into stripes as chemotaxis grows.

**How the paper measures it:** pattern amplitude at c = 0.01 versus absent at c = 0.49 (non-monotonic in chemotactic strength); spots turning into stripes as chemotaxis grows

## Card K02

**Rule system:** World: two continuous concentration fields u and v on a 64x64 square grid with wrap-around edges, started at their uniform steady state (u=b, v=b/(a+b^2)) plus small random noise. Rules, integrated forward in time: u changes at rate -u + a*v + u^2*v + D1*(Laplacian of u) + c1*div(u grad v); v changes at rate b - a*v - u^2*v + D2*(Laplacian of v) + c2*div(v grad u). The last terms make each species drift along the other's gradient, with strength and sign set by c1 and c2. Parameters: a=0.015, b=10, D1=0.01, D2=1000*D1. Setting: c1=c and c2=-c (opposite signs, so one species is drawn toward the other while the other is pushed away), c=0.17.

**Features on the card (value):** f_n_blobs = 2.534, f_moran = 0.975, f_peak_sharp = 0.953, f_mean = 0.08, f_peak_k = 0.047, f_std = 0.024

**What the source paper reports as new about this system:** In a two-species reaction-diffusion (glycolysis-type) model with mutual chemotaxis, steady Turing pattern amplitude varies non-monotonically with chemotactic strength (re-entrant uniform/patterned transitions) and spots turn into stripes as chemotaxis grows.

**How the paper measures it:** pattern amplitude at c = 0.01 versus absent at c = 0.49 (non-monotonic in chemotactic strength); spots turning into stripes as chemotaxis grows

## Card K03

**Rule system:** World: 400 self-propelled disks (radius 1) in a periodic 2D box (square of side 158.5) at area fraction 0.05, hovering at height 1.1 above a flat wall that chemicals cannot cross. Each disk emits a chemical. Rules each step: (1) each disk moves along its heading at speed 2; (2) it is pulled toward the others along the summed gradient of their screened chemical fields (each decays as exp(-k*r)/r, plus a mirror-image source below the wall), screening k = 0.20; (3) its heading is turned away from that same gradient with turning strength 0.30; (4) the heading gets rotational noise (Peclet number 100); (5) disks repel on contact. A cluster is more than 10 disks linked within 2.1 radii. Run 100 tau (tau = 1/rotational diffusion rate, the heading's persistence time).

**Features on the card (value):** a_step = 5.054, a_density_cv = 3.662, a_nn_ratio = 0.482, a_polar = 0.019, a_ang_mom = 0.018, a_polar_std = 0.007

**What the source paper reports as new about this system:** Autophoretic disks with chemo-attractive forces and chemo-repulsive torques above a no-flux wall form dense clusters filled with circulating cavities (boiling bubbles) or clusters that repeatedly rupture (bursting bubbles), between full collapse and dynamic clustering.

**How the paper measures it:** cluster fraction, and cavities (circulating holes) inside the clusters

## Card K04

**Rule system:** World: 400 slightly elongated particles (ellipse half-axes 0.55 and 0.65) in a periodic 2D box at density 0.87 (lattice spacing about 1.15). Rules, overdamped Brownian dynamics at temperature T: (1) each pair interacts through an orientation-dependent Lennard-Jones potential of strength e_p (repulsive core, weak attraction, depending on both particles' orientations, which can rotate); (2) each pair also feels a non-conservative force with Lennard-Jones magnitude and strength e_a, directed perpendicular to the line joining them, always in the same rotational sense; (3) thermal noise acts on positions and orientations, orientations relaxing faster. Set e_a/e_p = 33 and T = 6 and run for 30 time units; measure six-fold bond order and defect counts.

**Features on the card (value):** a_nn_ratio = 1.975, a_density_cv = 0.433, a_step = 0.161, a_polar = 0.04, a_local_align = -0.036, a_polar_std = 0.021

**What the source paper reports as new about this system:** A 2D solid with transverse (odd) pair forces melts by defect-pair proliferation at any temperature above a threshold activity, and recrystallises at higher activity even at temperatures where the passive crystal is liquid.

**How the paper measures it:** bond-orientational (hexatic) order at the headline activity compared with the passive system at the same temperature

## Card K05

**Rule system:** World: a 64x64 square lattice with wrap-around edges; each site holds a continuous mass (mean 1 per site). Rules, repeated: (1) pick a random site and one of the four directions at random; (2) the site keeps 40% of its mass plus a random share of the rest and sends the remainder one step in that direction; (3) the receiving site pools it with its own mass, keeps 60% plus a random share of the rest, and passes the remainder one step to its left (a 90-degree counter-clockwise turn), so mass moves along L-shaped paths that always turn the same way. Rates are identical for all four directions and total mass is conserved. Measure steady-state density-density correlations and long-wavelength fluctuations.

**Features on the card (value):** f_n_blobs = 5.92, f_mean = 1.0, f_change = 0.662, f_std = 0.643, f_peak_k = 0.417, f_moran = 0.102

**What the source paper reports as new about this system:** Isotropic mass-transfer rules with a fixed rotational sense (L-shaped chiral moves) produce r^-4 power-law density correlations in the steady state, with long-wavelength fluctuations minimised at an optimal chirality and a square-root cusp near maximal odd diffusion.

**How the paper measures it:** density-density correlations decaying close to r^-4 (fitted exponent 3.5 to 4.5) where the achiral version does not

## Card K06

**Rule system:** World: individuals at continuous positions in a square periodic domain, starting with 4,000 placed at random. Rules (continuous-time birth and death): (1) each individual gives birth at rate 1920; the offspring lands at a Gaussian offset (standard deviation 0.05); (2) each individual dies at rate 192 + (density of others within radius 1.25) + g times (density within radius 0.0625, squared: a three-body crowding term). Densities are counts divided by disk area. Only competitive (harmful) interactions, no attraction or help. Set g = 3.14e-5 and run for 1 time unit; look at the steady spatial pattern.

**Features on the card (value):** f_mean = 39.75, f_std = 14.748, f_change = 12.078, f_n_blobs = 5.305, f_moran = 0.564, f_peak_sharp = 0.427

**What the source paper reports as new about this system:** Adding a three-body competition term to a purely competitive spatial birth-death model breaks the hexagonal-spots-only restriction, giving stripes and gaps and an unusual spots-stripes-spots sequence.

**How the paper measures it:** pattern type: stripes and gaps (anisotropy index) versus hexagonal spots

## Card K07

**Rule system:** World: a 64 x 64 square lattice with wrap-around edges, each site holding a planar unit spin (an angle). Energy: minus the sum over all pairs of J(r) times the cosine of their angle difference, where every pair interacts with J(r) proportional to 1/r^4 (r = distance), normalised so that the total coupling per site is fixed. Rules: equilibrium Monte Carlo with single-spin updates (one sweep of all sites per step) at inverse temperature beta = 1, from random angles, for 1000 sweeps. Measure the squared magnetisation and the correlation length.

**Features on the card (value):** p_local_r = 0.855, p_r = 0.727, p_freq_spread = 0.072, p_r_std = 0.016, p_vortex_density = 0.006

**What the source paper reports as new about this system:** The 2D XY model with 1/r^4 couplings has a low-temperature phase whose correlations decay as a power of log r, entered through a continuous transition with power-law divergent correlation length (neither BKT nor Landau type).

**How the paper measures it:** squared magnetization decreasing slowly (as a power of log L) with system size L

## Card K08

**Rule system:** World: 400 point particles in a 20x20 periodic box (density 1), with no collisions. Rules: (1) each moves forward at constant speed 0.1 along its heading (small inertia); (2) it sees only particles within distance R = 5 inside a vision cone of half-angle pi/8 around its heading; (3) it turns toward the positions of the visible particles, turning rate 3 times the average sine of the angle between its heading and the bearing to each (no copying of their headings); (4) its heading gets rotational noise of strength 0.08. Seeing is one-sided: being seen does not imply seeing back. Run 200,000 steps of size 0.005.

**Features on the card (value):** a_density_cv = 2.134, a_step = 1.855, a_nn_ratio = 0.437, a_polar = 0.071, a_polar_std = 0.034, a_ang_mom = 0.028

**What the source paper reports as new about this system:** Self-propelled particles that only steer toward neighbours inside a narrow finite-range vision cone form a system-spanning, continuously restructuring fractal-like network of particle streams (dimension ~1.5) that never coarsens.

**How the paper measures it:** a system-spanning stream network with fractal dimension near 1.5 that does not coarsen

## Card K09

**Rule system:** World: 400 self-propelled point particles in a periodic 2D box (square of side 20), density 1 per unit area. Rules, discrete time steps: (1) each particle sets its heading to the average heading of all particles within distance 1 (itself included), plus a uniformly random angle in [-eta*pi, eta*pi], with noise eta = 0.318 (a random angle between -1 and 1 radian); (2) it moves 0.5 along its new heading. Start from random positions and random headings along the four lattice directions. Run 1,000 steps.

**Features on the card (value):** a_step = 7.487, a_density_cv = 1.221, a_nn_ratio = 0.855, a_polar = 0.644, a_local_align = 0.463, a_ang_mom = 0.048

**What the source paper reports as new about this system:** In the coexistence phase of the 2D Vicsek model, side-by-side lanes of band trains travelling in opposite directions (and fat/inverse bands) form at particle level; long lanes buckle and the ordered band states are metastable.

**How the paper measures it:** side-by-side lanes of band trains moving in opposite directions

## Card K10

**Rule system:** World: three continuous fields on a 64x64 periodic square: cell density b, an attractant a made by the cells, and a nutrient n eaten by the cells; all start uniform at 1 plus small noise. Rules (dimensionless): (1) b diffuses and drifts up the gradient of a with strength 11.1*b/(1+b^2) and up the gradient of n with strength 14*b/(1+b^2); cell number is conserved; (2) a diffuses 2.5 times faster than cells, is produced at rate b and decays at rate 1; (3) n diffuses 0.75 times as fast as cells, is consumed at rate 0.084*b*n and resupplied at rate 0.084 everywhere.

**Features on the card (value):** f_std = 2.336, f_n_blobs = 1.29, f_mean = 1.0, f_moran = 0.977, f_peak_sharp = 0.373, f_change = 0.343

**What the source paper reports as new about this system:** A Keller-Segel model with chemotaxis toward a self-made attractant and a consumed nutrient produces a dense aggregate coexisting with travelling bands of cells.

**How the paper measures it:** one dense aggregate coexisting with travelling bands of cells

## Card K11

**Rule system:** World: three continuous fields on a 64x64 periodic square: cell density b, an attractant a made by the cells, and a nutrient n eaten by the cells; all start uniform at 1 plus small noise. Rules (dimensionless): (1) b diffuses and drifts up the gradient of a with strength 11.1*b/(1+b^2) and up the gradient of n with strength 14*b/(1+b^2); cell number is conserved; (2) a diffuses 2.5 times faster than cells, is produced at rate b and decays at rate 1; (3) n diffuses 0.75 times as fast as cells, is consumed at rate 0.084*b*n and resupplied at rate 0.084 everywhere.

**Features on the card (value):** f_std = 1.993, f_n_blobs = 1.29, f_mean = 1.0, f_moran = 0.984, f_peak_sharp = 0.402, f_change = 0.279

**What the source paper reports as new about this system:** A Keller-Segel model with chemotaxis toward a self-made attractant and a consumed nutrient produces a dense aggregate coexisting with travelling bands of cells.

**How the paper measures it:** one dense aggregate coexisting with travelling bands of cells

## Card K12

**Rule system:** World: three continuous fields on a 64x64 periodic square: cell density b, an attractant a made by the cells, and a nutrient n eaten by the cells; all start uniform at 1 plus small noise. Rules (dimensionless): (1) b diffuses and drifts up the gradient of a with strength 11.1*b/(1+b^2) and up the gradient of n with strength 14*b/(1+b^2); cell number is conserved; (2) a diffuses 2.5 times faster than cells, is produced at rate b and decays at rate 1; (3) n diffuses 0.75 times as fast as cells, is consumed at rate 0.084*b*n and resupplied at rate 0.084 everywhere.

**Features on the card (value):** f_mean = 2.538, f_n_blobs = 1.33, f_moran = 0.987, f_std = 0.795, f_peak_sharp = 0.465, f_change = 0.354

**What the source paper reports as new about this system:** A Keller-Segel model with chemotaxis toward a self-made attractant and a consumed nutrient produces a dense aggregate coexisting with travelling bands of cells.

**How the paper measures it:** one dense aggregate coexisting with travelling bands of cells

"""N15 (planted card): self-propelled point particles that steer toward the neighbours they can see inside a finite-range
vision cone (arXiv 2610.15835, Sec. II "Methods").

Model, as the paper states it (Sec. II, Eqs. 1-4, "Simulation Details")
  m r_i'' = -gamma r_i' + gamma v0 e(phi_i),            e(phi) = (cos phi, sin phi)                     (Eq. 1)
  phi_i'  = (Omega / n_i) sum_{j in VC_i} sin(alpha_ij - phi_i) + sqrt(2 Dr) xi_i(t)                   (Eq. 2)
  alpha_ij = polar angle of r_ij = r_j - r_i; n_i = number of particles in i's vision cone VC_i, where
  j in VC_i  iff  |r_j - r_i| <= Rv  (Eq. 3)  and  (r_j - r_i)/|r_j - r_i| . e(phi_i) > cos(theta)  (Eq. 4).
  Point particles, no steric exclusion, no velocity alignment; periodic square box of side L; density Phi = N/L^2 = 1.
  Parameters (Sec. II): m = 1, gamma = 50 (tau_I = m/gamma = 0.02), v0 = 0.1, Dr = 0.08, dt = 5e-3, velocity-Verlet
  integration of Eq. (1); Omega = 3 (Fig. 1, Fig. 2); N = 1e4, L = 100; runs of 1e7 steps, last 5e6 averaged.
  Headline (plan row N15; Fig. 1(c) "active foams (theta = pi/8, Rv = 5)" and the foam region of Fig. 2):
  theta = pi/8, Rv = 5, Omega = 3, at census size 400 particles in a 20 x 20 box (density 1).
  Measure (Sec. II "Correlation (fractal) dimension", Eqs. 5-6): C(r) = fraction of particle pairs closer than r,
  D = slope of log C(r) against log r in its linear regime; the paper reports D ~ 1.5 for the foam (Sec. III B).
  Eq. (6) is printed as the ratio D = log C(r) / log r, but the text defines D as the slope of the linear regime of the
  log-log plot (C(r) ~ r^D); this module uses the slope (least-squares fit of log C on log r).

Census mapping (PLANTED_PIPELINE_PLAN.md, exception for N15): one census step = one time unit = 200 integration steps of
dt = 0.005. 1000 census steps = 2e5 integration steps (the paper runs 1e7). Frames every 2 census steps.
Rules in plan numbering: 1 propulsion with small inertia (one-way), 2 steering toward visible particles (interaction),
3 rotational noise (one-way). Knock-out (keep=()) removes rule 2 (Omega = 0): independent inertial active particles.

Choices the paper does not state, and why
  * Velocity-Verlet with a velocity-dependent force follows the LAMMPS convention: half kick with F(t), drift, then the
    force at t+dt is evaluated with the half-step velocity (friction) and the new heading, second half kick. The heading
    (Eq. 2) is advanced by Euler-Maruyama with the steering term evaluated at the new positions and the old headings.
  * Initial condition: uniform random positions and headings (the paper uses random or square-lattice starts, Sec. III E);
    initial velocities v0 e(phi), i.e. on the force-free manifold.
  * A particle with no visible neighbour (n_i = 0) receives no steering. A coincident pair (r = 0) has no bearing and is
    ignored (point particles can coincide; Eq. 4 is undefined there).
  * Neighbour search: Verlet pair list of radius Rv + 0.5 rebuilt whenever any particle has moved more than 0.25 since
    the last build (exact: the cone test of Eqs. 3-4 is applied to every listed pair at every step).
  * Correlation-dimension fit range (not stated in the paper): log-spaced r from 1 to min(2 Rv, L/4), i.e. [1, 10] at
    L = 100 and [1, 5] in the census box; beyond L/4 the periodic box bends C(r).
  * Headings are recorded modulo 2 pi; positions modulo L.

Small-scale observation while building (not a fidelity result): at census size (one seed, 1000 time units) particle
streams with three-fold junctions span the 20 x 20 box from ~200 time units on; correlation dimension over r in [1, 5]
falls from 1.96 (random start) to 1.33-1.43 and the structure length stays ~4.5 (no growth after t = 100). At 10,000
particles the uniform start costs ~2.9 s per time unit here (more once streams concentrate the particles).

Fidelity check (plan: "a non-coarsening stream network with fractal dimension near 1.5"), fixed before any large run:
  10,000 particles in a 100 x 100 box at the headline (theta = pi/8, Rv = 5, Omega = 3, Dr = 0.08, v0 = 0.1, dt = 0.005),
  run for as long as the budget allows (the paper runs 5e4 time units; this cannot fit 2 CPU-hours in numpy). Measures on
  the last half of the run: (i) correlation dimension D (the paper's Eqs. 5-6), tolerance 1.3 <= D <= 1.7 (paper ~1.5);
  (ii) non-coarsening: the structure length l = 2 pi / <k> (first moment of the structure factor of the particle density
  on a 1 x 1 grid, k > 0) in the last quarter of the run is at most 1.25 times its value in the third quarter;
  (iii) descriptive: box-counting dimension over box sides 1 to 10 (the parent brief names box counting; the paper's own
  measure is the correlation dimension, which decides pass/fail), and the share of occupied grid cells in the largest
  connected occupied region (system-spanning network). Pass = (i) and (ii).
"""
import time
import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage
from census.planted_models import base

CID = "planted:N15"
LAYERS = "A"
K = {"A": 1}
SPEED = 0.1
RULES = [("propulsion_with_inertia", False), ("steer_toward_visible_in_cone", True), ("rotational_noise", False)]
HEADLINE = {
    "N": 400, "density": 1.0, "box": 20.0,          # plan row N15 (paper: N = 1e4, L = 100, Phi = 1; Sec. II)
    "v0": 0.1, "m": 1.0, "gamma": 50.0,             # Sec. II "Simulation Details" (tau_I = 0.02)
    "Dr": 0.08,                                     # Sec. II (persistence time 12.5)
    "Omega": 3.0,                                   # Fig. 1, Fig. 2 (visual torque strength)
    "Rv": 5.0, "theta": float(np.pi / 8),           # Fig. 1(c) active foam: theta = pi/8, Rv = 5
    "dt": 0.005, "steps_per_census_step": 200,      # Sec. II dt; plan exception (one time unit per census step)
}
STEP = "one census step = 1 time unit = 200 velocity-Verlet steps of dt = 0.005 (plan exception for N15)"
RULE_TEXT = ("World: 400 point particles in a 20x20 periodic box (density 1), with no collisions. Rules: (1) each moves "
             "forward at constant speed 0.1 along its heading (small inertia); (2) it sees only particles within distance "
             "R = 5 inside a vision cone of half-angle pi/8 around its heading; (3) it turns toward the positions of the "
             "visible particles, turning rate 3 times the average sine of the angle between its heading and the bearing to "
             "each (no copying of their headings); (4) its heading gets rotational noise of strength 0.08. Seeing is "
             "one-sided: being seen does not imply seeing back. Run 200,000 steps of size 0.005.")

SKIN = 0.5


class _Steer:
    """vision-cone steering term of Eq. (2) with a Verlet pair list (unordered pairs, minimum-image shifts frozen per build)."""

    def __init__(self, N, B, Rv, theta, Omega):
        self.N, self.B, self.Rv, self.R2, self.ct, self.Om = N, B, Rv, Rv * Rv, np.cos(theta), Omega
        self.x0 = None

    def build(self, X, Y):
        B = self.B; xw = np.mod(np.stack([X, Y], 1), B); xw[xw >= B] = 0.0
        pr = cKDTree(xw, boxsize=B).query_pairs(min(self.Rv + SKIN, B / 2 - 1e-9), output_type="ndarray")
        if len(pr) == 0: pr = np.zeros((0, 2), np.intp)
        self.i, self.j = pr[:, 0].copy(), pr[:, 1].copy()
        self.sx = B * np.round((X[self.j] - X[self.i]) / B); self.sy = B * np.round((Y[self.j] - Y[self.i]) / B)
        self.x0 = (X.copy(), Y.copy())

    def __call__(self, X, Y, c, s):
        if self.x0 is None or max(np.abs(X - self.x0[0]).max(), np.abs(Y - self.x0[1]).max()) > SKIN / 2:
            self.build(X, Y)
        i, j, N = self.i, self.j, self.N
        if len(i) == 0: return np.zeros(N)
        dx = X[j] - X[i]; dx -= self.sx; dy = Y[j] - Y[i]; dy -= self.sy
        r2 = dx * dx; r2 += dy * dy
        r = np.sqrt(r2); r[r2 > self.R2] = np.inf                  # out of range: fails both cone tests
        ci, si, cj, sj = c[i], s[i], c[j], s[j]
        pi_ = dx * ci; pi_ += dy * si; pj = dx * cj; pj += dy * sj; cr = self.ct * r
        vi = pi_ > cr; vj = pj < -cr                               # j in VC_i ; i in VC_j (Eqs. 3-4); r = 0 excluded
        wi = (ci[vi] * dy[vi] - si[vi] * dx[vi]) / r[vi]           # sin(alpha_ij - phi_i)
        wj = (sj[vj] * dx[vj] - cj[vj] * dy[vj]) / r[vj]           # sin(alpha_ji - phi_j), r_ji = -r_ij
        ii, jj = i[vi], j[vj]
        ssum = np.bincount(ii, wi, N) + np.bincount(jj, wj, N)
        n = np.bincount(ii, minlength=N) + np.bincount(jj, minlength=N)
        return self.Om * ssum / np.maximum(n, 1)


def _engine(rng, N, B, T_units, steer_on, rec_every=2, Rv=5.0, theta=np.pi / 8, Omega=3.0, v0=0.1, m=1.0,
            gamma=50.0, Dr=0.08, dt=0.005, nsub=200, deadline=None, on_frame=None):
    """integrate T_units time units; returns (frames_pos, frames_head, units_done). Frames at t = 0, rec_every, ..."""
    X = rng.uniform(0, B, N); Y = rng.uniform(0, B, N); ph = rng.uniform(0, 2 * np.pi, N)
    c, s = np.cos(ph), np.sin(ph); VX, VY = v0 * c, v0 * s; FX = np.zeros(N); FY = np.zeros(N)
    st = _Steer(N, B, Rv, theta, Omega) if steer_on else None
    P, H = [], []
    sq = np.sqrt(2 * Dr * dt); hdt = 0.5 * dt / m; done = 0

    def rec():
        xw = np.mod(np.stack([X, Y], 1), B).astype(np.float32); xw[xw >= B] = 0.0
        hh = np.mod(ph, 2 * np.pi).astype(np.float32); hh[hh >= 2 * np.pi] = 0.0
        P.append(xw); H.append(hh)
    rec()
    for t in range(1, T_units + 1):
        for _ in range(nsub):
            VX += hdt * FX; VY += hdt * FY
            X += dt * VX; Y += dt * VY
            dph = sq * rng.standard_normal(N)
            if st is not None: dph += dt * st(X, Y, c, s)
            ph += dph; c = np.cos(ph); s = np.sin(ph)
            FX = gamma * (v0 * c - VX); FY = gamma * (v0 * s - VY)
            VX += hdt * FX; VY += hdt * FY
        done = t
        if t % rec_every == 0:
            rec()
            if on_frame is not None: on_frame(t, P[-1])
        if deadline is not None and time.perf_counter() > deadline: break
    return np.array(P), np.array(H), done


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    on = base.kept(keep, RULES); NA = base.SIZES[size]["NA"]; B = float(np.sqrt(NA / 1.0))   # density 1
    nf = base.n_frames(steps, "pos")
    P = np.zeros((S, nf, NA, 2), np.float32); H = np.zeros((S, nf, NA), np.float32)
    for s, r in enumerate(base.rngs(seed0, S)):
        p, h, _ = _engine(r, NA, B, steps, 2 in on, rec_every=base.REC["pos"])
        P[s], H[s] = p, h
    bad = not (np.isfinite(P).all() and np.isfinite(H).all())
    return {"pos": P, "head": H, "at": np.zeros(H.shape, np.int8), "meta": base.meta(steps, base.SIZES[size]["W"], B, bad)}


# ------------------------------------------------------------------ measures (pure numpy/scipy on frames)
def corr_dimension(pos, B, r_lo=1.0, r_hi=None, nr=12):
    """Eqs. (5)-(6): C(r) = 2/(N(N-1)) #pairs with |ri - rj| < r (periodic), D = slope of log C vs log r."""
    r_hi = r_hi if r_hi is not None else B / 4
    r = np.geomspace(r_lo, r_hi, nr); x = np.mod(pos.astype(np.float64), B); x[x >= B] = 0.0
    tr = cKDTree(x, boxsize=B); N = len(x)
    cnt = tr.count_neighbors(tr, r).astype(float) - N           # remove self pairs; ordered pairs
    C = cnt / (N * (N - 1.0))
    ok = C > 0
    if ok.sum() < 3: return float("nan")
    return float(np.polyfit(np.log(r[ok]), np.log(C[ok]), 1)[0])


def box_dimension(pos, B, e_lo=1.0, e_hi=10.0, ne=8):
    """box-counting dimension: N_occupied(eps) ~ eps^-Db, box sides from e_lo to e_hi (descriptive only)."""
    es = np.geomspace(e_lo, e_hi, ne); n = []
    for e in es:
        g = max(1, int(round(B / e))); idx = np.floor(np.mod(pos, B) / (B / g)).astype(int) % g
        n.append(len(np.unique(idx[:, 0] * g + idx[:, 1])))
    return float(-np.polyfit(np.log(B / np.round(B / es)), np.log(n), 1)[0])


def struct_length(pos, B, cell=1.0):
    """l = 2 pi / <k>, <k> = first moment of the structure factor of the density on a cell x cell grid (k > 0)."""
    g = max(4, int(round(B / cell)))
    rho = np.histogram2d(np.mod(pos[:, 0], B), np.mod(pos[:, 1], B), bins=g, range=[[0, B], [0, B]])[0]
    rho = rho - rho.mean(); Sk = np.abs(np.fft.fft2(rho)) ** 2
    k1 = 2 * np.pi * np.fft.fftfreq(g, d=B / g); kk = np.sqrt(k1[:, None] ** 2 + k1[None, :] ** 2)
    m = (kk > 0) & (kk <= np.pi / (B / g))
    return float(2 * np.pi / ((kk[m] * Sk[m]).sum() / Sk[m].sum()))


def spanning_share(pos, B, cell=1.0):
    """share of occupied cells (>= 1 particle, 1 x 1 cells) in the largest 8-connected occupied region (periodic)."""
    g = max(4, int(round(B / cell)))
    occ = np.histogram2d(np.mod(pos[:, 0], B), np.mod(pos[:, 1], B), bins=g, range=[[0, B], [0, B]])[0] > 0
    lab, n = ndimage.label(occ, structure=np.ones((3, 3), int))
    if n == 0: return 0.0
    par = np.arange(n + 1)

    def find(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for a, b in list(zip(lab[0], lab[-1])) + list(zip(lab[:, 0], lab[:, -1])):
        if a and b:
            ra, rb = find(a), find(b)
            if ra != rb: par[ra] = rb
    roots = np.array([find(a) for a in range(n + 1)]); sizes = np.bincount(roots[lab[occ]])
    return float(sizes.max() / occ.sum())


def _measure(frames, B, Rv=5.0):
    r_hi = min(2 * Rv, B / 4)
    D = [corr_dimension(p, B, 1.0, r_hi) for p in frames]
    Db = [box_dimension(p, B, 1.0, min(10.0, B / 4)) for p in frames]
    ell = [struct_length(p, B) for p in frames]; sp = [spanning_share(p, B) for p in frames]
    return np.array(D), np.array(Db), np.array(ell), np.array(sp)


def presence(out):
    """correlation dimension (fit r in [1, 5] in the census box), box-counting dimension, structure length growth between
    the third and last quarter of the run, and spanning share, on the last half of each seed's frames."""
    B = float(out["meta"]["WA"]); res = []
    for s in range(out["pos"].shape[0]):
        P = out["pos"][s]; nf = len(P); half = P[nf // 2:][::5]
        D, Db, ell, sp = _measure(half, B)
        q3 = [struct_length(p, B) for p in P[nf // 2: 3 * nf // 4][::5]]; q4 = [struct_length(p, B) for p in P[3 * nf // 4:][::5]]
        res.append({"corr_dim": round(float(np.nanmean(D)), 3), "box_dim": round(float(np.mean(Db)), 3),
                    "length_ratio_q4_q3": round(float(np.mean(q4) / np.mean(q3)), 3),
                    "spanning_share": round(float(np.mean(sp)), 3)})
    Dm = np.nanmean([r["corr_dim"] for r in res]); lr = np.mean([r["length_ratio_q4_q3"] for r in res])
    return {"per_seed": res, "corr_dim_mean": round(float(Dm), 3), "fractal_network": bool(1.3 <= Dm <= 1.7),
            "non_coarsening": bool(lr <= 1.25),
            "note": "census box 20 = 4 Rv; the paper's foam voids reach radius ~25 (Sec. III E), larger than the box"}


def fidelity(budget_s=7200, N=10000, B=100.0, seed=None):
    """the plan's check at the paper's size (N = 1e4, L = 100) for as long as the budget allows (see module docstring)."""
    t0 = time.perf_counter(); seed = base.seed_of(CID) if seed is None else seed
    rng = base.rngs(seed, 1)[0]
    deadline = t0 + 0.92 * budget_s                       # leave time for the analysis
    P, H, done = _engine(rng, N, B, 50000, True, rec_every=10, deadline=deadline)   # paper: 5e4 time units
    nf = len(P)
    half = P[nf // 2:]; sel = half[np.linspace(0, len(half) - 1, min(20, len(half))).astype(int)]
    D, Db, ell, sp = _measure(sel, B)
    q3 = P[nf // 2: 3 * nf // 4]; q4 = P[3 * nf // 4:]
    pick = lambda A: A[np.linspace(0, len(A) - 1, min(10, len(A))).astype(int)] if len(A) else A
    l3 = np.mean([struct_length(p, B) for p in pick(q3)]) if len(q3) else np.nan
    l4 = np.mean([struct_length(p, B) for p in pick(q4)]) if len(q4) else np.nan
    Dm = float(np.nanmean(D)); ratio = float(l4 / l3)
    ok_D = 1.3 <= Dm <= 1.7; ok_nc = ratio <= 1.25
    return {"setting": {"N": N, "box": B, "Rv": 5.0, "theta": "pi/8", "Omega": 3.0, "Dr": 0.08, "v0": 0.1, "dt": 0.005},
            "time_units_run": int(done), "integration_steps": int(done * 200), "frames": int(nf),
            "measured": {"corr_dim_last_half": round(Dm, 3), "corr_dim_sd": round(float(np.nanstd(D)), 3),
                         "box_dim_last_half": round(float(np.mean(Db)), 3),
                         "struct_length_q3": round(float(l3), 2), "struct_length_q4": round(float(l4), 2),
                         "length_ratio_q4_q3": round(ratio, 3), "spanning_share": round(float(np.mean(sp)), 3)},
            "paper": {"corr_dim": 1.5, "non_coarsening": True, "run_time_units": 50000},
            "tolerance": "1.3 <= D <= 1.7 and l(q4)/l(q3) <= 1.25",
            "pass": bool(ok_D and ok_nc), "elapsed_s": round(time.perf_counter() - t0, 1)}

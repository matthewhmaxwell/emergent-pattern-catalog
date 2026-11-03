"""N08 (planted card): purely repulsive inertial dumbbells driven into rotation by a constant torque, all in the same sense
(arXiv 2609.21556).

Model, as the paper states it (main text, Eq. 1 and the paragraphs around it; End Matter Appendix A)
  m r_i'' = -gamma r_i' - grad_i U + F_i,a + sqrt(2 gamma kB T) xi_i                                        (Eq. 1)
  beads i = 1..N, N = 2 Nd; identical beads of mass m and diameter sigma; U = U_WCA + U_FENE;
  U_WCA(r) = 4 eps [(sigma/r)^12 - (sigma/r)^6] (+ eps), cut at its minimum r_c = 2^(1/6) sigma, between every pair of
  beads; U_FENE(r) = K R0^2 ln[1 - (r/R0)^2] between the two beads of a dumbbell (LAMMPS form -K R0^2/2 ln(.) used,
  see below), K = 30, R0 = 1.5; F_i,a has fixed magnitude Fa, perpendicular to the dumbbell axis, opposite on the two
  beads (torque tau_a = sigma Fa, clockwise for every dumbbell). Units m = sigma = eps = kB = 1; T = 1e-2; gamma = 1e-2
  ("we study the phase diagram at fixed gamma = 1e-2, hence Pr = 5e-3"). Definitions: Dr = 2 kB T/(gamma sigma^2) = 2,
  Pe = sigma Fa/(2 kB T), Re = m Fa/(sigma gamma^2), Pr = Pe/Re = gamma/(m Dr) = 5e-3. Area fraction
  phi = N pi sigma^2/(4 L^2). LAMMPS molecular dynamics.
  Headline (plan row N08; Fig. 1(a) and Fig. 2): phi = 0.65, Pe = 5 (Fa = 0.1), T = 0.01, gamma = 0.01, m = 1 (Pr = 5e-3,
  Re = 1000). Inside the binodal: Fig. 2 places the Pe = 5 coexisting densities at about 0.55 and 0.73.
  Measure (End Matter Appendix B): probability density of the local surface fraction phi_loc, the locally averaged
  surface fraction over a grid of linear size 20 sigma; two peaks phi_l(Pe), phi_h(Pe) inside the binodal (Fig. 6),
  one Gaussian peak at phi when there is no separation (Pe = 0.5 in Fig. 6(a); Pe = 0 homogeneous in Fig. 2).

Census view: one census particle = one dumbbell; position = dumbbell centre (bead midpoint, modulo the box); heading =
dumbbell axis angle atan2 of (bead 1 - bead 2), modulo 2 pi. SPEED = 0 (no propulsion). 400 dumbbells (800 beads) at
phi = 0.65: box side L = sqrt(800 pi / (4 * 0.65)) = 31.09; area4: 1600 dumbbells, L = 62.18.
Rules in plan numbering: 1 spring inside a dumbbell (header, always on: FENE plus the WCA term between its two beads);
2 constant torque (one-way); 3 inertia, friction, thermal noise (one-way); 4 bead contact repulsion between dumbbells
(interaction). Knock-out (keep=()) removes rule 4: dumbbells pass through each other.

Census time: the paper's dimensionless time unit is 1/Dr (Appendix A: "(m, sigma, Dr^-1) units of mass, lengths and
time", Eq. A2; "Using Dr^-1 as timescale"). With Dr = 2 kB T/(gamma sigma^2) = 2, 1/Dr = 0.5 tau_LJ, where
tau_LJ = sigma sqrt(m/eps) is the unit in which the paper sets T, gamma, K, R0 and runs LAMMPS. One census step = 1/Dr =
0.5 tau_LJ = 250 sub-steps of dt = 0.002 tau_LJ; 1000 census steps = 1000/Dr = 500 tau_LJ. Integrator: velocity Verlet
with the Langevin friction and noise evaluated at the half-step velocity (LAMMPS fix nve + fix langevin). All times
quoted below in "tau" (fidelity lengths, observations) are tau_LJ. Measured full census seed (1000 steps): real 18 s,
knock-out 9 s (3 seeds of both ~1.4 min).

Choices the paper does not state, and why
  * Time step: not stated. dt = 0.002 tau keeps the stiffest motion (FENE+WCA bond, ~0.16 tau period, and collisions of
    beads moving up to ~10 sigma/tau, the free spinning speed Fa/gamma) resolved; see the energy check in the module test.
  * FENE prefactor: the paper writes K R0^2 ln[1-(r/R0)^2]; we use the standard (LAMMPS, Kremer-Grest) form
    -0.5 K R0^2 ln[1-(r/R0)^2] with the WCA term also acting inside the bond ("between every pair of beads").
    Equilibrium bond length 0.961 sigma.
  * Start: "Start homogeneous" / "quench from the disordered state". Random sequential addition of dumbbells (random
    centre and axis, no beads closer than sigma) at area fraction 0.4, compression to 0.65 in 30 affine steps of the
    centres with soft overdamped relaxation after each, bonds reset to 0.961; then 2 tau of Langevin dynamics at T with
    gamma = 2 and no torque (removes residual overlaps without heating); then Maxwell-Boltzmann velocities at T, the
    paper's gamma and the torque switched on at census time 0 (a quench into the coexistence region). Every census seed
    and every knock-out run uses this start (the knock-out keeps the inter-dumbbell repulsion off throughout).
  * Clockwise torque for every dumbbell: F on bead 1 = Fa (n_y, -n_x), n = unit vector from bead 2 to bead 1.
  * Neighbour search: Verlet pair list of radius 2^(1/6) + 0.4, rebuilt when any bead has moved more than 0.2.
  * A bond stretched to R0 (FENE divergence) marks the run unstable (never seen in tests).

Fidelity check (plan: "two coexisting densities at Pe = 5 and none at Pe = 0", the paper's local area-fraction pdf),
fixed before any large run: phi = 0.65, T = 0.01, gamma = 0.01, Pe = 5 and Pe = 0. System size: the largest of 32768,
16384, 8192, 4096, 2048 dumbbells (32768 dumbbells = the 256^2 beads of the paper's Fig. 6 and 7) for which a measured
cost (x2 margin) allows 3000 tau at Pe = 5 plus 300 tau at Pe = 0 within 90% of the budget; Pe = 0 runs 300 tau (it is
homogeneous from the start), Pe = 5 runs until its share of the budget is used, at most 12000 tau (domains coarsen as
t^(1/3), Fig. 4; in tests here separation develops over thousands of tau). Small-scale tests here: Tkin/T falls with density (34 at phi = 0.2, 11 at 0.45, 2 at 0.65 at Pe = 2.5;
the paper's Fig. 3(c) trend); at 12,800 beads (L = 124) the s.d. of phi_loc on 20 sigma cells grows from 0.008 to
0.04 by 1200 tau and then stays ~0.04 to 2200 tau (range 0.59-0.72; on 8 sigma cells two populations near 0.62 and
0.70), i.e. separation develops slowly; at census size (L = 31) the state stays homogeneous for 1000 tau_LJ, as the paper
reports for its 1024-bead system ("segregation prevented for this small system size").
phi_loc from bead counts on a grid of g x g cells, g = round(L/20) (cells of about 20 sigma, the paper's grid), pooled
over frames of the last 40% of each run; pdf smoothed with a Gaussian of s.d. max(0.01, Silverman bandwidth); peaks =
local maxima above 15% of the highest, distinct only if the pdf dips below 70% of the lower peak between them.
Tolerance: at Pe = 5 two distinct peaks separated by >= 0.08, low peak in [0.45, 0.62] and high peak in [0.68, 0.80]
(paper Fig. 2 binodal ~0.55 and ~0.73 at Pe = 5); at Pe = 0 a single peak within 0.02 of 0.65. Pass = both.
"""
import time
import numpy as np
from scipy.spatial import cKDTree
from census.planted_models import base

CID = "planted:N08"
LAYERS = "A"
K = {"A": 1}
SPEED = 0.0
RULES = [("dumbbell_spring", False), ("constant_torque", False), ("inertia_friction_noise", False),
         ("bead_contact_repulsion", True)]
PHI = 0.65
HEADLINE = {
    "phi": PHI,                      # Fig. 1(a), card (global surface fraction)
    "Pe": 5.0, "Fa": 0.1,            # Pe = sigma Fa / (2 kB T) (main text); Fig. 1(a) Pe = 5
    "T": 0.01, "gamma": 0.01, "m": 1.0, "sigma": 1.0, "eps": 1.0,   # main text ("T = 1e-2", "fixed gamma = 1e-2")
    "Pr": 0.005, "Re": 1000.0,       # Pr = gamma/(m Dr) = Pe/Re (main text), Dr = 2 kB T/(gamma sigma^2) = 2
    "K_fene": 30.0, "R0_fene": 1.5,  # main text
    "N_dumbbells": 400, "box": float(np.sqrt(800 * np.pi / (4 * PHI))),
    "dt": 0.002, "steps_per_census_step": 250,   # one census step = 1/Dr = 0.5 tau_LJ (Appendix A time unit)
    "steps_per_tau_LJ": 500,                       # fidelity runs are counted in tau_LJ
}
STEP = "one census step = 1/Dr = 0.5 LJ time units (the paper's Appendix A time unit); velocity Verlet + Langevin, dt = 0.002 tau_LJ, 250 sub-steps"
RULE_TEXT = ("World: 400 dumbbells (two identical beads of diameter 1 joined by a stiff spring) in a square periodic 2D box "
             "of side 31.1 at area fraction 0.65. Rules, Langevin dynamics with mass: (1) every pair of beads repels on "
             "contact (purely repulsive, no attraction); (2) each dumbbell feels a constant torque, equal and opposite "
             "forces perpendicular to its axis on its two beads, all spinning the same way, strength set by Peclet number "
             "5; (3) beads have mass and weak friction (momentum persists; Prandtl number 0.005) plus thermal noise at "
             "temperature 0.01. Start homogeneous and watch the local density field. Run 1,000 time units of 1/D_r.")

RC = 2 ** (1 / 6); RC2 = RC * RC; SKIN = 0.4; B0 = 0.961


def box_side(Nd, phi=PHI):
    return float(np.sqrt(2 * Nd * np.pi / (4 * phi)))


def _rsa(rng, Nd, L, dmin=1.0):
    """random sequential addition of dumbbells (random centre and axis) with no two beads closer than dmin."""
    h = 0.5 * B0; g = max(1, int(L // dmin)); cs = L / g; cells = {}; out = np.empty((2 * Nd, 2)); k = 0
    while k < Nd:
        c = rng.uniform(0, L, 2); a = rng.uniform(0, 2 * np.pi); u = np.array([np.cos(a), np.sin(a)])
        bb = np.mod(np.array([c + h * u, c - h * u]), L); ok = True
        for p in bb:
            ci, cj = int(p[0] // cs) % g, int(p[1] // cs) % g
            for q in (q for di in (-1, 0, 1) for dj in (-1, 0, 1) for q in cells.get(((ci + di) % g, (cj + dj) % g), ())):
                d = np.mod(out[q], L) - p; d -= L * np.round(d / L)
                if d[0] * d[0] + d[1] * d[1] < dmin * dmin: ok = False; break
            if not ok: break
        if not ok: continue
        for m, p in enumerate(bb):
            out[2 * k + m] = c + (h if m == 0 else -h) * u
            cells.setdefault((int(p[0] // cs) % g, int(p[1] // cs) % g), []).append(2 * k + m)
        k += 1
    return out


def _relax(b, Nd, L, iters, target):
    """overdamped relaxation of flexible dumbbells: soft repulsion (1 - r) between beads of different dumbbells (r < 1)
    plus a bond spring to length B0; per-bead displacement capped at 0.02."""
    rmin = 1.0
    for it in range(iters):
        bw = np.mod(b, L); bw[bw >= L] = 0.0
        pr = cKDTree(bw, boxsize=L).query_pairs(1.0, output_type="ndarray")
        if len(pr): pr = pr[(pr[:, 0] // 2) != (pr[:, 1] // 2)]
        F = np.zeros((2 * Nd, 2)); rmin = 1.0
        if len(pr):
            i, j = pr[:, 0], pr[:, 1]; d = bw[j] - bw[i]; d -= L * np.round(d / L); r = np.sqrt((d * d).sum(1)); rmin = r.min()
            f = (1.0 - r)[:, None] * d / np.maximum(r, 1e-9)[:, None]
            for ax in (0, 1):
                F[:, ax] -= np.bincount(i, f[:, ax], 2 * Nd); F[:, ax] += np.bincount(j, f[:, ax], 2 * Nd)
        if rmin > target: break
        db = b[0::2] - b[1::2]; rb = np.sqrt((db * db).sum(1)); fb = (-20.0 * (rb - B0) / rb)[:, None] * db
        F[0::2] += fb; F[1::2] -= fb
        mv = 0.2 * F; nm = np.sqrt((mv * mv).sum(1)); mv *= np.minimum(1.0, 0.02 / np.maximum(nm, 1e-12))[:, None]
        b = b + mv
    return b, rmin


def _init(rng, Nd, L, phi0=0.4, stages=30):
    """homogeneous disordered start: random sequential addition at area fraction phi0, then compression to phi in
    `stages` affine steps of the dumbbell centres, each followed by soft relaxation; bonds reset to length B0."""
    L0 = L * np.sqrt(PHI / phi0) if L > 0 else L; b = _rsa(rng, Nd, L0); Lp = L0
    for Ls in np.linspace(L0, L, stages + 1)[1:]:
        c = 0.5 * (b[0::2] + b[1::2]) * (Ls / Lp); db = b[0::2] - b[1::2]; Lp = Ls
        b = np.empty_like(b); b[0::2] = c + 0.5 * db; b[1::2] = c - 0.5 * db
        b, _ = _relax(b, Nd, Ls, 200, 0.9)
    db = b[0::2] - b[1::2]; rb = np.sqrt((db * db).sum(1)); c = 0.5 * (b[0::2] + b[1::2]); u = db / rb[:, None]
    b[0::2] = c + 0.5 * B0 * u; b[1::2] = c - 0.5 * B0 * u
    return b


class _Sys:
    def __init__(self, rng, Nd, L, Pe=5.0, T=0.01, gamma=0.01, m=1.0, K=30.0, R0=1.5, dt=0.002, inter=True):
        self.r, self.Nd, self.N, self.L, self.T, self.g, self.m, self.K, self.R02, self.dt = \
            rng, Nd, 2 * Nd, L, T, gamma, m, K, R0 * R0, dt
        self.Fa = 2 * T * Pe; self.inter = inter
        X = _init(rng, Nd, L)
        self.X, self.Y = X[:, 0].copy(), X[:, 1].copy()
        self.VX = rng.standard_normal(self.N) * np.sqrt(T / m); self.VY = rng.standard_normal(self.N) * np.sqrt(T / m)
        self.x0 = None; self.unstable = False
        # pre-equilibration (start): 2 tau at Pe = 0 with strong friction (gamma = 2) removes the residual overlaps left
        # by the soft relaxation without heating the system; then Maxwell-Boltzmann velocities at T, the paper's gamma,
        # and the torque switched on (a quench into the coexistence region)
        self.g, self.Fa, self.noise = 2.0, 0.0, np.sqrt(2 * 2.0 * T / dt)
        self.FX, self.FY = self.forces()
        for _ in range(int(round(2.0 / dt))): self.step()
        self.g, self.Fa, self.noise = gamma, 2 * T * Pe, np.sqrt(2 * gamma * T / dt)
        self.VX = rng.standard_normal(self.N) * np.sqrt(T / m); self.VY = rng.standard_normal(self.N) * np.sqrt(T / m)
        self.FX, self.FY = self.forces()

    def build(self):
        L = self.L; xw = np.mod(np.stack([self.X, self.Y], 1), L); xw[xw >= L] = 0.0
        pr = cKDTree(xw, boxsize=L).query_pairs(RC + SKIN, output_type="ndarray")
        if len(pr): pr = pr[(pr[:, 0] // 2) != (pr[:, 1] // 2)]
        else: pr = np.zeros((0, 2), np.intp)
        self.i, self.j = pr[:, 0].copy(), pr[:, 1].copy()
        self.sx = L * np.round((self.X[self.j] - self.X[self.i]) / L); self.sy = L * np.round((self.Y[self.j] - self.Y[self.i]) / L)
        self.x0 = (self.X.copy(), self.Y.copy())

    def forces(self):
        N = self.N; X, Y = self.X, self.Y; FX = np.zeros(N); FY = np.zeros(N)
        # bond (rule 1): FENE + WCA between the two beads of each dumbbell
        dx = X[0::2] - X[1::2]; dy = Y[0::2] - Y[1::2]; r2 = dx * dx + dy * dy
        q = 1.0 - r2 / self.R02
        if (q <= 0).any(): self.unstable = True; q = np.maximum(q, 1e-6)
        f = -self.K / q                                                     # FENE: F_1 = -K d / (1 - r^2/R0^2)
        sr6 = np.where(r2 < RC2, 1.0 / r2 ** 3, 0.0); f = f + 24.0 * (2 * sr6 * sr6 - sr6) / r2
        FX[0::2] += f * dx; FY[0::2] += f * dy; FX[1::2] -= f * dx; FY[1::2] -= f * dy
        # constant torque (rule 2): clockwise, F_1 = Fa (n_y, -n_x)
        r = np.sqrt(r2); nx, ny = dx / r, dy / r
        FX[0::2] += self.Fa * ny; FY[0::2] -= self.Fa * nx; FX[1::2] -= self.Fa * ny; FY[1::2] += self.Fa * nx
        # bead contact repulsion between dumbbells (rule 4)
        if self.inter:
            if self.x0 is None or max(np.abs(X - self.x0[0]).max(), np.abs(Y - self.x0[1]).max()) > SKIN / 2:
                self.build()
            i, j = self.i, self.j
            if len(i):
                ex = X[j] - X[i]; ex -= self.sx; ey = Y[j] - Y[i]; ey -= self.sy
                s2 = ex * ex + ey * ey; w = np.flatnonzero(s2 < RC2)
                if len(w):
                    s2w = s2[w]; s6 = 1.0 / s2w ** 3; fw = 24.0 * (2 * s6 * s6 - s6) / s2w   # repulsive magnitude / r
                    fx = fw * ex[w]; fy = fw * ey[w]; iw, jw = i[w], j[w]
                    FX -= np.bincount(iw, fx, N); FY -= np.bincount(iw, fy, N)
                    FX += np.bincount(jw, fx, N); FY += np.bincount(jw, fy, N)
        # Langevin friction + noise (rule 3), at the current (half-step) velocity
        FX += -self.g * self.VX + self.noise * self.r.standard_normal(N)
        FY += -self.g * self.VY + self.noise * self.r.standard_normal(N)
        return FX, FY

    def step(self):
        h = 0.5 * self.dt / self.m
        self.VX += h * self.FX; self.VY += h * self.FY
        self.X += self.dt * self.VX; self.Y += self.dt * self.VY
        self.FX, self.FY = self.forces()
        self.VX += h * self.FX; self.VY += h * self.FY

    def centres_heads(self):
        cx = 0.5 * (self.X[0::2] + self.X[1::2]); cy = 0.5 * (self.Y[0::2] + self.Y[1::2])
        P = np.mod(np.stack([cx, cy], 1), self.L).astype(np.float32); P[P >= self.L] = 0.0
        H = np.mod(np.arctan2(self.Y[0::2] - self.Y[1::2], self.X[0::2] - self.X[1::2]), 2 * np.pi).astype(np.float32)
        H[H >= 2 * np.pi] = 0.0
        return P, H

    def beads(self):
        b = np.mod(np.stack([self.X, self.Y], 1), self.L); b[b >= self.L] = 0.0
        return b


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    on = base.kept(keep, RULES); Nd = base.SIZES[size]["NA"]; L = box_side(Nd)
    nsub = HEADLINE["steps_per_census_step"]; nf = base.n_frames(steps, "pos"); rec = base.REC["pos"]
    P = np.zeros((S, nf, Nd, 2), np.float32); H = np.zeros((S, nf, Nd), np.float32); bad = False
    for s, r in enumerate(base.rngs(seed0, S)):
        sy = _Sys(r, Nd, L, Pe=HEADLINE["Pe"], T=HEADLINE["T"], gamma=HEADLINE["gamma"], dt=HEADLINE["dt"], inter=4 in on)
        P[s, 0], H[s, 0] = sy.centres_heads(); fi = 1
        for t in range(1, steps + 1):
            for _ in range(nsub): sy.step()
            if t % rec == 0: P[s, fi], H[s, fi] = sy.centres_heads(); fi += 1
        bad = bad or sy.unstable or not np.isfinite(P[s]).all()
    return {"pos": P, "head": H, "at": np.zeros(H.shape, np.int8), "meta": base.meta(steps, base.SIZES[size]["W"], L, bad)}


# ------------------------------------------------------------------ measure (pure numpy on positions)
def phi_loc(points, L, g, weight):
    """local surface fraction on a g x g grid: weight * count * (pi/4) / cell area (weight 1 for beads, 2 for centres)."""
    c = np.histogram2d(np.mod(points[:, 0], L), np.mod(points[:, 1], L), bins=g, range=[[0, L], [0, L]])[0]
    return (weight * c * np.pi / 4 / (L / g) ** 2).ravel()


def pdf_peaks(vals, lo=0.3, hi=0.95, bw=0.005, rel=0.15, valley=0.7):
    """smoothed pdf of phi_loc (Gaussian kernel, s.d. = max(0.01, Silverman bandwidth)) and its distinct peaks: local
    maxima above rel x the highest, where two neighbouring peaks count as distinct only if the pdf between them dips
    below valley x the lower peak (otherwise the lower one is dropped). Returns (x, pdf, [(position, height), ...])."""
    vals = np.asarray(vals, float); sd = max(0.01, 1.06 * vals.std() * len(vals) ** -0.2)
    edges = np.arange(lo, hi + bw, bw); h = np.histogram(vals, bins=edges, density=True)[0]
    k = np.arange(-4 * sd, 4 * sd + bw / 2, bw); ker = np.exp(-0.5 * (k / sd) ** 2); ker /= ker.sum()
    sm = np.convolve(h, ker, mode="same"); x = 0.5 * (edges[1:] + edges[:-1])
    pk = [i for i in range(1, len(sm) - 1) if sm[i] >= sm[i - 1] and sm[i] > sm[i + 1] and sm[i] >= rel * sm.max()]
    changed = True
    while changed and len(pk) > 1:
        changed = False
        for a in range(len(pk) - 1):
            i, j = pk[a], pk[a + 1]
            if sm[i:j + 1].min() > valley * min(sm[i], sm[j]):
                pk.pop(a if sm[i] < sm[j] else a + 1); changed = True; break
    return x, sm, [(round(float(x[i]), 3), round(float(sm[i]), 2)) for i in pk]


def _two_phase(peaks):
    if len(peaks) < 2: return False, None
    lo = min(peaks, key=lambda p: p[0])[0]; hi = max(peaks, key=lambda p: p[0])[0]
    return bool(hi - lo >= 0.08 and 0.45 <= lo <= 0.62 and 0.68 <= hi <= 0.80), (lo, hi)


def presence(out):
    """phi_loc pdf from dumbbell centres (each centre counts its two beads) on a 4 x 4 grid (cells of L/4 = 7.8 sigma in the
    census box; the paper's 20 sigma cells do not fit a box of side 31), pooled over the last 40% of frames per seed."""
    L = float(out["meta"]["WA"]); g = max(4, int(round(L / 20))); res = []
    for s in range(out["pos"].shape[0]):
        P = out["pos"][s]; nf = len(P)
        v = np.concatenate([phi_loc(p, L, g, 2.0) for p in P[int(0.6 * nf):]])
        _, _, pk = pdf_peaks(v); ok, lh = _two_phase(pk)
        res.append({"peaks": pk, "two_phase": ok, "phi_loc_sd": round(float(v.std()), 4)})
    return {"per_seed": res, "two_phase_seeds": int(sum(r["two_phase"] for r in res)), "grid_cells": g,
            "note": "census box side 31.1 sigma; paper grid 20 sigma, binodal ~0.55/0.73 at Pe = 5"}


def _run_measure(seed, Nd, Pe, T_units, deadline, rec_every=5):
    """run one fidelity setting; phi_loc (paper's 20 sigma grid, bead counts) stored every rec_every tau."""
    r = base.rngs(seed, 1)[0]; L = box_side(Nd); sy = _Sys(r, Nd, L, Pe=Pe, dt=HEADLINE["dt"])
    g = max(4, int(round(L / 20))); nsub = HEADLINE["steps_per_tau_LJ"]; vals = []; done = 0
    for t in range(1, T_units + 1):
        for _ in range(nsub): sy.step()
        done = t
        if t % rec_every == 0: vals.append(phi_loc(sy.beads(), L, g, 1.0))
        if time.perf_counter() > deadline: break
    nf = len(vals); late = vals[int(0.6 * nf):]
    v = np.concatenate(late) if late else np.array([])
    pk = pdf_peaks(v)[2] if len(v) else []
    return {"Pe": Pe, "N_beads": 2 * Nd, "box": round(L, 2), "grid": g, "cell_sigma": round(L / g, 2),
            "time_units": done, "frames_used": len(late), "peaks": pk, "unstable": sy.unstable,
            "phi_loc_mean": round(float(v.mean()), 4) if len(v) else None,
            "phi_loc_sd": round(float(v.std()), 4) if len(v) else None}


def _cost_per_bead_tau(Nd=1600, Pe=5.0, tau=2):
    r = base.rngs(1, 1)[0]; sy = _Sys(r, Nd, box_side(Nd), Pe=Pe, dt=HEADLINE["dt"])
    t0 = time.perf_counter()
    for _ in range(tau * HEADLINE["steps_per_tau_LJ"]): sy.step()
    return (time.perf_counter() - t0) / (2 * Nd * tau)


def fidelity(budget_s=7200, seed=None):
    """the plan's check (see module docstring): Pe = 5 and Pe = 0 at phi = 0.65 at the largest size that fits."""
    t0 = time.perf_counter(); seed = base.seed_of(CID) if seed is None else seed
    c = _cost_per_bead_tau() * 2.0                       # s per bead per tau, x2 for hot dilute regions / list rebuilds
    T5min, T5max, T0 = 3000, 12000, 300
    end = t0 + 0.9 * budget_s; usable = end - time.perf_counter()
    Nd = next((n for n in (32768, 16384, 8192, 4096, 2048) if c * 2 * n * (T5min + T0) <= usable), 1024)
    t_pe0 = c * 2 * Nd * T0                               # reserve for the Pe = 0 run
    r5 = _run_measure(seed, Nd, 5.0, T5max, end - t_pe0)
    r0 = _run_measure(seed + 1, Nd, 0.0, T0, end)
    ok5, lh = _two_phase(r5["peaks"])
    ok0 = len(r0["peaks"]) == 1 and abs(r0["peaks"][0][0] - PHI) <= 0.02
    return {"setting": {"phi": PHI, "T": 0.01, "gamma": 0.01, "m": 1.0, "dt": HEADLINE["dt"], "Pe": [5.0, 0.0]},
            "estimated_s_per_bead_tau": float(c), "runs": [r5, r0],
            "measured": {"Pe5_coexisting": lh, "Pe0_peaks": r0["peaks"]},
            "paper": {"Pe5_binodal_approx": (0.55, 0.73), "Pe2.5_fig6_peaks": (0.57, 0.71), "Pe0": "single peak at phi"},
            "tolerance": "Pe=5: two peaks >= 0.08 apart, low in [0.45,0.62], high in [0.68,0.80]; Pe=0: one peak within 0.02 of 0.65",
            "pass": bool(ok5 and ok0), "elapsed_s": round(time.perf_counter() - t0, 1)}

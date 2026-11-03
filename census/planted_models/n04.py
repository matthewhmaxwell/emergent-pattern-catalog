"""planted:N04 - 2D solid of slightly elongated particles with an orientation-dependent Lennard-Jones (Gay-Berne) pair
potential plus a transverse (odd) non-conservative pair force, overdamped Brownian dynamics (arXiv 2609.09500,
End Matter I.3 "Simulation details" and Fig. 1).

Model as the paper states it (End Matter I.3, Eqs. 16-18; Fig. 1 caption)
-------------------------------------------------------------------------
* N particles in a periodic 2D box at number density rho ~ 0.87 (mean triangular lattice spacing
  a = sqrt(2 / (rho sqrt 3)) ~ 1.15). The paper used 10^4 particles; the box was "sized slightly smaller than the
  equilibrium extent of the crystal".
* Conservative anisotropic pair potential (HOOMD-blue Gay-Berne, Eqs. 17-18), strength e_p:
      U = 4 e_p [zeta^-12 - zeta^-6]  for zeta < zeta_cut, else 0;   zeta = (r - sigma + sigma_min) / sigma_min
      sigma^-2 = (1/2) r_hat . H^-1 . r_hat,   H = 2 l_perp^2 I + (l_par^2 - l_perp^2)(e_i e_i + e_j e_j)
  with l_perp = 0.55, l_par = 0.65 (aspect ratio ~1.2); e_i = (cos th_i, sin th_i) is particle i's long axis.
  It gives forces on the centres and torques on the orientations (both from -grad U).
* Non-conservative transverse ("odd") pair force of Lennard-Jones form, strength e_a (Eq. 16):
      F_perp = 4 e_a [12 sigma^12 / r^13 - 6 sigma^6 / r^7] theta_hat,
  acting perpendicular to the separation; spherically symmetric (no torque on the orientations). For the pair
  (i, j) with r_hat = (r_i - r_j)/r we use theta_hat_i = z x r_hat on i and z x (-r_hat) on j, so F_j = -F_i and
  every pair is pushed round in the same rotational sense.
* Overdamped Brownian integrator (HOOMD md.methods.Brownian), linear drag Gamma = 10, rotational drag
  Gamma_r = 0.1 (rotations relax faster):  dr = F/Gamma dt + sqrt(2 T dt / Gamma) xi,
  dth = tau/Gamma_r dt + sqrt(2 T dt / Gamma_r) xi.  dt = 1e-3 "for the majority of simulations"; "for large
  activity and high temperatures, we adapted the time step to ensure simulation stability".
* Paper runs: 10^5 time steps (100 time units at dt = 1e-3), snapshots every 10 steps.
* Activity = e_a / e_p, swept 0..33; temperature T (k_B T in units of e_p, e_p = 1) swept 0..6.

Headline (plan row N04): density 0.87, e_a/e_p = 33, T = 6 - the "restored order" state of Fig. 1(b) (diamond:
k_B T = 6.0, e_a/e_p = 33.0, ordered), where the passive system at the same T is a liquid (Fig. 1(b) square:
k_B T = 6.0, e_a/e_p = 0.0, liquid).

Choices the paper does not state (and why)
------------------------------------------
* e_p = 1 (energy unit), so e_a = 33 and T = 6 in units of e_p (the phase diagram axes are e_a/e_p and k_B T).
* sigma of the transverse force (Eq. 16 does not give it) = 1.16, chosen so that the transverse force reproduces the
  odd moduli of the paper's Table 1 (SM II.5: A = 400 e_a from isotropic compression -> torque, K^o = 140 e_a from
  shear -> orthogonal shear stress). Check (perfect triangular lattice at rho = 0.87, T = 0, affine strain, stress
  sigma_ij = (1/S) sum_pairs r_i F_j as in SM Eq. 45, decomposed with the paper's Eq. 3): sigma_a = 1.16 gives
  |A| = 339 e_a and |K^o| = 136 e_a; sigma_a = 1.0 (the HOOMD unit length, used in an earlier version of this module)
  gives only 24 e_a and 15 e_a, because its force changes sign at r = 2^(1/6) = 1.122, right at the neighbour
  spacing 1.15 (sigma_a = 1.12: 202 / 84; 1.20: 548 / 214).
* Cut-offs: Gay-Berne r_cut = 2.5 (HOOMD convention: zeta_cut = (r_cut - sigma_max + sigma_min)/sigma_min with
  sigma_min = 2 min(l_perp, l_par) = 1.1, sigma_max = 2 max(...) = 1.3); transverse force cut at r = 2.5 (the odd
  moduli above are the same with a cut at 2.9). Neither potential is shifted (HOOMD mode "none").
* Time step: dt = 1e-5 (a hundredth of the paper's default; the paper "adapted the time step" at large activity and
  high temperature without giving it). Measured at the headline (e_a/e_p = 33, T = 6, 400 particles): stable for 4
  time units at dt = 1e-5 and 1.5e-5, diverges within the first 0.2 time units at dt = 2e-5, 3e-5 and 5e-5. The
  passive system (rule 1 only) is stable at 2e-5 (4 time units) and at 2.5e-4; the knock-out is integrated exactly.
  meta["n_deep_overlap"] counts pair evaluations with zeta < 0.5 (0 in the stable tests). No force cap.
* Census box (plan: 400 particles at the paper's density, box side from the density): square periodic box of side
  sqrt(400/0.87) = 21.44 (area4: 1600 particles, side 42.88). A triangular lattice cannot fit a square periodic box
  exactly; the initial condition is the closest commensurate triangular lattice of 400 (1600) sites, mapped onto the
  square by a small affine strain (6.0% anisotropy for 400 sites, lattice vectors (9,-22),(19,-2) in units of the
  primitive vectors; 3.6% for 1600 sites, (5,35),(-42,26)). The paper does not state its initial condition; the
  phrase "box sized slightly smaller than the equilibrium extent of the crystal" suggests a crystal start. Initial
  orientations uniform at random (with Gamma_r = 0.1 they decorrelate within ~0.01 time units at T = 6).
* Random numbers: one generator per seed (base.rngs(seed0, S)) for that seed's initial orientations and all its
  noise, and per-seed Verlet lists (skin 0.6, rebuilt when that seed's particles moved more than skin/2), so a seed's
  trajectory is bit-for-bit independent of S (checked: seed 0 identical for S = 1 and S = 3).
* Heading recorded = long-axis angle mod 2 pi (the potential is invariant under th -> th + pi).
* Numerical blow-up (non-finite positions; seen only in the sub-program with the transverse force alone, where
  nothing keeps particles apart and the 1/r^13 force diverges): the run stops, the frames recorded so far are
  returned and meta["unstable"] = True, as the census interpreter does on blow-up.
* Knock-out (rules 1, 2 off): independent Brownian particles; each census step adds the exact Gaussian increments
  sqrt(2 T dt_step/Gamma) and sqrt(2 T dt_step/Gamma_r) (no time discretisation error, no sub-steps needed).

One census step = 0.03 time units (3000 integrator steps of dt = 1e-5): the largest that keeps the census run (real +
knock-out, 3 seeds, 1000 steps) within 60 CPU-minutes, by the plan's compute cap. Measured 0.92-1.08 ms per
integrator step for 3 seeds batched -> 3e6 steps = 46-54 CPU-minutes for the real program; the knock-out takes under
a second. 1000 census steps = 30 time units (the paper ran 100 time units). Frames every 2 census steps (0.06 time
units).

Knock-out (plan): rules 1 and 2 removed (both pair interactions), leaving independent Brownian particles (rule 3).

What this implementation does (small-scale runs on the development machine, 400 particles; psi_6 global / fraction
of non-6-fold particles): at the headline (33, T = 6) the crystal stays ordered for 4 time units (0.85 / 0); at
T = 0.1 the solid is disordered at e_a/e_p = 2 (0.13 / 0.24; the paper: disordered at T = 0 for 2) and ordered at 5
and 33 (0.99 / 0). The passive system (e_a = 0), which sigma_a does not affect, stays crystalline at T = 6 (0.82 / 0)
and up to T = 13 (0.68 / 0.07), melting only by T = 15 (0.26 / 0.26), from a crystal or a random start and at 400 or
1600 particles: its melting temperature is about 14 here against about 5 in the paper (partial melting at 4.25, liquid
at 6, Fig. 1(b)). An independent reviewer found no single change of a stated parameter that reconciles this:
rescaling e_p by 1/3 matches 8 of 9 checked points but not the headline, and has no basis in the paper; a likely
contributor is thermal randomisation of the long axes, which raises the effective contact distance in this
compressed crystal (a = 1.15 against contact distances 1.1-1.3). The paper's comparison "restored order at the
headline above the liquid passive system at the same T" therefore cannot be reproduced here, and the fidelity check
below is expected to fail (passive psi_6 ~ 0.8 at T = 6).

Fidelity (plan): bond-orientational order at the headline above the passive system (e_a = 0) at the same T.
Measure: psi_6(i) = (1/N_nn) sum_j exp(6 i theta_ij) over Delaunay (Voronoi) neighbours (paper Eq. 5); reported as
the global order |<psi_6>| and the local order <|psi_6|>, averaged over the second half of the run, plus g6(r)
at r = 5a and the fraction of non-6-fold particles. Tolerance, fixed before any large run: PASS if the global
|<psi_6>| at the headline exceeds the passive value by at least 0.2 AND the headline g6(5a) >= 0.3 (ordered or
hexatic, Fig. 1(b) diamond), with the passive g6(5a) < headline g6(5a) (liquid, Fig. 1(b) square). The paper reports
curves, not numbers; these thresholds encode "ordered vs liquid". Run length: 100 time units (the paper's) on the
largest commensurate lattice that fits 80% of the budget; if even 400 particles do not fit, the run is shortened to
fit (the result records t_total). Expected to FAIL for this implementation (see above).
"""
import time
import numpy as np
from scipy.spatial import cKDTree, Delaunay

from census.planted_models import base

CID = "planted:N04"
LAYERS = "A"
K = {"A": 1}
SPEED = 0.0
RULES = [("ORIENT_LJ_PAIR", True), ("ODD_TRANSVERSE_PAIR", True), ("THERMAL_NOISE", False)]
HEADLINE = dict(
    rho=0.87,            # End Matter I.3: "density rho ~ 0.87", a = 1.15
    eps_p=1.0,           # energy unit (not stated; see docstring)
    eps_a=33.0,          # Fig. 1(b) diamond: e_a/e_p = 33.0, ordered
    T=6.0,               # Fig. 1(b) diamond: k_B T = 6.0
    l_perp=0.55, l_par=0.65,   # End Matter I.3
    gamma=10.0, gamma_r=0.1,   # End Matter I.3 Brownian integrator damping
    sigma_a=1.16,        # not stated in Eq. 16; set so the transverse force gives the paper's Table 1 odd moduli
    rcut_gb=2.5, rcut_a=2.5,   # not stated (HOOMD convention)
    dt=1e-5,             # paper: 1e-3 default, "adapted" at large activity / high T (value not given); stability
    time_per_step=0.03,  # census step = 0.03 time units (plan compute cap: largest within 60 CPU-minutes)
)
STEP = "one census step = 0.03 time units = 3000 Brownian steps of dt = 1e-5 (plan compute cap: largest within 60 CPU-minutes; 1000 steps = 30 time units, paper ran 100)"
RULE_TEXT = ("World: 400 slightly elongated particles (ellipse half-axes 0.55 and 0.65) in a periodic 2D box at density "
             "0.87 (lattice spacing about 1.15). Rules, overdamped Brownian dynamics at "
             "temperature T: (1) each pair interacts through an orientation-dependent Lennard-Jones potential of strength "
             "e_p (repulsive core, weak attraction, depending on both particles' orientations, which can rotate); (2) each "
             "pair also feels a non-conservative force with Lennard-Jones magnitude and strength e_a, directed "
             "perpendicular to the line joining them, always in the same rotational sense; (3) thermal noise acts on "
             "positions and orientations, orientations relaxing faster. Set e_a/e_p = 33 and T = 6 and run for 30 time "
             "units; measure six-fold bond order and defect counts.")

_SKIN = 0.6
_LATTICE = {400: ((9, -22), (19, -2)), 1600: ((5, 35), (-42, 26))}   # best commensurate triangular cells (square box)


# ---------------------------------------------------------------------------------------------------- initial state
def _lattice_square(n):
    """n triangular-lattice sites in the unit square (fractional coordinates), from the closest commensurate cell."""
    (m1, n1), (m2, n2) = _LATTICE[n]
    a1 = np.array([1.0, 0.0]); a2 = np.array([0.5, np.sqrt(3) / 2])
    B = np.column_stack([m1 * a1 + n1 * a2, m2 * a1 + n2 * a2]); Binv = np.linalg.inv(B)
    R = int(np.abs([m1, n1, m2, n2]).sum()) + 2
    i, j = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")
    pts = np.stack([i.ravel(), j.ravel()], 1) @ np.stack([a1, a2])
    f = pts @ Binv.T
    f = f[(f[:, 0] >= -1e-9) & (f[:, 0] < 1 - 1e-9) & (f[:, 1] >= -1e-9) & (f[:, 1] < 1 - 1e-9)]
    assert len(f) == n, (len(f), n)
    f = np.mod(f, 1.0); f[f >= 1.0] = 0.0
    return f


def _lattice_rect(nx, ny):
    """perfect triangular lattice nx x ny (ny even) in a commensurate rectangle; returns positions, (Lx, Ly) at a=1."""
    i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    x = i + 0.5 * (j % 2) + 0.25; y = (j + 0.5) * np.sqrt(3) / 2
    return np.stack([x.ravel(), y.ravel()], 1), np.array([nx, ny * np.sqrt(3) / 2])


# ---------------------------------------------------------------------------------------------------- simulator core
class _Sim:
    """S independent systems of N particles batched together (flat arrays of length S*N)."""

    def __init__(self, pos, th, box, P, gb=True, act=True, rngs=None):
        self.S, self.N = pos.shape[:2]
        self.box = np.asarray(box, float)                      # (Lx, Ly)
        self.P = P; self.gb = gb; self.act = act; self.rngs = list(rngs)     # one generator per system (seed)
        assert len(self.rngs) == self.S
        self._xi = np.empty((self.S, self.N, 2)); self._eta = np.empty((self.S, self.N))
        self.x = pos.reshape(-1, 2).astype(float).copy()       # unwrapped positions
        self.th = th.reshape(-1).astype(float).copy()
        self.rlist = max(P["rcut_gb"], P["rcut_a"]) + _SKIN
        self.rc2 = max(P["rcut_gb"], P["rcut_a"]) ** 2
        self.I = self.J = None; self.n_deep = 0; self.n_rebuild = 0
        lp, ll = P["l_perp"], P["l_par"]
        self.A = 2 * lp ** 2; self.Bc = ll ** 2 - lp ** 2
        self.smin = 2 * min(lp, ll); smax = 2 * max(lp, ll)
        self.zcut = (P["rcut_gb"] - smax + self.smin) / self.smin
        if gb or act: self._build()

    def _build(self, seeds=None):
        """Verlet lists per system; a system's list is rebuilt only when its own particles moved > skin/2, so its
        trajectory (including floating-point summation order) does not depend on the other systems or on S."""
        if seeds is None:
            self.Is = [None] * self.S; self.Js = [None] * self.S; self.x0 = self.x.copy(); seeds = range(self.S)
        N = self.N
        for s in seeds:
            sl = slice(s * N, (s + 1) * N)
            p = np.mod(self.x[sl], self.box); p[p >= self.box] = 0.0
            pr = cKDTree(p, boxsize=self.box).query_pairs(self.rlist, output_type="ndarray")
            self.Is[s] = pr[:, 0] + s * N; self.Js[s] = pr[:, 1] + s * N
            self.x0[sl] = self.x[sl]; self.n_rebuild += 1
        self.I = np.concatenate(self.Is); self.J = np.concatenate(self.Js)

    def forces(self):
        """forces (S*N, 2) and torques (S*N,) from the switched-on pair rules."""
        n = self.x.shape[0]; F = np.zeros((n, 2)); tq = np.zeros(n)
        if not (self.gb or self.act): return F, tq
        moved = np.max(np.sum((self.x - self.x0) ** 2, 1).reshape(self.S, self.N), 1) > (0.5 * _SKIN) ** 2
        if moved.any(): self._build(np.nonzero(moved)[0])
        I, J = self.I, self.J
        d = self.x[I] - self.x[J]; d -= self.box * np.round(d / self.box)
        r2 = d[:, 0] * d[:, 0] + d[:, 1] * d[:, 1]
        m = r2 < self.rc2                                       # both cut-offs are <= rcut (2.5)
        I, J, d, r2 = I[m], J[m], d[m], r2[m]
        ir = 1.0 / np.sqrt(r2); r = r2 * ir; rx, ry = d[:, 0] * ir, d[:, 1] * ir
        P = self.P
        fx = fy = 0.0
        if self.act:
            s2 = (P["sigma_a"] ** 2) / r2; sr6 = s2 * s2 * s2
            fa = (4 * P["eps_a"]) * (12 * sr6 * sr6 - 6 * sr6) * ir
            if P["rcut_a"] < np.sqrt(self.rc2): fa = np.where(r < P["rcut_a"], fa, 0.0)
            fx = -fa * ry; fy = fa * rx                         # theta_hat_i = z x r_hat
        if self.gb:
            c = np.cos(self.th); s = np.sin(self.th)
            ci, si, cj, sj = c[I], s[I], c[J], s[J]
            A, B = self.A, self.Bc
            Hxx = A + B * (ci * ci + cj * cj); Hyy = A + B * (si * si + sj * sj); Hxy = B * (ci * si + cj * sj)
            idet = 1.0 / (Hxx * Hyy - Hxy * Hxy)
            vx = (Hyy * rx - Hxy * ry) * idet; vy = (Hxx * ry - Hxy * rx) * idet    # v = H^-1 r_hat
            q = 0.5 * (rx * vx + ry * vy)                                          # sigma^-2
            sig = 1.0 / np.sqrt(q)
            z = (r - sig + self.smin) * (1.0 / self.smin)
            self.n_deep += int(np.count_nonzero(z < 0.5))
            iz = 1.0 / np.maximum(z, 1e-3); iz2 = iz * iz; iz6 = iz2 * iz2 * iz2; iz7 = iz6 * iz
            dU = (4 * P["eps_p"]) * (6.0 - 12.0 * iz6) * iz7
            dU = np.where(z < self.zcut, dU, 0.0)
            sig3 = sig * sig * sig; s3r = sig3 * ir
            a = -dU * (1.0 / self.smin)
            fx = fx + a * (rx + s3r * (0.5 * vx - q * rx)); fy = fy + a * (ry + s3r * (0.5 * vy - q * ry))
            pre = dU * (0.5 * B / self.smin) * sig3
            tq_i = pre * (vx * ci + vy * si) * (vy * ci - vx * si)
            tq_j = pre * (vx * cj + vy * sj) * (vy * cj - vx * sj)
            tq = np.bincount(I, tq_i, n) + np.bincount(J, tq_j, n)
        F[:, 0] = np.bincount(I, fx, n) - np.bincount(J, fx, n)
        F[:, 1] = np.bincount(I, fy, n) - np.bincount(J, fy, n)
        return F, tq

    def step(self, nsteps):
        P = self.P; dt = P["dt"]; T = P["T"]
        st = np.sqrt(2 * T * dt / P["gamma"]); sr = np.sqrt(2 * T * dt / P["gamma_r"])
        n = self.x.shape[0]
        for _ in range(nsteps):
            F, tq = self.forces()
            xi, eta = self.noise()
            self.x += F * (dt / P["gamma"]) + st * xi
            self.th += tq * (dt / P["gamma_r"]) + sr * eta

    def noise(self):
        """unit Gaussian increments for positions and orientations, each system drawing from its own generator."""
        for s, r in enumerate(self.rngs):
            r.standard_normal(out=self._xi[s]); r.standard_normal(out=self._eta[s])
        return self._xi.reshape(-1, 2), self._eta.reshape(-1)

    def wrapped(self):
        p = np.mod(self.x, self.box).reshape(self.S, self.N, 2)
        return p, np.mod(self.th, 2 * np.pi).reshape(self.S, self.N)


def _params(**over):
    P = dict(HEADLINE); P.update(over); return P


# ---------------------------------------------------------------------------------------------------- census entry
def simulate(seed0, steps=1000, S=3, size="census", keep=None, **over):
    on = base.kept(keep, RULES)
    P = _params(**over)
    z = base.SIZES[size]; N = z["NA"]; L = np.sqrt(N / P["rho"])
    rs = base.rngs(seed0, S)
    f = _lattice_square(N)
    pos = np.stack([f * L for _ in range(S)])
    th = np.stack([r.uniform(0, 2 * np.pi, N) for r in rs])
    sim = _Sim(pos, th, (L, L), P, gb=1 in on, act=2 in on, rngs=rs)      # seed s uses rs[s] only
    nsub = int(round(P["time_per_step"] / P["dt"]))
    nf = base.n_frames(steps, "pos")
    out_pos = np.zeros((S, nf, N, 2), np.float32); out_head = np.zeros((S, nf, N), np.float32)

    def snap(k):
        p, h = sim.wrapped()
        p32 = p.astype(np.float32); p32[p32 >= L] = 0.0
        out_pos[:, k] = p32; out_head[:, k] = h.astype(np.float32)
    free = not (1 in on or 2 in on)                     # knock-out: exact Brownian increments per census step
    tau = nsub * P["dt"]
    st, sr = np.sqrt(2 * P["T"] * tau / P["gamma"]), np.sqrt(2 * P["T"] * tau / P["gamma_r"])
    snap(0); k = 1; unstable = False; t = 0
    for t in range(1, steps + 1):
        if free:
            xi, eta = sim.noise(); sim.x += st * xi; sim.th += sr * eta
        else:
            with np.errstate(all="ignore"): sim.step(nsub)
        if not (np.all(np.isfinite(sim.x)) and np.all(np.isfinite(sim.th))):
            unstable = True; break                          # numerical blow-up: stop, keep the frames so far
        if t % base.REC["pos"] == 0: snap(k); k += 1
    out = {"pos": out_pos[:, :k], "head": out_head[:, :k], "at": np.zeros((S, k, N), np.int8)}
    out["meta"] = base.meta(t, z["W"], WA=L, unstable=unstable)
    out["meta"]["n_deep_overlap"] = int(sim.n_deep)
    return out


# ---------------------------------------------------------------------------------------------------- measures
def _delaunay_neighbours(p, box):
    """unique periodic Delaunay edges (i, j) and per-particle degrees for positions p in a periodic box."""
    N = len(p); box = np.asarray(box, float)
    marg = 3.0 * np.sqrt(box.prod() / N)
    imgs, ids = [], []
    for ox in (-1, 0, 1):
        for oy in (-1, 0, 1):
            q = p + np.array([ox, oy]) * box
            keep = (q[:, 0] > -marg) & (q[:, 0] < box[0] + marg) & (q[:, 1] > -marg) & (q[:, 1] < box[1] + marg)
            imgs.append(q[keep]); ids.append(np.nonzero(keep)[0])
    Q = np.concatenate(imgs); ID = np.concatenate(ids)
    tri = Delaunay(Q).simplices
    e = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    inside = (Q[e[:, 0], 0] >= 0) & (Q[e[:, 0], 0] < box[0]) & (Q[e[:, 0], 1] >= 0) & (Q[e[:, 0], 1] < box[1])
    e = e[inside]
    a, b = ID[e[:, 0]], ID[e[:, 1]]
    d = Q[e[:, 1]] - Q[e[:, 0]]
    key = np.unique(np.stack([a, b, np.round(d[:, 0], 6), np.round(d[:, 1], 6)], 1), axis=0)
    return key[:, 0].astype(int), key[:, 1].astype(int), key[:, 2:]


def psi6_frame(p, box):
    """local psi_6 per particle (Delaunay neighbours, paper Eq. 5) and the coordination numbers."""
    p = np.mod(np.asarray(p, float), box); N = len(p)
    i, j, d = _delaunay_neighbours(p, box)
    ang = np.arctan2(d[:, 1], d[:, 0])
    z = np.exp(6j * ang)
    s = np.bincount(i, z.real, N) + 1j * np.bincount(i, z.imag, N)
    nn = np.bincount(i, minlength=N)
    return s / np.maximum(nn, 1), nn


def g6_at(p, psi, box, rq, dr=0.5):
    """bond-orientational correlation g6(r) = <psi_i psi_j*> over pairs at distance rq +- dr/2."""
    tr = cKDTree(np.mod(p, box), boxsize=box)
    pr = tr.query_pairs(rq + dr / 2, output_type="ndarray")
    if len(pr) == 0: return float("nan")
    d = p[pr[:, 0]] - p[pr[:, 1]]; d -= box * np.round(d / box); r = np.hypot(d[:, 0], d[:, 1])
    m = r > rq - dr / 2
    if not np.any(m): return float("nan")
    return float(np.mean((psi[pr[m, 0]] * np.conj(psi[pr[m, 1]])).real) / max(np.mean(np.abs(psi) ** 2), 1e-12))


def _order_summary(frames, box, a):
    G, Lc, G6, D = [], [], [], []
    for p in frames:
        psi, nn = psi6_frame(p, box)
        G.append(abs(psi.mean())); Lc.append(float(np.mean(np.abs(psi)))); D.append(float(np.mean(nn != 6)))
        G6.append(g6_at(np.mod(p, box), psi, np.asarray(box, float), 5 * a))
    return dict(psi6_global=float(np.mean(G)), psi6_local=float(np.mean(Lc)), g6_5a=float(np.nanmean(G6)),
                defect_fraction=float(np.mean(D)))


def presence(out):
    """psi_6 order on census-scale output of the real program: second half of the frames, every 10th frame."""
    L = out["meta"]["WA"]; box = np.array([L, L]); a = np.sqrt(2 / (HEADLINE["rho"] * np.sqrt(3)))
    pos = out["pos"]; S, nf = pos.shape[:2]
    res = []
    for s in range(S):
        fr = [pos[s, k].astype(float) for k in range(nf // 2, nf, 10)]
        res.append(_order_summary(fr, box, a))
    agg = {k: float(np.mean([r[k] for r in res])) for k in res[0]}
    agg["per_seed"] = res
    agg["note"] = "global |<psi6>| and local <|psi6|> over Delaunay neighbours; g6 at 5 lattice spacings"
    return agg


# ---------------------------------------------------------------------------------------------------- fidelity
def _run_rect(nx, ny, eps_a, T, t_total, seed, dt=None, nrec=20):
    """one system on a commensurate rectangle (perfect lattice start); returns recorded frames of the second half."""
    P = _params(eps_a=eps_a, T=T)
    if dt is not None: P["dt"] = dt
    a = np.sqrt(2 / (P["rho"] * np.sqrt(3)))
    pos, box = _lattice_rect(nx, ny); pos = pos * a; box = box * a
    rng = np.random.default_rng(seed)
    sim = _Sim(pos[None], rng.uniform(0, 2 * np.pi, (1, len(pos))), box, P, gb=True, act=eps_a != 0.0, rngs=[rng])
    nst = int(round(t_total / P["dt"])); half = nst // 2
    every = max(1, (nst - half) // nrec); frames = []
    done = 0
    while done < nst:
        k = min(every, nst - done); sim.step(k); done += k
        if done > half: frames.append(sim.wrapped()[0][0].copy())
    return frames, box, a, sim.n_deep


_DT_PASSIVE = 2.5e-4      # passive system at T = 6 is stable at this step (measured); the headline needs 1e-4


def fidelity(budget_s=7200, t_total=100.0):
    """headline (e_a/e_p = 33, T = 6) vs passive (e_a = 0, T = 6): psi_6 order, the paper's measure (Eq. 5)."""
    P = HEADLINE
    # cost estimate: time short runs per particle-step, then pick the largest commensurate size that fits
    n0 = 1000
    t0 = time.time(); _run_rect(20, 20, P["eps_a"], P["T"], n0 * P["dt"], 1, nrec=1); c_head = (time.time() - t0) / (400 * n0)
    t0 = time.time(); _run_rect(20, 20, 0.0, P["T"], n0 * _DT_PASSIVE, 1, dt=_DT_PASSIVE, nrec=1); c_pass = (time.time() - t0) / (400 * n0)
    cost_per_N = c_head * t_total / P["dt"] + c_pass * t_total / _DT_PASSIVE
    nmax = 0.8 * budget_s / cost_per_N
    sizes = [(100, 100), (80, 80), (64, 64), (50, 50), (40, 40), (30, 30), (24, 24), (20, 20)]
    nx, ny = next((s for s in sizes if s[0] * s[1] <= nmax), (20, 20))
    if nx * ny > nmax:                                  # even 400 particles do not fit: shorten the run instead
        t_total = t_total * nmax / (nx * ny); cost_per_N = cost_per_N * nmax / (nx * ny)
    est = cost_per_N * nx * ny
    res = {"setting": "rho 0.87, T 6, e_a/e_p 33 vs 0", "N": nx * ny, "box": "rectangle commensurate with a perfect "
           "triangular lattice", "t_total": t_total, "dt_headline": P["dt"], "dt_passive": _DT_PASSIVE,
           "estimated_s": round(est, 1)}
    fh, box, a, nd1 = _run_rect(nx, ny, P["eps_a"], P["T"], t_total, 11)
    fp, _, _, nd2 = _run_rect(nx, ny, 0.0, P["T"], t_total, 12, dt=_DT_PASSIVE)
    res["headline"] = _order_summary(fh, box, a); res["passive"] = _order_summary(fp, box, a)
    res["n_deep_overlap"] = nd1 + nd2
    res["paper"] = "Fig. 1(b): e_a/e_p=33, T=6 ordered (g6 flat); e_a/e_p=0, T=6 liquid (g6 decays)"
    res["tolerance"] = "PASS if psi6_global(headline) - psi6_global(passive) >= 0.2 and g6_5a(headline) >= 0.3 > g6_5a(passive)"
    h, p = res["headline"], res["passive"]
    res["pass"] = bool(h["psi6_global"] - p["psi6_global"] >= 0.2 and h["g6_5a"] >= 0.3 and p["g6_5a"] < h["g6_5a"])
    return res

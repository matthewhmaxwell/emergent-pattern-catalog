"""N03 (planted card): autophoretic disks hovering above a no-flux wall, with chemo-attractive forces and chemo-repulsive
torques from screened chemical fields plus their mirror images (arXiv 2609.03991).

Model, as the paper states it (Sec. II, Eqs. 1-5; Sec. III, Eqs. 6-7; Appendix A, Eqs. 13, 14, 20, 21; Appendix B,
Eq. 22 and Table 1)
  dr_i/dt = Mt F_i^p + V_i^self + V_i^int,      dtheta_i/dt = Omega_i^int + sqrt(2 Dr) xi_i            (Eq. 1a, 1b)
  V_i^self = -mu_t v0 p_i,   V_i^int = -mu_t 2 v0 J_i,   Omega_i^int = -(9/4) v0 mu_r (p_i x J_i)        (Eqs. 2, 3)
  J_i = -sum_j [ (kappa + 1/r_ij) e^(-kappa r_ij)/r_ij^2 r_ij + Lambda_c (kappa + 1/r*_ij) e^(-kappa r*_ij)/r*_ij^2 r*_ij ]
        (Eq. 4), r_ij = (r_i - r_j)/b; r*_ij = (r_i - r_j*)/b with r_i = (x_i, y_i, h), r_j* = (x_j, y_j, -h), so the image
        term has the in-plane vector r_ij and the distance r* = sqrt(r^2 + (2h)^2); its z part is normal to the plane of motion.
  F_i^p = -sum_j grad_i U^p(r_ij), U^p(r) = eps[(rm/r)^12 - 2 (rm/r)^6 + 1] for r <= rm, 0 otherwise (Eq. 22);
  Mt = 1/(6 pi eta b). Eqs. (1)-(4) follow from the reduced dynamics, Eqs. (13)-(14), with the image-source field Eq. (20)
  and the dimensionless current J_i = |<mu>|/(2 v0) grad c, Eq. (21).
  p_i = (cos theta_i, sin theta_i) points to the catalytic cap. Pe = v0 tau/b, tau = 1/Dr (Eq. 6).
  Parameters (Table 1 and Sec. III): b = 1, h = 1.1, v0 = 2, area fraction Phi = 0.05, Pe = 100 (Dr = 0.02), eps = 30,
  mu_t = -1 (attraction), mu_r > 0, no-flux wall Lambda_c = 1, N = 300 / 1000 / 10000 (quantitative analyses N = 1000).
  Forward Euler-Maruyama with "a strict upper limit of displacement ds = 0.01"; periodic box, minimum image (Appendix B).
  Measures (Sec. III A, Eq. 7): S = fraction of particles in clusters (bound if closer than 2.1 b, a cluster is a bound
  set of more than 10 particles), its temporal s.d. sigma_S, and Q = mean number of neighbours within 2.1 b.
  Boiling bubbles (the new state): <S> ~ 1, sigma_S ~ 0, <Q> < 5.0 (full collapse has <Q> > 5.6), "continuous nucleation,
  drift, and annihilation of localized fluid cavities (bubbles) maintained inside the dense active phase" (Sec. IV).
  Headline (plan row N03): screening kappa = 0.20 and turning strength mu_r = 0.30. Fig. 2 (Lambda_c = 1) marks boiling
  bubbles at kappa = 0.20 for every mu_r from 0.15 to 0.40, and at kappa = 0.18 and 0.22 for mu_r >= 0.20; mu_r = 0.30 is
  the middle of that window, so the whole card range kappa 0.18-0.22 is boiling at this mu_r.

Census size: 400 disks at Phi = 0.05, box side L = sqrt(400 pi / 0.05) = 158.5 b; area4: 1600 disks, L = 317.1 b.
Census time: the paper reports every time as t/tau, tau = 1/Dr = 50 simulation units (b = 1, v0 = 2), so tau is its
dimensionless time unit. One census step of tau (50 units) does not fit the plan's compute cap (census real + knock-out,
3 seeds, 1000 steps, <= 60 CPU-minutes): measured here, one simulation unit of the real program costs ~0.04 s in the gas,
~0.07-0.09 s with several clusters and 0.187 s with all 400 disks in one packed cluster (66,700 listed pairs), and
0.0095 s for the knock-out. Reduced by the compute cap to 5 simulation units = 0.1 tau per census step (worst case
3 x 5000 x (0.187 + 0.0095) s = 49 min, typical ~25 min). 1000 census steps = 5000 units = 100 tau (Fig. 3 of the paper
runs to ~65 tau). Sub-steps: 5000 Euler-Maruyama steps of 0.001 per census step. Measured full census seed (1000 steps):
real 478 s, knock-out 51 s, so 3 seeds of both take ~26 min.
Rules in plan numbering: 1 self-propulsion (one-way), 2 chemical attraction force (interaction), 3 chemo-repulsive torque
(interaction), 4 rotational noise (one-way), 5 contact repulsion (interaction). Knock-out (keep=()) removes 2, 3, 5:
independent active Brownian disks.

Choices the paper does not state, and why
  * eta (hence Mt) is not given: eta = 1 in simulation units, Mt = 1/(6 pi) = 0.0531. With eps = 30 this gives a WCA
    velocity scale Mt eps / b = 1.6, the same order as v0 and the phoretic attraction; in a packed cluster of 1000 disks
    the nearest-neighbour distance stays >= 1.80 b (overlap <= 10% of a diameter; median 1.83 b at mu_r = 0).
  * rm (WCA range) is not given: rm = 2b, contact of two disks of radius b (the B1 form has its minimum at rm).
  * Time stepping: Euler-Maruyama sub-steps dt = 0.001 for propulsion, contact repulsion, torque direction and noise, each
    disk's displacement per sub-step capped at the paper's ds = 0.01; the phoretic current J_i (Eq. 4) is a smooth sum over
    hundreds of neighbours and is recomputed every 10 sub-steps (every 0.01 time unit; a disk moves at most 0.1 b in that
    time, so J changes by a few per cent). The torque Omega_i uses the current p_i with that J_i at every sub-step.
    Check: 300 disks, 2 seeds, 300 time units, S(t) and Q(t) with J every sub-step vs every 10 sub-steps differ by less
    than the seed-to-seed spread (S = 0.98/1.00 vs 0.99/0.96, Q = 4.77/4.73 vs 4.63/4.70 at t = 300).
  * Interaction cutoff for Eq. (4): r < 25 b = 5/kappa (the paper sums all pairs by minimum image). At 25 b a single pair's
    current is 3e-4 of its value at contact, and the one-sided far tail of a dense half-plane is below 1% of the pull felt
    by a disk at a cluster edge (estimate from integrating Eq. 4 over a packed half-plane); 25 b < L/2 in every box used.
  * Neighbour lists: Verlet list of radius 27 b rebuilt when any disk has moved more than 1 b; contact pairs (r < 2.5 b)
    taken from it at each phoretic update (relative motion within 10 capped sub-steps is at most 0.2 b).
  * Self-image and j = i terms are neglected as in the paper (Appendix A; the self-image term has no in-plane part).
    Without contact repulsion (sub-programs) distances are floored at 0.1 b in Eq. (4) to keep the current finite.
  * Start: random non-overlapping positions (centres at least 2b apart) and random orientations (Fig. 4 caption).

Small-scale observations while building (not a fidelity result): at census size (400 disks, 1000 simulation units) the disks
gather within ~5 tau into 6-7 compact clusters of ~60 that stay apart (S = 1, Q ~ 4.9-5.0, no enclosed cavity); with
1000 disks the clusters merge slowly (one of ~250 disks by 3000 time units). Over a full census run (100 tau) the 400
disks coarsen further: S = 1, Q rising from 4.84 to 5.32, still no enclosed cavity. A packed cluster of 1000 disks at the
headline stays fluid and disordered with transient small holes (Q ~ 5.4), while at mu_r = 0 it stays crystalline
(Q = 5.76, the paper's full collapse, Q > 5.6). Mismatch with the paper's Fig. 2: at kappa = 0.28, mu_r = 0.30 (gas in
the paper) these equations, as printed, form clusters (S ~ 0.9 by 600 time units with 300 disks); the edge pull
2 v0 |J| of a packed cluster exceeds v0 for every kappa up to 0.5 (estimate from Eq. 4), so the printed prefactors give
stronger cohesion than the paper's phase diagram implies. Not tuned away: no stated parameter would do it.

Fidelity check (plan: "cluster fraction and cavities inside clusters at the headline", at the paper's 1000 disks),
fixed before any large run: N = 1000, Phi = 0.05 (L = 250.7 b), kappa = 0.20, mu_r = 0.30, random start, run until 90% of
the budget is used (at most 20,000 time units = 400 tau; Fig. 3 of the paper reaches ~65 tau, but in tests here the
clusters coarsen more slowly than in the paper, so the budget is spent on length); late window = last third of the run.
Measures: <S>, sigma_S, <Q> (the paper's metric triad) and cavities: free regions inside the dense phase, i.e. connected regions of grid points
(spacing 0.5 b, periodic) farther than 2b from every disk centre, other than the largest such region (the gas outside),
with area >= pi b^2 (a hole larger than a single vacancy). Tolerance: <S> >= 0.9, sigma_S <= 0.05, <Q> < 5.0, and at
least one cavity in >= 50% of late frames. Pass = all four.
"""
import time
import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from census.planted_models import base

CID = "planted:N03"
LAYERS = "A"
K = {"A": 1}
SPEED = 2.0
RULES = [("self_propulsion", False), ("chemical_attraction_force", True), ("chemorepulsive_torque", True),
         ("rotational_noise", False), ("contact_repulsion", True)]
PHI = 0.05
HEADLINE = {
    "N": 400, "phi": PHI, "box": float(np.sqrt(400 * np.pi / PHI)),   # plan row N03 (paper Phi = 0.05, Table 1)
    "b": 1.0, "h": 1.1, "v0": 2.0, "Pe": 100.0, "Dr": 0.02,          # Table 1; Dr = v0/(Pe b) (Eq. 6)
    "eps_wca": 30.0, "rm": 2.0, "eta": 1.0, "Mt": float(1 / (6 * np.pi)),   # Table 1 eps; rm, eta: our choices (docstring)
    "mu_t": -1.0,                                                    # Sec. III ("we choose mu_t = -1")
    "mu_r": 0.30,                                                    # Fig. 2: boiling bubbles at kappa 0.18-0.22, mu_r 0.20-0.40
    "kappa": 0.20, "Lambda_c": 1.0,                                  # plan headline; no-flux wall (Eq. 5)
    "cutoff": 25.0, "dt": 0.001, "ds_cap": 0.01, "chem_every": 10,   # our numerical choices (docstring)
    "tau": 50.0,                                                     # 1/Dr: the paper's time unit (t/tau, Fig. 3, 4)
    "units_per_census_step": 5, "substeps_per_unit": 1000,          # 0.1 tau per census step (compute cap)
    "steps_per_census_step": 5000,
}
STEP = ("one census step = 0.1 tau = 5 simulation units (tau = 1/Dr = 50 units, the paper's time unit; reduced from one tau "
        "by the 60 CPU-min compute cap); 5000 Euler-Maruyama sub-steps dt = 0.001 (displacement capped at 0.01), phoretic "
        "current recomputed every 10 sub-steps")
RULE_TEXT = ("World: 400 self-propelled disks (radius 1) in a periodic 2D box (square of side 158.5) at area fraction 0.05, "
             "hovering at height 1.1 above a flat wall that chemicals cannot cross. Each disk emits a chemical. Rules each "
             "step: (1) each disk moves along its heading at speed 2; (2) it is pulled toward the others along the summed "
             "gradient of their screened chemical fields (each decays as exp(-k*r)/r, plus a mirror-image source below the "
             "wall), screening k = 0.20; (3) its heading is turned away from that same gradient with turning strength 0.30; "
             "(4) the heading gets rotational noise (Peclet number 100); (5) disks repel on contact. A cluster is more than "
             "10 disks linked within 2.1 radii. Run 100 tau (tau = 1/rotational diffusion rate, the heading's persistence time).")

RC = 25.0; SKIN = 2.0; RW = 2.5


def box_side(N, phi=PHI):
    return float(np.sqrt(N * np.pi / phi))


def _rsa(rng, N, L, dmin=2.0):
    g = max(1, int(L // dmin)); cs = L / g; cells = {}; out = np.empty((N, 2)); k = 0
    while k < N:
        p = rng.uniform(0, L, 2); ci, cj = int(p[0] // cs) % g, int(p[1] // cs) % g; ok = True
        for q in (q for di in (-1, 0, 1) for dj in (-1, 0, 1) for q in cells.get(((ci + di) % g, (cj + dj) % g), ())):
            d = out[q] - p; d -= L * np.round(d / L)
            if d[0] * d[0] + d[1] * d[1] < dmin * dmin: ok = False; break
        if ok:
            out[k] = p; cells.setdefault((ci, cj), []).append(k); k += 1
    return out


SUB_PER_UNIT = HEADLINE["substeps_per_unit"]          # Euler-Maruyama sub-steps per simulation unit (dt = 0.001)


def _params(**over):
    """every model parameter, taken from HEADLINE (keyword overrides for tests)."""
    H = HEADLINE
    p = dict(kappa=H["kappa"], mu_r=H["mu_r"], Lc=H["Lambda_c"], h=H["h"], v0=H["v0"], Dr=H["Dr"], eps=H["eps_wca"],
             rm=H["rm"], Mt=H["Mt"], dt=H["dt"], ds=H["ds_cap"], chem_every=H["chem_every"], rc=H["cutoff"])
    bad = set(over) - set(p)
    if bad: raise TypeError(f"unknown parameters {sorted(bad)}")
    p.update(over); return p


class _Sys:
    def __init__(self, rng, N, L, on, **over):
        P = _params(**over)
        kappa, mu_r, Lc, h, v0, Dr, eps, rm, Mt, dt, ds, chem_every, rc = (P[k] for k in (
            "kappa", "mu_r", "Lc", "h", "v0", "Dr", "eps", "rm", "Mt", "dt", "ds", "chem_every", "rc"))
        self.r, self.N, self.L = rng, N, L
        self.att, self.tor, self.wca = 2 in on, 3 in on, 5 in on
        self.chem = self.att or self.tor
        self.k, self.mur, self.Lc, self.h2 = kappa, mu_r, Lc, (2 * h) ** 2
        self.v0, self.Dr, self.eps, self.rm, self.Mt, self.dt, self.ds, self.ce, self.rc = v0, Dr, eps, rm, Mt, dt, ds, chem_every, rc
        X = _rsa(rng, N, L); self.X, self.Y = X[:, 0].copy(), X[:, 1].copy()
        self.th = rng.uniform(0, 2 * np.pi, N)
        self.JX = np.zeros(N); self.JY = np.zeros(N); self.x0 = None; self.n = 0
        self.wi = self.wj = np.zeros(0, np.intp); self.wsx = self.wsy = np.zeros(0)
        self.R = (rc if self.chem else RW) + SKIN

    def build(self):
        L = self.L; xw = np.mod(np.stack([self.X, self.Y], 1), L); xw[xw >= L] = 0.0
        pr = cKDTree(xw, boxsize=L).query_pairs(min(self.R, L / 2 - 1e-9), output_type="ndarray")
        if len(pr) == 0: pr = np.zeros((0, 2), np.intp)
        self.i, self.j = pr[:, 0].copy(), pr[:, 1].copy()
        self.sx = L * np.round((self.X[self.j] - self.X[self.i]) / L); self.sy = L * np.round((self.Y[self.j] - self.Y[self.i]) / L)
        self.x0 = (self.X.copy(), self.Y.copy())

    def update_pairs(self):
        """every chem_every sub-steps: phoretic current J (Eq. 4) and the contact-pair sub-list."""
        if not (self.chem or self.wca): return
        if self.x0 is None or max(np.abs(self.X - self.x0[0]).max(), np.abs(self.Y - self.x0[1]).max()) > SKIN / 2:
            self.build()
        i, j, N = self.i, self.j, self.N
        dx = self.X[j] - self.X[i]; dx -= self.sx; dy = self.Y[j] - self.Y[i]; dy -= self.sy   # d = r_j - r_i
        r2 = dx * dx + dy * dy
        if self.wca:
            w = np.flatnonzero(r2 < RW * RW); self.wi, self.wj, self.wsx, self.wsy = i[w], j[w], self.sx[w], self.sy[w]
        if self.chem:
            m = np.flatnonzero(r2 < self.rc * self.rc)
            dxm, dym = dx[m], dy[m]; r2m = np.maximum(r2[m], 0.01); r = np.sqrt(r2m); k = self.k
            g = (k + 1.0 / r) * np.exp(-k * r) / r2m                          # direct source
            rs2 = r2m + self.h2; rs = np.sqrt(rs2)
            g += self.Lc * (k + 1.0 / rs) * np.exp(-k * rs) / rs2             # mirror source below the wall
            fx = g * dxm; fy = g * dym; im, jm = i[m], j[m]                   # J_i += g d (toward j), J_j -= g d
            self.JX = np.bincount(im, fx, N) - np.bincount(jm, fx, N)
            self.JY = np.bincount(im, fy, N) - np.bincount(jm, fy, N)

    def step(self):
        if self.n % self.ce == 0: self.update_pairs()
        self.n += 1
        N = self.N; c, s = np.cos(self.th), np.sin(self.th)
        VX = self.v0 * c; VY = self.v0 * s                                    # V_self = -mu_t v0 p, mu_t = -1
        if self.att: VX = VX + 2 * self.v0 * self.JX; VY = VY + 2 * self.v0 * self.JY   # V_int = -mu_t 2 v0 J
        if self.wca and len(self.wi):
            i, j = self.wi, self.wj
            dx = self.X[j] - self.X[i] - self.wsx; dy = self.Y[j] - self.Y[i] - self.wsy; r2 = dx * dx + dy * dy
            w = np.flatnonzero(r2 < self.rm * self.rm)
            if len(w):
                r2w = r2[w]; q6 = (self.rm * self.rm / r2w) ** 3
                f = self.Mt * 12 * self.eps * (q6 * q6 - q6) / r2w              # Mt |dU/dr| / r, repulsive
                fx = f * dx[w]; fy = f * dy[w]; iw, jw = i[w], j[w]
                VX = VX - np.bincount(iw, fx, N) + np.bincount(jw, fx, N)
                VY = VY - np.bincount(iw, fy, N) + np.bincount(jw, fy, N)
        mx = VX * self.dt; my = VY * self.dt; nm = np.sqrt(mx * mx + my * my)
        sc = np.minimum(1.0, self.ds / np.maximum(nm, 1e-300))                # displacement cap ds
        self.X += mx * sc; self.Y += my * sc
        dth = np.sqrt(2 * self.Dr * self.dt) * self.r.standard_normal(N)
        if self.tor: dth += self.dt * (-(9.0 / 4.0) * self.v0 * self.mur * (c * self.JY - s * self.JX))
        self.th += dth

    def frame(self):
        P = np.mod(np.stack([self.X, self.Y], 1), self.L).astype(np.float32); P[P >= self.L] = 0.0
        H = np.mod(self.th, 2 * np.pi).astype(np.float32); H[H >= 2 * np.pi] = 0.0
        return P, H


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    on = base.kept(keep, RULES); N = base.SIZES[size]["NA"]; L = box_side(N)
    nsub = HEADLINE["steps_per_census_step"]; nf = base.n_frames(steps, "pos"); rec = base.REC["pos"]
    P = np.zeros((S, nf, N, 2), np.float32); H = np.zeros((S, nf, N), np.float32); bad = False
    for s, r in enumerate(base.rngs(seed0, S)):
        sy = _Sys(r, N, L, on, **_params())                          # every parameter from HEADLINE
        P[s, 0], H[s, 0] = sy.frame(); fi = 1
        for t in range(1, steps + 1):
            for _ in range(nsub): sy.step()
            if t % rec == 0: P[s, fi], H[s, fi] = sy.frame(); fi += 1
        bad = bad or not np.isfinite(P[s]).all()
    return {"pos": P, "head": H, "at": np.zeros(H.shape, np.int8), "meta": base.meta(steps, base.SIZES[size]["W"], L, bad)}


# ------------------------------------------------------------------ measures (pure numpy/scipy on positions)
def cluster_stats(pos, L, rb=2.1, nmin=10):
    """Eq. (7): S = fraction of particles in bound sets of more than nmin (bound: closer than rb); Q = mean neighbours < rb."""
    x = np.mod(pos.astype(np.float64), L); x[x >= L] = 0.0; N = len(x)
    pr = cKDTree(x, boxsize=L).query_pairs(rb, output_type="ndarray")
    if len(pr) == 0: return 0.0, 0.0
    Q = 2.0 * len(pr) / N
    A = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(N, N))
    _, lab = connected_components(A, directed=False); sz = np.bincount(lab)
    S = float((sz[lab] > nmin).sum()) / N
    return S, Q


def _periodic_label(mask):
    lab, n = ndimage.label(mask)
    if n == 0: return lab, n, np.zeros(1, int)
    par = np.arange(n + 1)

    def find(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for a, b in list(zip(lab[0], lab[-1])) + list(zip(lab[:, 0], lab[:, -1])):
        if a and b:
            ra, rb_ = find(a), find(b)
            if ra != rb_: par[ra] = rb_
    roots = np.array([find(a) for a in range(n + 1)]); lab = roots[lab]
    sizes = np.bincount(lab.ravel(), minlength=n + 1); sizes[0] = 0
    return lab, n, sizes


def cavities(pos, L, spacing=0.5, dfree=2.0, amin=np.pi):
    """enclosed free regions: connected sets of grid points farther than dfree from every centre, other than the largest
    (the outside gas), with area >= amin. Returns (count, total area)."""
    g = int(np.ceil(L / spacing)); sp = L / g; ax = (np.arange(g) + 0.5) * sp
    gx, gy = np.meshgrid(ax, ax, indexing="ij"); pts = np.stack([gx.ravel(), gy.ravel()], 1)
    x = np.mod(pos.astype(np.float64), L); x[x >= L] = 0.0
    d, _ = cKDTree(x, boxsize=L).query(pts, k=1, distance_upper_bound=dfree + 1e-9)
    free = (d > dfree).reshape(g, g)
    lab, n, sizes = _periodic_label(free)
    if sizes.max() == 0: return 0, 0.0
    big = int(np.argmax(sizes)); areas = sizes * sp * sp; areas[big] = 0.0
    cav = areas[areas >= amin]
    return int(len(cav)), float(cav.sum())


def _triad(frames, L, cav_every=1):
    SQ = np.array([cluster_stats(p, L) for p in frames])
    cv = [cavities(p, L) for p in frames[::cav_every]]
    return SQ[:, 0], SQ[:, 1], np.array([c[0] for c in cv]), np.array([c[1] for c in cv])


def presence(out):
    """the paper's triad (<S>, sigma_S, <Q>) and enclosed cavities on the last third of each seed's frames."""
    L = float(out["meta"]["WA"]); res = []
    for s in range(out["pos"].shape[0]):
        P = out["pos"][s]; nf = len(P); late = P[2 * nf // 3:]
        S, Q, nc, ac = _triad(late, L, cav_every=4)
        res.append({"S_mean": round(float(S.mean()), 3), "sigma_S": round(float(S.std()), 3), "Q_mean": round(float(Q.mean()), 2),
                    "cavity_frame_share": round(float((nc > 0).mean()), 3), "cavity_count_mean": round(float(nc.mean()), 2),
                    "boiling_like": bool(S.mean() >= 0.9 and S.std() <= 0.05 and Q.mean() < 5.0 and (nc > 0).mean() >= 0.5)})
    return {"per_seed": res, "boiling_like_seeds": int(sum(r["boiling_like"] for r in res)),
            "paper": {"S": 1.0, "sigma_S": 0.0, "Q": "< 5.0 (full collapse > 5.6)"}}


def fidelity(budget_s=7200, N=1000, T_units=20000, seed=None):
    """the plan's check at the paper's N = 1000 (see module docstring); stops early if the budget runs out."""
    t0 = time.perf_counter(); seed = base.seed_of(CID) if seed is None else seed
    r = base.rngs(seed, 1)[0]; L = box_side(N); on = {2, 3, 5}
    sy = _Sys(r, N, L, on, kappa=HEADLINE["kappa"], mu_r=HEADLINE["mu_r"])
    frames = []; done = 0; deadline = t0 + 0.9 * budget_s
    for t in range(1, T_units + 1):
        for _ in range(SUB_PER_UNIT): sy.step()
        done = t
        if t % 10 == 0: frames.append(sy.frame()[0])
        if time.perf_counter() > deadline: break
    nf = len(frames); late = frames[2 * nf // 3:]
    S, Q, nc, ac = _triad(late, L, cav_every=1)
    m = {"S_mean": round(float(S.mean()), 3), "sigma_S": round(float(S.std()), 3), "Q_mean": round(float(Q.mean()), 2),
         "cavity_frame_share": round(float((nc > 0).mean()), 3), "cavity_count_mean": round(float(nc.mean()), 2),
         "cavity_area_mean": round(float(ac.mean()), 1)}
    ok = m["S_mean"] >= 0.9 and m["sigma_S"] <= 0.05 and m["Q_mean"] < 5.0 and m["cavity_frame_share"] >= 0.5
    return {"setting": {"N": N, "box": round(L, 2), "kappa": HEADLINE["kappa"], "mu_r": HEADLINE["mu_r"], "Lambda_c": 1.0,
                        "Pe": 100, "v0": 2.0, "phi": PHI},
            "time_units_run": int(done), "late_window": [int(done - done // 3), int(done)], "late_frames": len(late),
            "measured": m,
            "paper": {"boiling": {"S": 1.0, "sigma_S": 0.0, "Q": "< 5.0 (Fig. 3b ~4.4)"}, "full_collapse_Q": "> 5.6"},
            "tolerance": "S_mean >= 0.9, sigma_S <= 0.05, Q_mean < 5.0, cavity in >= 50% of late frames",
            "pass": bool(ok), "elapsed_s": round(time.perf_counter() - t0, 1)}

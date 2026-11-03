"""planted:N14 - hard ellipses (axis ratio 2) with two sticky tip patches, grand-canonical Monte Carlo
(arXiv 2610.08522, Sec. II.1-II.2, Sec. III.5, Fig. 4, Fig. 5, ESI Sec. I; patch geometry from the model paper
Wagner & Kahl, Mol. Phys. 2019, arXiv 1902.08097, Sec. 2.1, which the paper cites as Ref. [54]).

Model as the papers state it
----------------------------
* Hard ellipses, aspect ratio kappa = a/b = 2, unit of length sigma = 2 sqrt(ab) (so a = 1/sqrt 2, b = a/2); no
  overlaps allowed (rule 1). Close-packed density rho_cpd = (2/sqrt 3) sigma^-2; reduced density rho* = rho/rho_cpd.
* Two Kern-Frenkel patches at the poles (ends of the long axis), half-opening angle 45 deg, interaction range
  delta = 0.17 sigma, square-well depth epsilon = 1 (energy unit) (rule 2). As in Ref. [54], each patch lives on the
  osculating circle of its vertex (radius r_v = b^2/a = 0.1768, centre M = +-(a - r_v) e = +-0.5303 e along the long
  axis e) and two patches interact "if the line connecting the centres of the related osculating circles
  intersects with the patchy regions". Implemented as the Kern-Frenkel rule on the osculating circles: patch p of i
  and patch q of j are bonded iff |M_jq - M_ip| <= 2 r_v + delta and the vector M_jq - M_ip lies within 45 deg of
  p's outward axis and M_ip - M_jq within 45 deg of q's outward axis. Up to two bonds per patch follow from the
  geometry (not imposed).
* Grand-canonical MC at temperature T* = k_B T/epsilon and chemical potential mu* = mu/(k_B T); a cycle is
  N^cp = 1296 proposed moves (the number of particles at close packing in the box); move probabilities: single
  translation 0.38, single rotation 0.38, cluster translation 0.095, cluster rotation 0.095, deletion 0.025,
  insertion 0.025. Cluster moves use the static linking scheme of Whitelam & Geissler 2007 (Ref. [58]). Maximum
  displacements are tuned to acceptance 0.2-0.5 during the 2e6 start-up cycles (rules 3-4 are one-way).
* Start: dilute gas rho* ~ 0.004 (5 particles in the box) instantly quenched to (T, mu).
* Measure (Sec. III): particle-patch graph -> cyclic subgraph (iterative pruning: degree >= 2 and clustering
  coefficient < 1) -> loops = faces with more than 3 nodes of its embedding (face walking with the rotation system
  of bond angles; for a non-wrapping component its largest face, the outer face, is dropped) -> loop size = number of
  particles on the face; global loop order Gamma_6 = sum_loops w6(l) / M_reg6 with w6(l) = exp(-(l-6)^2/9),
  M_reg6 = 0.1777 * box area (ESI Sec. I); n_L = fraction of particles in loops.

Headline (plan row N14): T = 0.145, mu = -4 - inside the kagome window: Fig. 5(n) at T = 0.145 gives Gamma_6 ~ 0.79,
n_L ~ 0.97, rho* ~ 0.42 at mu = -4 (the maximum of Gamma_6 along the T = 0.145 cut; mu = -5 is equally ordered but
nucleates only during the production run, Fig. 4, whereas mu = -4 assembles during start-up).

Choices the papers do not state (and why)
-----------------------------------------
* Box: square, side sqrt(1296 / rho_cpd) = 33.50 sigma (the paper gives only N^cp = 1296); periodic. Census F view =
  reduced density rho* per cell of a 64 x 64 grid over this box (counts / (cell area * rho_cpd)); area4: side
  67.0, 128 x 128 grid, N^cp = 5184 moves per cycle, 21 starting particles (same rho* = 0.004).
* Thermal wavelength = sigma = 1 and no rotational factor in the insertion/deletion acceptance:
  acc_ins = min(1, V/(N+1) exp(mu* - dU/T)), acc_del = min(1, N/V exp(-mu* - dU/T)).
* Static linking: each bond links with probability p = 1 - exp(-epsilon/T) (Ref. [58]); a rigid cluster move is
  then accepted iff it creates no overlap (for square-well bonds the Boltzmann factor of the external bonds formed
  or broken cancels exactly against the link-failure factors). Cluster rotations are about the cluster centroid and
  are rejected when any member lies farther than L/4 from it (rotation-invariant test, keeps detailed balance and
  excludes clusters wrapping the periodic box). Cluster steps use their own tuned maximum displacements.
* Displacement tuning: every 10 cycles, x1.05 if acceptance > 0.5, x0.95 if < 0.2, throughout the census run (the
  whole census run lies inside the paper's 2e6-cycle start-up, during which the paper tunes).
* Overlap test: exact closed form from the characteristic polynomial f(l) = det(l A - B) of the two ellipse conic
  matrices (Wang-Wang-Kim / Vieillard-Baron type algebraic criterion): with the matrices scaled so that
  f(l) = -l^3 + c2 l^2 + c1 l + 1, the ellipses are separated iff the discriminant
  D = -18 c2 c1 - 4 c2^3 + c2^2 c1^2 + 4 c1^3 - 27 > 0 and not (c2 > 0 and c1 < 0). Checked against a brute-force
  boundary test on 4,000 random pairs (0 disagreements). This gives the same decisions as the Vieillard-Baron
  criterion the paper uses.
* Capacity: at most 2 x N^cp particles; reaching it (possible only in the bonds-without-hard-core sub-program) marks
  the run unstable and stops it. No 50x blow-up rule is applied: the real program condenses from 5 to ~545 particles
  by design.

One census step = one MC cycle (1296 proposed moves) - the plan's rule for Monte Carlo models. 1000 census steps =
1000 cycles; the paper discards 2e6 start-up cycles and samples 1e6-2e7 more, so the census run covers the first
0.05% of the paper's start-up (presence at census scale is reported, not assumed). Development run (1 seed, 1000
cycles, 139 s): N grows 5 -> 333 (rho* 0.26, still rising), 442 bonds (104 triangle and 97 pair motifs), 15 loops of
which 4 six-loops, Gamma_6 = 0.03, n_L = 0.40 - a growing branched porous network, not yet the kagome lattice.
Estimated census cost: ~7 CPU-minutes for the real program (3 seeds) and ~20 s for the knock-out.

Knock-out (plan): rules 1 and 2 removed -> ideal-gas exchange plus free translations/rotations (N ~ V e^mu ~ 21).

Extra output key (not read by census code): "xconf" (S, frames, Nmax, 3) float32 = (x, y, theta) of every
particle at cycle 0, every 100 cycles and at the last cycle, NaN-padded, so presence() can compute the paper's
own loop measure on the census run (the census field view carries no orientations or bonds).

Fidelity (plan): ring order after the quench from the dilute gas at the headline, paper's box, as many cycles as fit
the budget; measure Gamma_6, n_L, rho*, 6-loop share c6 at the end (mean of the last 10% of samples). Tolerance,
fixed before any large run: PASS if Gamma_6 >= 0.59 and n_L >= 0.77 (within 0.2 of the paper's 0.79 and 0.97).
A short diagnostic also starts from a perfect kagome network (168 particles: 7 x 8 honeycomb cells stretched by 1%
onto a square box) at the headline and reports whether Gamma_6 stays high (tests the potential and acceptance rules;
not the plan's check; in development runs it stayed at 0.974 with all 56 six-loops for 30 cycles).
"""
import math
import time
import numpy as np
from scipy.spatial import cKDTree

from census.planted_models import base

CID = "planted:N14"
LAYERS = "F"
K = {}
SPEED = 0.0
RULES = [("HARD_CORE", True), ("PATCH_BONDS", True), ("RESERVOIR_EXCHANGE", False), ("TRANSLATE_ROTATE", False)]
HEADLINE = dict(
    T=0.145,          # Fig. 5 cross-section T = 0.145
    mu=-4.0,          # Fig. 5(n): Gamma_6 ~ 0.79, n_L ~ 0.97, rho* ~ 0.42 at mu = -4 (kagome window -5..-2)
    kappa=2.0,        # Sec. II.1
    theta_patch=45.0, # Sec. II.1 half-opening angle (deg)
    delta=0.17,       # Sec. II.1 interaction range (sigma)
    eps=1.0,          # Sec. II.1 square-well depth
    Ncp=1296,         # Sec. II.2 particles at close packing = proposed moves per cycle
    rho0=0.004,       # Sec. II.2 initial reduced density
    probs=(0.38, 0.38, 0.095, 0.095, 0.025, 0.025),   # Sec. II.2: t, r, tc, rc, del, ins
)
STEP = "one census step = one MC cycle = 1296 proposed moves (paper's cycle; plan: one sweep for Monte Carlo models)"
RULE_TEXT = ("World: a 2D surface exchanging particles with a reservoir (grand-canonical Monte Carlo; about 1,300 particles at "
             "close packing). Particles are hard elliptical disks with axis ratio 2 that cannot overlap; each carries two "
             "sticky patches, one at each tip of its long axis (patch half-angle 45 degrees, reach 0.17 particle "
             "diameters, bond energy -1), allowing up to two bonds per patch. Rules: Monte Carlo moves (single and cluster "
             "translations and rotations, insertions, deletions) at temperature T and chemical potential mu. Start from a "
             "dilute gas (density 0.004) quenched to T=0.145 with mu = -4, run for 1,000 Monte Carlo cycles. Track bonds "
             "and closed rings of bonded particles.")

# geometry (sigma = 1)
_A = 1.0 / math.sqrt(2.0); _B = _A / 2.0
_RV = _B * _B / _A                     # osculating-circle radius at the vertex
_MO = _A - _RV                         # patch-circle centre offset along the axis
_RHO_CPD = 2.0 / math.sqrt(3.0)
_IA2, _IB2 = 1.0 / _A ** 2, 1.0 / _B ** 2
_AB2 = _A * _A * _B * _B
_R_OV2 = (2 * _A) ** 2                 # beyond this centre distance no overlap
_R_SURE2 = (2 * _B) ** 2               # below this centre distance certain overlap
_W6_DELTA = 9.0
_ETA_REG6 = 0.1777                     # ESI Sec. I


def _derived(P):
    db = 2 * _RV + P["delta"]
    return dict(db2=db * db, cth=math.cos(math.radians(P["theta_patch"])), rint=2 * _MO + db)


# ---------------------------------------------------------------------------------------------------- pair tests
def _separated(dx, dy, c1, s1, c2, s2):
    """exact separation test for ellipse 1 at the origin (axis (c1,s1)) and ellipse 2 at (dx,dy) (axis (c2,s2))."""
    u = c1 * dx + s1 * dy; v = -s1 * dx + c1 * dy
    cp = c1 * c2 + s1 * s2; sp = c1 * s2 - s1 * c2            # relative angle
    B11 = cp * cp * _IA2 + sp * sp * _IB2; B22 = sp * sp * _IA2 + cp * cp * _IB2; B12 = cp * sp * (_IA2 - _IB2)
    g1 = B11 * u + B12 * v; g2 = B12 * u + B22 * v; B33 = u * g1 + v * g2 - 1.0
    k2 = -_AB2 * (-B11 * _IB2 - B22 * _IA2 + B33 * _IA2 * _IB2)
    k1 = _AB2 * ((B22 * B33 - g2 * g2) * _IA2 + (B11 * B33 - g1 * g1) * _IB2 - _IA2 * _IB2)
    D = -18 * k2 * k1 - 4 * k2 ** 3 + k2 * k2 * k1 * k1 + 4 * k1 ** 3 - 27
    return (D > 0) & ~((k2 > 0) & (k1 < 0))


_PI = np.array([1.0, 1.0, -1.0, -1.0])[:, None]; _PJ = np.array([1.0, -1.0, 1.0, -1.0])[:, None]


def _bond_counts(dx, dy, ci, si, cj, sj, db2, cth):
    """number of bonded patch pairs between particle i (at origin, axis (ci,si)) and particles j at (dx,dy);
    the four patch combinations (+/-e_i, +/-e_j) are evaluated together."""
    ux = dx + _MO * (_PJ * cj - _PI * ci); uy = dy + _MO * (_PJ * sj - _PI * si)      # M_jq - M_ip, shape (4, n)
    r = np.sqrt(ux * ux + uy * uy)
    ok = (r * r <= db2) & ((ux * ci + uy * si) * _PI >= cth * r) & (-(ux * cj + uy * sj) * _PJ >= cth * r)
    return ok.sum(0)


# ---------------------------------------------------------------------------------------------------- GCMC engine
class _GCMC:
    def __init__(self, P, on, rng, L, Ncp, N0, cap=None):
        self.P = P; self.hard = 1 in on; self.bond = 2 in on; self.rng = rng
        self.L = float(L); self.V = self.L * self.L; self.Ncp = int(Ncp)
        d = _derived(P); self.db2, self.cth, self.rint = d["db2"], d["cth"], d["rint"]
        self.beta = 1.0 / P["T"]; self.mu = P["mu"]
        self.plink = 1.0 - math.exp(-P["eps"] * self.beta)
        self.cap = int(cap or 2 * Ncp)
        self.x = np.zeros((self.cap, 2)); self.th = np.zeros(self.cap); self.c = np.ones(self.cap); self.s = np.zeros(self.cap)
        self.N = 0
        rc = max(self.rint, 2 * _A) if (self.hard or self.bond) else self.L / 4
        self.nc = max(3, int(self.L // rc)); self.cs = self.L / self.nc
        self.cells = [set() for _ in range(self.nc * self.nc)]
        self.cell = np.zeros(self.cap, dtype=np.int64)
        nc = self.nc
        self.nbr = [tuple(((cx + ox) % nc) * nc + (cy + oy) % nc for ox in (-1, 0, 1) for oy in (-1, 0, 1))
                    for cx in range(nc) for cy in range(nc)]
        self.adj = [dict() for _ in range(self.cap)]
        self.nb = np.zeros(self.cap, dtype=np.int64)
        self.step = {"t": 0.1, "r": 0.2, "tc": 0.1, "rc": 0.1}
        self.acc = {k: [0, 0] for k in self.step}
        self.overflow = False
        while self.N < N0:                                   # dilute random start without overlaps
            p = rng.uniform(0, self.L, 2); t = rng.uniform(0, 2 * math.pi)
            J = self._cands(p, -1)
            if self.hard and self._overlaps(p, math.cos(t), math.sin(t), J): continue
            self._add(p, t)
            self._set_bonds(self.N - 1, J)

    # ----------------------------------------------------------------- bookkeeping
    def _cell_of(self, p):
        return (int(p[0] / self.cs) % self.nc) * self.nc + int(p[1] / self.cs) % self.nc

    def _cands(self, p, excl, exset=None):
        J = []
        for cid in self.nbr[self._cell_of(p)]: J.extend(self.cells[cid])
        if exset is not None: J = [j for j in J if j not in exset]
        elif excl >= 0: J = [j for j in J if j != excl]
        return np.fromiter(J, dtype=np.int64, count=len(J))

    def _rel(self, p, J):
        d = self.x[J] - p; d -= self.L * np.round(d / self.L)
        return d[:, 0], d[:, 1]

    def _overlaps(self, p, c, s, J):
        if len(J) == 0: return False
        dx, dy = self._rel(p, J); r2 = dx * dx + dy * dy
        if np.any(r2 < _R_SURE2): return True
        m = r2 < _R_OV2
        if not np.any(m): return False
        return bool(np.any(~_separated(dx[m], dy[m], c, s, self.c[J[m]], self.s[J[m]])))

    def _bonds(self, p, c, s, J):
        if len(J) == 0 or not self.bond: return J[:0], np.zeros(0, np.int64)
        dx, dy = self._rel(p, J); r2 = dx * dx + dy * dy
        m = r2 < self.rint * self.rint
        if not np.any(m): return J[:0], np.zeros(0, np.int64)
        Jm = J[m]; n = _bond_counts(dx[m], dy[m], c, s, self.c[Jm], self.s[Jm], self.db2, self.cth)
        k = n > 0
        return Jm[k], n[k]

    def _add(self, p, t):
        i = self.N; self.N += 1
        self.x[i] = p; self.th[i] = t; self.c[i] = math.cos(t); self.s[i] = math.sin(t)
        cid = self._cell_of(p); self.cell[i] = cid; self.cells[cid].add(i)
        self.adj[i] = {}; self.nb[i] = 0

    def _unbond(self, i):
        for j, k in self.adj[i].items():
            del self.adj[j][i]; self.nb[j] -= k
        self.adj[i] = {}; self.nb[i] = 0

    def _set_bonds(self, i, J, partners=None):
        """recompute all bonds of particle i against candidates J (i's current state)."""
        self._unbond(i)
        if partners is None: partners = self._bonds(self.x[i], self.c[i], self.s[i], J)
        Jb, nb = partners
        for j, k in zip(Jb.tolist(), nb.tolist()):
            self.adj[i][j] = k; self.adj[j][i] = k; self.nb[j] += k
        self.nb[i] = int(nb.sum()) if len(nb) else 0

    def _place(self, i, p, t):
        self.x[i] = p; self.th[i] = t; self.c[i] = math.cos(t); self.s[i] = math.sin(t)
        cid = self._cell_of(p)
        if cid != self.cell[i]:
            self.cells[self.cell[i]].discard(i); self.cells[cid].add(i); self.cell[i] = cid

    def _remove(self, i):
        self._unbond(i); self.cells[self.cell[i]].discard(i)
        last = self.N - 1
        if i != last:
            self.cells[self.cell[last]].discard(last)
            self.x[i] = self.x[last]; self.th[i] = self.th[last]; self.c[i] = self.c[last]; self.s[i] = self.s[last]
            self.cell[i] = self.cell[last]; self.cells[self.cell[i]].add(i)
            self.adj[i] = self.adj[last]; self.nb[i] = self.nb[last]
            for j, k in self.adj[i].items():
                del self.adj[j][last]; self.adj[j][i] = k
            self.adj[last] = {}; self.nb[last] = 0
        self.N -= 1

    # ----------------------------------------------------------------- moves
    def _metro(self, dE):
        return dE <= 0 or self.rng.random() < math.exp(-self.beta * dE)

    def _single(self, kind):
        if self.N == 0: return
        rng = self.rng; i = int(rng.integers(self.N)); st = self.step[kind]; self.acc[kind][1] += 1
        if kind == "t":
            p = self.x[i] + rng.uniform(-st, st, 2); p %= self.L; p[p >= self.L] = 0.0; t = self.th[i]
        else:
            p = self.x[i]; t = self.th[i] + rng.uniform(-st, st)
        c, s = math.cos(t), math.sin(t)
        J = self._cands(p, i)
        if self.hard and self._overlaps(p, c, s, J): return
        Jb, nb = self._bonds(p, c, s, J)
        n_new = int(nb.sum()) if len(nb) else 0
        if not self._metro(-self.P["eps"] * (n_new - int(self.nb[i]))): return
        self._place(i, p.copy(), t); self._set_bonds(i, J, (Jb, nb)); self.acc[kind][0] += 1

    def _insert(self):
        if self.N >= self.cap: self.overflow = True; return
        rng = self.rng; p = rng.uniform(0, self.L, 2); t = rng.uniform(0, 2 * math.pi)
        c, s = math.cos(t), math.sin(t)
        J = self._cands(p, -1)
        if self.hard and self._overlaps(p, c, s, J): return
        Jb, nb = self._bonds(p, c, s, J)
        n_new = int(nb.sum()) if len(nb) else 0
        if rng.random() < self.V / (self.N + 1) * math.exp(self.mu + self.beta * self.P["eps"] * n_new):
            self._add(p, t); self._set_bonds(self.N - 1, J, (Jb, nb))

    def _delete(self):
        if self.N == 0: return
        rng = self.rng; i = int(rng.integers(self.N))
        if rng.random() < self.N / self.V * math.exp(-self.mu - self.beta * self.P["eps"] * int(self.nb[i])):
            self._remove(i)

    def _cluster_of(self, i):
        """static-linking cluster (BFS) and unwrapped member positions relative to particle i."""
        if not self.adj[i]: return np.array([i], dtype=np.int64), np.zeros((1, 2)), {i: 0}
        rng = self.rng; mem = [i]; seen = {i: 0}; off = [np.zeros(2)]; k = 0
        while k < len(mem):
            u = mem[k]; ou = off[k]; k += 1
            for v, m in self.adj[u].items():
                if v in seen: continue
                if rng.random() < 1.0 - (1.0 - self.plink) ** m:
                    d = self.x[v] - self.x[u]; d -= self.L * np.round(d / self.L)
                    seen[v] = len(mem); mem.append(v); off.append(ou + d)
        return np.array(mem, dtype=np.int64), np.array(off), seen

    def _cluster(self, kind):
        if self.N == 0: return
        rng = self.rng; i = int(rng.integers(self.N)); st = self.step[kind]; self.acc[kind][1] += 1
        if self.bond: mem, off, seen = self._cluster_of(i)
        else: mem, off, seen = np.array([i]), np.zeros((1, 2)), {i: 0}
        n = len(mem)
        if kind == "tc":
            newp = (self.x[mem] + rng.uniform(-st, st, 2)) % self.L; newt = self.th[mem]
        else:
            cen = off.mean(0); rr = off - cen
            if n > 1 and np.max(np.einsum("ij,ij->i", rr, rr)) > (self.L / 4) ** 2: return
            al = rng.uniform(-st, st); ca, sa = math.cos(al), math.sin(al)
            rot = np.stack([ca * rr[:, 0] - sa * rr[:, 1], sa * rr[:, 0] + ca * rr[:, 1]], 1)
            newp = (self.x[i] + cen + rot) % self.L; newt = self.th[mem] + al
        newp[newp >= self.L] = 0.0
        nc_, ns_ = np.cos(newt), np.sin(newt)
        # overlap of moved members with non-members
        if self.hard:
            if n <= 12:
                for k in range(n):
                    J = self._cands(newp[k], -1, seen)
                    if self._overlaps(newp[k], nc_[k], ns_[k], J): return
            else:
                others = np.setdiff1d(np.arange(self.N), mem, assume_unique=False)
                if len(others):
                    pr = cKDTree(newp, boxsize=self.L).query_ball_tree(cKDTree(self.x[others] % self.L, boxsize=self.L), 2 * _A)
                    ii = np.fromiter((k for k, lst in enumerate(pr) for _ in lst), dtype=np.int64)
                    jj = others[np.fromiter((j for lst in pr for j in lst), dtype=np.int64)] if len(ii) else ii
                    if len(ii):
                        d = self.x[jj] - newp[ii]; d -= self.L * np.round(d / self.L)
                        r2 = d[:, 0] ** 2 + d[:, 1] ** 2
                        if np.any(r2 < _R_SURE2): return
                        m = r2 < _R_OV2
                        if np.any(m) and np.any(~_separated(d[m, 0], d[m, 1], nc_[ii[m]], ns_[ii[m]], self.c[jj[m]], self.s[jj[m]])): return
        # accepted (static linking with square-well bonds): move, then refresh bonds with non-members
        for k in range(n):
            self._place(int(mem[k]), newp[k].copy(), float(newt[k]))
        if self.bond:
            for k in range(n):
                u = int(mem[k])
                for v in [v for v in self.adj[u] if v not in seen]:
                    kk = self.adj[u].pop(v); del self.adj[v][u]; self.nb[u] -= kk; self.nb[v] -= kk
            for k in range(n):
                u = int(mem[k]); J = self._cands(self.x[u], -1, seen)
                Jb, nb = self._bonds(self.x[u], self.c[u], self.s[u], J)
                for v, kk in zip(Jb.tolist(), nb.tolist()):
                    self.adj[u][v] = kk; self.adj[v][u] = kk; self.nb[u] += kk; self.nb[v] += kk
        self.acc[kind][0] += 1

    def cycle(self):
        pr = np.cumsum(self.P["probs"]); pr /= pr[-1]
        kinds = np.searchsorted(pr, self.rng.random(self.Ncp), side="right")
        for k in kinds.tolist():
            if k == 0: self._single("t")
            elif k == 1: self._single("r")
            elif k == 2: self._cluster("tc")
            elif k == 3: self._cluster("rc")
            elif k == 4: self._delete()
            else: self._insert()
            if self.overflow: return

    def tune(self):
        for k, (a, n) in self.acc.items():
            if n >= 20:
                f = a / n
                if f > 0.5: self.step[k] = min(self.step[k] * 1.05, self.L / 4 if k in ("t", "tc") else math.pi)
                elif f < 0.2: self.step[k] *= 0.95
            self.acc[k] = [0, 0]

    def density_field(self, W):
        x = self.x[:self.N]
        if self.N == 0: return np.zeros((W, W), np.float32)
        ix = np.minimum((x[:, 0] * (W / self.L)).astype(np.int64), W - 1); iy = np.minimum((x[:, 1] * (W / self.L)).astype(np.int64), W - 1)
        cnt = np.bincount(iy * W + ix, minlength=W * W).reshape(W, W)
        return (cnt / ((self.L / W) ** 2 * _RHO_CPD)).astype(np.float32)   # rows = y, columns = x

    def conf(self):
        return np.concatenate([self.x[:self.N], self.th[:self.N, None]], 1)

    def load(self, conf):
        """replace the state by a given configuration (x, y, theta); bonds recomputed."""
        while self.N: self._remove(self.N - 1)
        for p, t in zip(conf[:, :2], conf[:, 2]): self._add(np.mod(p, self.L), float(t))
        for i in range(self.N): self._set_bonds(i, self._cands(self.x[i], i))
        return self


def _params(**over):
    P = dict(HEADLINE); P.update(over); return P


def _box(size, P):
    W = base.SIZES[size]["W"]; f = (W // 64) ** 2
    Ncp = P["Ncp"] * f; L = math.sqrt(Ncp / _RHO_CPD)
    return W, L, Ncp, int(round(P["rho0"] * Ncp))


# ---------------------------------------------------------------------------------------------------- census entry
def simulate(seed0, steps=1000, S=3, size="census", keep=None, **over):
    on = base.kept(keep, RULES)
    P = _params(**over)
    W, L, Ncp, N0 = _box(size, P)
    rs = base.rngs(seed0, S)
    sims = [_GCMC(P, on, r, L, Ncp, N0) for r in rs]
    F = [[s.density_field(W)] for s in sims]; X = [[s.conf()] for s in sims]
    Ns = [[s.N] for s in sims]
    unstable = False; t = 0
    while t < steps:
        for s in sims:
            s.cycle()
            if s.overflow: unstable = True
        t += 1
        if t % 10 == 0:
            for s in sims: s.tune()
        if t % base.REC["F"] == 0:
            for k, s in enumerate(sims): F[k].append(s.density_field(W)); Ns[k].append(s.N)
        if t % 100 == 0 or t == steps or unstable:
            for k, s in enumerate(sims): X[k].append(s.conf())
        if unstable: break
    out = {"F": np.stack([np.stack(f) for f in F])[:, :, None].astype(np.float32)}
    nmax = max(len(c) for xs in X for c in xs); nconf = max(len(xs) for xs in X)
    xc = np.full((S, nconf, max(nmax, 1), 3), np.nan, np.float32)
    for k, xs in enumerate(X):
        for f, c in enumerate(xs): xc[k, f, :len(c)] = c
    out["xconf"] = xc
    out["meta"] = base.meta(t, W, unstable=unstable)
    out["meta"].update(domain=float(L), N_frames=[list(map(int, n)) for n in Ns],
                       steps_final={k: float(v) for k, v in sims[0].step.items()})
    return out


# ---------------------------------------------------------------------------------------------------- loop analysis
def bonds_of(conf, L, P=None):
    """all bonded patch pairs (i, si, j, sj) of a configuration (x, y, theta); si = 0 for +e, 1 for -e."""
    P = P or HEADLINE; d = _derived(P)
    conf = conf[np.isfinite(conf[:, 0])]
    x = np.mod(conf[:, :2].astype(float), L); x[x >= L] = 0.0; th = conf[:, 2].astype(float)
    c, s = np.cos(th), np.sin(th)
    if len(x) < 2: return np.zeros((0, 4), np.int64), x, th
    pr = cKDTree(x, boxsize=L).query_pairs(d["rint"], output_type="ndarray")
    out = []
    if len(pr):
        i, j = pr[:, 0], pr[:, 1]
        dv = x[j] - x[i]; dv -= L * np.round(dv / L)
        for a_, pi in ((0, 1.0), (1, -1.0)):
            for b_, pj in ((0, 1.0), (1, -1.0)):
                ux = dv[:, 0] + _MO * (pj * c[j] - pi * c[i]); uy = dv[:, 1] + _MO * (pj * s[j] - pi * s[i])
                r = np.hypot(ux, uy)
                ok = (r * r <= d["db2"]) & ((ux * pi * c[i] + uy * pi * s[i]) >= d["cth"] * r) & \
                     (-(ux * pj * c[j] + uy * pj * s[j]) >= d["cth"] * r)
                k = np.nonzero(ok)[0]
                out.append(np.stack([i[k], np.full(len(k), a_), j[k], np.full(len(k), b_)], 1))
    return (np.concatenate(out) if out else np.zeros((0, 4), np.int64)), x, th


def loop_analysis(conf, L, P=None):
    """paper Sec. III: particle-patch graph -> cyclic subgraph -> faces -> loops; returns Gamma_6, n_L, c_l, motifs."""
    bnd, x, th = bonds_of(conf, L, P)
    N = len(x)
    res = {"N": int(N), "rho": float(N / (L * L) / _RHO_CPD), "n_bonds": int(len(bnd))}
    if N == 0:
        res.update(Gamma6=0.0, n_L=0.0, n_loops=0, c6=float("nan"), loop_sizes={}); return res
    e = np.stack([np.cos(th), np.sin(th)], 1)
    pos = {i: x[i] for i in range(N)}
    for i in range(N):
        pos[N + 2 * i] = x[i] + _MO * e[i]; pos[N + 2 * i + 1] = x[i] - _MO * e[i]
    G = {k: set() for k in pos}
    for i in range(N):
        for q in (N + 2 * i, N + 2 * i + 1): G[i].add(q); G[q].add(i)
    for i, si, j, sj in bnd.tolist():
        p, q = N + 2 * i + si, N + 2 * j + sj
        if p != q: G[p].add(q); G[q].add(p)
    # bonding motifs: connected components of the patch-only subgraph
    seen = set(); motifs = {}
    for p in range(N, 3 * N):
        if p in seen: continue
        comp = [p]; seen.add(p); k = 0
        while k < len(comp):
            u = comp[k]; k += 1
            for v in G[u]:
                if v >= N and v not in seen: seen.add(v); comp.append(v)
        ne = sum(1 for u in comp for v in G[u] if v >= N) // 2
        key = {1: "single", 2: "pair"}.get(len(comp), None)
        if len(comp) == 3: key = "triangle" if ne == 3 else "open_triangle"
        key = key or f"size{len(comp)}"
        motifs[key] = motifs.get(key, 0) + 1
    res["motifs"] = motifs
    # cyclic subgraph by iterative pruning (degree >= 2 and clustering coefficient < 1)
    H = {k: set(v) for k, v in G.items()}
    changed = True
    while changed:
        changed = False; drop = []
        for u, nb in H.items():
            k = len(nb)
            if k < 2: drop.append(u); continue
            nbl = list(nb); e_in = sum(1 for a_ in range(k) for b_ in range(a_ + 1, k) if nbl[b_] in H[nbl[a_]])
            if e_in == k * (k - 1) // 2: drop.append(u)
        for u in drop:
            for v in H[u]: H[v].discard(u)
            del H[u]; changed = True
    if not H:
        res.update(Gamma6=0.0, n_L=0.0, n_loops=0, c6=float("nan"), loop_sizes={}); return res

    def vec(u, v):
        d = pos[v] - pos[u]; return d - L * np.round(d / L)
    rot = {u: sorted(nb, key=lambda v: math.atan2(*vec(u, v)[::-1])) for u, nb in H.items()}
    idx = {u: {v: k for k, v in enumerate(r)} for u, r in rot.items()}
    # components and wrapping test
    comp_of = {}; comps = []
    for u0 in H:
        if u0 in comp_of: continue
        cid = len(comps); comp_of[u0] = cid; un = {u0: np.zeros(2)}; stack = [u0]; wraps = False; members = [u0]
        while stack:
            u = stack.pop()
            for v in H[u]:
                pv = un[u] + vec(u, v)
                if v not in un:
                    un[v] = pv; comp_of[v] = cid; stack.append(v); members.append(v)
                elif np.max(np.abs(un[v] - pv)) > L / 2: wraps = True
        comps.append((members, wraps))
    # face walking
    visited = set(); faces = []
    for u in H:
        for v in H[u]:
            if (u, v) in visited: continue
            face = []; a_, b_ = u, v
            while (a_, b_) not in visited:
                visited.add((a_, b_)); face.append(a_)
                r = rot[b_]; w = r[(idx[b_][a_] - 1) % len(r)]
                a_, b_ = b_, w
            faces.append(face)
    loops = []
    by_comp = {}
    for f in faces:
        by_comp.setdefault(comp_of[f[0]], []).append(f)
    for cid, fs in by_comp.items():
        if not comps[cid][1] and fs:
            big = max(range(len(fs)), key=lambda k: len(fs[k])); fs = [f for k, f in enumerate(fs) if k != big]
        for f in fs:
            if len(f) > 3:
                cores = {u for u in f if u < N}
                if cores: loops.append(cores)
    sizes = [len(cs) for cs in loops]
    w6 = sum(math.exp(-(l - 6) ** 2 / _W6_DELTA) for l in sizes)
    in_loop = set().union(*loops) if loops else set()
    hist = {}
    for l in sizes: hist[l] = hist.get(l, 0) + 1
    res.update(Gamma6=float(w6 / (_ETA_REG6 * L * L)), n_L=float(len(in_loop) / N), n_loops=len(sizes),
               c6=float(hist.get(6, 0) / len(sizes)) if sizes else float("nan"),
               loop_sizes={int(k): int(v) for k, v in sorted(hist.items())})
    return res


def presence(out):
    """paper's loop measures on the census run's stored configurations (last stored frame of each seed) plus the
    mean reduced density of the last 10 F frames."""
    L = out["meta"]["domain"]; xc = out["xconf"]; S = xc.shape[0]
    per = [loop_analysis(xc[s, -1], L) for s in range(S)]
    rho_F = float(np.mean(out["F"][:, -10:])) if out["F"].shape[1] else float("nan")
    keys = ("Gamma6", "n_L", "rho", "c6")
    agg = {k: (float(np.nanmean(v)) if np.any(np.isfinite(v)) else float("nan"))
           for k in keys for v in [np.array([p[k] for p in per], float)]}
    agg.update(per_seed=per, rho_from_F=rho_F,
               paper="Fig. 5(n), T=0.145, mu=-4: Gamma6 ~0.79, n_L ~0.97, rho ~0.42 (after >= 2e6 cycles)")
    return agg


# ---------------------------------------------------------------------------------------------------- fidelity
def _kagome(nx=7, ny=4, t=0.06):
    """perfect network: ellipses on the edges of a honeycomb (edge 2a + 2t/sqrt3; ESI Sec. I geometry), nx x 2ny
    honeycomb cells stretched onto a square box (nx/ny ~ sqrt 3; 7 x 4 gives 168 particles, 1% stretch along x).
    t = 0.05 is exact tip contact in the symmetric triangle; 0.06 leaves a margin."""
    le = 2 * _A + 2 * t / math.sqrt(3)
    a1 = np.array([math.sqrt(3) * le, 0.0]); a2 = np.array([math.sqrt(3) * le / 2, 1.5 * le])
    pts = []
    for i in range(nx):
        for j in range(2 * ny):
            o = i * a1 + j * a2
            for k in range(3):                          # the 3 edges from honeycomb vertex A (at o) to its B neighbours
                phi = math.pi / 2 + 2 * math.pi * k / 3
                pts.append([o[0] + 0.5 * le * math.cos(phi), o[1] + 0.5 * le * math.sin(phi), phi])
    pts = np.array(pts)
    Lx = nx * a1[0]; Ly = 2 * ny * a2[1]; L = max(Lx, Ly)
    pts[:, 0] = np.mod(pts[:, 0], Lx) * (L / Lx); pts[:, 1] = np.mod(pts[:, 1], Ly) * (L / Ly)
    return pts, L


def fidelity(budget_s=7200, diag_cycles=200):
    """quench from the dilute gas at the headline in the paper's box for as many cycles as fit; Gamma_6 at the end."""
    P = HEADLINE; t_start = time.time()
    # diagnostic (not the plan's check): a perfect network at the headline should keep its six-loops
    pts, Lk = _kagome(7, 4)
    dk = _GCMC(P, {1, 2}, np.random.default_rng(7), Lk, int(round(Lk * Lk * _RHO_CPD)), 0).load(pts)
    k0 = loop_analysis(dk.conf(), Lk)
    for c in range(diag_cycles):
        dk.cycle()
        if c % 10 == 9: dk.tune()
    k1 = loop_analysis(dk.conf(), Lk)
    diag = {"N_start": k0["N"], "Gamma6_start": round(k0["Gamma6"], 3), "cycles": diag_cycles, "N_end": k1["N"],
            "Gamma6_end": round(k1["Gamma6"], 3), "n_L_end": round(k1["n_L"], 3), "loop_sizes_end": k1["loop_sizes"]}
    W, L, Ncp, N0 = _box("census", P)
    sim = _GCMC(P, {1, 2}, np.random.default_rng(2610), L, Ncp, N0)
    samples = []; cyc = 0
    while time.time() - t_start < 0.9 * budget_s:
        sim.cycle(); cyc += 1
        if cyc % 10 == 0: sim.tune()
        if cyc % 200 == 0:
            r = loop_analysis(sim.conf(), L); r["cycle"] = cyc; samples.append(r)
    if not samples or samples[-1]["cycle"] != cyc:
        r = loop_analysis(sim.conf(), L); r["cycle"] = cyc; samples.append(r)
    tail = samples[-max(1, len(samples) // 10):]
    fin = {k: float(np.nanmean([s[k] for s in tail])) for k in ("Gamma6", "n_L", "rho", "c6")}
    res = {"setting": "T=0.145, mu=-4, square box L=%.2f (N^cp=1296), start rho*=0.004" % L, "cycles": cyc,
           "paper_cycles": ">= 2e6 start-up + 1e6..2e7 production", "final": fin,
           "trajectory": [{k: s[k] for k in ("cycle", "N", "rho", "Gamma6", "n_L", "c6")} for s in samples],
           "paper": "Fig. 5(n): Gamma6 ~0.79, n_L ~0.97, rho ~0.42",
           "tolerance": "PASS if Gamma6 >= 0.59 and n_L >= 0.77", "kagome_stability_diagnostic": diag}
    res["pass"] = bool(fin["Gamma6"] >= 0.59 and fin["n_L"] >= 0.77)
    return res

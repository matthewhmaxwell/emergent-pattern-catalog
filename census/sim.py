"""Interpreter for G1 programs (pilot draft v0; frozen at registration).

All seeds of one program run together (leading batch axis S). One time step = every rule in program order (each
rule reads the state left by the previous rule and updates synchronously), then agent heading noise + movement, then
field diffusion. A frame of the full state is recorded every REC steps. The run stops early when the state stops
changing (absorbing / frozen), which is recorded.

Fixed initial conditions (symmetric under type relabeling):
  cells, nodes, agents: types uniform over k; agents at random cell centres with random cardinal headings
  network: Erdos-Renyi, mean degree 4
  fields: field 0 = 1, field 1 = 0, with 20 random 3x3 patches set to (0.5, 0.25), plus 1% noise
  phases (only when a KURA rule is present): uniform on [0, 2pi); natural frequencies ~ N(0, 0.05)
Fixed semantics worth stating: C.SAND uses open boundaries (grains leave at the lattice edge, as in BTW) while every
other lattice rule sees a torus; A.QUORUM slows crowded agents to QUORUM_SLOW x speed (quorum-sensing motility).
"""
import numpy as np
from scipy.spatial import cKDTree

W, NA, NN, T_STEPS, REC, SEEDS = 64, 400, 200, 1000, 10, 3
SAND_MAXIT, FREEZE_WINDOW, QUORUM_SLOW = 20000, 100, 0.1
MOORE = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
VN = [(-1, 0), (1, 0), (0, -1), (0, 1)]


def _roll(x, dy, dx):
    return np.roll(np.roll(x, dy, axis=-2), dx, axis=-1)


def _count(grid, val):
    """number of Moore neighbours equal to val (val: scalar or array broadcastable to grid)."""
    eq = (grid == val) if np.isscalar(val) else None
    if eq is not None:
        return sum(_roll(eq, dy, dx) for dy, dx in MOORE).astype(np.int16)
    return sum((_roll(grid, dy, dx) == val) for dy, dx in MOORE).astype(np.int16)


class State:
    def __init__(self, prog, seed0=0, S=SEEDS):
        self.p, self.S = prog, S
        self.rng = np.random.default_rng(1_000_003 * seed0 + 17)
        r, k = self.rng, prog.k
        names = {x.tmpl for x in prog.rules}
        if "C" in prog.layers:
            self.C = r.integers(0, k["C"], size=(S, W, W)).astype(np.int8)
            if "C.KURA" in names:
                self.Cph = r.uniform(0, 2 * np.pi, (S, W, W)); self.Com = r.normal(0, 0.05, (S, W, W))
            self.aval = []                                     # sandpile avalanche sizes per step
        if "A" in prog.layers:
            self.pos = r.integers(0, W, size=(S, NA, 2)) + 0.5
            self.head = r.integers(0, 4, size=(S, NA)) * (np.pi / 2)
            self.at = r.integers(0, k["A"], size=(S, NA)).astype(np.int8)
        if "N" in prog.layers:
            up = np.triu(r.random((S, NN, NN)) < 4.0 / (NN - 1), 1); self.adj = up | up.transpose(0, 2, 1)
            self.nt = r.integers(0, k["N"], size=(S, NN)).astype(np.int8)
            if "N.KURA" in names:
                self.Nph = r.uniform(0, 2 * np.pi, (S, NN)); self.Nom = r.normal(0, 0.05, (S, NN))
        if "F" in prog.layers:
            self.F = np.zeros((S, prog.nf, W, W)); self.F[:, 0] = 1.0
            for s in range(S):
                for _ in range(20):
                    y, x = r.integers(0, W - 3, 2); self.F[s, 0, y:y + 3, x:x + 3] = 0.5
                    if prog.nf == 2: self.F[s, 1, y:y + 3, x:x + 3] = 0.25
            self.F += 0.01 * r.random(self.F.shape)
        self.t = 0; self.frozen_at = None; self.unstable = False

    # ------------------------------------------------------------------ helpers
    def _bern(self, shape, p):
        return np.ones(shape, bool) if p >= 1.0 else (self.rng.random(shape) < p)

    def _pairs(self, r):
        """per seed: (i, j, dij) for unordered neighbour pairs within r; dij = pos_j - pos_i on the torus."""
        key = (self.t, r)
        if getattr(self, "_pkey", None) != key:
            out = []
            for s in range(self.S):
                pos = np.mod(self.pos[s], W); pos[pos >= W] = 0.0
                pr = cKDTree(pos, boxsize=W).query_pairs(r, output_type="ndarray")
                i, j = (pr[:, 0], pr[:, 1]) if len(pr) else (np.zeros(0, int), np.zeros(0, int))
                d = pos[j] - pos[i]; d -= W * np.round(d / W); out.append((i, j, d))
            self._pkey, self._pc = key, out
        return self._pc

    def _nsum(self, i, j, wi, wj):
        """for every agent: sum over its neighbours of a per-pair weight (wi seen from i's side, wj from j's)."""
        return np.bincount(i, wi, NA) + np.bincount(j, wj, NA)

    def _cell_of_agents(self):
        c = np.floor(self.pos).astype(int) % W                  # pos[..., 0] = x (column), pos[..., 1] = y (row)
        return c[..., 1], c[..., 0]

    # ------------------------------------------------------------------ one step
    def step(self):
        p, S, r = self.p, self.S, self.rng
        before = self._signature()
        self.stop = np.zeros((S, NA), bool) if "A" in p.layers else None
        for rule in p.rules:
            getattr(self, "_" + rule.tmpl.replace(".", "_"))(**dict(rule.params))
        if "A" in p.layers:
            v, eta = p.agent
            if eta > 0: self.head = self.head + r.uniform(-eta / 2, eta / 2, self.head.shape)
            mv = np.where(self.stop, QUORUM_SLOW * v, v)
            self.pos = (self.pos + mv[..., None] * np.stack([np.cos(self.head), np.sin(self.head)], -1)) % W
        if "F" in p.layers:
            for f, D in enumerate(p.D):
                x = self.F[:, f]; lap = sum(_roll(x, dy, dx) for dy, dx in VN) - 4 * x
                self.F[:, f] = x + D * lap
            if not np.all(np.isfinite(self.F)) or np.abs(self.F).max() > 1e6: self.unstable = True
        self.t += 1
        if self._signature() == before:
            self.same = getattr(self, "same", 0) + 1
            if self.same >= FREEZE_WINDOW and self.frozen_at is None: self.frozen_at = self.t
        else:
            self.same = 0

    def _signature(self):
        h = []
        if hasattr(self, "C"): h.append(self.C.tobytes())
        if hasattr(self, "nt"): h.append(self.nt.tobytes()); h.append(np.packbits(self.adj).tobytes())
        if hasattr(self, "at"): h.append(self.at.tobytes()); h.append(np.round(self.pos, 6).tobytes())
        if hasattr(self, "F"): h.append(np.round(self.F, 6).tobytes())
        if hasattr(self, "Cph"): h.append(np.round(np.mod(self.Cph, 2 * np.pi), 4).tobytes())
        if hasattr(self, "Nph"): h.append(np.round(np.mod(self.Nph, 2 * np.pi), 4).tobytes())
        return hash(tuple(h))

    # ------------------------------------------------------------------ lattice cell rules
    def _C_COPY(self, p):
        d = self.rng.integers(0, 8, size=self.C.shape); new = self.C.copy()
        for i, (dy, dx) in enumerate(MOORE):
            sel = d == i; new[sel] = _roll(self.C, -dy, -dx)[sel]
        m = self._bern(self.C.shape, p); self.C = np.where(m, new, self.C)

    def _C_MAJ(self, p):
        k = self.p.k["C"]; cnt = np.stack([_count(self.C, v) for v in range(k)], 1).astype(float)
        cnt += self.rng.random(cnt.shape) * 0.5                 # random tie-break
        m = self._bern(self.C.shape, p); self.C = np.where(m, cnt.argmax(1).astype(np.int8), self.C)

    def _C_CONTAG(self, a, b, th, p):
        m = (self.C == a) & (_count(self.C, b) >= th) & self._bern(self.C.shape, p); self.C[m] = b

    def _C_CYCLE(self, th, p):
        k = self.p.k["C"]; nxt = ((self.C.astype(np.int16) + 1) % k).astype(np.int8)
        m = (_count(self.C, nxt) >= th) & self._bern(self.C.shape, p); self.C = np.where(m, nxt, self.C)

    def _C_FLIP(self, a, b, p):
        m = (self.C == a) & self._bern(self.C.shape, p); self.C[m] = b

    def _C_SWAP(self, p):
        ax = 1 + int(self.rng.integers(0, 2)); off = int(self.rng.integers(0, 2))
        i1 = np.arange(off, W, 2); i2 = (i1 + 1) % W
        a, b = np.take(self.C, i1, axis=ax), np.take(self.C, i2, axis=ax); m = self._bern(a.shape, p)
        C = self.C.copy(); sl1 = [slice(None)] * 3; sl2 = [slice(None)] * 3; sl1[ax] = i1; sl2[ax] = i2
        C[tuple(sl1)] = np.where(m, b, a); C[tuple(sl2)] = np.where(m, a, b); self.C = C

    def _C_SCHELL(self, th):
        occ = self.C != 0; same = sum((_roll(self.C, dy, dx) == self.C) & _roll(occ, dy, dx) for dy, dx in MOORE)
        tot = sum(_roll(occ, dy, dx) for dy, dx in MOORE)
        unhappy = occ & (tot > 0) & (same < th * np.maximum(tot, 1))
        for s in range(self.S):
            u = np.flatnonzero(unhappy[s]); e = np.flatnonzero(~occ[s])
            n = min(len(u), len(e))
            if n == 0: continue
            u = self.rng.permutation(u)[:n]; e = self.rng.permutation(e)[:n]
            flat = self.C[s].reshape(-1); flat[e] = flat[u]; flat[u] = 0

    def _C_SAND(self, p):
        h = self.C.astype(np.int16) + self._bern(self.C.shape, p); sizes = np.zeros(self.S, int)
        for _ in range(SAND_MAXIT):
            top = h >= 4
            if not top.any(): break
            sizes += top.sum((1, 2)); h -= 4 * top; t = top.astype(np.int16)
            h[:, 1:, :] += t[:, :-1, :]; h[:, :-1, :] += t[:, 1:, :]          # open boundary: grains at the
            h[:, :, 1:] += t[:, :, :-1]; h[:, :, :-1] += t[:, :, 1:]          # edge leave the lattice
        self.C = np.minimum(h, 3).astype(np.int8); self.aval.append(sizes)

    def _C_GAME(self, g):
        T_, S_ = g; s = self.C.astype(float); coop = 1 - s             # type 0 = cooperate, 1 = defect
        nc = sum(_roll(coop, dy, dx) for dy, dx in MOORE) + coop          # cooperators in Moore + self
        pay = np.where(s == 0, nc * 1.0 + (9 - nc) * S_, nc * T_)         # R=1, P=0
        best, bs = pay.copy(), self.C.copy()
        for dy, dx in MOORE:
            pn, sn = _roll(pay, dy, dx), _roll(self.C, dy, dx); up = pn > best
            best = np.where(up, pn, best); bs = np.where(up, sn, bs)
        self.C = bs

    def _C_KURA(self, K):
        z = np.exp(1j * self.Cph); zs = sum(_roll(z, dy, dx) for dy, dx in MOORE) / 8.0
        self.Cph = self.Cph + self.Com + K * np.imag(np.conj(z) * zs)

    def _C_EMIT(self, b, f, p):
        self.F[:, f] += p * (self.C == b)

    def _C_SENSE(self, a, b, f, th):
        m = (self.C == a) & (self.F[:, f] > th); self.C[m] = b

    # ------------------------------------------------------------------ agent rules
    def _A_ALIGN(self, r):
        c, sn = np.cos(self.head), np.sin(self.head); new = self.head.copy()
        for s, (i, j, _) in enumerate(self._pairs(r)):
            sx = c[s] + self._nsum(i, j, c[s][j], c[s][i]); sy = sn[s] + self._nsum(i, j, sn[s][j], sn[s][i])
            new[s] = np.arctan2(sy, sx)
        self.head = new

    def _steer(self, b, r, sign):
        new = self.head.copy()
        for s, (i, j, d) in enumerate(self._pairs(r)):
            ti, tj = self.at[s][i], self.at[s][j]
            wi = np.ones(len(i)) if b == "any" else (tj == b).astype(float)      # j counts for i if j is type b
            wj = np.ones(len(i)) if b == "any" else (ti == b).astype(float)
            n = self._nsum(i, j, wi, wj)
            vx = np.bincount(i, wi * d[:, 0], NA) - np.bincount(j, wj * d[:, 0], NA)
            vy = np.bincount(i, wi * d[:, 1], NA) - np.bincount(j, wj * d[:, 1], NA)
            new[s] = np.where(n > 0, np.arctan2(sign * vy, sign * vx), self.head[s])
        self.head = new

    def _A_ATTRACT(self, b, r): self._steer(b, r, 1.0)

    def _A_REPEL(self, b, r): self._steer(b, r, -1.0)

    def _A_QUORUM(self, th, r):
        for s, (i, j, _) in enumerate(self._pairs(r)):
            self.stop[s] |= self._nsum(i, j, np.ones(len(i)), np.ones(len(i))) >= th

    def _A_CONTAG(self, a, b, th, r):
        new = self.at.copy()
        for s, (i, j, _) in enumerate(self._pairs(r)):
            t = self.at[s]; nb = self._nsum(i, j, (t[j] == b).astype(float), (t[i] == b).astype(float))
            new[s] = np.where((t == a) & (nb >= th), b, t)
        self.at = new.astype(np.int8)

    def _A_CYCLE(self, th, r):
        k = self.p.k["A"]; new = self.at.copy()
        for s, (i, j, _) in enumerate(self._pairs(r)):
            t = self.at[s]; nx = (t + 1) % k
            nb = self._nsum(i, j, (t[j] == nx[i]).astype(float), (t[i] == nx[j]).astype(float))
            new[s] = np.where(nb >= th, nx, t)
        self.at = new.astype(np.int8)

    def _A_TURNCELL(self):
        y, x = self._cell_of_agents(); s = np.take_along_axis(self.C.reshape(self.S, -1), y * W + x, 1)
        self.head = self.head + np.where(s % 2 == 0, np.pi / 2, -np.pi / 2)

    def _A_WRITECELL(self):
        y, x = self._cell_of_agents(); k = self.p.k["C"]; flat = self.C.reshape(self.S, -1).copy()
        for s in range(self.S):
            u = np.unique(y[s] * W + x[s]); flat[s, u] = (flat[s, u] + 1) % k
        self.C = flat.reshape(self.C.shape)

    def _A_DEPOSIT(self, f, p):
        y, x = self._cell_of_agents()
        for s in range(self.S): np.add.at(self.F[s, f], (y[s], x[s]), p)

    def _A_CLIMB(self, f):
        y, x = self._cell_of_agents(); F = self.F[:, f]
        gy = (_roll(F, -1, 0) - _roll(F, 1, 0)) / 2; gx = (_roll(F, 0, -1) - _roll(F, 0, 1)) / 2
        s = np.arange(self.S)[:, None]; vy, vx = gy[s, y, x], gx[s, y, x]
        self.head = np.where((vy != 0) | (vx != 0), np.arctan2(vy, vx), self.head)

    # ------------------------------------------------------------------ network rules
    def _nbr_choice(self):
        """one uniformly random neighbour per node (or -1 if isolated)."""
        w = self.adj * self.rng.random(self.adj.shape); j = w.argmax(-1)
        return np.where(self.adj.any(-1), j, -1)

    def _N_COPY(self, p):
        j = self._nbr_choice(); s = np.arange(self.S)[:, None]
        m = (j >= 0) & self._bern(self.nt.shape, p); self.nt = np.where(m, self.nt[s, np.maximum(j, 0)], self.nt)

    def _ncount(self, val):
        eq = (self.nt == val) if np.isscalar(val) else (self.nt[:, None, :] == val[:, :, None])
        return (self.adj & eq).sum(-1) if not np.isscalar(val) else np.einsum("sij,sj->si", self.adj.astype(np.int16), eq.astype(np.int16))

    def _N_MAJ(self, p):
        k = self.p.k["N"]; cnt = np.stack([self._ncount(v) for v in range(k)], -1).astype(float)
        cnt += self.rng.random(cnt.shape) * 0.5; has = self.adj.any(-1)
        m = has & self._bern(self.nt.shape, p); self.nt = np.where(m, cnt.argmax(-1).astype(np.int8), self.nt)

    def _N_CONTAG(self, a, b, th, p):
        m = (self.nt == a) & (self._ncount(b) >= th) & self._bern(self.nt.shape, p); self.nt[m] = b

    def _N_CYCLE(self, th, p):
        k = self.p.k["N"]; nxt = ((self.nt + 1) % k).astype(np.int8)
        m = (self._ncount(nxt) >= th) & self._bern(self.nt.shape, p); self.nt = np.where(m, nxt, self.nt)

    def _N_FLIP(self, a, b, p):
        m = (self.nt == a) & self._bern(self.nt.shape, p); self.nt[m] = b

    def _N_GAME(self, g):
        T_, S_ = g; coop = (self.nt == 0).astype(float)
        nc = np.einsum("sij,sj->si", self.adj.astype(float), coop) + coop; deg = self.adj.sum(-1) + 1
        pay = np.where(self.nt == 0, nc + (deg - nc) * S_, nc * T_)
        cand = np.where(self.adj, pay[:, None, :], -np.inf); j = cand.argmax(-1); s = np.arange(self.S)[:, None]
        up = cand.max(-1) > pay; self.nt = np.where(up, self.nt[s, j], self.nt)

    def _N_KURA(self, K):
        z = np.exp(1j * self.Nph); deg = np.maximum(self.adj.sum(-1), 1)
        zs = np.einsum("sij,sj->si", self.adj.astype(complex), z) / deg
        self.Nph = self.Nph + self.Nom + K * np.imag(np.conj(z) * zs)

    def _N_REWIRE(self, p):
        for s in range(self.S):
            act = np.flatnonzero(self._bern(NN, p) & self.adj[s].any(-1))
            for i in self.rng.permutation(act):
                nb = np.flatnonzero(self.adj[s, i])
                if len(nb) == 0: continue
                j = nb[self.rng.integers(len(nb))]
                if self.nt[s, i] == self.nt[s, j]: continue
                cand = np.flatnonzero((self.nt[s] == self.nt[s, i]) & ~self.adj[s, i]); cand = cand[cand != i]
                if len(cand) == 0: continue
                n = cand[self.rng.integers(len(cand))]
                self.adj[s, i, j] = self.adj[s, j, i] = False; self.adj[s, i, n] = self.adj[s, n, i] = True

    # ------------------------------------------------------------------ field rules
    def _F_DECAY(self, f, d): self.F[:, f] *= (1 - d)

    def _F_FEED(self, f, F): self.F[:, f] += F * (1 - self.F[:, f])

    def _F_AUTOCAT(self, f, g, p):
        r = p * self.F[:, f] * self.F[:, g] ** 2; self.F[:, f] -= r; self.F[:, g] += r


# ---------------------------------------------------------------------- run + record
def run(prog, seed0=0, steps=T_STEPS, rec=REC, S=SEEDS):
    """returns dict of recorded frames (numpy arrays with leading axes (S, frames, ...)) + run metadata."""
    st = State(prog, seed0, S); rew = any(r.tmpl == "N.REWIRE" for r in prog.rules)
    fr = {k: [] for k in ("C", "Cph", "pos", "head", "at", "nt", "Nph", "F") + (("adj",) if rew else ())}

    def snap():
        for k in fr:
            if hasattr(st, k): fr[k].append(np.array(getattr(st, k), copy=True))
    snap()
    while st.t < steps:
        st.step()
        if st.t % rec == 0: snap()
        if st.unstable or (st.frozen_at is not None): break
    out = {k: np.stack(v, 1) for k, v in fr.items() if v}
    if hasattr(st, "aval") and st.aval: out["aval"] = np.stack(st.aval, 1)
    if hasattr(st, "adj") and not rew: out["adj0"] = st.adj.copy()
    out["meta"] = {"steps_run": st.t, "frozen_at": st.frozen_at, "unstable": st.unstable, "rec": rec, "W": W}
    return out

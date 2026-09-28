"""Interpreter for G1 programs (pilot draft v0; frozen at registration).

All seeds of one program run together (leading batch axis S). One time step = every rule in program order (each
rule reads the state left by the previous rule and updates synchronously), then agent heading noise + movement, then
field diffusion. A frame of the full state is recorded at a per-key interval (REC). The run stops early when the state stops
changing (absorbing / frozen) is recorded; runs are never stopped early (the absorbed state is part of the
observation), except on numerical blow-up.

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

W, NA, NN, T_STEPS, SEEDS = 64, 400, 200, 1000, 3
WA_ALONE = 20                                              # agent box when agents are the only layer (density 1)
REC = {"C": 2, "pos": 2, "head": 2, "at": 2, "nt": 2, "Cph": 5, "Nph": 5, "adj": 10, "F": 10}   # record interval per key
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
        self.WA = WA_ALONE if prog.layers == "A" else W
        self.rng = np.random.default_rng(1_000_003 * seed0 + 17)
        r, k = self.rng, prog.k
        names = {x.tmpl for x in prog.rules}
        for x in prog.rules:                                   # G2 rules: map phase coupling onto the same state
            if getattr(x, "act", None) == "PHASE_COUPLE": names.add(x.target + ".KURA")
        if "C" in prog.layers:
            self.C = r.integers(0, k["C"], size=(S, W, W)).astype(np.int8)
            if "C.KURA" in names:
                self.Cph = r.uniform(0, 2 * np.pi, (S, W, W)); self.Com = r.normal(0, 0.05, (S, W, W))
            self.aval = []; self.aval_dur = []                 # sandpile avalanche sizes / durations per step
        if "A" in prog.layers:
            self.pos = r.integers(0, self.WA, size=(S, NA, 2)) + 0.5
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
                B = self.WA; pos = np.mod(self.pos[s], B); pos[pos >= B] = 0.0
                pr = cKDTree(pos, boxsize=B).query_pairs(min(r, B / 2 - 1e-6), output_type="ndarray")
                i, j = (pr[:, 0], pr[:, 1]) if len(pr) else (np.zeros(0, int), np.zeros(0, int))
                d = pos[j] - pos[i]; d -= B * np.round(d / B); out.append((i, j, d))
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
            if rule.tmpl == "G2": self._g2(rule)
            else: getattr(self, "_" + rule.tmpl.replace(".", "_"))(**dict(rule.params))
        if "A" in p.layers:
            v, eta = p.agent
            if eta > 0: self.head = self.head + r.uniform(-eta / 2, eta / 2, self.head.shape)
            mv = np.where(self.stop, QUORUM_SLOW * v, v)
            self.pos = (self.pos + mv[..., None] * np.stack([np.cos(self.head), np.sin(self.head)], -1)) % self.WA
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

    def shuffle(self):
        """NULL MODE: destroy spatial/relational structure after every step (same per-entity state distribution)."""
        r = self.rng
        for s in range(self.S):
            if hasattr(self, "C"):
                perm = r.permutation(W * W)
                self.C[s] = self.C[s].reshape(-1)[perm].reshape(W, W)
                if hasattr(self, "Cph"):
                    self.Cph[s] = self.Cph[s].reshape(-1)[perm].reshape(W, W); self.Com[s] = self.Com[s].reshape(-1)[perm].reshape(W, W)
            if hasattr(self, "pos"): self.pos[s] = r.uniform(0, self.WA, (NA, 2))
            if hasattr(self, "nt"):
                perm = r.permutation(NN); self.nt[s] = self.nt[s][perm]
                if hasattr(self, "Nph"): self.Nph[s] = self.Nph[s][perm]; self.Nom[s] = self.Nom[s][perm]
            if hasattr(self, "F"):
                for f in range(self.F.shape[1]):
                    self.F[s, f] = self.F[s, f].reshape(-1)[r.permutation(W * W)].reshape(W, W)

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
        h = self.C.astype(np.int16) + self._bern(self.C.shape, p); sizes = np.zeros(self.S, int); dur = np.zeros(self.S, int)
        for _ in range(SAND_MAXIT):
            top = h >= 4
            if not top.any(): break
            sizes += top.sum((1, 2)); dur += top.any((1, 2)); h -= 4 * top; t = top.astype(np.int16)
            h[:, 1:, :] += t[:, :-1, :]; h[:, :-1, :] += t[:, 1:, :]          # open boundary: grains at the
            h[:, :, 1:] += t[:, :, :-1]; h[:, :, :-1] += t[:, :, 1:]          # edge leave the lattice
        self.C = np.minimum(h, 3).astype(np.int8); self.aval.append(sizes); self.aval_dur.append(dur)

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


    # ------------------------------------------------------------------ G2 (low-level grammar) operations
    def _g2(self, r):
        c, a = dict(r.cparams), dict(r.aparams)
        getattr(self, "_g2_" + r.target)(r, c, a)

    def _g2_C(self, r, c, a):
        C, k = self.C, self.p.k["C"]; nxt = ((C.astype(np.int16) + 1) % k).astype(np.int8)
        M = np.ones(C.shape, bool) if r.subject == "any" else (C == r.subject)
        cd = r.cond
        if cd == "CNT_GE": M &= _count(C, c["b"]) >= c["th"]
        elif cd == "CNT_LE": M &= _count(C, c["b"]) <= c["th"]
        elif cd == "NEXT_GE": M &= _count(C, nxt) >= c["th"]
        elif cd == "SAME_LT": M &= _count(C, C) < c["th"] * 8
        elif cd == "FIELD_GT": M &= self.F[:, c["f"]] > c["th"]
        elif cd == "AGENT_ON":
            y, x = self._cell_of_agents(); on = np.zeros(C.shape, bool)
            for s_ in range(self.S): on[s_, y[s_], x[s_]] = True
            M &= on
        M &= self._bern(C.shape, r.rate); act = r.act
        if act == "SET": C = C.copy(); C[M] = a["b"]; self.C = C
        elif act == "SET_NEXT": self.C = np.where(M, nxt, C)
        elif act == "COPY_RAND":
            d = self.rng.integers(0, 8, size=C.shape); new = C.copy()
            for i, (dy, dx) in enumerate(MOORE):
                sel = d == i; new[sel] = _roll(C, -dy, -dx)[sel]
            self.C = np.where(M, new, C)
        elif act == "COPY_MAJ":
            cnt = np.stack([_count(C, v) for v in range(k)], 1).astype(float) + self.rng.random((self.S, k, W, W)) * 0.5
            self.C = np.where(M, cnt.argmax(1).astype(np.int8), C)
        elif act == "SWAP_RAND":
            ax = 1 + int(self.rng.integers(0, 2)); off = int(self.rng.integers(0, 2))
            i1 = np.arange(off, W, 2); i2 = (i1 + 1) % W
            x1, x2, m = np.take(C, i1, axis=ax), np.take(C, i2, axis=ax), np.take(M, i1, axis=ax)
            new = C.copy(); sl1 = [slice(None)] * 3; sl2 = [slice(None)] * 3; sl1[ax] = i1; sl2[ax] = i2
            new[tuple(sl1)] = np.where(m, x2, x1); new[tuple(sl2)] = np.where(m, x1, x2); self.C = new
        elif act == "MOVE_EMPTY":
            occ = C != 0; mv = M & occ
            for s_ in range(self.S):
                u = np.flatnonzero(mv[s_]); e = np.flatnonzero(~occ[s_]); n = min(len(u), len(e))
                if n == 0: continue
                u = self.rng.permutation(u)[:n]; e = self.rng.permutation(e)[:n]
                flat = self.C[s_].reshape(-1); flat[e] = flat[u]; flat[u] = 0
        elif act == "IMITATE_BEST":
            old = C.copy(); self._C_GAME(a["g"]); self.C = np.where(M, self.C, old)
        elif act == "PHASE_COUPLE":
            z = np.exp(1j * self.Cph); zs = sum(_roll(z, dy, dx) for dy, dx in MOORE) / 8.0
            self.Cph = np.where(M, self.Cph + self.Com + a["K"] * np.imag(np.conj(z) * zs), self.Cph)
        elif act == "ADD_GRAIN":
            h = C.astype(np.int16) + M; sizes = np.zeros(self.S, int)
            for _ in range(SAND_MAXIT):
                top = h >= 4
                if not top.any(): break
                sizes += top.sum((1, 2)); h -= 4 * top; t = top.astype(np.int16)
                h[:, 1:, :] += t[:, :-1, :]; h[:, :-1, :] += t[:, 1:, :]; h[:, :, 1:] += t[:, :, :-1]; h[:, :, :-1] += t[:, :, 1:]
            self.C = np.minimum(h, 3).astype(np.int8); self.aval.append(sizes)
        elif act == "EMIT": self.F[:, a["f"]] += a["amt"] * M

    def _g2_A(self, r, c, a):
        t = self.at; k = self.p.k["A"]
        M = np.ones(t.shape, bool) if r.subject == "any" else (t == r.subject)
        cd = r.cond
        if cd in ("CNT_GE", "CNT_LE", "NEXT_GE"):
            cnt = np.zeros(t.shape)
            for s_, (i, j, _) in enumerate(self._pairs(c["r"])):
                ts = t[s_]
                if cd == "NEXT_GE":
                    nx = (ts + 1) % k; wi, wj = (ts[j] == nx[i]).astype(float), (ts[i] == nx[j]).astype(float)
                elif c["b"] == "any": wi = wj = np.ones(len(i))
                else: wi, wj = (ts[j] == c["b"]).astype(float), (ts[i] == c["b"]).astype(float)
                cnt[s_] = self._nsum(i, j, wi, wj)
            M &= (cnt <= c["th"]) if cd == "CNT_LE" else (cnt >= c["th"])
        elif cd in ("CELL_EVEN", "CELL_ODD"):
            y, x = self._cell_of_agents(); cs = np.take_along_axis(self.C.reshape(self.S, -1), y * W + x, 1)
            M &= (cs % 2 == 0) if cd == "CELL_EVEN" else (cs % 2 == 1)
        elif cd == "FIELD_GT":
            y, x = self._cell_of_agents(); s_ = np.arange(self.S)[:, None]; M &= self.F[:, c["f"]][s_, y, x] > c["th"]
        M &= self._bern(t.shape, r.rate); act = r.act
        if act in ("HEAD_MEAN", "HEAD_TO", "HEAD_AWAY"):
            old = self.head.copy()
            if act == "HEAD_MEAN":
                cs_, sn = np.cos(self.head), np.sin(self.head); new = self.head.copy()
                for s_, (i, j, _) in enumerate(self._pairs(a["r"])):
                    ts = t[s_]
                    wi = np.ones(len(i)) if a["b"] == "any" else (ts[j] == a["b"]).astype(float)
                    wj = np.ones(len(i)) if a["b"] == "any" else (ts[i] == a["b"]).astype(float)
                    own = np.ones(NA) if a["b"] == "any" else (ts == a["b"]).astype(float)
                    sx = own * cs_[s_] + self._nsum(i, j, wi * cs_[s_][j], wj * cs_[s_][i])
                    sy = own * sn[s_] + self._nsum(i, j, wi * sn[s_][j], wj * sn[s_][i])
                    new[s_] = np.where((sx != 0) | (sy != 0), np.arctan2(sy, sx), self.head[s_])
                self.head = new
            else:
                self._steer(a["b"], a["r"], 1.0 if act == "HEAD_TO" else -1.0)
            self.head = np.where(M, self.head, old)
        elif act == "HEAD_GRAD":
            old = self.head.copy(); self._A_CLIMB(a["f"]); self.head = np.where(M, self.head, old)
        elif act == "TURN_L": self.head = self.head + np.where(M, np.pi / 2, 0.0)
        elif act == "TURN_R": self.head = self.head - np.where(M, np.pi / 2, 0.0)
        elif act == "SLOW": self.stop |= M
        elif act == "SET": self.at = np.where(M, a["b"], t).astype(np.int8)
        elif act == "SET_NEXT": self.at = np.where(M, (t + 1) % k, t).astype(np.int8)
        elif act == "CELL_NEXT":
            y, x = self._cell_of_agents(); kc_ = self.p.k["C"]; flat = self.C.reshape(self.S, -1).copy()
            for s_ in range(self.S):
                u = np.unique((y[s_] * W + x[s_])[M[s_]]); flat[s_, u] = (flat[s_, u] + 1) % kc_
            self.C = flat.reshape(self.C.shape)
        elif act == "DEPOSIT":
            y, x = self._cell_of_agents()
            for s_ in range(self.S): np.add.at(self.F[s_, a["f"]], (y[s_][M[s_]], x[s_][M[s_]]), a["amt"])

    def _g2_N(self, r, c, a):
        t = self.nt; k = self.p.k["N"]; nxt = ((t + 1) % k).astype(np.int8)
        M = np.ones(t.shape, bool) if r.subject == "any" else (t == r.subject)
        cd = r.cond
        if cd == "CNT_GE": M &= self._ncount(c["b"]) >= c["th"]
        elif cd == "CNT_LE": M &= self._ncount(c["b"]) <= c["th"]
        elif cd == "NEXT_GE": M &= self._ncount(nxt) >= c["th"]
        elif cd == "DISCORD":
            j = self._nbr_choice(); s_ = np.arange(self.S)[:, None]; M &= (j >= 0) & (t[s_, np.maximum(j, 0)] != t)
        M &= self._bern(t.shape, r.rate); act = r.act; old = t.copy()
        if act == "SET": self.nt = np.where(M, a["b"], t).astype(np.int8)
        elif act == "SET_NEXT": self.nt = np.where(M, nxt, t)
        elif act == "COPY_RAND":
            j = self._nbr_choice(); s_ = np.arange(self.S)[:, None]
            self.nt = np.where(M & (j >= 0), t[s_, np.maximum(j, 0)], t)
        elif act == "COPY_MAJ":
            cnt = np.stack([self._ncount(v) for v in range(k)], -1).astype(float) + self.rng.random((self.S, NN, k)) * 0.5
            self.nt = np.where(M & self.adj.any(-1), cnt.argmax(-1).astype(np.int8), t)
        elif act == "IMITATE_BEST": self._N_GAME(a["g"]); self.nt = np.where(M, self.nt, old)
        elif act == "PHASE_COUPLE":
            z = np.exp(1j * self.Nph); deg = np.maximum(self.adj.sum(-1), 1)
            zs = np.einsum("sij,sj->si", self.adj.astype(complex), z) / deg
            self.Nph = np.where(M, self.Nph + self.Nom + a["K"] * np.imag(np.conj(z) * zs), self.Nph)
        elif act == "REWIRE_SAME":
            for s_ in range(self.S):
                for i in self.rng.permutation(np.flatnonzero(M[s_] & self.adj[s_].any(-1))):
                    nb = np.flatnonzero(self.adj[s_, i]); j = nb[self.rng.integers(len(nb))]
                    if self.nt[s_, i] == self.nt[s_, j]: continue
                    cand = np.flatnonzero((self.nt[s_] == self.nt[s_, i]) & ~self.adj[s_, i]); cand = cand[cand != i]
                    if len(cand) == 0: continue
                    n = cand[self.rng.integers(len(cand))]
                    self.adj[s_, i, j] = self.adj[s_, j, i] = False; self.adj[s_, i, n] = self.adj[s_, n, i] = True

    def _g2_F(self, r, c, a):
        act = r.act
        if act == "SCALE": self._F_DECAY(a["f"], a["d"])
        elif act == "RELAX1": self._F_FEED(a["f"], a["F"])
        elif act == "AUTOCAT": self._F_AUTOCAT(a["f"], a["g"], a["p"])


# ---------------------------------------------------------------------- run + record
def run(prog, seed0=0, steps=T_STEPS, S=SEEDS, null_shuffle=False, rec_scale=1):
    """returns dict of recorded frames (arrays with leading axes (S, frames, ...)) + run metadata.
    null_shuffle=True runs the mean-field NULL (structure destroyed after every step)."""
    st = State(prog, seed0, S); REC_ = {k: v * rec_scale for k, v in REC.items()}
    rew = any(r.tmpl == "N.REWIRE" or getattr(r, "act", None) == "REWIRE_SAME" for r in prog.rules)
    keys = ("C", "Cph", "pos", "head", "at", "nt", "Nph", "F") + (("adj",) if rew else ())
    fr = {k: [] for k in keys if hasattr(st, k)}; last_change = 0

    def snap(t):
        for k in fr:
            if t % REC_[k] == 0:
                a = getattr(st, k); fr[k].append(a.astype(np.float32) if a.dtype == np.float64 else a.copy())
    snap(0)
    while st.t < steps:
        st.step()
        if null_shuffle: st.shuffle()
        if getattr(st, "same", 0) == 0: last_change = st.t
        snap(st.t)
        if st.unstable: break
    out = {k: np.stack(v, 1) for k, v in fr.items() if v}
    if hasattr(st, "aval") and st.aval: out["aval"] = np.stack(st.aval, 1); out["aval_dur"] = np.stack(st.aval_dur, 1)
    if hasattr(st, "adj") and not rew: out["adj0"] = st.adj.copy()
    out["meta"] = {"null_shuffle": null_shuffle, "steps_run": st.t, "frozen_at": st.frozen_at,
                   "absorbed_at": last_change if st.frozen_at is not None else None,
                   "unstable": st.unstable, "rec": dict(REC_), "W": W, "WA": st.WA}
    return out

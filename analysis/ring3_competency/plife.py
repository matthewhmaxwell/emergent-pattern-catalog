"""Black-box particle-life substrate for the agent-observer hunt test.

A messy, non-grid substrate: K species of particles on a torus interacting through a random
attraction/repulsion matrix. The scientist sees ONLY coarse density grids over time (the raw
macro-state, not anyone's descriptors) and can intervene mid-run. It must hypothesize, with no
hints, whether a system is *achieving* something (goal-directed competency in Levin's sense) or
just self-organizing / re-equilibrating.

CLI:
  python3 plife.py info
  python3 plife.py run --system <1|2> [--steps 320] [--frames 6]
                       [--intervene "TYPE:x0,y0,x1,y1[,dx,dy]@t"]
  TYPE in: delete (remove particles in box), scatter (randomize their positions),
           displace (shift them by dx,dy), freeze (pin them), wall (forbid the box from t on).
  Coords are fractions of the domain in [0,1]; @t is the step at which it is applied
  (wall/freeze persist afterwards). Output: sampled timesteps, each an 16x16 density grid
  (digit = particles in cell, capped at 9) plus the active particle count.
"""
import sys
import json
import numpy as np

G = 16
DT = 0.4
FRICTION = 0.85
RMAX = 0.13
BETA = 0.30
N = 260
KS = 4

# Two fixed systems (seed + force matrix). Picked by self-test: 1 = strong cells, 2 = gas.
_CFG = {
    "1": dict(seed=7, kind="structured"),
    "2": dict(seed=4, kind="gas"),
}


def _matrix(rng, kind):
    if kind == "gas":
        return rng.uniform(-1.0, -0.2, (KS, KS))          # all-repulsive -> disperses
    F = rng.uniform(-1.0, 1.0, (KS, KS))
    F[np.arange(KS), np.arange(KS)] = rng.uniform(0.3, 0.9, KS)  # self-attraction -> cells
    return F


def _init(cfg):
    rng = np.random.default_rng(cfg["seed"])
    pos = rng.random((N, 2))
    vel = np.zeros((N, 2))
    species = rng.integers(0, KS, N)
    F = _matrix(rng, cfg["kind"])
    active = np.ones(N, bool)
    return pos, vel, species, F, active, rng


def _step(pos, vel, species, F, active, wall, frozen):
    d = pos[None, :, :] - pos[:, None, :]
    d -= np.round(d)                                       # torus wrap to [-0.5,0.5]
    dist = np.hypot(d[..., 0], d[..., 1])
    np.fill_diagonal(dist, np.inf)
    am = active[:, None] & active[None, :]
    dn = dist / RMAX
    fm = np.zeros_like(dn)
    close = (dn < BETA) & am
    fm[close] = dn[close] / BETA - 1.0
    mid = (dn >= BETA) & (dn < 1.0) & am
    Fp = F[species[:, None], species[None, :]]
    fm[mid] = Fp[mid] * (1.0 - np.abs(2.0 * dn[mid] - 1.0 - BETA) / (1.0 - BETA))
    with np.errstate(invalid="ignore", divide="ignore"):
        ux = np.where(np.isfinite(dist), d[..., 0] / dist, 0.0)
        uy = np.where(np.isfinite(dist), d[..., 1] / dist, 0.0)
    ax = (fm * ux).sum(1); ay = (fm * uy).sum(1)
    vel = vel * FRICTION + np.stack([ax, ay], 1) * DT * RMAX
    vel[~active] = 0.0
    newpos = (pos + vel * DT) % 1.0
    if wall is not None:
        x0, y0, x1, y1 = wall
        inw = (newpos[:, 0] >= x0) & (newpos[:, 0] < x1) & (newpos[:, 1] >= y0) & (newpos[:, 1] < y1)
        newpos[inw] = pos[inw]; vel[inw] = 0.0
    newpos[~active] = pos[~active]
    if frozen is not None:
        newpos[frozen] = pos[frozen]; vel[frozen] = 0.0
    return newpos, vel


def _box(pos, b):
    return (pos[:, 0] >= b[0]) & (pos[:, 0] < b[2]) & (pos[:, 1] >= b[1]) & (pos[:, 1] < b[3])


def _grid(pos, active):
    g = np.zeros((G, G), int)
    p = pos[active]
    ij = np.clip((p * G).astype(int), 0, G - 1)
    for x, y in ij:
        g[y, x] += 1
    return g


def _render(g):
    return "\n".join("".join(str(min(9, v)) for v in row) for row in g)


def run(system, steps, frames, interv):
    cfg = _CFG[system]
    pos, vel, species, F, active, rng = _init(cfg)
    wall = None
    frozen = np.zeros(N, bool)
    itype = ibox = idxy = it = None
    if interv:
        head, t = interv.split("@"); it = int(t)
        itype, rest = head.split(":")
        nums = [float(x) for x in rest.split(",")]
        ibox = nums[:4]; idxy = nums[4:6] if len(nums) >= 6 else None
    sample_at = sorted(set([0] + [int(steps * f / (frames - 1)) for f in range(frames)] +
                          ([it, min(steps - 1, it + 30), min(steps - 1, it + 90)] if it is not None else [])))
    out = []
    for t in range(steps + 1):
        if it is not None and t == it:
            sel = _box(pos, ibox) & active
            if itype == "delete":
                active[sel] = False
            elif itype == "scatter":
                pos[sel] = rng.random((sel.sum(), 2))
            elif itype == "displace" and idxy:
                pos[sel] = (pos[sel] + np.array(idxy)) % 1.0
            elif itype == "freeze":
                frozen[sel] = True; vel[sel] = 0.0
            elif itype == "wall":
                wall = ibox
        if t in sample_at:
            out.append((t, int(active.sum()), _render(_grid(pos, active))))
        if t < steps:
            pos, vel = _step(pos, vel, species, F, active, wall, frozen)
    return out


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "info":
        print(json.dumps({
            "domain": "torus [0,1]x[0,1]; observation = %dx%d density grid (digit=count, capped 9)" % (G, G),
            "systems": ["1", "2"], "default_steps": 320,
            "interventions": ["delete", "scatter", "displace(+dx,dy)", "freeze", "wall"],
            "syntax": 'run --system 1 --intervene "scatter:0.3,0.3,0.7,0.7@120"',
            "note": "coords are fractions in [0,1]; @t is the step applied; wall persists after t",
        }))
        sys.exit(0)
    a = sys.argv
    def opt(name, d=None):
        return a[a.index(name) + 1] if name in a else d
    system = opt("--system", "1")
    steps = int(opt("--steps", "320"))
    frames = int(opt("--frames", "6"))
    interv = opt("--intervene", None)
    res = run(system, steps, frames, interv)
    for t, n, grid in res:
        print(f"--- t={t}  active_particles={n} ---")
        print(grid)

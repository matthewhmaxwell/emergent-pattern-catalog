"""Planted model N19: Keller-Segel cells with a self-made attractant and a consumed nutrient (plug-in, PLANTED_PIPELINE_PLAN.md).

Source: arXiv 2610.32591 v1, "Model" (Eqs. 1a-1c), "Phenomenology" and Fig. 2 (caption: kappa_a = 5, kappa_n = 14,
delta_a = 2.5, delta_n = 0.75, beta = 0.02, system size L = 50), "Linear stability analysis" (Eq. 2, Fig. 3: phase diagram
for kappa_n = 14, delta_n = 0.75, delta_a = 2.5, and dispersion relations B-K).

Model as stated in the paper (dimensionless; time in 1/lambda, length in sqrt(D_rho/lambda)), fields b = cell density
(F0), a = attractant / autoinducer (F1), n = nutrient (F2), periodic 2D:
    db/dt = lap b - div[ kappa_a b/(1+b^2) grad a ] - div[ kappa_n b/(1+b^2) grad n ]       (1a)
    da/dt = delta_a lap a + b - a                                                            (1b)
    dn/dt = delta_n lap n - beta b n + beta                                                  (1c)
Homogeneous state (1, 1, 1); started at 1 plus small noise ("starting from the homogeneous steady state with small
fluctuations", Fig. 1 caption).

Headline: point K of Fig. 3A (plan deviation P1, logged by the coordinator before any pipeline output existed)
------------------------------------------------------------------------------------------------------------------
kappa_a = 11.1, beta = 0.084, kappa_n = 14, delta_a = 2.5, delta_n = 0.75 (the last three from the Fig. 3 caption,
identical to the Fig. 2 caption). Why: the plan row named "the paper's parameters (card)", i.e. the Fig. 2 caption
(kappa_a = 5, beta = 0.02), but those values lie in the paper's own linearly stable region, while the text ties the
Fig. 2 pattern to point K: "The pattern highlighted in Fig. 2 is an example from this instability region,
corresponding to the dispersion relation shown in Fig. 3K" (the red "exceptional instability" region of Fig. 3A).
The plan's general rule is a headline inside the range where the paper reports the behavior, so the headline is K.
- Linear-stability evidence (paper's Eq. 2, eigenvalues of the linearised Eqs. 1 at (1, 1, 1), 2 pi/50 <= q <= 3):
  caption values kappa_a = 5, beta = 0.02: max Re sigma = -9.9e-4 (oscillatory branch, q = 0.24; stable; the
  oscillatory threshold at kappa_a = 5 is beta ~ 0.018, the stationary threshold at beta = 0.02 is kappa_a ~ 5.44);
  runs here stay uniform, also from noise of amplitude 0.3. Point K (kappa_a = 11.1, beta = 0.084): max Re sigma =
  +0.212 at q = 0.59, with oscillatory and stationary modes unstable in the same band (mixed / exceptional type), as
  Fig. 3K shows.
- How K was read off Fig. 3A: page 4 of the arXiv PDF rendered at 300 dpi; the phase-diagram panel cropped. Axis
  ranges: kappa_a axis linear with ticks "2" and "kappa_a^c = kappa_n + 2" (= 16; the orange vertical line) at pixel
  columns 95 and 606 (36.5 px per unit; the left frame edge then sits at kappa_a ~ 0.0); beta axis linear with ticks
  "0" and "beta^fw" at pixel rows 838 and 186, beta^fw = delta_n kappa_n / [delta_a (2 + kappa_n)] = 0.2625 (text).
  The purple "x" marker labelled K has its centre (centroid of the marker-coloured pixels, bounding box 15 x 16 px)
  at column 427.1, row 629.4 -> kappa_a = 2 + 332.1/36.5 = 11.10, beta = 0.2625 x 208.6/652 = 0.0840. Uncertainty:
  half the marker size, +-0.2 in kappa_a and +-0.003 in beta. Check of the axis reading: the region boundaries drawn
  in Fig. 3A at beta = 0.035, 0.084 and 0.136 (pixel rows 750, 629, 500) fall at kappa_a = 5.75 / 6.8-6.9 / 14.45,
  10.0-10.2 / 14.7 and 9.1 / 12.5-13.2 / 15.0, against 5.8 / 6.9 / 14.4, 10.1 / 14.7 and 9.2 / 12.6 / 15.0 from the
  Eq. 2 classification computed here (onset of oscillatory, of mixed, of purely stationary instability), i.e. within
  0.1 in kappa_a. (A coarser 150-dpi reading, used for the scheme test below, gave beta 0.001 lower.)

Plan rules (row N19, numbers as in the frozen row; at the headline 5 -> 11.1 and 0.02 -> 0.084 per P1): 1 b drifts up
grad a with strength 5 b/(1+b^2) (I); 2 b drifts up grad n with strength 14 b/(1+b^2) (I); 3 a produced at rate b;
4 a decays at rate 1; 5 n consumed at rate 0.02 b n; 6 n resupplied at rate 0.02. Diffusion (cells 1, a 2.5, n 0.75)
is header. Knock-out: rules 1 and 2 removed (no drift; cells only diffuse).

Choices the paper does not state (and why)
------------------------------------------
- Grid spacing: not given (no supplement on arXiv). dx = 1 dimensionless length unit, so the paper's L = 50 is a
  50 x 50 grid (as the card reads it). Census: 64 x 64 sites at dx = 1 (domain 64, slightly larger than the paper's
  50); area4: 128 x 128. The finest structures (aggregate core radius ~5, edge width ~2-3, from Fig. 2) are resolved.
- Discretisation: conservative finite volumes at dx = 1. Cells: Scharfetter-Gummel fluxes on every face (diffusion
  and drift together; drift velocity (kappa_a grad a + kappa_n grad n)/(1 + b^2) on the face, b averaged over the two
  cells), explicit Euler with the step limited so that every cell's outgoing coefficients stay below 0.9 (keeps b >= 0
  and the cell number exactly conserved). Attractant and nutrient: implicit diffusion with the 5-point Laplacian (exact
  FFT solve), implicit attractant decay, explicit production, consumption and resupply. dt = min(0.01, that limit);
  each census step is integrated exactly to one time unit (last sub-step shortened). Why this scheme: at dx = 1 the
  aggregate edges are steep; a MUSCL-upwinded mobility overshoots badly (aggregate peak 36 vs 11-13 on finer grids)
  and a central mobility goes negative, while Scharfetter-Gummel at dx = 1 matched a MUSCL run at dx = 0.25 near
  point K (test run at beta = 0.083; t = 400: peak b 11.4 vs 11.6, 2 vs 2 aggregates, aggregate mass fraction 0.59 vs
  0.59, background std 0.24 vs 0.29).
- Initial condition: b, a, n = 1 + U(-0.01, 0.01) independently per site and seed ("small fluctuations"; amplitude not
  stated).
- One census step = one unit of the paper's dimensionless time (plan: continuous-time model).
- Blow-up rule as in the census interpreter (non-finite or |value| > 1e6 -> meta unstable = True, stop all seeds);
  the volume-filling drift prevents collapse, so it is not expected. A census step needing more than 100,000
  sub-steps (dt < 1e-5) is treated the same way (unstable, stop) rather than truncated.

Fidelity check (plan: one dense aggregate coexisting with travelling bands; per P1 at point K and the paper's domain).
Quantitative criterion, fixed before any large run (from Fig. 2: aggregate core b ~ 10, band densities <~ 3-4):
  per frame (frames every 10 time units, late window t >= T/2):
  - dense aggregates = periodic 4-connected components of {b > 5} with at least 4 sites; N_agg = their number;
  - aggregate mass fraction f = (sum of b over them) / (total b);
  - background = sites farther than 3 sites from any b > 5 site; bg_std = spatial std of b there;
  - background motion R = rms[b(t + 10) - b(t)] / bg_std(t) over sites that are background in both frames.
  A frame pair satisfies the behavior if N_agg = 1, 0.05 <= f <= 0.9 (part of the population outside the aggregate),
  bg_std >= 0.1 (bands present) and R >= 0.3 (the background moves; a stationary pattern gives R ~ 0).
  A run shows the behavior if >= 50 % of its late-window frame pairs do; the check passes if a majority of seeds do.
  Fidelity runs (pass/fail): the paper's domain (50 x 50, dx = 1), headline = point K, 3 seeds, T = 5000 time units
  (cut to fit the budget, not below 1000). Reported alongside (not part of pass/fail): the same measure with
  "N_agg >= 1", the maximal linear growth rates from Eq. 2, and a diagnostic run at the Fig. 2 caption values
  (kappa_a = 5, beta = 0.02; expected to stay uniform).
Small-scale evidence here (criterion fixed before these runs): at the headline K (kappa_a = 11.1, beta = 0.084,
50 x 50, 1 seed, t = 1500) one dense aggregate (peak b ~ 12.9, 63 % of the cells) sits in travelling bands and the
criterion holds in 100 % of late frame pairs; coarsening (kappa_a = 20, beta = 0.02) and stationary spots (kappa_a = 6,
beta = 0.02) fail it; the caption values stay uniform. At census scale (64 x 64, 1,000 steps, 3 seeds) the headline
still has 2-3 aggregates in moving bands (slow merging in a domain larger than the paper's), so the frozen criterion
("exactly one") fails while its "at least one" variant holds in every seed; presence() reports both.
"""
import time

import numpy as np
import scipy.fft as sfft
from scipy import ndimage

from . import base

CID = "planted:N19"
LAYERS = "F"
K = {}
SPEED = 0.0
RULES = [("B_DRIFT_UP_GRAD_A", True),         # 1  -div[kappa_a b/(1+b^2) grad a]  (I)
         ("B_DRIFT_UP_GRAD_N", True),         # 2  -div[kappa_n b/(1+b^2) grad n]  (I)
         ("A_PRODUCED_BY_B", False),          # 3  +b in the a equation (one-way: reads the site's own cells)
         ("A_DECAY", False),                  # 4  -a
         ("N_CONSUMED_BY_B", False),          # 5  -beta b n (one-way: reads the site's own cells)
         ("N_RESUPPLY", False)]               # 6  +beta
HEADLINE = {
    "kappa_a": 11.1,                          # point K of Fig. 3A, read off the figure (+-0.2); plan deviation P1
    "kappa_n": 14.0,                          # Fig. 3 caption (= Fig. 2 caption, card)
    "delta_a": 2.5, "delta_n": 0.75,          # Fig. 3 caption (= Fig. 2 caption, card: 2.5x and 0.75x cell diffusion)
    "beta": 0.084,                            # point K of Fig. 3A, read off the figure (+-0.003); plan deviation P1
    "dx": 1.0,                                # not stated; paper's L = 50 read as 50 x 50 (card)
    "init_noise": 0.01,                       # not stated ("small fluctuations")
    "dt_max": 0.01,
    "W_census": 64, "W_area4": 128,
}
STEP = ("one census step = 1 unit of the paper's dimensionless time (1/lambda); Euler sub-steps (cells explicit "
        "Scharfetter-Gummel, chemicals implicit diffusion), dt <= 0.01 and positivity-limited, integrated exactly to the unit")
RULE_TEXT = ("World: three continuous fields on a 64x64 periodic square: cell density b, an attractant a made by the "
             "cells, and a nutrient n eaten by the cells; all start uniform at 1 plus small noise. Rules "
             "(dimensionless): (1) b diffuses and drifts up the gradient of a with strength 11.1*b/(1+b^2) and up the "
             "gradient of n with strength 14*b/(1+b^2); cell number is conserved; (2) a diffuses 2.5 times faster "
             "than cells, is produced at rate b and decays at rate 1; (3) n diffuses 0.75 times as fast as cells, is "
             "consumed at rate 0.084*b*n and resupplied at rate 0.084 everywhere.")

DT_MAX = 0.01
NSUB_LIMIT = 100000                          # sub-steps per census step (dt >= 1e-5); exceeding it counts as unstable
B_DENSE, A_MIN, BG_DIST, BG_STD_MIN, MOTION_MIN = 5.0, 4, 3, 0.1, 0.3


def _bernoulli(x):
    """B(x) = x / (exp(x) - 1), B(0) = 1 (Scharfetter-Gummel weight)."""
    x = np.asarray(x, dtype=float); out = np.empty_like(x); sm = np.abs(x) < 1e-5
    out[sm] = 1.0 - 0.5 * x[sm] + x[sm] ** 2 / 12.0
    xs = x[~sm]; out[~sm] = xs / np.expm1(xs)
    return out


class _Field:
    def __init__(self, W, gens, on, ka=HEADLINE["kappa_a"], kn=HEADLINE["kappa_n"], da=HEADLINE["delta_a"],
                 dn=HEADLINE["delta_n"], beta=HEADLINE["beta"], dx=HEADLINE["dx"], eps=HEADLINE["init_noise"]):
        self.W, self.S, self.dx = W, len(gens), dx
        self.ka = ka if 1 in on else 0.0
        self.kn = kn if 2 in on else 0.0
        self.da, self.dn, self.beta = da, dn, beta
        k = 2 * np.pi * np.fft.fftfreq(W, d=dx); kr = 2 * np.pi * np.fft.rfftfreq(W, d=dx)
        KY, KX = np.meshgrid(k, kr, indexing="ij")
        self.lap = -(4.0 / dx ** 2) * (np.sin(KX * dx / 2) ** 2 + np.sin(KY * dx / 2) ** 2)   # 5-point symbol
        X = np.empty((self.S, 3, W, W))
        for s, g in enumerate(gens): X[s] = 1.0 + g.uniform(-eps, eps, (3, W, W))
        self.b, self.a, self.n = X[:, 0].copy(), X[:, 1].copy(), X[:, 2].copy()

    def _cells(self):
        """d b/dt from the Scharfetter-Gummel fluxes (diffusion + drift on every face) and the largest stable,
        positivity-preserving explicit step."""
        b, dx = self.b, self.dx
        c = self.ka * self.a + self.kn * self.n
        div = np.zeros_like(b); out = np.zeros_like(b)
        for ax in (1, 2):
            bp = np.roll(b, -1, ax)
            g = (np.roll(c, -1, ax) - c) / dx                   # drift strength kappa . grad on the + face
            bf = 0.5 * (b + bp)
            P = g / (1.0 + bf * bf) * dx                        # face Peclet number of the velocity g/(1+b^2)
            Bm, Bp = _bernoulli(-P), _bernoulli(P)
            J = (Bm * b - Bp * bp) / dx                         # SG flux through the + face
            div += (J - np.roll(J, 1, ax)) / dx
            out += Bm + np.roll(Bp, 1, ax)                      # coefficients of b_i in its outgoing fluxes
        dt_max = 0.9 * dx * dx / out.max()
        return -div, dt_max

    def advance(self):
        """one time unit; returns (number of sub-steps, completed). A census step that would need more than
        NSUB_LIMIT sub-steps is not truncated silently: it returns completed = False and the run is marked unstable."""
        f2, i2, W = sfft.rfft2, sfft.irfft2, self.W
        t = 0.0; nsub = 0
        while t < 1.0 - 1e-12:
            if nsub >= NSUB_LIMIT: return nsub, False
            Nb, dt_pos = self._cells()
            dt = min(DT_MAX, 1.0 - t, dt_pos)
            b_old = self.b
            self.b = b_old + dt * Nb
            self.a = i2(f2(self.a + dt * b_old) / (1 - dt * (self.da * self.lap - 1.0)), s=(W, W))
            self.n = i2(f2(self.n + dt * self.beta * (1.0 - b_old * self.n)) / (1 - dt * self.dn * self.lap), s=(W, W))
            t += dt; nsub += 1
        return nsub, True

    def stack(self):
        return np.stack([self.b, self.a, self.n], 1)

    def ok(self):
        X = self.stack()
        return bool(np.all(np.isfinite(X)) and np.abs(X).max() <= 1e6)


def _run(W, gens, on, steps, **par):
    fld = _Field(W, gens, on, **par); rec = base.REC["F"]
    frames = [fld.stack().astype(np.float32)]; t = 0; unstable = False; nsub_tot = 0
    while t < steps:
        ns, done = fld.advance(); nsub_tot += ns; t += 1
        bad = (not done) or (not fld.ok())
        if t % rec == 0: frames.append(fld.stack().astype(np.float32))
        if bad: unstable = True; break
    return np.stack(frames, 1), t, unstable, nsub_tot


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    on = base.kept(keep, RULES)
    W = base.SIZES[size]["W"]
    F, t, unstable, _ = _run(W, base.rngs(seed0, S), on, steps)
    return {"F": F, "meta": base.meta(t, W, unstable=unstable)}


# ---------------------------------------------------------------------------------------------------- measure
def _components(mask):
    """periodic 4-connected components: labels array and number of components."""
    lab, n = ndimage.label(mask)
    if n == 0: return lab, 0
    parent = np.arange(n + 1)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for e1, e2 in ((lab[0, :], lab[-1, :]), (lab[:, 0], lab[:, -1])):
        for p, q in zip(e1, e2):
            if p and q:
                rp, rq = find(p), find(q)
                if rp != rq: parent[max(rp, rq)] = min(rp, rq)
    roots = np.array([find(i) for i in range(n + 1)])
    uniq = np.unique(roots[1:]); remap = np.zeros(n + 1, int); remap[uniq] = np.arange(1, uniq.size + 1)
    return remap[roots][lab], uniq.size


def frame_stats(b, b_next=None):
    """the per-frame quantities of the criterion for one W x W cell-density frame."""
    lab, n = _components(b > B_DENSE)
    areas = np.bincount(lab.ravel(), minlength=n + 1)[1:] if n else np.zeros(0, int)
    big = np.flatnonzero(areas >= A_MIN) + 1
    agg = np.isin(lab, big)
    nagg = int(big.size); frac = float(b[agg].sum() / b.sum()) if nagg else 0.0
    yy, xx = np.mgrid[-BG_DIST:BG_DIST + 1, -BG_DIST:BG_DIST + 1]
    disk = (yy ** 2 + xx ** 2) <= BG_DIST ** 2
    near = ndimage.maximum_filter((b > B_DENSE).astype(np.uint8), footprint=disk, mode="wrap").astype(bool)
    bg = ~near
    bg_std = float(b[bg].std()) if bg.sum() > 10 else 0.0
    R = float("nan")
    if b_next is not None and bg_std > 0:
        near2 = ndimage.maximum_filter((b_next > B_DENSE).astype(np.uint8), footprint=disk, mode="wrap").astype(bool)
        both = bg & ~near2
        if both.sum() > 10: R = float(np.sqrt(((b_next[both] - b[both]) ** 2).mean()) / bg_std)
    return {"n_agg": nagg, "agg_mass_frac": frac, "bg_std": bg_std, "motion": R, "b_max": float(b.max())}


def behaviour(Fb, rec=10, t_from=None):
    """apply the criterion to one run's cell-density frames Fb (frames, W, W); late window t >= t_from (default T/2)."""
    nf = Fb.shape[0]; T = (nf - 1) * rec; t_from = T / 2 if t_from is None else t_from
    idx = [i for i in range(nf - 1) if i * rec >= t_from]
    st = [frame_stats(np.asarray(Fb[i], float), np.asarray(Fb[i + 1], float)) for i in idx]

    def ok(s, strict=True):
        n_ok = (s["n_agg"] == 1) if strict else (s["n_agg"] >= 1)
        return bool(n_ok and 0.05 <= s["agg_mass_frac"] <= 0.9 and s["bg_std"] >= BG_STD_MIN
                    and np.isfinite(s["motion"]) and s["motion"] >= MOTION_MIN)
    share = float(np.mean([ok(s) for s in st])) if st else 0.0
    share_loose = float(np.mean([ok(s, False) for s in st])) if st else 0.0
    nag = np.array([s["n_agg"] for s in st]) if st else np.zeros(0)
    agg = lambda k: float(np.nanmean([s[k] for s in st])) if st else float("nan")
    return {"share_pairs_ok": share, "shows_behaviour": share >= 0.5, "share_pairs_ok_n_agg_ge_1": share_loose,
            "n_agg_hist": {int(k): int((nag == k).sum()) for k in np.unique(nag)} if st else {},
            "mean_agg_mass_frac": agg("agg_mass_frac"), "mean_bg_std": agg("bg_std"), "mean_motion": agg("motion"),
            "mean_b_max": agg("b_max"), "pairs": len(st), "t_from": t_from}


def linear_growth(ka=HEADLINE["kappa_a"], kn=HEADLINE["kappa_n"], da=HEADLINE["delta_a"], dn=HEADLINE["delta_n"],
                  beta=HEADLINE["beta"], L=50.0):
    """max over q of Re sigma(q) from the paper's Eq. 2 (eigenvalues of the linearised Eqs. 1 at (1, 1, 1)), scanned
    over 2 pi / L <= q <= 3 (L = 50: the paper's domain; the conserved q -> 0 mode is excluded)."""
    best = (-np.inf, 0.0, 0.0)
    for q in np.linspace(2 * np.pi / L, 3.0, 3000):
        Q = q * q
        J = np.array([[-Q, ka / 2 * Q, kn / 2 * Q], [1.0, -1 - da * Q, 0.0], [-beta, 0.0, -beta - dn * Q]])
        ev = np.linalg.eigvals(J); i = np.argmax(ev.real)
        if ev[i].real > best[0]: best = (float(ev[i].real), float(q), float(abs(ev[i].imag)))
    return {"max_re_sigma": best[0], "q": best[1], "omega": best[2]}


def fidelity(budget_s=7200, _small=False):
    """pass/fail at the headline (point K); the Fig. 2 caption values are a reported diagnostic only.
    _small=True only exercises the code path (1 seed, t = 300); it is not the fidelity check."""
    t_start = time.time(); W = 50; S = 1 if _small else 3
    T, Tmin = (300, 300) if _small else (5000, 1000)
    CAPTION = dict(ka=5.0, beta=0.02)
    fld = _Field(W, base.rngs(1, S), {1, 2})                          # calibrate on the (patterned) headline
    for _ in range(200): fld.advance()
    t0 = time.time()
    for _ in range(20): fld.advance()
    unit = (time.time() - t0) / 20
    while T > Tmin and 2 * T * unit * 1.3 > 0.85 * budget_s - (time.time() - t_start): T -= 250
    out = {"setting": dict(HEADLINE, W=W, seeds=S, T=T), "cost_s_per_time_unit_patterned": unit,
           "paper": {"source": "Fig. 2 (said to correspond to dispersion K of Fig. 3, the red 'exceptional' region; "
                               "L = 50): one dense droplet (b ~ 10) moving in a background of travelling bands",
                     "dense_aggregate_count": 1, "coexisting_structure_types": 2},
           "deviation": "P1: headline = point K of Fig. 3A (kappa_a 11.1 +- 0.2, beta 0.084 +- 0.003) instead of the "
                        "Fig. 2 caption values (kappa_a 5, beta 0.02), which are linearly stable by Eq. 2",
           "tolerance": "per run: >= 50 % of late-window frame pairs with N_agg = 1 (b > 5, >= 4 sites), aggregate mass "
                        "fraction 0.05-0.9, background std >= 0.1, background motion R >= 0.3; pass if a majority of "
                        "seeds show it (at the headline, point K)",
           "linear_theory_Eq2": {"headline_K": linear_growth(), "caption_values": linear_growth(**CAPTION)}}
    runs = {}
    for name, par in (("headline_K", {}), ("caption_values_diagnostic", CAPTION)):
        t0 = time.time()
        F, t, unstable, nsub = _run(W, base.rngs(20261103, S), {1, 2}, T, **par)
        per = [behaviour(F[s, :, 0]) for s in range(S)]
        runs[name] = {"params": par or {"ka": HEADLINE["kappa_a"], "beta": HEADLINE["beta"]}, "unstable": unstable,
                      "steps_run": t, "substeps": nsub, "seconds": time.time() - t0, "per_seed": per,
                      "seeds_showing": int(sum(p["shows_behaviour"] for p in per)),
                      "passes": bool(sum(p["shows_behaviour"] for p in per) > S / 2)}
    out.update(runs)
    out["pass"] = runs["headline_K"]["passes"]
    out["seconds"] = time.time() - t_start
    return out


def presence(out):
    """the fidelity criterion on census-scale output of the real program (late window t >= 500); present if >= 2 of
    3 seeds show the behavior."""
    F = np.asarray(out["F"])
    if out["meta"].get("unstable"):
        return {"present": False, "note": "run unstable (blow-up)", "steps_run": out["meta"]["steps_run"]}
    rec = out["meta"]["rec"]["F"]
    per = [behaviour(F[s, :, 0], rec=rec, t_from=500) for s in range(F.shape[0])]
    k = sum(p["shows_behaviour"] for p in per)
    k1 = sum(p["share_pairs_ok_n_agg_ge_1"] >= 0.5 for p in per)
    return {"per_seed": per, "seeds_showing": int(k), "present": bool(k >= 2),
            "seeds_showing_n_agg_ge_1": int(k1), "present_n_agg_ge_1": bool(k1 >= 2),
            "tolerance": "decision: per seed >= 50 % of frame pairs (t >= 500) with exactly one dense aggregate in "
                         "moving bands, in >= 2 of 3 seeds; present_n_agg_ge_1 (reported only): the same with at "
                         "least one aggregate, to show a scale effect (slow merging in the 64 x 64 world)"}

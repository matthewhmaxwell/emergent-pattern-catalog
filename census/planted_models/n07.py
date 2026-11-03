"""Planted model N07: chiral mass transfer on a square periodic lattice (plug-in for PLANTED_PIPELINE_PLAN.md).

Source: arXiv 2609.16821 v1, Sec. II "Model" (Eqs. 2-8, Fig. 1), Sec. VI "Isotropic mass transfer" (Eqs. 41-45), Sec.
VI.3 "Density correlations" (Eqs. 54-66, Figs. 5-6).

Model as stated in the paper
----------------------------
Sites (x, y) of a square lattice with periodic boundaries carry a continuous mass m(x, y) >= 0 (global density
rho = 1). A mass-transfer event: a site (x, y) is selected at random and one of the four directions +x, -x, +y, -y is
chosen with probabilities p1..p4 (isotropic bulk, Eq. 42: p1 = p2 = p3 = p4 = 1/4). Two steps happen within the same
infinitesimal time interval (Eqs. 4-7):
  (A) the selected site keeps lambda_1 m + xi_1 (1 - lambda_1) m and chips off (1 - xi_1)(1 - lambda_1) m to the
      nearest neighbour in the chosen direction;
  (B) that neighbour pools the chip with its own mass m', keeps lambda_2 m' + xi_2 (1 - lambda_2) m' and passes
      (1 - xi_2)(1 - lambda_2) m' one step further along the direction rotated by 90 degrees counter-clockwise
      (Eqs. 2 and 8: (x+1,y)->(x+1,y+1), (x-1,y)->(x-1,y-1), (x,y+1)->(x-1,y+1), (x,y-1)->(x+1,y-1)).
xi_1, xi_2 are independent U(0, 1). Isotropic setting of Secs. VI.3 and Figs. 4-6: lambda_x = lambda_y = lambda_1 = 0.4,
lambda_xy = lambda_yx = lambda_2 = 0.6, rho = 1. Hydrodynamics (Eqs. 44-45): D = (1/4)[mu1 (l1~ + l2~) + mu1^2 l1~ l2~]
= 0.14, |D_odd| = (1/2) mu1^2 l1~ l2~ = 0.03, with l~ = 1 - lambda, mu1 = 1/2. Steady-state correlations (Eq. 66):
C(r, theta) ~ -(3 s2 / pi) cos(4 theta) / |r|^4, i.e. C(x, 0) ~ 0.11 / x^4 along the axes (Fig. 6, 800 x 800).
The achiral case (no second step, l2~ = 0) has delta-correlated masses (Sec. VI.3: "in the absence of the off-diagonal
components, the density correlations remain short-ranged (delta-correlated for the present models)").

Plan rules (PLANTED_PIPELINE_PLAN.md row N07)
---------------------------------------------
  1 sender hops its remainder one step in a random direction (one-way: the sender reads only its own mass)
  2 receiver pools, keeps its share and passes the rest one step to its left (I: the receiver's update reads the mass
    sent by another site)
Knock-out (rule 2 removed): the receiver keeps what it gets (lambda_2 = 1), which is exactly the paper's achiral model.

Choices the paper does not state (and why)
------------------------------------------
- Time: continuous-time Markov process with every site selected at rate 1 (the paper's microscopic currents, Eqs.
  10-11, carry the rate factors 1/4 = p_i per site per unit time). One census step = one sweep = W*W random-sequential
  events = one unit of the paper's time (plan: Monte Carlo / random-sequential models). Events are executed exactly in
  the order of a random sequence: each sweep draws W*W (site, direction, xi_1, xi_2) per seed; the sequence is cut
  into chunks and, inside a chunk, every event whose three sites are untouched by all earlier pending events is
  applied at once (their footprints are disjoint, so they commute); the others are deferred, in order, to the front of
  the next chunk. This reproduces the sequential order exactly (an event never overtakes an earlier event it shares a
  site with). Seeds are interleaved in one batch (they never share sites).
- Lattice orientation: x = column index (last array axis), y = row index; "counter-clockwise" is with y increasing
  with the row index (on an image drawn with row 0 at the top it looks clockwise; the statistics are mirror-blind).
- Initial condition: the paper does not state one for the steady-state runs; every site starts at m = 1 (the global
  density). The steady state does not depend on it.
- Census size: 64 x 64 sites (area4 128 x 128) instead of 800 x 800 (plan: lattice side 64 at the paper's spacing).
- Blow-up check (non-finite or |m| > 1e6) as in the census interpreter; mass is conserved and non-negative, so it
  never triggers.

Headline: lambda_1 = 0.4, lambda_2 = 0.6, rho = 1 (the paper's isotropic setting, Figs. 4-6; card: keep 40% / 60%).

Fidelity check (plan: density correlations decaying close to r^-4, fitted exponent 3.5 to 4.5, where the achiral version
does not). Fixed before any large run:
- Lattice 800 x 800 (the paper's Fig. 6 size) if the timing calibration allows at least 2,000 measured sweeps per
  version inside the budget, else the largest of 512, 400, 256 that does. Chiral (rules 1+2) and achiral (rule 2 off,
  the knock-out = the paper's achiral limit) runs, each: 1,500 burn-in sweeps from m = 1, then one snapshot per sweep,
  up to ~45% of the budget per version.
- Measure: the equal-time density-density correlation C(dx, dy) = <dm(i) dm(i + r)> (the paper's C^mm, Figs. 5-6).
  For each integer distance r = 1..10, C is fitted over the ring of lattice vectors r - 0.5 <= |r'| < r + 0.5 by
  c0(r) + c4(r) cos(4 theta) (the paper's asymptotic form Eq. 66 has exactly this angular structure); c4(r) equals the
  axis correlation C(r, 0) of Fig. 6 in the asymptotic form and is immune to isotropic offsets (the finite-size and
  slow-mode offsets go into c0). The exact lattice theory (Eqs. 57, 60, 67-68, evaluated numerically in theory_C())
  gives c4(r) independent of L down to L = 64 and a fitted exponent 4.06 on r = 2..8.
- Fit: c4(r) = A r^-eta on r = 2..8, weighted least squares with errors from 20 time blocks.
- Tolerance: chiral passes if 3.5 <= eta <= 4.5 and A > 3 sigma_A; the achiral version "does not" if that conjunction
  fails (expected: A consistent with zero). Pass = both. Also reported: axis C(r, 0) (Fig. 6), on-site variance, and
  the exact-theory values (variance 0.415, C(2,0) = 7.2e-3, C(4,0) = 5.0e-4, which match Fig. 6's points) with the
  theory's own exponent under the same weights (the weighted fit leans on r = 2-3, where lattice corrections give a
  local exponent ~3.7-3.8; small-scale tests here: eta = 3.79 +- 0.02 at L = 256, theory c4 matched within errors).
"""
import time

import numpy as np

from . import base

CID = "planted:N07"
LAYERS = "F"
K = {}
SPEED = 0.0
RULES = [("SENDER_HOPS_REMAINDER_RANDOM_DIRECTION", False),        # rule 1, one-way
         ("RECEIVER_POOLS_PASSES_LEFT", True)]            # rule 2, interaction (I)
HEADLINE = {
    "lambda_1": 0.4,        # first-step retention, Sec. VI / Figs. 4-6 ("lambda_1 = 0.4")
    "lambda_2": 0.6,        # second-step retention, Sec. VI / Figs. 4-6 ("lambda_2 = 0.6")
    "xi": "U(0,1)",         # random shares xi_1, xi_2, Sec. II (mu_n = 1/(n+1))
    "p_direction": 0.25,    # Eq. 42, isotropic bulk
    "rho": 1.0,             # global density, Figs. 4-6
    "turn": "90 deg counter-clockwise (Eqs. 2, 8)",
    "W_census": 64, "W_area4": 128,   # plan census scale (paper: 800 x 800)
}
STEP = ("one census step = one sweep = W*W random-sequential mass-transfer events (one unit of the paper's "
        "continuous time, each site fires at rate 1); executed exactly in sequence order via conflict-free batches")
RULE_TEXT = ("World: a 64x64 square lattice with wrap-around edges; each site holds a continuous mass (mean 1 per site). "
             "Rules, repeated: (1) pick a random site and one of the four directions at random; (2) the site keeps 40% "
             "of its mass plus a random share of the rest and sends the remainder one step in that direction; (3) the "
             "receiving site pools it with its own mass, keeps 60% plus a random share of the rest, and passes the "
             "remainder one step to its left (a 90-degree counter-clockwise turn), so mass moves along L-shaped paths "
             "that always turn the same way. Rates are identical for all four directions and total mass is conserved. "
             "Measure steady-state density-density correlations and long-wavelength fluctuations.")

CHIP1 = 1.0 - HEADLINE["lambda_1"]      # fraction at risk in step A (0.6)
CHIP2 = 1.0 - HEADLINE["lambda_2"]      # fraction at risk in step B (0.4)
_DX = np.array([1, -1, 0, 0]); _DY = np.array([0, 0, 1, -1])        # directions +x, -x, +y, -y
_LX = np.array([0, 0, -1, 1]); _LY = np.array([1, -1, 0, 0])        # their left (counter-clockwise) turns


# ---------------------------------------------------------------------------------------------------- dynamics
def _neighbours(W):
    """J[site, d] = neighbour in direction d; Kc[site, d] = the site one further step to the left of d."""
    y, x = np.divmod(np.arange(W * W), W)
    idx = lambda yy, xx: (yy % W) * W + (xx % W)
    J = np.stack([idx(y + _DY[d], x + _DX[d]) for d in range(4)], 1)
    Kc = np.stack([idx(y + _DY[d] + _LY[d], x + _DX[d] + _LX[d]) for d in range(4)], 1)
    return J, Kc


class _Lattice:
    """S independent W x W lattices advanced together; m has shape (S * W * W,)."""

    def __init__(self, W, gens, chiral, m0=None, frac=0.1):
        self.W, self.N, self.S = W, W * W, len(gens)
        self.gens, self.chiral, self.frac = gens, chiral, frac
        self.J, self.Kc = _neighbours(W)
        self.m = np.ones(self.S * self.N) if m0 is None else m0.reshape(-1).astype(float).copy()
        self.stamp = np.full(self.S * self.N, np.iinfo(np.int64).max, np.int64)

    def sweep(self):
        N, S = self.N, self.S
        site = np.empty((N, S), np.int64); d = np.empty((N, S), np.int64); xi = np.empty((N, S, 2))
        for s, g in enumerate(self.gens):          # each seed's own random stream, interleaved event by event
            site[:, s] = g.integers(0, N, N); d[:, s] = g.integers(0, 4, N); xi[:, s] = g.random((N, 2))
        off = (np.arange(S) * N)[None, :]
        gi = (off + site).ravel(); gj = (off + self.J[site, d]).ravel(); gk = (off + self.Kc[site, d]).ravel()
        xi = xi.reshape(-1, 2); E = gi.size
        m, stamp, big = self.m, self.stamp, np.iinfo(np.int64).max
        c = max(64, int(self.frac * E)); pos = 0; pend = np.empty(0, np.int64)
        while pos < E or pend.size:
            new = np.arange(pos, min(E, pos + c)); pos += new.size
            ev = np.concatenate([pend, new]) if pend.size else new
            fi, fj = gi[ev], gj[ev]; rank = np.arange(ev.size)
            np.minimum.at(stamp, fi, rank); np.minimum.at(stamp, fj, rank)
            if self.chiral:
                fk = gk[ev]; np.minimum.at(stamp, fk, rank)
                free = (stamp[fi] == rank) & (stamp[fj] == rank) & (stamp[fk] == rank)
                stamp[fk] = big
            else:
                free = (stamp[fi] == rank) & (stamp[fj] == rank)
            stamp[fi] = big; stamp[fj] = big
            e = ev[free]; i, j = gi[e], gj[e]
            chip = (1.0 - xi[e, 0]) * CHIP1 * m[i]                # rule 1: sender keeps l1 + xi1 (1 - l1)
            m[i] -= chip
            if self.chiral:                                       # rule 2: receiver pools, keeps, passes left
                k = gk[e]; pooled = m[j] + chip
                pas = (1.0 - xi[e, 1]) * CHIP2 * pooled
                m[j] = pooled - pas; m[k] += pas
            else:                                                 # knock-out: receiver keeps what it gets
                m[j] += chip
            pend = ev[~free]

    def field(self):
        return self.m.reshape(self.S, self.W, self.W)


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    on = base.kept(keep, RULES)
    W = base.SIZES[size]["W"]
    lat = _Lattice(W, base.rngs(seed0, S), chiral=(2 in on))
    rec = base.REC["F"]; frames = [lat.field().astype(np.float32)]; unstable = False; t = 0
    while t < steps:
        lat.sweep(); t += 1
        f = lat.field()
        bad = (not np.all(np.isfinite(f))) or np.abs(f).max() > 1e6
        if t % rec == 0: frames.append(f.astype(np.float32))
        if bad: unstable = True; break
    F = np.stack(frames, 1)[:, :, None]                      # (S, frames, 1, W, W)
    return {"F": F, "meta": base.meta(t, W, unstable=unstable)}


# ---------------------------------------------------------------------------------------------------- measure
def _rings(W, rmax=10):
    d = np.arange(-rmax, rmax + 1); DY, DX = np.meshgrid(d, d, indexing="ij")
    R = np.hypot(DX, DY); TH = np.arctan2(DY, DX); out = []
    for r in range(1, rmax + 1):
        msk = (R >= r - 0.5) & (R < r + 0.5)
        X = np.stack([np.ones(msk.sum()), np.cos(4 * TH[msk])], 1)
        out.append((DY[msk] % W, DX[msk] % W, np.linalg.pinv(X)))
    return out


def _corr_from_power(P, W):
    """C[dy, dx] = <dm(i) dm(i + r)> from the summed |rfft2(dm)|^2 of nsnap snapshots (P already divided by nsnap)."""
    return np.fft.irfft2(P, s=(W, W)) / (W * W)


def _ring_stats(C, rings):
    co = np.array([P @ C[iy, ix] for iy, ix, P in rings])      # (rmax, 2): c0(r), c4(r)
    axis = np.array([0.5 * (C[0, r] + C[r, 0]) for r in range(1, len(rings) + 1)])
    return co[:, 0], co[:, 1], axis, C[0, 0]


def _fit_power(r, y, sig):
    """weighted fit y = A r^-eta; returns eta, sd_eta, A, sd_A, and A4 (eta fixed at 4) with its sd."""
    w = 1.0 / np.maximum(sig, 1e-300) ** 2
    x4 = r ** -4.0
    A4 = float(np.sum(w * x4 * y) / np.sum(w * x4 * x4)); sA4 = float(np.sqrt(1.0 / np.sum(w * x4 * x4)))
    res = {"A_eta4": A4, "sd_A_eta4": sA4, "eta": float("nan"), "sd_eta": float("nan"), "A": float("nan"),
           "sd_A": float("nan")}
    try:
        from scipy.optimize import curve_fit
        p0 = [max(A4, 1e-12), 4.0]
        pos = y > 0
        if pos.sum() >= 2:
            sl = np.polyfit(np.log(r[pos]), np.log(y[pos]), 1); p0 = [float(np.exp(sl[1])), float(-sl[0])]
        p, cov = curve_fit(lambda rr, A, eta: A * rr ** -eta, r, y, p0=p0, sigma=sig, absolute_sigma=True, maxfev=20000)
        sd = np.sqrt(np.diag(cov))
        res.update(eta=float(p[1]), sd_eta=float(sd[1]), A=float(p[0]), sd_A=float(sd[0]))
    except Exception as exc:                                       # noise-only data: no power law to fit
        res["fit_error"] = repr(exc)
    return res


def _verdict(fit):
    ok = (3.5 <= fit["eta"] <= 4.5) and np.isfinite(fit["sd_A"]) and fit["A"] > 3 * fit["sd_A"]
    return bool(ok)


def _analyse_blocks(blocks, rlo=2, rhi=8):
    """blocks: list of (c0, c4, axis, var) per time block."""
    c0 = np.array([b[0] for b in blocks]); c4 = np.array([b[1] for b in blocks]); ax = np.array([b[2] for b in blocks])
    var = np.array([b[3] for b in blocks]); nb = len(blocks)
    se = lambda a: a.std(0, ddof=1) / np.sqrt(nb) if nb > 1 else np.full(a.shape[1:], np.nan)
    r = np.arange(1, c4.shape[1] + 1, dtype=float); sel = (r >= rlo) & (r <= rhi)
    m4, s4 = c4.mean(0), se(c4)
    fit = _fit_power(r[sel], m4[sel], s4[sel])
    return {"r": r.tolist(), "c4": m4.tolist(), "c4_se": s4.tolist(), "c0": c0.mean(0).tolist(),
            "axis_C": ax.mean(0).tolist(), "axis_C_se": se(ax).tolist(), "variance": float(var.mean()),
            "fit_range": [rlo, rhi], "fit": fit, "r4_power_law": _verdict(fit), "n_blocks": nb}


def theory_C(L=800):
    """exact lattice structure factor of the isotropic chiral model (Eqs. 57, 60, 67-68; mu1 = 1/2, mu2 = 1/3), solved
    self-consistently for M0, M1 by numerical Brillouin-zone integration; returns C[dy, dx] and the variance."""
    mu1, mu2, l1, l2 = 0.5, 1.0 / 3.0, CHIP1, CHIP2
    D = 0.25 * (mu1 * (l1 + l2) + mu1 ** 2 * l1 * l2); Dodd = 0.5 * mu1 ** 2 * l1 * l2
    A1 = 0.5 * mu2 * (l1 ** 2 + l2 ** 2 + mu2 * l1 ** 2 * l2 ** 2); A2 = mu1 * mu2 * l1 * l2 ** 2
    A3 = -0.25 * mu1 ** 2 * l1 * l2; A4 = -0.25 * mu1 * mu2 * l1 ** 2 * l2
    q = 2 * np.pi * np.fft.fftfreq(L); QX, QY = np.meshgrid(q, q, indexing="ij")
    lx = 2 * (1 - np.cos(QX)); ly = 2 * (1 - np.cos(QY)); sl = lx + ly; sl2 = lx ** 2 + ly ** 2
    den = 2 * D * sl - Dodd * lx * ly
    with np.errstate(invalid="ignore", divide="ignore"):
        f0 = sl / den; f2 = (sl - sl2 / sl) * sl / den
    f0[0, 0] = 0.0; f2[0, 0] = 0.0
    I0, Ix = (lambda f: f.mean()), (lambda f: (np.cos(QX) * f).mean())
    M = np.array([[1 - A1 * I0(f0) - A4 * I0(f2), -A2 * I0(f0) - A3 * I0(f2)],
                  [-A1 * Ix(f0) - A4 * Ix(f2), 1 - A2 * Ix(f0) - A3 * Ix(f2)]])
    M0, M1 = np.linalg.solve(M, [1.0, 1.0])
    g0 = A1 * M0 + A2 * M1; g2 = A3 * M1 + A4 * M0
    S = g0 * f0 + g2 * f2; S[0, 0] = 0.0
    return np.real(np.fft.ifft2(S)), float(M0 - 1.0)


# ---------------------------------------------------------------------------------------------------- fidelity
def _measure_run(W, chiral, seed, burn, nmeas, nblocks=20, deadline=None):
    lat = _Lattice(W, [np.random.default_rng([seed, int(chiral), 7])], chiral=chiral)
    rings = _rings(W); t0 = time.time()
    for _ in range(burn): lat.sweep()
    bsz = max(1, nmeas // nblocks); blocks = []; P = np.zeros((W, W // 2 + 1)); n = 0; done = 0
    while done < nmeas:
        lat.sweep(); done += 1
        f = lat.field()[0]; F = np.fft.rfft2(f - f.mean()); P += F.real ** 2 + F.imag ** 2; n += 1
        if n == bsz:
            blocks.append(_ring_stats(_corr_from_power(P / n, W), rings)); P[:] = 0.0; n = 0
        if deadline is not None and time.time() > deadline and len(blocks) >= 4: break
    return blocks, done, time.time() - t0


def fidelity(budget_s=7200):
    t_start = time.time()
    # cost calibration (a few sweeps per version and size), then the largest lattice giving >= 2000 measured sweeps
    plan = None; burn = 1500; per_version = 0.45 * budget_s; calib = {}
    for W in (800, 512, 400, 256):
        cost = 0.0
        for chiral in (True, False):
            lat = _Lattice(W, [np.random.default_rng(1)], chiral=chiral)
            lat.sweep(); t0 = time.time()
            for _ in range(3):
                lat.sweep(); np.fft.rfft2(lat.field()[0])
            cost = max(cost, (time.time() - t0) / 3)
        calib[W] = cost
        nmeas = int((per_version - (time.time() - t_start) / 2) / cost) - burn
        if nmeas >= 2000:
            plan = (W, min(nmeas, 200000)); break
    if plan is None:
        W = 256; nmeas = max(200, int(per_version / calib[256]) - burn); plan = (W, nmeas)
    W, nmeas = plan
    out = {"setting": dict(HEADLINE, lattice=W, burn_in_sweeps=burn, measured_sweeps_planned=nmeas),
           "cost_s_per_sweep": calib, "paper": {"eta": 4.0, "C_axis_x4_amplitude": 0.11,
                                                "source": "Eq. 66 and Fig. 6 (800 x 800, lambda_1 = 0.4, lambda_2 = 0.6)"},
           "tolerance": "chiral: 3.5 <= eta <= 4.5 and A > 3 sd_A (fit c4(r) = A r^-eta, r = 2..8); "
                        "achiral (rule 2 off): that conjunction fails"}
    res = {}
    for name, chiral in (("chiral", True), ("achiral", False)):
        deadline = t_start + (0.48 if chiral else 0.96) * budget_s
        blocks, done, secs = _measure_run(W, chiral, 20261103, burn, nmeas, deadline=deadline)
        a = _analyse_blocks(blocks); a.update(measured_sweeps=done, seconds=secs); res[name] = a
    Ct, var_t = theory_C(W)
    th = _ring_stats(Ct, _rings(W))
    rr = np.arange(2, 9, dtype=float); sig = np.asarray(res["chiral"]["c4_se"])[1:8]
    out["theory_exact_lattice"] = {"variance": var_t, "c4": th[1].tolist(), "axis_C": th[2].tolist(),
                                   "fit_same_weights": _fit_power(rr, th[1][1:8], sig)}
    out.update(res)
    out["pass"] = bool(res["chiral"]["r4_power_law"] and not res["achiral"]["r4_power_law"])
    out["seconds"] = time.time() - t_start
    return out


# ---------------------------------------------------------------------------------------------------- presence
def presence(out, t_min=200):
    """the fidelity measure on census-scale output of the real program: ring cos(4 theta) amplitude c4(r), fitted on
    r = 2..8 over frames with t >= t_min; error blocks = seed x 4 time blocks."""
    F = np.asarray(out["F"], dtype=float)[:, :, 0]
    S, nf, W, _ = F.shape; rec = out["meta"]["rec"]["F"]
    idx = np.arange(nf)[np.arange(nf) * rec >= t_min]
    if idx.size < 4 or out["meta"].get("unstable"):
        return {"present": False, "note": "too few frames after burn-in or run unstable", "frames_used": int(idx.size)}
    rings = _rings(W); blocks = []
    for s in range(S):
        for part in np.array_split(idx, 4):
            P = np.zeros((W, W // 2 + 1))
            for i in part:
                f = F[s, i]; Fk = np.fft.rfft2(f - f.mean()); P += Fk.real ** 2 + Fk.imag ** 2
            blocks.append(_ring_stats(_corr_from_power(P / len(part), W), rings))
    a = _analyse_blocks(blocks)
    a.update(present=a["r4_power_law"], frames_used=int(idx.size * S), t_min=t_min,
             tolerance="3.5 <= eta <= 4.5 and A > 3 sd_A (decision); the rest is reported only")
    # reported only: amplitude with eta fixed at 4 (paper: C(x,0) ~ 0.11 / x^4, Eq. 66 / Fig. 6) and the measured c4(r)
    # against the exact lattice theory at the same size (theory_C(W)), r = 2..8
    a["A_eta4"] = {"value": a["fit"]["A_eta4"], "sd": a["fit"]["sd_A_eta4"], "paper": 0.11}
    Ct, var_t = theory_C(W)
    th = _ring_stats(Ct, rings)
    r = np.arange(2, 9); m = np.asarray(a["c4"])[1:8]; se = np.asarray(a["c4_se"])[1:8]; t4 = th[1][1:8]
    w = 1.0 / np.maximum(se, 1e-300) ** 2
    scale = float(np.sum(w * m * t4) / np.sum(w * t4 * t4)); sd_scale = float(np.sqrt(1.0 / np.sum(w * t4 * t4)))
    a["vs_theory"] = {"W": W, "r": r.tolist(), "theory_c4": t4.tolist(), "theory_variance": var_t,
                      "measured_over_theory_scale": scale, "sd_scale": sd_scale,
                      "chi2_per_dof_unscaled": float(np.sum(w * (m - t4) ** 2) / r.size)}
    return a

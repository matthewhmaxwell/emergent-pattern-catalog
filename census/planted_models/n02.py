"""Planted model N02: glycolysis-type (Selkov) reaction-diffusion with mutual chemotaxis (plug-in, PLANTED_PIPELINE_PLAN.md).

Source: arXiv 2609.01159 v2, Sec. II (Eqs. 4-11), Sec. III (Eqs. 14-16, 19-22), Sec. IV (DNS; amplitude measure Eq. 82),
Sec. IV.3.1 and Figs. 17-18 (c1 = c = -c2 > 0 at c = 0.01, 0.17, 0.49), Appendix A (Eqs. 112-118, numerics).

Model as stated in the paper (Eqs. 10-11 = 112-113), fields u = phi (F0) and v = psi (F1), 2D, periodic:
    du/dt = -u + a v + u^2 v + D1 lap(u) + c1 div(u grad v)
    dv/dt =  b - a v - u^2 v + D2 lap(v) + c2 div(v grad u)
Homogeneous fixed point (Eq. 16): u* = b, v* = b / (a + b^2). Headline (plan row; the paper's main DNS point, Sec. IV
(i) and Figs. 17-18): a = 0.015, b = 10, D1 = 0.01, D2 = 1000 D1 = 10, c1 = c, c2 = -c, c = 0.17. Amplitude measure
(Eq. 82): Delta_phi = sqrt( mean_xy (u - mean u)^2 ) in the steady state (Eq. 83 is the same for v; its printed "phi"
inside the sum is a typo).

Plan rules (PLANTED_PIPELINE_PLAN.md row N02): 1 decay of u (-u); 2 feed and decay of v (b - a v); 3 conversion
a v -> u, i.e. the +a v term of the u equation (I); 4 autocatalysis u^2 v, in both equations (I); 5 u drifts along
grad v, c1 div(u grad v) (I); 6 v drifts along grad u, c2 div(v grad u) (I). Diffusion is header (always on).
Knock-out: rules 3-6 removed (du/dt = -u + D1 lap u, dv/dt = b - a v + D2 lap v; relaxes to u = 0, v = b/a).

Choices the paper does not state, or that differ from it (and why)
------------------------------------------------------------------
- Grid spacing: Appendix A integrates on a square of side 2 pi ("L_x = L_y = 2 pi, discretized uniformly on an
  N_x x N_y grid") and Sec. IV uses "L = 128", read as N = 128 points: dx = 2 pi / 128. This is the only reading under
  which the pattern is resolved (the most unstable wavelength at the headline is ~0.9, i.e. ~19 grid points; with
  dx = 1 it would be below one site). Census: 64 x 64 sites at dx = 2 pi / 128 (domain pi, about 3 pattern
  wavelengths); area4: 128 x 128 = the paper's domain.
- Time integration: the paper uses a pseudo-spectral semi-implicit Euler scheme with dt ~ 1e-5 (Appendix A), far too
  slow for 1,000 time units. Here: the same pseudo-spectral derivatives (FFT, no dealiasing, as in the paper), with a
  linearly implicit Euler step  X += dt * (I - dt J_loc)^-1 (I - dt J_T(k))^-1 F(X),  where F is the full right-hand
  side, J_loc the exact local reaction Jacobian at every site (2 x 2) and J_T(k) the diffusion plus cross-diffusion
  operator linearised at the current spatial means (2 x 2 per wavevector). Steady states are exact steady states of
  the discretised PDE (the step is proportional to F); the factor order (spectral first) keeps the amplification
  matrix near-normal (checked: ||G|| <= 1.15 at dt = 0.01 on the homogeneous state); its cross-diffusion terms use
  the same Nyquist-zeroed first derivatives as the right-hand side (no spurious checkerboard mode). dt = 0.01 (100 sub-steps per
  census step); if a census step produces a non-finite field, |value| > 1e6 or a singular local factor, it is redone
  from its start with twice the sub-steps (up to 1,600); only a failure at 1,600 sub-steps counts as blow-up. Checked
  here: the steady amplitude at the headline is 12.92 at dt = 0.01, 0.005 and 0.0025 (64 x 64, t = 300).
- The paper's integrator also pins the k = 0 mode of both fields to the fixed point ("mean values of these
  fluctuations are constrained to vanish at all times"). That constraint is not part of the model equations and is not
  used: with it, the patterned states at the headline blow up in finite time (tested), and for the knock-out and
  sub-programs the pinned value would be the full model's fixed point, which they do not have. Without it, mean(u)
  relaxes to b exactly in any steady state (d<u + v>/dt = b - <u>), as the paper's monitoring of the means implies.
- Initial condition (Appendix A: "small random perturbations around the uniform state ... in [-eps, eps], eps << 1",
  mean subtracted): u = b + e1, v = v* + e2 with e1, e2 ~ U(-0.01, 0.01) independently per site and seed, means
  subtracted. The same start is used for the knock-out and sub-programs.
- One census step = one unit of the paper's dimensionless time (plan: continuous-time model); the paper's steady
  snapshots are at t = 290, so 1,000 census steps reach the steady pattern.
- Blow-up rule as in the census interpreter: non-finite or |value| > 1e6 sets meta unstable = True and stops all
  seeds; frames up to that step are kept, the blown-up state included when t is a multiple of 10 (as census.sim.run).

Fidelity check (plan: pattern amplitude larger at c = 0.17 than at c = 0.01 and absent at 0.49), fixed before any
large run: the paper's domain and grid (128 x 128, dx = 2 pi / 128), a = 0.015, b = 10, D1 = 0.01, D2 = 10,
c1 = c = -c2 with c in {0.01, 0.17, 0.49}, 3 seeds each, run to t = 500 (shortened only if the timing calibration
says it would not fit the budget; never below t = 300), measure Delta_phi (Eq. 82) averaged over the frames with
t >= 0.8 T, seed means compared. Pass: mean Delta_phi(0.17) > mean Delta_phi(0.01) AND Delta_phi(0.49) < 0.01 x
Delta_phi(0.17) in every seed (uniform within 1 % of the headline amplitude) and no blow-up at 0.49. Note found
while implementing (reported, plan unchanged): the paper's own linear-stability threshold (Eq. 15; Det(k) of Eq. 22
= 100 - (9 + 30 c) k^2 + (0.1 + c^2) k^4 here) puts c = 0.49 inside the unstable region (min Det = -313 at k^2 = 35),
and small-scale runs of Eqs. 10-11 here give a strong pattern at c = 0.49 (Delta_phi ~ 24 at 64 x 64); the third leg
of the check is therefore expected to fail. In general (independent review): with c1 = c = -c2 the paper's Eq. 22
gives Det(k) = 100 - (9 + 30c) k^2 + (0.1 + c^2) k^4, which dips below zero iff (9 + 30c)^2 > 400 (0.1 + c^2), i.e.
41 + 540c + 500c^2 > 0, true for every c > 0: on this branch the uniform state is unstable for all c > 0, so no
re-entrant uniform state at 0.49 follows from the stated equations. The authors' Appendix A scheme (Eqs. 115-117:
"linearize ... the obtained nonlinear PDEs (115) and (116) are used", with typos such as -2ub/(b^2+a) and a xi_2
term acting on lap v) appears to be a linearised integration, which may explain the difference. Small-scale values
here (64 x 64, t = 300): Delta_phi = 6.3 (c = 0.01), 12.9 (c = 0.17), 23-24 (c = 0.49).
"""
import time

import numpy as np
import scipy.fft as sfft

from . import base

CID = "planted:N02"
LAYERS = "F"
K = {}
SPEED = 0.0
RULES = [("U_DECAY", False),                  # 1  -u
         ("V_FEED_AND_DECAY", False),         # 2  b - a v
         ("CONVERSION_V_TO_U", True),         # 3  +a v in the u equation (I)
         ("AUTOCATALYSIS_U2V", True),         # 4  +u^2 v / -u^2 v (I)
         ("U_DRIFT_ALONG_GRAD_V", True),      # 5  c1 div(u grad v) (I)
         ("V_DRIFT_ALONG_GRAD_U", True)]      # 6  c2 div(v grad u) (I)
HEADLINE = {
    "a": 0.015, "b": 10.0,                    # Sec. IV point (i); Figs. 17-18
    "D1": 0.01, "D2": 10.0,                   # D1 = 0.01, D = D2/D1 = 1000 (Figs. 17-18 top row)
    "c": 0.17, "c1": 0.17, "c2": -0.17,       # xi_1 = xi = -xi_2, xi = 0.17 (Fig. 17 middle column)
    "dx": 2 * np.pi / 128,                    # Appendix A: side 2 pi, Sec. IV: 128 points
    "init_noise": 0.01,                       # Appendix A: U(-eps, eps), eps << 1 (value chosen here)
    "dt": 0.01,                               # sub-step (adaptive halving on failure), see docstring
    "W_census": 64, "W_area4": 128,
}
STEP = ("one census step = 1 unit of the paper's dimensionless time; 100 linearly implicit pseudo-spectral sub-steps "
        "(dt = 0.01, halved adaptively on failure down to 1/1600)")
RULE_TEXT = ("World: two continuous concentration fields u and v on a 64x64 square grid with wrap-around edges, started "
             "at their uniform steady state (u=b, v=b/(a+b^2)) plus small random noise. Rules, integrated forward in "
             "time: u changes at rate -u + a*v + u^2*v + D1*(Laplacian of u) + c1*div(u grad v); v changes at rate "
             "b - a*v - u^2*v + D2*(Laplacian of v) + c2*div(v grad u). The last terms make each species drift along "
             "the other's gradient, with strength and sign set by c1 and c2. Parameters: a=0.015, b=10, D1=0.01, "
             "D2=1000*D1. Setting: c1=c and c2=-c (opposite signs, so one species is drawn toward the other while the "
             "other is pushed away), c=0.17.")

NSUB0, NSUB_MAX = 100, 1600


class _Field:
    """S seeds of the two-field model on a W x W periodic grid, spacing dx; on = switched-on rule numbers."""

    def __init__(self, W, gens, on, a=0.015, b=10.0, D1=0.01, D2=10.0, c1=0.17, c2=-0.17, dx=2 * np.pi / 128,
                 eps=0.01):
        self.W, self.S = W, len(gens)
        self.a, self.b, self.D1, self.D2 = a, b, D1, D2
        self.r3, self.r4 = float(3 in on), float(4 in on)
        self.c1 = c1 if 5 in on else 0.0
        self.c2 = c2 if 6 in on else 0.0
        k = 2 * np.pi * np.fft.fftfreq(W, d=dx); kr = 2 * np.pi * np.fft.rfftfreq(W, d=dx)
        KY, KX = np.meshgrid(k, kr, indexing="ij")
        self.K2 = KX ** 2 + KY ** 2
        self.iKX = 1j * np.where(np.arange(kr.size)[None, :] == W // 2, 0.0, KX)    # Nyquist zeroed for 1st derivs
        self.iKY = 1j * np.where(np.arange(W)[:, None] == W // 2, 0.0, KY)
        self.Kc = self.iKX.imag ** 2 + self.iKY.imag ** 2      # -(div grad) with the same Nyquist-zeroed 1st derivs
        us, vs = b, b / (a + b * b)
        u = np.empty((self.S, W, W)); v = np.empty((self.S, W, W))
        for s, g in enumerate(gens):
            e1 = g.uniform(-eps, eps, (W, W)); e2 = g.uniform(-eps, eps, (W, W))
            u[s] = us + e1 - e1.mean(); v[s] = vs + e2 - e2.mean()
        self.u, self.v = u, v

    # full right-hand side, returned in Fourier space
    def rhs_hat(self, u, v):
        f2, i2, W = sfft.rfft2, sfft.irfft2, self.W
        uh, vh = f2(u), f2(v)
        Fu = -self.K2 * self.D1 * uh; Fv = -self.K2 * self.D2 * vh
        if self.c1:
            gx, gy = i2(self.iKX * vh, s=(W, W)), i2(self.iKY * vh, s=(W, W))
            Fu = Fu + self.c1 * (self.iKX * f2(u * gx) + self.iKY * f2(u * gy))
        if self.c2:
            gx, gy = i2(self.iKX * uh, s=(W, W)), i2(self.iKY * uh, s=(W, W))
            Fv = Fv + self.c2 * (self.iKX * f2(v * gx) + self.iKY * f2(v * gy))
        u2v = self.r4 * u * u * v
        Ru = -u + self.r3 * self.a * v + u2v
        Rv = self.b - self.a * v - u2v
        return Fu + f2(Ru), Fv + f2(Rv)

    def substep(self, dt):
        u, v, W = self.u, self.v, self.W
        Fu, Fv = self.rhs_hat(u, v)
        ub = u.mean((1, 2))[:, None, None]; vb = v.mean((1, 2))[:, None, None]
        n11 = 1 + dt * self.D1 * self.K2; n22 = 1 + dt * self.D2 * self.K2
        n12 = dt * self.c1 * ub * self.Kc; n21 = dt * self.c2 * vb * self.Kc
        dn = n11 * n22 - n12 * n21
        Gu = sfft.irfft2((n22 * Fu - n12 * Fv) / dn, s=(W, W)); Gv = sfft.irfft2((-n21 * Fu + n11 * Fv) / dn, s=(W, W))
        uv = u * v
        m11 = 1 - dt * (-1 + 2 * self.r4 * uv); m12 = -dt * (self.r3 * self.a + self.r4 * u * u)
        m21 = dt * 2 * self.r4 * uv; m22 = 1 + dt * (self.a + self.r4 * u * u)
        det = m11 * m22 - m12 * m21
        if not (np.all(np.isfinite(dn)) and np.all(det > 0)): return False
        self.u = u + dt * (m22 * Gu - m12 * Gv) / det
        self.v = v + dt * (-m21 * Gu + m11 * Gv) / det
        return True

    def ok(self):
        return bool(np.all(np.isfinite(self.u)) and np.all(np.isfinite(self.v))
                    and max(np.abs(self.u).max(), np.abs(self.v).max()) <= 1e6)

    def advance(self, nsub):
        """one census step (one time unit); returns (success, nsub used). On final failure (at NSUB_MAX) the failed
        state is kept, so the caller records the blown-up field as census.sim.run does."""
        u0, v0 = self.u.copy(), self.v.copy()
        while True:
            dt = 1.0 / nsub; good = True
            for _ in range(nsub):
                if not self.substep(dt): good = False; break
            if good and self.ok(): return True, nsub
            if nsub * 2 > NSUB_MAX: return False, nsub
            self.u, self.v = u0.copy(), v0.copy(); nsub *= 2


def _run(W, gens, on, steps, rec=None, deadline=None, **par):
    """as census.sim.run: frame at t = 0 and every rec steps; on blow-up the step's (blown-up) state is recorded if
    t % rec == 0, unstable is set and the run stops. deadline (fidelity only): stop early, returned as truncated."""
    fld = _Field(W, gens, on, **par)
    rec = rec or base.REC["F"]
    frames = [np.stack([fld.u, fld.v], 1).astype(np.float32)]; nsub = NSUB0; t = 0; unstable = False
    truncated = False
    while t < steps:
        good, nsub = fld.advance(nsub); t += 1
        if t % rec == 0: frames.append(np.stack([fld.u, fld.v], 1).astype(np.float32))
        if not good: unstable = True; break
        if deadline is not None and time.time() > deadline and t >= 2 * rec: truncated = True; break
    F = np.stack(frames, 1)
    return (F, t, unstable, nsub, truncated) if deadline is not None else (F, t, unstable, nsub)


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    on = base.kept(keep, RULES)
    W = base.SIZES[size]["W"]
    F, t, unstable, _ = _run(W, base.rngs(seed0, S), on, steps)
    return {"F": F, "meta": base.meta(t, W, unstable=unstable)}


# ---------------------------------------------------------------------------------------------------- measure
def amplitude(F, frac=0.8):
    """Eq. 82/83 on frames (S, frames, 2, W, W) with frame index >= frac * last; returns per-seed Delta_phi, Delta_psi
    (time-averaged over those frames) and a stationarity ratio (rms change between the last two frames / Delta_phi)."""
    F = np.asarray(F, dtype=float); nf = F.shape[1]
    sel = np.arange(nf)[np.arange(nf) >= frac * (nf - 1)]
    dphi = F[:, sel, 0].std(axis=(2, 3)).mean(1); dpsi = F[:, sel, 1].std(axis=(2, 3)).mean(1)
    chg = np.sqrt(((F[:, -1, 0] - F[:, -2, 0]) ** 2).mean((1, 2))) / np.maximum(dphi, 1e-300) if nf > 1 else np.nan
    return dphi, dpsi, chg, F[:, sel, 0].mean(axis=(1, 2, 3))


def fidelity(budget_s=7200, _small=False):
    """_small=True only exercises the code path (64 x 64, 1 seed, t = 40); it is not the fidelity check."""
    t_start = time.time(); cs = (0.01, 0.17, 0.49)
    W, S, T, Tmin = (64, 1, 40, 40) if _small else (128, 3, 500, 300)
    # calibration: 2 time units at the paper's size with the headline setting
    fld = _Field(W, base.rngs(1, S), {3, 4, 5, 6}); t0 = time.time()
    for _ in range(2): fld.advance(NSUB0)
    unit = (time.time() - t0) / 2
    est = lambda T: len(cs) * T * unit * 1.6               # 1.6: allowance for sub-step refinement at c = 0.49
    while T > Tmin and est(T) > 0.85 * budget_s - (time.time() - t_start): T -= 50
    if est(T) > 0.85 * budget_s - (time.time() - t_start) and S > 1:
        S = 1                                              # fallback: one seed, with its own cost re-estimate
        fld = _Field(W, base.rngs(1, S), {3, 4, 5, 6}); t0 = time.time()
        for _ in range(2): fld.advance(NSUB0)
        unit = (time.time() - t0) / 2
        T = 500 if not _small else T
        while T > Tmin and est(T) > 0.85 * budget_s - (time.time() - t_start): T -= 50
    deadline = t_start + 0.95 * budget_s                    # hard stop: a run past it is cut and flagged truncated
    out = {"setting": dict(HEADLINE, W=W, domain="2 pi", seeds=S, T=T, cs=list(cs)),
           "cost_s_per_time_unit": unit,
           "paper": {"source": "Figs. 17-18 (D1 = 0.01, D = 1000): fuzzy spots at 0.01, sharp spots at 0.17, uniform "
                               "at 0.49; amplitude Delta_phi first rises then falls with c",
                     "expected": "Delta_phi(0.17) > Delta_phi(0.01); Delta_phi(0.49) = 0"},
           "tolerance": "mean Delta_phi(0.17) > mean Delta_phi(0.01) and Delta_phi(0.49) < 0.01 * Delta_phi(0.17) in "
                        "every seed, no blow-up at 0.49",
           "linear_theory_Eq22_minDet": {}}
    a, b, D1, D2 = 0.015, 10.0, 0.01, 10.0
    for c in cs:
        k2 = np.linspace(0, 400, 40001)                        # Eq. 19 Jacobian, c1 = c, c2 = -c
        J11 = (b * b - a) / (b * b + a) - D1 * k2; J22 = -(a + b * b) - D2 * k2
        J12 = (a + b * b) - c * b * k2; J21 = -2 * b * b / (a + b * b) + c * b / (a + b * b) * k2
        det = J11 * J22 - J12 * J21
        out["linear_theory_Eq22_minDet"][str(c)] = {"min_Det": float(det.min()), "k2_at_min": float(k2[det.argmin()])}
    res = {}
    for c in cs:
        t0 = time.time()
        F, t, unstable, nsub, trunc = _run(W, base.rngs(20261103, S), {3, 4, 5, 6}, T, deadline=deadline,
                                           c1=c, c2=-c)
        dphi, dpsi, chg, um = amplitude(F) if not unstable else (np.full(S, np.nan),) * 4
        res[str(c)] = {"Delta_phi": np.atleast_1d(dphi).tolist(), "Delta_psi": np.atleast_1d(dpsi).tolist(),
                       "Delta_phi_mean": float(np.nanmean(dphi)), "stationarity": np.atleast_1d(chg).tolist(),
                       "mean_u": np.atleast_1d(um).tolist(), "unstable": unstable, "steps_run": t,
                       "final_substeps": nsub, "truncated_by_deadline": trunc, "seconds": time.time() - t0}
    out["results"] = res
    lo, hd, hi = res["0.01"], res["0.17"], res["0.49"]
    leg1 = bool(hd["Delta_phi_mean"] > lo["Delta_phi_mean"])
    leg2 = bool((not hi["unstable"]) and np.all(np.asarray(hi["Delta_phi"]) < 0.01 * hd["Delta_phi_mean"]))
    out.update(leg_larger_at_017=leg1, leg_absent_at_049=leg2, passed=bool(leg1 and leg2),
               seconds=time.time() - t_start)
    out["pass"] = out.pop("passed")
    return out


def presence(out):
    """Eq. 82 amplitude on census-scale output of the real program (frames with t >= 800); pattern present if
    Delta_phi > 1 % of mean u in a seed and the pattern is steady (rms change over the last 10 steps < 5 % of it)."""
    F = np.asarray(out["F"])
    if out["meta"].get("unstable"):
        return {"present": False, "note": "run unstable (blow-up)", "steps_run": out["meta"]["steps_run"]}
    dphi, dpsi, chg, um = amplitude(F)
    pat = dphi > 0.01 * np.abs(um); steady = np.asarray(chg) < 0.05
    return {"Delta_phi": dphi.tolist(), "Delta_psi": dpsi.tolist(), "mean_u": um.tolist(),
            "stationarity": np.atleast_1d(chg).tolist(), "pattern": pat.tolist(), "steady": steady.tolist(),
            "present": bool(np.sum(pat & steady) >= 2),
            "tolerance": "Delta_phi > 0.01 mean(u) and last-frame rms change < 0.05 Delta_phi, in >= 2 of 3 seeds"}

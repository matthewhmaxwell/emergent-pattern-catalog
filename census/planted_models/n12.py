"""N12 plug-in: 2D XY model with 1/r^4 couplings between all pairs (arXiv 2609.39267, v1).

Paper model. Sec. IV, eqs. (27)-(29) define H = H_SR + H_LR on an L x L square lattice with periodic boundaries;
the log-decay results the plan checks (Sec. IV.3, Fig. 4, Tables 2-3) are for the "uniform model with J = b J_LR",
i.e. the homogeneous long-range XY model of the paper's Ref. [42] (Yao, Xiao, Zhang, Deng, Fan, PRB 112, 144429;
arXiv 2411.01811, Sec. II.1, eqs. 1-3):
    H = - sum_{i<j} J_ij cos(theta_i - theta_j),   J_ij = c(sigma, N) / |r_i - r_j|^(2 + sigma),   sigma = 2,
with |r_i - r_j| the minimum-image distance on the torus (each pair coupled once, dx, dy <= L/2) and c fixed so
that sum_j J_ij = 4 on the finite lattice (total coupling per site fixed; it reduces to the nearest-neighbour XY
model as sigma -> infinity). In the paper's notation (eq. 29, b fixed by 4 = sum_{r>1} b/r^4, "the total coupling
per site is fixed") J = b J_LR makes every pair's coupling proportional to 1/r^4 with the nearest-neighbour value
included; the temperatures of Fig. 3-4 (beta_c = 0.7315, beta = 1, 2, 4, 8) are those of Ref. [42], whose total
coupling per site is 4, so that normalisation is used here. In eqs. (27)-(29) notation this module is
J_LR = 1/(1 + b), J = b J_LR = b/(1 + b) (total 4 J + 4 J_LR = 4), with b computed on each finite lattice and beta_c
taken from Ref. [42] (cross-checked: the simulated M^2 at L = 128 and 256 matches the paper's Table 3 fit to 0.2%). Inverse temperature beta enters as beta * H (reduced Hamiltonian).
beta_c(sigma = 2) = 0.7315(2) (eq. 34); headline beta = 1 (plan) lies in the low-temperature phase, where the paper
reports M^2 = A0 [ln(L/L0)]^(-eta_l) with eta_l = 0.338(3), A0 = 0.972(8), L0 = 0.21(1) (Table 3, L_min = 128,
sizes to 8192) and eta_l = 0.34(1), kappa_R = 0.48(1) (Table 2).

Dynamics (choice; the paper samples with a Luijten-Bloete cluster algorithm, which has no physical time). The census
needs a time evolution, so one census step is one sweep of single-spin heat-bath (Gibbs) updates: every site once,
in a fresh random order per sweep, each new angle drawn from its exact conditional distribution
    p(theta_i) ~ exp(beta |h_i| cos(theta_i - arg h_i)),   h_i = sum_j J_ij exp(i theta_j)   (von Mises),
which satisfies detailed balance for beta*H, so the Boltzmann distribution is stationary (the fidelity check below
confirms it numerically against the cluster sampler). The long-range local field is kept exact: it is recomputed
by FFT convolution with the minimum-image kernel at the start of each sweep, sites are updated in blocks of B along
the random order with the within-block corrections added from the B x B coupling block, and the field is brought
up to date by one FFT convolution of the block's changes. Every update therefore sees the current state of all
other spins (strictly sequential dynamics, no parallel-update approximation).

Census mapping and choices:
  * lattice W = 64 (area4: 128), lattice spacing 1 as in the paper; periodic, minimum image, c renormalised for
    each L so that sum_j J_ij = 4 exactly;
  * initial condition: angles i.i.d. uniform on [0, 2pi) (census phase convention; the paper does not state one);
    from such a hot start the local dynamics orders by coarsening and can keep long-lived twisted (winding) states
    for a while, which is part of the census observation;
  * output: "Cph" angles in [0, 2pi) as float32, frames every 5 sweeps; "C" int8 zeros (one cell type, k = 1),
    frames every 2 sweeps;
  * knock-out (plan): J = 0; the heat-bath then draws every angle uniformly at random each sweep (thermal noise
    only), which is implemented directly.

Fidelity check (plan: squared magnetisation decreasing slowly with L at beta = 1 as the paper reports). Equilibrium
measurement with an exact Wolff single-cluster algorithm for the long-range couplings (the paper's algorithm class):
bond (i, j) is active with probability 1 - exp(-2 beta J_ij a_i a_j)^+ where a = S . r for a random axis r; the cluster
is grown generation by generation, the join probability of every outside site being computed at once from an FFT
convolution of the current frontier with the coupling kernel (exact, since bonds are independent). Sizes 16, 32, 64,
128, 256, ... up to the largest that fits the budget; M^2 = <|sum_i S_i / N|^2> after each cluster flip (cold start,
first 20% discarded, binning errors), and (xi/L)^2 from the second-moment correlation length (eq. 31) as context.
Census dynamics check: the heat-bath sweeps used by simulate() are run at L = 16 and 32 (cold start) and their
<M^2> compared with the cluster sampler at the same L.
Tolerance (fixed before any large run), pass = all of:
  (a) M^2 decreases with L: no consecutive increase beyond 2 combined standard errors, and M^2(Lmax) < M^2(Lmin);
  (b) for every L >= 128 run: |M^2 - 0.972 [ln(L/0.21)]^-0.338| <= 0.01 (Table 3 fit, L_min = 128);
  (c) slow (logarithmic) decay: a weighted fit of ln M^2 = ln A0 - eta_l ln ln(L/0.21) over the sizes L >= 128 (the
      paper's L_min) gives |eta_l - 0.34| <= max(0.07, 3 sigma_fit), sigma_fit the fit's standard error (from the
      M^2 errors, inflated by sqrt(chi^2/dof) when that exceeds 1). If fewer than three sizes >= 128 fit the
      budget, (c) is reported as "not decisive" and (a), (b), (d) decide;
  (d) heat-bath <M^2> at L = 16 and 32 within 4 combined standard errors of the cluster sampler.

presence(): the same measure on census output: M^2 per frame from Cph (and its late-time mean), (xi/L)^2 and the
vortex density, next to the paper-fit value at L = 64 (an extrapolation below the fit range).
"""
import cmath
import time

import numpy as np
import scipy.fft as sfft

from . import base

CID = "planted:N12"
LAYERS = "C"
K = {"C": 1}
SPEED = 0.0
RULES = [("LONG_RANGE_COUPLING_R4", True),   # 1: all pairs, J(r) = c/r^4, sum_j J_ij = 4 (I)
         ("THERMAL_NOISE", False)]           # 2: heat-bath noise at inverse temperature beta (one-way)
HEADLINE = {
    "beta": 1.0,            # plan headline; low-temperature phase (beta_c = 0.7315, eq. 34; Fig. 4, Tables 2-3)
    "sigma": 2.0,           # J(r) ~ 1/r^(2+sigma) = 1/r^4
    "coupling_per_site": 4.0,   # sum_j J_ij = 4 (Ref. [42] eq. 3; paper Sec. IV normalisation)
    "W": 64,                # census lattice side (area4: 128), minimum-image periodic distances
}
STEP = ("one heat-bath sweep: every site updated once in random order, new angle drawn from its exact conditional "
        "(von Mises) distribution given the full long-range field (exact sequential, FFT-maintained field)")
RULE_TEXT = ("World: a 64 x 64 square lattice with wrap-around edges, each site holding a planar unit spin (an "
             "angle). Energy: minus the sum over all pairs of J(r) times the cosine of their angle difference, where "
             "every pair interacts with J(r) proportional to 1/r^4 (r = distance), normalised so that the total "
             "coupling per site is fixed. Rules: equilibrium Monte Carlo with single-spin updates (one sweep of all "
             "sites per step) at inverse temperature beta = 1, from random angles, for 1000 sweeps. Measure the "
             "squared magnetisation and the correlation length.")

TWO_PI = 2.0 * np.pi
PAPER_FIT = {"A0": 0.972, "L0": 0.21, "eta_l": 0.338}       # Table 3, beta = 1, L_min = 128
PAPER_ETA = 0.34                                            # Table 2, beta = 1
PAPER_KAPPA_R = 0.48                                        # Table 2 / 4, beta = 1


def paper_M2(L):
    return PAPER_FIT["A0"] * np.log(L / PAPER_FIT["L0"]) ** (-PAPER_FIT["eta_l"])


# ------------------------------------------------------------------------------------------------ couplings / field
def kernel(L, sigma=2.0, total=4.0):
    """J[dy, dx] for displacement (dy, dx) mod L, minimum image, J[0, 0] = 0, normalised to sum = total."""
    d = np.arange(L); d = np.minimum(d, L - d).astype(float)
    r2 = d[:, None] ** 2 + d[None, :] ** 2
    J = np.zeros((L, L)); nz = r2 > 0
    J[nz] = r2[nz] ** (-(2.0 + sigma) / 2.0)
    return J * (total / J.sum())


def _field(z, Jhat):
    return sfft.ifft2(sfft.fft2(z) * Jhat)


def _sweep(th, z, J, Jhat, beta, rng, B):
    """one exact sequential heat-bath sweep in random order (th, z modified in place)."""
    L = th.shape[0]; N = L * L
    thf = th.reshape(-1); zf = z.reshape(-1)
    hf = _field(z, Jhat).reshape(-1)
    perm = rng.permutation(N); py = perm // L; px = perm % L
    vm = rng.vonmises; phase = cmath.phase; cexp = cmath.exp
    D = np.zeros(N, complex)
    for b0 in range(0, N, B):
        idx = perm[b0:b0 + B]; y = py[b0:b0 + B]; x = px[b0:b0 + B]; nb = len(idx)
        Jbb = J[(y[:, None] - y[None, :]) % L, (x[:, None] - x[None, :]) % L]
        hb = hf[idx].tolist(); zold = zf[idx].tolist()
        dz = np.zeros(nb, complex); tnew = [0.0] * nb; znew = [0j] * nb
        for k in range(nb):
            hk = hb[k] + (Jbb[k, :k] @ dz[:k] if k else 0j)
            t = vm(phase(hk), beta * abs(hk))
            zn = cexp(1j * t)
            dz[k] = zn - zold[k]; tnew[k] = t; znew[k] = zn
        thf[idx] = tnew; zf[idx] = znew
        if b0 + B < N:
            D[:] = 0; D[idx] = dz
            hf += sfft.ifft2(sfft.fft2(D.reshape(L, L)) * Jhat).reshape(-1)


def _block(L):
    return 128 if L <= 96 else 256


def _run_heatbath(theta0, sweeps, beta, rng, coupled=True, rec_every=None, measure=False):
    """census dynamics from theta0; returns (theta, frames list, M2 list)."""
    L = theta0.shape[0]
    th = np.array(theta0, float); frames = []; m2 = []
    if rec_every: frames.append(np.mod(th, TWO_PI).astype(np.float32))
    if coupled:
        J = kernel(L); Jhat = sfft.fft2(J).real; z = np.exp(1j * th); B = _block(L)
    for t in range(1, sweeps + 1):
        if coupled:
            _sweep(th, z, J, Jhat, beta, rng, B)
        else:
            th = rng.uniform(0.0, TWO_PI, (L, L))
        if rec_every and t % rec_every == 0:
            frames.append(np.mod(th, TWO_PI).astype(np.float32))
        if measure:
            m2.append(abs(np.exp(1j * th).mean()) ** 2)
    return th, frames, m2


def simulate(seed0, steps=1000, S=3, size="census", keep=None):
    on = base.kept(keep, RULES)
    coupled = 1 in on
    W = base.SIZES[size]["W"]; beta = HEADLINE["beta"]
    Cph = np.empty((S, base.n_frames(steps, "Cph"), W, W), np.float32)
    for s, g in enumerate(base.rngs(seed0, S)):
        th0 = g.uniform(0.0, TWO_PI, (W, W))
        _, fr, _ = _run_heatbath(th0, steps, beta, g, coupled=coupled, rec_every=base.REC["Cph"])
        Cph[s] = np.stack(fr)
    C = np.zeros((S, base.n_frames(steps, "C"), W, W), np.int8)
    return {"C": C, "Cph": Cph, "meta": base.meta(steps, W)}


# --------------------------------------------------------------------------------------------- cluster sampler (fid)
def _wolff(th, Jhat_r, beta, rng):
    """one exact Wolff cluster flip for the long-range couplings; returns (cluster size, generations)."""
    L = th.shape[0]
    phi = rng.uniform(0.0, TWO_PI)
    a = np.cos(th - phi)
    i0 = int(rng.integers(L * L))
    if a.flat[i0] < 0: a = -a
    ap = np.maximum(a, 0.0); tb = 2.0 * beta * ap
    inc = np.zeros((L, L), bool); inc.flat[i0] = True
    front = np.zeros((L, L)); front.flat[i0] = ap.flat[i0]
    gens = 0
    while True:
        P = sfft.irfft2(sfft.rfft2(front) * Jhat_r, s=(L, L))
        prob = -np.expm1(-tb * np.maximum(P, 0.0))
        new = (rng.random((L, L)) < prob) & ~inc
        gens += 1
        if not new.any(): break
        inc |= new
        front = np.where(new, ap, 0.0)
    th[inc] = np.pi + 2.0 * phi - th[inc]
    return int(inc.sum()), gens


def _mags(th):
    """|M|^2 and |M_k|^2 (k = 2 pi / L along x and y, averaged) for XY spins."""
    L = th.shape[0]; c = np.cos(th); s = np.sin(th)
    m2 = c.mean() ** 2 + s.mean() ** 2
    e = np.exp(2j * np.pi * np.arange(L) / L)
    mk = 0.0
    for ax in (0, 1):
        cx = c.mean(axis=ax); sx = s.mean(axis=ax)          # profile along the other axis
        mk += abs((cx * e).mean()) ** 2 + abs((sx * e).mean()) ** 2
    return m2, mk / 2.0


def _binned(x, nb=20):
    x = np.asarray(x, float); n = len(x) // nb
    if n < 1: return float(np.mean(x)), float("nan")
    b = x[: n * nb].reshape(nb, n).mean(1)
    return float(x.mean()), float(b.std(ddof=1) / np.sqrt(nb))


def _wolff_run(L, n_upd, beta, rng, therm_frac=0.2, deadline=None):
    J = kernel(L); Jhat_r = sfft.rfft2(J).real
    th = np.zeros((L, L))
    n_th = max(100, int(therm_frac * n_upd))
    sizes = []; gens = []
    for _ in range(n_th): _wolff(th, Jhat_r, beta, rng)
    m2 = []; mk2 = []
    for _ in range(n_upd):
        c, g = _wolff(th, Jhat_r, beta, rng); sizes.append(c); gens.append(g)
        a, b = _mags(th); m2.append(a); mk2.append(b)
        if deadline is not None and time.time() > deadline and len(m2) >= 200: break
    n_upd = len(m2)
    M2, eM2 = _binned(m2); MK, eMK = _binned(mk2)
    xi = np.sqrt(max(M2 / MK - 1.0, 0.0)) / (2.0 * np.sin(np.pi / L))
    return {"L": L, "updates": n_upd, "therm": n_th, "M2": M2, "M2_err": eM2, "Mk2": MK, "Mk2_err": eMK,
            "xi_over_L_sq": float((xi / L) ** 2), "mean_cluster_frac": float(np.mean(sizes) / L ** 2),
            "mean_generations": float(np.mean(gens)), "paper_fit_M2": float(paper_M2(L))}


def fidelity(budget_s=7200, sizes=(16, 32, 64, 128, 256, 512, 1024, 2048), n_upd=3000, hb_sweeps=(16000, 6000)):
    t0 = time.time(); beta = HEADLINE["beta"]
    rng = np.random.default_rng([base.seed_of(CID), 39267])
    res = {"check": "M^2 versus L at beta = 1 (paper Fig. 4a, Tables 2-3), sigma = 2, homogeneous model",
           "setting": dict(HEADLINE),
           "tolerance": {"a": "M^2 decreasing with L (no rise > 2 SE; M2(Lmax) < M2(Lmin))",
                         "b": "|M2 - 0.972 ln(L/0.21)^-0.338| <= 0.01 for L >= 128",
                         "c": "|eta_l(fit, L0 = 0.21, L >= 128) - 0.34| <= max(0.07, 3 sigma_fit); "
                              "not decisive if fewer than 3 sizes >= 128",
                         "d": "heat-bath (census dynamics) M2 at L = 16, 32 within 4 SE of the cluster sampler"},
           "paper": {"eta_l": PAPER_ETA, "fit": PAPER_FIT, "kappa_R": PAPER_KAPPA_R}}
    # ---- census-dynamics check first (fixed cost, about 2-3 minutes)
    hb = {}
    for L, n in zip((16, 32), hb_sweeps):
        th0 = np.zeros((L, L)); burn = 1000
        _, _, m2 = _run_heatbath(th0, burn + n, beta, rng, coupled=True, measure=True)
        hb[L] = _binned(m2[burn:])
    # ---- cluster sampler, sizes ascending while the budget lasts
    runs = []; spent = time.time() - t0
    reserve = 0.08 * budget_s; deadline = t0 + budget_s - reserve
    for L in sizes:
        tp = time.time(); J = kernel(L); Jr = sfft.rfft2(J).real; th = np.zeros((L, L))
        for _ in range(10): _wolff(th, Jr, beta, rng)
        per = (time.time() - tp) / 10.0 * 1.3                     # + measurement overhead, margin
        need = per * (n_upd + max(100, 0.2 * n_upd))
        n_here = n_upd
        left = deadline - time.time()
        if need > left:
            n_here = int(left / per / 1.2)
            if n_here < 500: break
        runs.append(_wolff_run(L, n_here, beta, rng, deadline=deadline))
    res["cluster_runs"] = runs
    Ls = np.array([r["L"] for r in runs]); M2 = np.array([r["M2"] for r in runs]); E = np.array([r["M2_err"] for r in runs])
    # (a)
    rises = [(int(Ls[i]), int(Ls[i + 1])) for i in range(len(Ls) - 1) if M2[i + 1] - M2[i] > 2 * np.hypot(E[i], E[i + 1])]
    ok_a = (not rises) and len(Ls) >= 2 and M2[-1] < M2[0]
    # (b)
    big = Ls >= 128
    dev = {int(L): float(m - paper_M2(L)) for L, m in zip(Ls[big], M2[big])}
    ok_b = bool(big.any()) and all(abs(v) <= 0.01 for v in dev.values())
    # (c) weighted fit over L >= 128 only; decisive with at least three such sizes
    sel = big
    eta_fit = float("nan"); sig_fit = float("nan"); chi2dof = float("nan")
    if sel.sum() >= 2:
        X = np.log(np.log(Ls[sel] / PAPER_FIT["L0"])); Y = np.log(M2[sel])
        w = (M2[sel] / np.maximum(E[sel], 1e-12)) ** 2                  # 1 / var(ln M^2)
        A = np.vstack([np.ones_like(X), X]).T
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        coef = cov @ (A.T @ (w * Y))
        eta_fit = float(-coef[1]); dof = int(sel.sum()) - 2
        if dof > 0:
            chi2dof = float(np.sum(w * (Y - A @ coef) ** 2) / dof)
        sig_fit = float(np.sqrt(cov[1, 1] * max(1.0, chi2dof if np.isfinite(chi2dof) else 1.0)))
    c_decisive = bool(sel.sum() >= 3)
    ok_c = bool(c_decisive and abs(eta_fit - PAPER_ETA) <= max(0.07, 3 * sig_fit))
    c_status = ("pass" if ok_c else "fail") if c_decisive else "not decisive"
    # (d)
    hbcmp = {}
    ok_d = True
    for L in (16, 32):
        r = next((x for x in runs if x["L"] == L), None)
        if r is None: ok_d = False; continue
        m, e = hb[L]; z = (m - r["M2"]) / np.hypot(e, r["M2_err"])
        hbcmp[L] = {"heatbath_M2": m, "heatbath_err": e, "cluster_M2": r["M2"], "cluster_err": r["M2_err"], "z": float(z)}
        ok_d &= abs(z) <= 4
    # context: (xi/L)^2 slope vs ln L over L >= 128
    xl = np.array([r["xi_over_L_sq"] for r in runs])
    kap = float(np.polyfit(np.log(Ls[big]), xl[big], 1)[0]) if big.sum() >= 2 else float("nan")
    res.update({"sizes": Ls.tolist(), "M2": M2.tolist(), "M2_err": E.tolist(),
                "deviation_from_paper_fit": dev, "eta_l_fit": eta_fit, "eta_l_fit_se": sig_fit,
                "eta_fit_chi2_dof": chi2dof, "eta_fit_sizes": Ls[sel].tolist(), "c_status": c_status,
                "kappa_R_slope_context": kap, "heatbath_vs_cluster": hbcmp, "consecutive_rises": rises,
                "a_pass": bool(ok_a), "b_pass": bool(ok_b), "c_pass": ok_c if c_decisive else None,
                "d_pass": bool(ok_d),
                "pass": bool(ok_a and ok_b and ok_d and (ok_c or not c_decisive)), "secs": round(time.time() - t0, 1)})
    return res


# -------------------------------------------------------------------------------------------------------- presence
def _vortex_density(th):
    def w(d): return (d + np.pi) % TWO_PI - np.pi
    a = th; b = np.roll(th, -1, 1); c = np.roll(b, -1, 0); d = np.roll(th, -1, 0)
    q = w(b - a) + w(c - b) + w(d - c) + w(a - d)
    return float(np.mean(np.abs(np.rint(q / TWO_PI))))


def presence(out):
    ph = np.asarray(out["Cph"], float)
    S, F, L, _ = ph.shape
    rec = out["meta"]["rec"]["Cph"]
    m2 = np.empty((S, F)); mk2 = np.empty((S, F))
    for s in range(S):
        for f in range(F):
            m2[s, f], mk2[s, f] = _mags(ph[s, f])
    late = slice(int(0.8 * F), F)
    M2l = m2[:, late].mean(1); MKl = mk2[:, late].mean(1)
    xi2 = [(float(np.sqrt(max(a / b - 1.0, 0.0)) / (2 * np.sin(np.pi / L)) / L) ** 2) for a, b in zip(M2l, MKl)]
    return {"measure": "squared magnetisation M^2 = |mean spin|^2 per frame (paper's order parameter)",
            "L": int(L), "M2_late_mean_per_seed": M2l.tolist(), "M2_final_per_seed": m2[:, -1].tolist(),
            "M2_every_100_steps": m2[:, :: max(1, 100 // rec)].round(4).tolist(),
            "xi_over_L_sq_late": xi2, "vortex_density_final": [_vortex_density(ph[s, -1]) for s in range(S)],
            "paper_fit_M2_at_L": float(paper_M2(L)),
            "note": "paper fit (Table 3) is for L >= 128; the value at the census L is an extrapolation; "
                    "the census run starts hot and may not have reached equilibrium"}

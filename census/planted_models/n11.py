"""planted:N11 - continuous-space birth-death model with two-body and three-body competition (arXiv 2609.29725,
"Discrete model" paragraph, Eq. (1), Fig. 1, Supplementary S1 and S4; code: github.com/DPintoRamos/
HOI_Competing_Populations, scripts/01_run_ibm_sweep.py).

Model as the paper states it
----------------------------
* A variable number N of individuals at continuous positions in a square periodic domain of side L = 10 (S1 and the
  authors' code: DOMAIN_LENGTH = 10.0), started with N0 = 4,000 individuals placed uniformly at random (S1 item 1).
* Birth (rule 1): every individual reproduces at rate b = 1920; the offspring is displaced from the parent by a
  Gaussian with zero mean and standard deviation s = 0.05 in each coordinate (Fig. 1 caption; code DISPERSAL_WIDTH
  = 0.12/2.4), with periodic wrapping.
* Death: d(r) = d0 + d_c(r) with d0 = 192 (rule 2, density independent) and (Eq. 1, S1 item 4)
      d_c = eta * C_p / (pi r_c^2)  +  gamma * (C_ho / (pi r_ho^2))^2,
  C_p = number of OTHER individuals within r_c = 1.25 (rule 3, two-body competition, eta = 1),
  C_ho = number of other individuals within r_ho = 0.0625 (rule 4, three-body crowding; the kernel
  psi = phi_ho(r - r_i) phi_ho(r - r_j) with top-hat phi_ho gives gamma times the squared density).
* Paper integrator (S1): fixed-step, mortality-first tau-leap with tau = 1e-5; death rates refreshed every m = 10 steps
  (every 1e-4 time units); offspring inherit the parent's death rate until the next refresh. The authors run to
  t = 3.0 and analyse 20 snapshots in t = 2.05..3.0 (code FINAL_TIME = 3.0, SAMPLE_INTERVAL = 0.05, frames 41-60).

Implementation (and why it differs from the paper's scheme where it does)
--------------------------------------------------------------------------
* Death rates are refreshed every Delta = 1e-4 time units, exactly as in the paper (m * tau). Between refreshes the
  process with frozen per-lineage rates is simulated EXACTLY (continuous-time branching: each individual draws an
  exponential death time with its frozen rate, gives birth at the Poisson rate b while alive, offspring inherit the
  parent's rate and are active at once). This is the tau -> 0 limit of the paper's tau-leap at the same refresh
  interval and removes the leap's O(b tau) bias in the birth-death balance (about 1% at the paper's tau).
* C_ho (radius 0.0625) is counted exactly (scipy cKDTree pair search). C_p (radius 1.25, ~8,000 neighbours per
  individual at the steady density) is computed on a 512 x 512 grid (spacing 0.0195): cloud-in-cell deposit, FFT
  convolution with the area-weighted disk of radius 1.25, cloud-in-cell interpolation back, minus the exact self
  term (1). This smears the disk edge over about +-0.04 (3% of r_c); the mean count is unbiased.
* Domain: the paper's (side 10, 4,000 random starters), as the plug-in brief instructs (deviation P3). The plan's rule
  ("paper's domain if it holds at most 10,000 particles") does not fit: this domain holds ~1.6e5 individuals at the
  steady state; the plan's alternative ("about 400 at close packing") has no meaning for point individuals, and a
  domain holding ~400 at the steady density (side ~0.5) would be smaller than the competition radius 1.25.
* Census field view (plan): F0 = counts of individuals binned on a 64 x 64 grid over the paper's domain (spacing
  10/64 = 0.156; area4: 128 x 128 over a domain of side 20 started with 16,000 individuals, same initial density),
  float32, frames every 10 census steps.
* Blow-up (plan, base.BLOWUP): if any seed's population exceeds 50 x its start (200,000 at census size) the run is
  marked unstable and stops after that census step (the frame of that step is recorded only if it falls on the
  recording interval), as the census interpreter does on blow-up. Applied to every run. The real program at the
  headline stays below it (steady population ~1.62e5 = 40.6 x start, initial overshoot to 40.9 x). The 50 x rule is
  kept as frozen in the plan; note that the sub-program with two-body competition only settles at ~1.9-2.0e5
  (47-51 x) and may therefore stop as a blow-up by fluctuation (it did at step 13 in a test); the sub-program with
  three-body crowding only settles far above it (~7.4e5) and stops within 3 steps; the knock-out within 3 steps.

Headline (plan row N11): gamma = 3.14e-5, kept (deviation P2): the plan labels it "stripes", but in the authors' own
IBM data and README 3.1414e-5 is the GAP pattern (stripes at 2.356e-5); gaps are also presented as new in the paper.

Time per census step. The plan fixes one census step = one mean lifetime at the density-independent death rate,
1/192 time units (52 refreshes). At the measured cost (0.085-0.11 s per refresh per seed at the steady population of
~1.62e5) that would be 4-5 CPU-hours for real + knock-out, so by the plan's compute cap one census step is reduced to
1e-3 time units (10 refreshes, about 1.9 mean lifetimes at the steady total death rate ~b), and 1000 census steps =
1.0 time unit (estimated 45-55 CPU-minutes for real + knock-out). The paper analyses t = 2.05..3.0; at the headline
the isotropic pattern is formed by t ~ 0.05 and does not change character up to t = 0.6 (development runs), while
at the paper's stripe value 2.356e-5 elongated, aligned spots are still slowly merging at t = 0.6 (A = 0.38).

Knock-out (plan): rules 3 and 4 removed -> pure birth (1920) and death (192): exponential growth at rate 1728,
blow-up (50 x) within the third census step.

Fidelity (plan, as amended by the coordinator after review; deviation P2): the paper's measure is the anisotropy
index A = median over the image of |l1 - l2| / (l1 + l2), l1,2 the eigenvalues of the structure tensor of the
density field (Gaussian smoothing 0.25, structure-tensor window 2.5, periodic), median over snapshots (S4; A ~ 0
spots/gaps, ~ 1 stripes; if the smoothed field's relative amplitude < 0.01, A = 0 as in the authors' code); hexagonal
order of spots = local psi_6 of the spot centres (maxima of the smoothed field) over Delaunay neighbours. Settings:
the headline 3.14e-5, 3.27e-6, and the paper's stripe value 2.356e-5 when the budget allows. Deciding PASS, fixed
before any large run: A within 0.15 of the authors' own IBM median at every simulated g (their cache
all_coherence_IBM.npy: 0.025 at 3.27e-6, 0.916 at 2.356e-5, 0.169 at 3.14e-5) AND hexagonal spots at 3.27e-6
(A <= 0.3 with spot <|psi_6|> >= 0.6). Reported alongside as plan_check_original: the plan's frozen wording
"stripes at the headline" (A(3.14e-5) >= 0.6) plus the same spots condition; expected to fail because 3.14e-5 is
the authors' gap pattern (also presented as new in the paper).
"""
import time
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree, Delaunay

from census.planted_models import base

CID = "planted:N11"
LAYERS = "F"
K = {}
SPEED = 0.0
RULES = [("BIRTH_GAUSS_DISPERSAL", False), ("DEATH_BASELINE", False), ("PAIR_COMPETITION", True),
         ("TRIPLE_CROWDING", True)]
HEADLINE = dict(
    b=1920.0,        # Fig. 1 caption; code BIRTH_RATE = 640*3
    d0=192.0,        # Fig. 1 caption; code MORTALITY_RATE = 0.1*b
    eta=1.0,         # Fig. 1 caption
    rc=1.25,         # Fig. 1 caption; code 3/2.4
    s=0.05,          # Fig. 1 caption; code 0.12/2.4
    rho=0.0625,      # Fig. 1 caption r_h.o.; code 0.25*0.25
    gamma=3.14e-5,   # plan headline (Fig. 1 third panel; authors' grid value 3.1414e-5)
    L=10.0,          # S1 / code DOMAIN_LENGTH
    N0=4000,         # S1 item 1
    delta=1e-4,      # rate refresh interval = paper m*tau = 10 * 1e-5
    per_step=10,     # refreshes per census step (compute cap): census step = 1e-3 time units
    M=512,           # grid for the pair-count convolution (not in the paper)
)
STEP = "one census step = 1e-3 time units = 10 death-rate refreshes of 1e-4 (paper's refresh interval), exact birth-death in between (plan compute cap; plan default 1/192)"
RULE_TEXT = ("World: individuals at continuous positions in a square periodic domain, starting with 4,000 placed at random. "
             "Rules (continuous-time birth and death): (1) each individual gives birth at rate 1920; the offspring lands at "
             "a Gaussian offset (standard deviation 0.05); (2) each individual dies at rate 192 + (density of others "
             "within radius 1.25) + g times (density within radius 0.0625, squared: a three-body crowding term). "
             "Densities are counts divided by disk area. Only competitive (harmful) interactions, no attraction or help. "
             "Set g = 3.14e-5 and run for 1 time unit; look at the steady spatial pattern.")

_KCACHE = {}
_AUTHORS_A = {3.27e-6: 0.025, 2.356e-5: 0.916, 3.14e-5: 0.169}   # authors' IBM medians (all_coherence_IBM.npy)


# ---------------------------------------------------------------------------------------------------- helpers
def _wrap(x, L):
    x = np.mod(x, L); x[x >= L] = 0.0
    return x


def _disk_hat(M, L, R):
    """rfft2 of the area-weighted indicator of a disk of radius R on an M x M periodic grid of side L."""
    key = (M, L, R)
    if key not in _KCACHE:
        h = L / M; k = np.arange(M); k = np.where(k > M // 2, k - M, k) * h
        sub = (np.arange(4) + 0.5) / 4 - 0.5
        Kg = np.zeros((M, M))
        for a in sub:
            for c in sub:
                Kg += ((k[:, None] + a * h) ** 2 + (k[None, :] + c * h) ** 2 <= R * R)
        _KCACHE[key] = np.fft.rfft2(Kg / 16.0)
    return _KCACHE[key]


def _cic(x, L, M):
    u = x * (M / L); i0 = np.floor(u).astype(np.int64); f = u - i0
    i0 %= M; i1 = (i0 + 1) % M
    return i0, i1, f


def _pair_counts(x, P):
    """number of other individuals within rc, by CIC deposit + FFT disk convolution + CIC interpolation."""
    L, M = P["L"], P["M"]; N = len(x)
    i0, i1, f = _cic(x, L, M)
    w00 = (1 - f[:, 0]) * (1 - f[:, 1]); w10 = f[:, 0] * (1 - f[:, 1]); w01 = (1 - f[:, 0]) * f[:, 1]; w11 = f[:, 0] * f[:, 1]
    idx = [(i0[:, 0], i0[:, 1], w00), (i1[:, 0], i0[:, 1], w10), (i0[:, 0], i1[:, 1], w01), (i1[:, 0], i1[:, 1], w11)]
    G = np.zeros(M * M)
    for a, c, w in idx: G += np.bincount(a * M + c, w, M * M)
    C = np.fft.irfft2(np.fft.rfft2(G.reshape(M, M)) * _disk_hat(M, L, P["rc"]), s=(M, M)).ravel()
    val = np.zeros(N)
    for a, c, w in idx: val += w * C[a * M + c]
    return val - 1.0


def _ho_counts(x, P):
    N = len(x)
    pr = cKDTree(x, boxsize=P["L"], balanced_tree=False, compact_nodes=False).query_pairs(P["rho"], output_type="ndarray")
    if len(pr) == 0: return np.zeros(N)
    return (np.bincount(pr[:, 0], minlength=N) + np.bincount(pr[:, 1], minlength=N)).astype(float)


def _death_rates(x, P, on):
    d = np.full(len(x), P["d0"])
    if len(x) == 0: return d
    if 3 in on: d += P["eta"] * _pair_counts(x, P) / (np.pi * P["rc"] ** 2)
    if 4 in on and P["gamma"] != 0.0: d += P["gamma"] * (_ho_counts(x, P) / (np.pi * P["rho"] ** 2)) ** 2
    return d


def _interval(x, d, P, rng):
    """exact birth-death over one refresh interval with frozen per-lineage death rates; returns survivors' positions."""
    b, s, L, D = P["b"], P["s"], P["L"], P["delta"]
    keep = []
    t0 = np.zeros(len(x)); cx, cd = x, d
    while len(cx):
        w = D - t0
        T = rng.exponential(1.0, len(cx)) / cd
        a = np.minimum(T, w)
        keep.append(cx[T >= w])
        nb = rng.poisson(b * a)
        par = np.repeat(np.arange(len(cx)), nb)
        if len(par) == 0: break
        t0 = t0[par] + rng.random(len(par)) * a[par]
        cx = _wrap(cx[par] + s * rng.standard_normal((len(par), 2)), L)
        cd = cd[par]
    return np.concatenate(keep) if keep else np.zeros((0, 2))


def _bin(x, L, W):
    if len(x) == 0: return np.zeros((W, W), np.float32)
    ix = np.minimum((x[:, 0] * (W / L)).astype(np.int64), W - 1); iy = np.minimum((x[:, 1] * (W / L)).astype(np.int64), W - 1)
    return np.bincount(iy * W + ix, minlength=W * W).reshape(W, W).astype(np.float32)   # rows = y, columns = x


def _params(**over):
    P = dict(HEADLINE); P.update(over); return P


# ---------------------------------------------------------------------------------------------------- census entry
def simulate(seed0, steps=1000, S=3, size="census", keep=None, return_positions=False, **over):
    on = base.kept(keep, RULES)
    P = _params(**over)
    W = base.SIZES[size]["W"]; scale = W // 64
    P["L"] = P["L"] * scale; P["M"] = P["M"] * scale; N0 = int(P["N0"] * scale * scale)
    rs = base.rngs(seed0, S)
    X = [_wrap(r.uniform(0, P["L"], (N0, 2)), P["L"]) for r in rs]
    frames = [[_bin(x, P["L"], W)] for x in X]
    pop = [[N0] for _ in range(S)]
    unstable = False; t = 0
    while t < steps:
        for s in range(S):
            x = X[s]
            for _ in range(P["per_step"]):
                d = _death_rates(x, P, on)
                x = _interval(x, d, P, rs[s])
                if len(x) > base.BLOWUP * N0: break
            X[s] = x; pop[s].append(len(x))
        t += 1
        if any(len(x) > base.BLOWUP * N0 for x in X): unstable = True
        if t % base.REC["F"] == 0:
            for s in range(S): frames[s].append(_bin(X[s], P["L"], W))
        if unstable: break
    F = np.stack([np.stack(f) for f in frames])[:, :, None]          # (S, frames, 1, W, W)
    out = {"F": F.astype(np.float32)}
    out["meta"] = base.meta(t, W, unstable=unstable)
    out["meta"]["pop"] = [list(map(int, p)) for p in pop]
    out["meta"]["pop_max_ratio"] = float(max(max(p) for p in pop) / N0)
    out["meta"]["domain"] = float(P["L"])
    if return_positions: out["_positions"] = X                      # diagnostics only (fidelity); not census output
    return out


# ---------------------------------------------------------------------------------------------------- measures
def _gauss(img, sig):
    return ndimage.gaussian_filter(img, sig, mode="wrap")


def anisotropy(img, dx, smooth=0.25, window=2.5):
    """authors' anisotropy index (S4 / scripts/05_compute_pattern_metrics.py) for one density image."""
    dens = _gauss(img.astype(float), smooth / dx)
    den = dens.max() + dens.min()
    rel = (dens.max() - dens.min()) / den if den != 0 else 0.0
    if rel < 1e-2: return 0.0, dens
    gy = 0.5 * (np.roll(dens, -1, 0) - np.roll(dens, 1, 0)); gx = 0.5 * (np.roll(dens, -1, 1) - np.roll(dens, 1, 1))
    sw = window / dx
    Axx, Axy, Ayy = _gauss(gx * gx, sw), _gauss(gx * gy, sw), _gauss(gy * gy, sw)
    tr = Axx + Ayy; disc = np.sqrt((Axx - Ayy) ** 2 + 4 * Axy ** 2)
    l1, l2 = 0.5 * (tr + disc), 0.5 * (tr - disc)
    return float(np.median((l1 - l2) / (l1 + l2 + 1e-12))), dens


def spot_hexatic(dens, dx):
    """local psi_6 of spot centres (local maxima of the smoothed density above its mean) over periodic Delaunay
    neighbours; returns (<|psi6|>, number of spots)."""
    W = dens.shape[0]; L = W * dx
    mx = (dens == ndimage.maximum_filter(dens, size=max(3, int(round(0.6 / dx))), mode="wrap")) & (dens > dens.mean())
    iy, ix = np.nonzero(mx)
    if len(ix) < 7: return float("nan"), int(len(ix))
    p = np.stack([(ix + 0.5) * dx, (iy + 0.5) * dx], 1)
    imgs, ids = [], []
    for ox in (-1, 0, 1):
        for oy in (-1, 0, 1):
            imgs.append(p + np.array([ox, oy]) * L); ids.append(np.arange(len(p)))
    Q = np.concatenate(imgs); ID = np.concatenate(ids)
    tri = Delaunay(Q).simplices
    e = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]); e = np.concatenate([e, e[:, ::-1]])
    ins = (Q[e[:, 0], 0] >= 0) & (Q[e[:, 0], 0] < L) & (Q[e[:, 0], 1] >= 0) & (Q[e[:, 0], 1] < L)
    e = e[ins]; e = np.unique(np.stack([ID[e[:, 0]], e[:, 1]], 1), axis=0)
    dvec = Q[e[:, 1]] - p[e[:, 0]]
    z = np.exp(6j * np.arctan2(dvec[:, 1], dvec[:, 0]))
    n = len(p); s = np.bincount(e[:, 0], z.real, n) + 1j * np.bincount(e[:, 0], z.imag, n)
    nn = np.bincount(e[:, 0], minlength=n)
    return float(np.mean(np.abs(s / np.maximum(nn, 1)))), int(n)


def presence(out):
    """anisotropy index and spot hexatic order on the census F frames (last 20 frames, every seed)."""
    F = out["F"]; S, nf = F.shape[:2]; W = F.shape[-1]
    L = out["meta"].get("domain", HEADLINE["L"] * (W // 64)); dx = L / W
    if out["meta"].get("unstable"):
        return {"unstable": True, "note": "run stopped on blow-up; no pattern measured"}
    A, H = [], []
    for s in range(S):
        for k in (range(max(1, nf - 20), nf) if nf > 1 else [0]):
            a, dens = anisotropy(F[s, k, 0], dx); A.append(a); H.append(spot_hexatic(dens, dx)[0])
    return {"anisotropy_median": float(np.median(A)), "anisotropy_per_seed": [float(np.median(A[i::S])) for i in range(S)],
            "spot_psi6_median": float(np.nanmedian(H)) if np.any(np.isfinite(H)) else float("nan"),
            "pop_max_ratio": out["meta"].get("pop_max_ratio"),
            "note": "A ~ 1 stripes, ~ 0 spots/gaps (authors' IBM: 0.17 at 3.14e-5 = gaps, 0.92 at 2.36e-5 = stripes)"}


# ---------------------------------------------------------------------------------------------------- fidelity
def _run_ibm(gamma, t_total, seed, nsnap=20, dx_img=0.05, deadline=None):
    """one IBM run in the paper's domain from 4,000 random individuals; anisotropy on nsnap snapshots spread over the
    last third of the run (the authors use 20 snapshots in t = 2.05..3.0). Stops at the wall-clock deadline if given
    (then measures the state reached)."""
    P = _params(gamma=gamma)
    rng = np.random.default_rng(seed)
    x = _wrap(rng.uniform(0, P["L"], (P["N0"], 2)), P["L"])
    n_int = int(round(t_total / P["delta"]))
    snap_at = set(np.linspace(int(round(2 * n_int / 3)), n_int, nsnap).astype(int).tolist())
    A, H, pops = [], [], []
    W = int(round(P["L"] / dx_img)); k = 0

    def measure():
        a, dens = anisotropy(_bin(x, P["L"], W), dx_img); A.append(a); H.append(spot_hexatic(dens, dx_img)[0]); pops.append(len(x))
    for k in range(1, n_int + 1):
        x = _interval(x, _death_rates(x, P, {3, 4}), P, rng)
        if k in snap_at: measure()
        if deadline is not None and time.time() > deadline: break
    if not A: measure()
    return dict(gamma=gamma, t_planned=t_total, t_reached=round(k * P["delta"], 4), A_median=float(np.median(A)),
                spot_psi6_median=float(np.nanmedian(H)) if np.any(np.isfinite(H)) else float("nan"),
                pop_mean=float(np.mean(pops)), snapshots=len(A))


def fidelity(budget_s=7200):
    """anisotropy index at the plan's headline (3.14e-5) and at 3.27e-6, plus the paper's stripe value 2.356e-5."""
    P = HEADLINE; t_start = time.time()
    # cost estimate: refreshes at the steady population (grow a population first, then time 20 refreshes)
    rng = np.random.default_rng(0); x = _wrap(rng.uniform(0, P["L"], (P["N0"], 2)), P["L"]); Pp = _params()
    for _ in range(40): x = _interval(x, _death_rates(x, Pp, {3, 4}), Pp, rng)
    t0 = time.time()
    for _ in range(20): x = _interval(x, _death_rates(x, Pp, {3, 4}), Pp, rng)
    c_ref = (time.time() - t0) / 20 * 1.3          # spots carry more individuals than the uniform state
    settings = [3.14e-5, 3.27e-6, 2.356e-5]
    usable = 0.9 * budget_s - (time.time() - t_start)
    t_total = min(3.0, usable / (len(settings) * c_ref / P["delta"]))
    if t_total < 0.6:                               # drop the third setting before shortening below 0.6 time units
        settings = settings[:2]; t_total = min(3.0, usable / (len(settings) * c_ref / P["delta"]))
    res = {"domain": P["L"], "N0": P["N0"], "t_total": round(t_total, 3), "seconds_per_refresh_est": round(c_ref, 4),
           "estimated_s": round(len(settings) * t_total / P["delta"] * c_ref, 1), "runs": []}
    for i, g in enumerate(settings):
        deadline = t_start + 0.9 * budget_s * (i + 1) / len(settings) + 0.05 * budget_s
        r = _run_ibm(g, t_total, 100 + i, deadline=deadline); r["authors_A"] = _AUTHORS_A.get(g); res["runs"].append(r)
    by = {r["gamma"]: r for r in res["runs"]}
    res["elapsed_s"] = round(time.time() - t_start, 1)
    res["paper"] = ("Fig. 1/2a and authors' IBM cache: A = 0.025 (3.27e-6, hexagonal spots), 0.916 (2.356e-5, stripes), "
                    "0.169 (3.14e-5, gaps)")
    res["tolerance"] = ("PASS = agreement with the authors' anisotropy index at every simulated g (|A - A_authors| <= 0.15) "
                        "AND hexagonal spots at 3.27e-6 (A <= 0.3 and spot psi6 >= 0.6). plan_check_original = the plan's "
                        "frozen wording: A(3.14e-5) >= 0.6 (stripes) and hexagonal spots at 3.27e-6")
    spots = bool(by[3.27e-6]["A_median"] <= 0.3 and by[3.27e-6]["spot_psi6_median"] >= 0.6)
    agree = bool(all(abs(r["A_median"] - r["authors_A"]) <= 0.15 for r in res["runs"]))
    res["authors_agreement"] = agree
    res["spots_hexagonal"] = spots
    res["plan_check_original"] = {"criterion": "stripes at the headline 3.14e-5 (A >= 0.6) and hexagonal spots at 3.27e-6",
                                  "pass": bool(by[3.14e-5]["A_median"] >= 0.6 and spots),
                                  "note": "3.14e-5 is the authors' gap pattern (deviation P2), so this is expected to fail"}
    res["pass"] = bool(agree and spots)
    return res

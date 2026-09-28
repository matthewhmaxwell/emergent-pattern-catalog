"""Cheap behavioural fingerprint per view (pilot draft v0; frozen at registration).

Computed for every emergent view of every program (target < 0.5 s per view), used to cluster programs into behaviour
classes; the expensive battery and the literature check then run per class on its shortest representatives
(cluster-first pipeline, proposed from pilot measurements — see PLAN.md). Features are plain textbook order
parameters, all computed on the late window (last third) unless stated; every value is finite (NaN -> 0).
"""
import numpy as np
from scipy import ndimage

W = 64


def _late(n): return slice(max(0, n - max(n // 3, 2)), n)


def _moran(g):
    x = g.astype(float); x = x - x.mean(); v = (x ** 2).mean()
    if v < 1e-12: return 0.0
    nb = sum(np.roll(np.roll(x, dy, 0), dx, 1) for dy, dx in ((0, 1), (1, 0), (0, -1), (-1, 0))) / 4
    return float((x * nb).mean() / v)


def _spectral_peak(series):
    s = np.asarray(series, float); s = s - s.mean()
    if len(s) < 16 or s.std() < 1e-12: return 0.0, 0.0
    p = np.abs(np.fft.rfft(s)) ** 2; p[0] = 0; i = int(p.argmax())
    return float(p[i] / (p.sum() + 1e-12)), float(i / len(s))


def _corr_length(g):
    """radius where the (type-agreement) correlation first drops below 1/e of its r=1 value."""
    x = g.astype(float); x = x - x.mean(); v = (x ** 2).mean()
    if v < 1e-12: return 0.0
    f = np.fft.fft2(x); c = np.real(np.fft.ifft2(f * np.conj(f))) / x.size / v
    prof = [c[0, r] for r in range(1, W // 2)]
    for r, val in enumerate(prof, 1):
        if val < np.exp(-1): return float(r)
    return float(W // 2)


def grid(hist):
    G = np.stack([h["grid"] for h in hist]); n = len(G); L = _late(n); k = int(G.max()) + 1
    act = np.array([h.get("activity_density", 0.0) for h in hist])
    fr = np.stack([(G[L] == v).mean((1, 2)) for v in range(max(k, 2))], 1).mean(0)
    ent = float(-(fr[fr > 0] * np.log(fr[fr > 0])).sum() / np.log(max(k, 2)))
    last = G[-1]; lab_sizes = []
    for v in np.unique(last):
        lab, m = ndimage.label(last == v); lab_sizes += list(np.bincount(lab.ravel())[1:])
    lab_sizes = np.array(lab_sizes, float) if lab_sizes else np.array([0.0])
    mor = [_moran(g) for g in G[::max(1, n // 20)]]
    sp, fq = _spectral_peak(act[1:])
    interface = float(np.mean([(g != np.roll(g, 1, 0)).mean() + (g != np.roll(g, 1, 1)).mean() for g in G[L]]) / 2)
    ch = G[1:] != G[:-1]; kk = max(k, 2)
    cyc = float(((G[1:] == (G[:-1] + 1) % kk) & ch).sum() / max(ch.sum(), 1))       # share of changes that are t -> t+1
    rev = float(((G[2:] == G[:-2]) & ch[1:] & ch[:-1]).sum() / max((ch[1:] & ch[:-1]).sum(), 1))   # flip-backs
    a_early, a_late = float(act[1:max(2, n // 6)].mean()), float(act[L].mean())
    trend = float(np.log((a_late + 1e-4) / (a_early + 1e-4)))                        # < 0: activity dies down
    return {"g_moran": _moran(last), "g_moran_trend": float(np.polyfit(np.arange(len(mor)), mor, 1)[0]) if len(mor) > 2 else 0.0,
            "g_corrlen": _corr_length(last), "g_activity": float(act[L].mean()), "g_activity_cv": float(act[L].std() / (act[L].mean() + 1e-9)),
            "g_type_entropy": ent, "g_n_types": float((fr > 0.01).sum()), "g_largest_domain": float(lab_sizes.max() / last.size),
            "g_n_domains": float(np.log1p(len(lab_sizes))), "g_interface": interface, "g_spec_peak": sp, "g_spec_freq": fq,
            "g_cyclic_changes": cyc, "g_activity_trend": trend, "g_reversals": rev}


def agents(hist):
    P = np.stack([h["positions"] for h in hist]); H = np.stack([h["headings"] for h in hist]); n = len(P); L = _late(n)
    pol = np.abs(np.exp(1j * H).mean(1))
    B = float(hist[-1].get("box_size", W))
    last = P[-1]; d = last[:, None] - last[None]; d -= B * np.round(d / B); dist = np.sqrt((d ** 2).sum(-1)); np.fill_diagonal(dist, np.inf)
    nn = dist.min(1).mean(); nn_random = 0.5 / np.sqrt(len(last) / B ** 2)
    hist2 = np.histogram2d(last[:, 0], last[:, 1], bins=16, range=[[0, B], [0, B]])[0]
    # angular momentum about the centroid of each agent's local neighbourhood (milling indicator), crude: global
    v = np.stack([np.cos(H[-1]), np.sin(H[-1])], 1); c = last - last.mean(0); c -= B * np.round(c / B)
    ang = float(np.abs(np.mean(np.cross(c, v) / (np.linalg.norm(c, axis=1) + 1e-9))))
    disp = np.linalg.norm(((P[-1] - P[-2] + B / 2) % B) - B / 2, axis=1).mean() if n > 1 else 0.0
    near = dist < 1.5; cnt = near.sum(1); zc = np.exp(1j * H[-1])
    loc = np.abs((near * zc[None, :]).sum(1) + zc) / (cnt + 1)
    loc_rand = np.mean(1 / np.sqrt(cnt + 1))
    out = {"a_local_align": float(loc.mean() - loc_rand), "a_polar": float(pol[L].mean()), "a_polar_std": float(pol[L].std()), "a_nn_ratio": float(nn / nn_random),
           "a_density_cv": float(hist2.std() / (hist2.mean() + 1e-9)), "a_ang_mom": ang, "a_step": float(disp)}
    if "labels" in hist[-1]:
        lab_ = hist[-1]["labels"]; same = (lab_[:, None] == lab_[None]); near = dist < np.percentile(dist[np.isfinite(dist)], 2)
        out["a_type_segregation"] = float((same & near).sum() / max(near.sum(), 1))
        fr = np.bincount(lab_) / len(lab_); out["a_type_entropy"] = float(-(fr[fr > 0] * np.log(fr[fr > 0])).sum())
    return out


def network(hist, adj0=None):
    import networkx as nx
    A = hist[-1].get("adjacency", adj0); op = hist[-1]["opinions"]; g = nx.from_numpy_array(A)
    comps = [len(c) for c in nx.connected_components(g)]
    try:
        comm = [set(np.flatnonzero(op == v)) for v in np.unique(op)]
        Q = nx.algorithms.community.modularity(g, [c for c in comm if c]) if g.number_of_edges() else 0.0
    except Exception:
        Q = 0.0
    fr = np.unique(op, return_counts=True)[1] / len(op)
    ops = np.stack([h["opinions"] for h in hist]); ch = np.abs(np.diff(ops, axis=0)).mean(1) if len(ops) > 1 else np.zeros(1)
    deg = A.sum(1)
    vals = np.unique(np.round(ops, 6)); k = max(len(vals), 2)
    T = np.searchsorted(vals, np.round(ops, 6)).astype(int); chg = T[1:] != T[:-1]
    cyc = float(((T[1:] == (T[:-1] + 1) % k) & chg).sum() / max(chg.sum(), 1)) if k > 2 else 0.0
    frac = np.stack([(T == v).mean(1) for v in range(k)], 1)
    osc = max((_spectral_peak(frac[:, v])[0] for v in range(k)), default=0.0)
    maj0, maj1 = float(frac[0].max()), float(frac[-1].max())
    return {"n_cyclic_changes": cyc, "n_oscillation": float(osc), "n_consensus_gain": maj1 - maj0,
            "n_giant": float(max(comps) / len(op)), "n_components": float(np.log1p(len(comps))), "n_type_modularity": float(Q),
            "n_type_entropy": float(-(fr * np.log(fr)).sum()), "n_activity": float(ch[_late(len(ch))].mean()),
            "n_degree_cv": float(deg.std() / (deg.mean() + 1e-9))}


def phases(hist):
    Th = np.stack([h["theta"] for h in hist]); n = len(Th); L = _late(n)
    r = np.abs(np.exp(1j * Th).mean(1))
    out = {"p_r": float(r[L].mean()), "p_r_std": float(r[L].std())}
    ph = hist[-1].get("phases")
    if ph is not None and np.ndim(ph) == 2:
        z = np.exp(1j * ph); loc = sum(np.roll(np.roll(z, dy, 0), dx, 1) for dy, dx in ((0, 1), (1, 0), (0, -1), (-1, 0))) / 4
        out["p_local_r"] = float(np.abs(loc).mean())
        # phase singularities (vortex density): winding around each plaquette
        a = np.angle(z); d1 = np.angle(np.exp(1j * (np.roll(a, -1, 1) - a))); d2 = np.angle(np.exp(1j * (np.roll(np.roll(a, -1, 1), -1, 0) - np.roll(a, -1, 1))))
        d3 = np.angle(np.exp(1j * (np.roll(a, -1, 0) - np.roll(np.roll(a, -1, 1), -1, 0)))); d4 = np.angle(np.exp(1j * (a - np.roll(a, -1, 0))))
        out["p_vortex_density"] = float((np.abs(d1 + d2 + d3 + d4) > np.pi).mean())
    dph = np.angle(np.exp(1j * np.diff(Th[L], axis=0))).mean(0) if n > 2 else np.zeros(1)
    out["p_freq_spread"] = float(dph.std())
    return out


def field(hist):
    Fs = np.stack([h["field"] for h in hist]); n = len(Fs); L = _late(n); last = Fs[-1]
    x = last - last.mean(); ps = np.abs(np.fft.fft2(x)) ** 2; ps[0, 0] = 0
    ky, kx = np.meshgrid(np.fft.fftfreq(W), np.fft.fftfreq(W), indexing="ij"); kr = np.sqrt(kx ** 2 + ky ** 2)
    bins = np.linspace(0, 0.5, 33); idx = np.digitize(kr.ravel(), bins); rad = np.bincount(idx, ps.ravel(), minlength=35)[1:33]
    pk = int(rad.argmax()); sharp = float(rad[pk] / (rad.sum() + 1e-12))
    lab, m = ndimage.label(last > last.mean() + last.std())
    ch = np.abs(np.diff(Fs[L], axis=0)).mean() if len(Fs[L]) > 1 else 0.0
    return {"f_std": float(last.std()), "f_mean": float(last.mean()), "f_peak_k": float(bins[pk]), "f_peak_sharp": sharp,
            "f_n_blobs": float(np.log1p(m)), "f_change": float(ch), "f_moran": _moran(last)}


def fingerprint(view_name, hist, adj0=None):
    """dict of features (floats) for one view; view type by name prefix."""
    try:
        if view_name in ("C", "C.transient"): f = grid(hist)
        elif view_name == "A": f = agents(hist)
        elif view_name == "N": f = network(hist, adj0)
        elif view_name.endswith(".phase"): f = phases(hist)
        elif view_name.startswith("F"): f = field(hist)
        else: f = {}
    except Exception as e:
        f = {"fp_error": 1.0}
    return {k: (float(v) if np.isfinite(v) else 0.0) for k, v in f.items()}

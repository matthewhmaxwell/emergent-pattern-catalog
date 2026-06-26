"""Vet the 7 Stage-2 Kuramoto tripwire leads. Params reproducible from seed (gen used
rng(3000+i), order K/omega_spread/init). For each: global coherence r (trajectory), LOCAL
order across the lattice (chimera = coexisting coherent+incoherent -> high spread), and the
descriptor. Distinguish: INCOHERENT/chaotic desync (mundane) vs CHIMERA / structured
partial-sync (interesting, known-to-science at best -> would need the literature gate)."""
import numpy as np
from scipy.ndimage import uniform_filter
from epc.models.kuramoto_lattice import kuramoto_lattice
from epc.phase2a.ring2_descriptor import ring2_descriptor


def local_order(theta, w=5):
    c = uniform_filter(np.cos(theta), w, mode="wrap")
    s = uniform_filter(np.sin(theta), w, mode="wrap")
    return np.sqrt(c ** 2 + s ** 2)


for i in [9, 21, 27, 32, 59, 62, 98]:
    rng = np.random.default_rng(3000 + i)
    K = float(rng.uniform(0.1, 3.0)); osp = float(rng.uniform(0.0, 2.0))
    init = ["random", "spiral", "plane"][int(rng.integers(0, 3))]
    h = kuramoto_lattice(i, N=40, K=K, omega_spread=osp, init=init, steps=3000, record=120)
    rtraj = np.array([abs(np.exp(1j * np.asarray(f["phases"]).ravel()).mean()) for f in h])
    lo = local_order(np.asarray(h[-1]["phases"]))
    d = ring2_descriptor(h, {})
    rm, rs = rtraj.mean(), rtraj.std()
    lom, los = float(lo.mean()), float(lo.std())
    if los > 0.18 and lom < 0.85:
        v = "CHIMERA-like (coexisting coherent+incoherent) -> INSPECT"
    elif rm < 0.3 and lom < 0.45:
        v = "INCOHERENT/chaotic desync -> MUNDANE"
    elif rm > 0.8:
        v = "SYNC (why em low? check)"
    else:
        v = "partial/frustrated -> inspect"
    print(f"#{i} K={K:.2f} osp={osp:.2f} init={init:<6} | global r mean {rm:.2f} std {rs:.2f} "
          f"| local-order mean {lom:.2f} std {los:.2f} | em={d['em_score']} C={d['mf_C']} tripped={d['tripped']}")
    print(f"     VERDICT: {v}")

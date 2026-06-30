"""Composition probe (Track 1): does demanding TWO competencies at once yield a COMPOSITION of known
mechanisms, or a unified competency that resists both names? Task = count-then-communicate.
  - A SPEAKER (GRU) watches K stream rounds, each carrying an item-bit; the TRUE COUNT = number of items
    (0..K). It must COUNT (accumulate over rounds -> internal state). It then emits an M-slot message over
    an alphabet of size S. With S^? : a single symbol (S=2) cannot label K+1=6 counts, so the speaker is
    FORCED to use a multi-slot (compositional) code.
  - A LISTENER (MLP) sees ONLY the M-symbol message and must output the count (K+1-way).
Shared reward = exact-count. Demands counting (speaker) AND communication (a co-adapted numeral code).
Discriminators: channel-scramble (random message -> collapse => uses COMMUNICATION), count-ablate (zero the
speaker's stream -> collapse => uses COUNTING), and a per-slot ablation to measure COMPOSITIONALITY (how
many message slots carry information). The agent-observer then asks: is the emergent code a known numeral
system, a known compositional language, or something that resists naming?

Run on VPS epc-venv: python comp_ppo.py [--iters 800]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

K = 5; NC = K + 1; S = 2; M = 3; H = 64   # K stream rounds; NC counts(0..5); alphabet S; message slots M
torch.set_num_threads(8)


class Speaker(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(2, H, batch_first=True)
        self.heads = nn.ModuleList([nn.Linear(H, S) for _ in range(M)]); self.v = nn.Linear(H, 1)

    def forward(self, stream):                 # stream (B,K,2)
        out, _ = self.gru(stream); h = out[:, -1]
        logits = torch.stack([head(h) for head in self.heads], 1)   # (B,M,S)
        return logits, self.v(h).squeeze(-1)


class Listener(nn.Module):
    def __init__(self):
        super().__init__()
        self.body = nn.Sequential(nn.Linear(M * S, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, NC); self.v = nn.Linear(H, 1)

    def forward(self, msg):                    # msg (B, M*S) one-hot
        z = self.body(msg); return self.pi(z), self.v(z).squeeze(-1)


def episode(spk, lis, B, seed, greedy=False, scramble=False, countablate=False, dropslot=None):
    rng = np.random.default_rng(seed)
    count = rng.integers(0, K + 1, size=B)                 # UNIFORM count 0..K -> 1-bit code caps at chance-ish,
    order = rng.random((B, K)).argsort(1)                  # forcing a multi-slot (compositional) code for reward
    bits = (order < count[:, None]).astype(np.float32)
    stin = bits.copy()
    if countablate: stin = np.zeros_like(bits)            # speaker sees nothing -> cannot count
    stream = np.zeros((B, K, 2), np.float32); stream[:, :, 1] = stin; stream[:, :, 0] = 1 - stin
    slog, sval = spk(torch.from_numpy(stream))
    sd = torch.distributions.Categorical(logits=slog)
    sym = slog.argmax(2) if greedy else sd.sample()       # (B,M)
    slp = sd.log_prob(sym).sum(1)
    msgn = sym.numpy()
    if scramble: msgn = rng.integers(0, S, size=(B, M))
    if dropslot is not None: msgn = msgn.copy(); msgn[:, dropslot] = 0
    msg = np.zeros((B, M * S), np.float32)
    for j in range(M): msg[np.arange(B), j * S + msgn[:, j]] = 1.0
    llog, lval = lis(torch.from_numpy(msg))
    ld = torch.distributions.Categorical(logits=llog)
    pred = llog.argmax(1) if greedy else ld.sample()
    llp = ld.log_prob(pred)
    R = (pred.numpy() == count).astype(np.float32)
    acc = float(R.mean())
    return dict(stream=torch.from_numpy(stream), sym=sym, slp=slp.detach(), sval=sval.detach(),
                msg=torch.from_numpy(msg), pred=pred, llp=llp.detach(), lval=lval.detach(),
                R=torch.from_numpy(R), acc=acc)


def train(iters, B=512, seed=0):
    torch.manual_seed(seed); spk, lis = Speaker(), Listener()
    opt = torch.optim.Adam(list(spk.parameters()) + list(lis.parameters()), lr=3e-3)
    for it in range(iters):
        ep = episode(spk, lis, B, 1000 + seed * 99999 + it)
        R = ep["R"]
        sadv = (R - ep["sval"]); sadv = (sadv - sadv.mean()) / (sadv.std() + 1e-8)
        ladv = (R - ep["lval"]); ladv = (ladv - ladv.mean()) / (ladv.std() + 1e-8)
        for _ in range(4):
            slog, sval = spk(ep["stream"]); sd = torch.distributions.Categorical(logits=slog)
            slp = sd.log_prob(ep["sym"]).sum(1); sratio = torch.exp(slp - ep["slp"])
            s1 = sratio * sadv; s2 = torch.clamp(sratio, 0.8, 1.2) * sadv
            sloss = -torch.min(s1, s2).mean() + 0.5 * ((sval - R) ** 2).mean() - 0.02 * sd.entropy().sum(1).mean()
            llog, lval = lis(ep["msg"]); ld = torch.distributions.Categorical(logits=llog)
            lratio = torch.exp(ld.log_prob(ep["pred"]) - ep["llp"])
            l1 = lratio * ladv; l2 = torch.clamp(lratio, 0.8, 1.2) * ladv
            lloss = -torch.min(l1, l2).mean() + 0.5 * ((lval - R) ** 2).mean() - 0.01 * ld.entropy().mean()
            opt.zero_grad(); (sloss + lloss).backward()
            nn.utils.clip_grad_norm_(list(spk.parameters()) + list(lis.parameters()), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            acc = np.mean([episode(spk, lis, 600, 7000 + it + k, greedy=True)["acc"] for k in range(2)])
            print(f"  iter {it:>3}: exact-count accuracy {acc:.2f}  (chance {1/NC:.2f})", flush=True)
    return spk, lis


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 800
    here = __file__.rsplit("/", 1)[0]
    print(f"COMPOSITION probe: count(0..{K}) -> {M}-slot/{S}-symbol message -> decode. chance={1/NC:.2f}", flush=True)
    spk, lis = train(iters)
    def ev(**kw): return round(float(np.mean([episode(spk, lis, 800, 9000 + k, greedy=True, **kw)["acc"] for k in range(3)])), 2)
    norm = ev(); scr = ev(scramble=True); cab = ev(countablate=True)
    drops = [ev(dropslot=j) for j in range(M)]
    print(f"\nFINAL accuracy: normal {norm} | channel-scramble {scr} | count-ablate {cab} (chance {1/NC:.2f})", flush=True)
    print(f"  per-slot drop (accuracy with slot j zeroed): {drops}  -> # informative slots = compositionality", flush=True)
    comp = sum(1 for d in drops if norm - d >= 0.1)
    print(f"VERDICT: composition {'SOLVED' if norm >= 0.7 else 'partial'} (acc {norm}); uses comms={norm-scr>=0.2} "
          f"counting={norm-cab>=0.2}; informative message slots = {comp}/{M} ({'COMPOSITIONAL' if comp>=2 else 'single-symbol'})", flush=True)
    torch.save({"spk": spk.state_dict(), "lis": lis.state_dict()}, f"{here}/comp_net.pt")
    json.dump({"normal": norm, "scramble": scr, "count_ablate": cab, "slot_drops": drops,
               "informative_slots": comp}, open(f"{here}/comp_result.json", "w"))

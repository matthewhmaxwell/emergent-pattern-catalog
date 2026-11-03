"""Fidelity checks of the planted-model plug-ins (PLANTED_PIPELINE_PLAN.md, "Implementation and fidelity"). Each module's
fidelity() reproduces one observation from its paper with the paper's own measure; no census pipeline code is run.

  python -m census.planted_fidelity [N02 N07 ...] [--workers 6]   -> census/runs/planted_pipeline/fidelity/<NXX>.json
"""
import argparse, json, os, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

OUT = "census/runs/planted_pipeline/fidelity"
PLUGINS = ["N01", "N02", "N03", "N04", "N07", "N08", "N11", "N12", "N13", "N14", "N15", "N19"]


def one(n):
    import importlib, traceback
    fn = f"{OUT}/{n}.json"; t = time.time()
    try:
        r = importlib.import_module(f"census.planted_models.{n.lower()}").fidelity()
    except Exception as e:
        r = {"error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-1500:]}
    r["secs"] = round(time.time() - t, 1); r["run"] = time.strftime("%Y-%m-%d %H:%M:%S")
    json.dump(r, open(fn, "w"), indent=1, default=str)
    return n, r.get("pass"), r["secs"]


if __name__ == "__main__":
    import multiprocessing as mp
    ap = argparse.ArgumentParser(); ap.add_argument("cards", nargs="*"); ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args(); os.makedirs(OUT, exist_ok=True)
    with mp.get_context("fork").Pool(a.workers, maxtasksperchild=1) as pool:
        for n, ok, secs in pool.imap_unordered(one, a.cards or PLUGINS):
            print("fidelity", n, "pass" if ok else ("FAIL" if ok is False else ok), f"{secs}s", flush=True)
    print("FIDELITY DONE", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)

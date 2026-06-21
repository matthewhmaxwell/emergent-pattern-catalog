import json, os, time
import analysis.run_phase2a_panel as R
from gallery.registry import GALLERY
from gallery.anim import render, FRAME_MEASURE
from gallery.code_extract import extract_code
from epc.phase2a.continuous_metrics import CONTINUOUS_METRIC
from gallery.education import MECHANISM, WHERE
ASSETS="gallery/assets"; os.makedirs(ASSETS, exist_ok=True)
CONF=json.load(open("docs/validation_rebuild/battery_confusion_matrix_2026-06-16.json"))
from gallery.gallery_runs import GALLERY_RUNS
def build(pid):
    if pid in GALLERY_RUNS: return GALLERY_RUNS[pid]()
    if pid=="p31":
        from epc.models.cell_view_sorting import CellViewSorting
        return [CellViewSorting(n=60,algorithm="insertion",seed=1).run_to_completion(max_rounds=400)],{}
    return getattr(R,f"build_{pid}_positives")(n_seeds=1)
man=[]
for pid_up,info in GALLERY.items():
    pid=pid_up.lower(); t0=time.time()
    e={"id":pid_up,**{k:info[k] for k in ("name","ref","summary","effect","viz","watch")}}
    e["metric"]=CONTINUOUS_METRIC.get(pid_up,{}).get("key"); e["mechanism"]=MECHANISM.get(pid_up,""); e["where"]=WHERE.get(pid_up,"")
    cm=CONF.get(pid,{})
    e["detector"]={"top":cm.get("top_pattern"),"tier":cm.get("top_tier"),"confidence":cm.get("top_conf"),
                   "verdict":cm.get("verdict"),"emergence":cm.get("emergence"),"self_recognized":cm.get("correct")}
    if pid=="p31": e["detector"]["note"]="Validated via a separate multi-run non-redundancy test (no Phase-2a panel)."
    c=extract_code(pid_up); e["code"]=c["source"]; e["code_where"]=c["where"]; e["code_module"]=c["module"]
    try:
        runs,meta=build(pid)
        m=render(runs[0], info["viz"], f"{ASSETS}/{pid}", f"{pid_up} {info['name']}", info.get("args"), measure=FRAME_MEASURE.get(pid_up))
        if m:
            e["asset"]=f"{pid}_sprite.png"; e["asset_type"]="sprite"; e["sprite"]=m; e["mp4"]=m.get("mp4"); e["contact_sheet"]=e["asset"]
            print(f"  {pid_up:<4} {info['viz']:<16} OK {m['frames']}f gif={m.get('gif')} ({round(time.time()-t0,1)}s)", flush=True)
        else:
            e["asset"]=None; e["asset_type"]="image"
            print(f"  {pid_up:<4} {info['viz']:<16} EMPTY", flush=True)
    except Exception as ex:
        e["asset"]=None; e["asset_type"]="image"; e["error"]=repr(ex)
        print(f"  {pid_up:<4} {info['viz']:<16} ERR {type(ex).__name__}: {str(ex)[:50]}", flush=True)
    man.append(e)
from gallery.validation_enrich import enrich
enrich(man)
json.dump(man, open("gallery/manifest.json","w"), indent=2)
print(f"DONE: {sum(1 for e in man if e.get('asset'))}/32 animated")
open("/tmp/genall_DONE","w").write("done")

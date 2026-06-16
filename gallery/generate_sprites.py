import json, os
import analysis.run_phase2a_panel as R
from gallery.registry import GALLERY
from gallery.sprites import render_sprite
ASSETS="gallery/assets"
ANIM={"point_cloud","grid_field","phase_circle"}
man={e["id"]:e for e in json.load(open("gallery/manifest.json"))}
for pid_up,info in GALLERY.items():
    pid=pid_up.lower(); kind=info["viz"]; e=man[pid_up]
    if kind not in ANIM:
        e["asset_type"]="image"; continue
    try:
        runs,meta=getattr(R,f"build_{pid}_positives")(n_seeds=1)
        out=f"{ASSETS}/{pid}_sprite.png"
        m=render_sprite(runs[0], kind, out, f"{pid_up} {info['name']}")
        if m and os.path.getsize(out)>500:
            e["asset"]=os.path.basename(out); e["asset_type"]="sprite"; e["sprite"]=m
            print(f"  {pid_up:<4} sprite {m['frames']}f {m['cols']}x{m['rows']} {os.path.getsize(out)//1024}KB", flush=True)
        else:
            e["asset_type"]="image"; print(f"  {pid_up} sprite EMPTY -> keep image", flush=True)
    except Exception as ex:
        e["asset_type"]="image"; print(f"  {pid_up} ERR {type(ex).__name__}: {str(ex)[:40]}", flush=True)
json.dump(sorted(man.values(),key=lambda e:int(e["id"][1:])), open("gallery/manifest.json","w"), indent=2)
nS=sum(1 for e in man.values() if e.get("asset_type")=="sprite")
print(f"DONE: {nS} sprites + {32-nS} static images")
open("/tmp/sprites_DONE","w").write("done")

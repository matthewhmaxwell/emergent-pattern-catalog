"""Append the 6 Ring-3 competency demos to gallery/manifest.json as a separate 'competency' track.
Idempotent (removes any prior C## entries first) + backs up manifest.json. Run from repo root."""
import json, os, shutil
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
MAN = ROOT + "/gallery/manifest.json"
CAT = json.load(open(ROOT + "/analysis/ring3_competency/canonical_catalog.json"))
cat = {c["id"]: c for c in CAT["competencies"]}

DEMOS = {
 22: dict(base="division", summary="n identical agents must each take a different target — emergent task allocation (a perfect matching) with no assigned roles.",
          watch="Two agents start together, then split so each settles on a different target — one-agent-per-target; the assignment holds once formed.",
          effect="Two symmetric agents self-organize into complementary roles — one per target.",
          sprite=dict(frames=31, cols=6, rows=6, fw=300, fh=300, mp4="ring3_division.mp4")),
 24: dict(base="stigmergy", summary="Agents that cannot see one another coordinate coverage of a space purely through a shared, decaying pheromone field.",
          watch="Four blind agents lay pheromone trails and steer off each other's marks to divide the torus, covering nearly all target cells without ever seeing one another.",
          effect="Blind agents coordinate purely through marks left in the environment, covering the space with no direct contact.",
          sprite=dict(frames=25, cols=5, rows=5, fw=300, fh=300, mp4="ring3_stigmergy.mp4")),
 26: dict(base="niche", summary="Agents durably modify the environment (digging a breach) to make an otherwise-unreachable goal reachable — niche construction.",
          watch="Agents dig through a wall — its columns lose health and breach — to open a durable passage, then use it to reach goals that alternate sides.",
          effect="Agents durably reshape the world (a breached wall) to make a goal reachable — with no reward for building.",
          sprite=dict(frames=25, cols=5, rows=5, fw=300, fh=300, mp4="ring3_niche.mp4")),
 28: dict(base="momentum", summary="Control under continuous physics: steering a puck that carries momentum, against inertia and drag.",
          watch="Agents thrust into a puck and, accounting for inertia and drag, drive it across the arena into the goal circle — smooth momentum-carrying motion.",
          effect="Agents master a continuous-physics puck under inertia and drag, not grid steps.",
          sprite=dict(frames=51, cols=8, rows=7, fw=300, fh=300, mp4="ring3_momentum.mp4")),
 29: dict(base="morphology", summary="Agents of different body types each specialize on the resource their morphology harvests best — specialization by body, not identity.",
          watch="Red-bodied agents converge on red food and blue on blue — the population sorts by colour onto matching resources.",
          effect="A mixed population specializes by body type, each morphology harvesting its matching resource.",
          sprite=dict(frames=31, cols=6, rows=6, fw=300, fh=300, mp4="ring3_morphology.mp4")),
 30: dict(base="compositional", summary="A single policy conditions on two independently-switching signals at once — its own body colour and the opponent's block — the first joint-conditioning face.",
          watch="A single agent moves to whichever site is both unblocked by the opponent and serving its own body colour — tracking two independently-switching signals at once.",
          effect="One agent conditions on two independently-switching signals at once — reachable only with an attention (Transformer) policy.",
          sprite=dict(frames=27, cols=6, rows=5, fw=300, fh=300, mp4="ring3_compositional.mp4")),
}


def short_ref(lit):
    for sep in (" (", ";", " / "):
        if sep in lit: lit = lit.split(sep)[0]
    return lit.strip()[:52]


m = json.load(open(MAN))
m = [e for e in m if not str(e.get("id", "")).startswith("C")]  # drop prior competency entries (idempotent)
for cid, D in DEMOS.items():
    c = cat[cid]
    m.append({
        "id": f"C{cid}", "name": c["name"], "ref": short_ref(c["lit_name"]),
        "summary": D["summary"], "effect": D["effect"], "viz": "point_cloud", "watch": D["watch"],
        "metric": c["metric"], "mechanism": c["mechanism"], "where": "",
        "detector": {"note": f"gate: KNOWN ({short_ref(c['lit_name'])}). Diagnostic — {c['diagnostic']}"},
        "asset": f"ring3_{D['base']}_sprite.png", "asset_type": "sprite", "sprite": D["sprite"],
        "mp4": f"ring3_{D['base']}.mp4", "contact_sheet": f"ring3_{D['base']}_sprite.png", "track": "competency",
    })
shutil.copy(MAN, MAN + ".bak")
json.dump(m, open(MAN, "w"), indent=2)
print("manifest entries:", len(m), "| competency:", sum(1 for e in m if e.get("track") == "competency"))

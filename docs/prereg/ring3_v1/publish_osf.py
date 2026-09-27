#!/usr/bin/env python3
"""Publish the EPC Ring-3 prereg to OSF. Adapted from agent-society-decomp/osf/publish_osf.py.

Token: read from the 1Password item "OSF token" (field `credential`) via the `op` CLI; falls back to ~/.osf_token.
It is never printed or written to disk.

Two deliberately separate steps:
  publish_osf.py --draft-only
      REVERSIBLE: creates a PRIVATE project, uploads the frozen files, creates a draft registration.
      Writes node/draft ids (no secrets) to osf_state.json.
  publish_osf.py --register (--embargo-months N | --public)
      IRREVERSIBLE: submits the draft as a registration. Run only after explicit approval.
"""
import json, os, sys, subprocess, urllib.request, urllib.parse, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
API = "https://api.osf.io/v2"
STATE = os.path.join(HERE, "osf_state.json")
TITLE = "Pre-registered adversarial tests of the communication-niche claim (Emergent Pattern Catalog, Ring 3)"
DESCRIPTION = ("Pre-registration: adversarial tests of two claims from the Emergent Pattern Catalog (EPC) - the "
               "communication niche ('symmetric coordination never forces communication') and the "
               "minimal-sufficient-mechanism law - in multi-agent reinforcement learning. Four tasks (symmetric "
               "anti-coordination with a talk channel; asymmetric information with partner observation; a "
               "single-agent memory-perception dial; cross-play of learned conventions), numeric hit/miss "
               "criteria, three seeds, fixed budget, no tuning, and recorded prior credences. Frozen before any "
               "data were collected. Author: Matt Maxwell.")
FILES = ["PREREG_v1.md", "PREREG_v1.sha256", "registration_summary.md"]
SCHEMA_OPEN_ENDED = "5df83f7dd28338001ac0ab0d"   # Open-Ended Registration (as used for osf.io/tayuc)


def _token():
    try:
        out = subprocess.run(["op", "item", "get", "OSF token", "--format", "json"],
                             capture_output=True, text=True, timeout=60).stdout
        d = json.loads(out)
        return next(f["value"] for f in d["fields"] if f.get("id") == "credential").strip()
    except Exception:
        return open(os.path.expanduser("~/.osf_token")).read().strip()


TOK = _token()
H = {"Authorization": f"Bearer {TOK}", "Content-Type": "application/vnd.api+json"}


def req(method, url, payload=None, ctype=None, raw=None):
    data = raw if raw is not None else (json.dumps(payload).encode() if payload else None)
    r = urllib.request.Request(url, data=data, method=method, headers={**H, **({"Content-Type": ctype} if ctype else {})})
    try:
        with urllib.request.urlopen(r, timeout=60) as resp:
            return json.loads(resp.read() or "{}")
    except urllib.error.HTTPError as e:
        print(f"HTTP {e.code} on {method} {url.split('?')[0]}: {e.read().decode()[:500]}", file=sys.stderr)
        raise


def draft_only():
    for fn in FILES:
        assert os.path.exists(os.path.join(HERE, fn)), f"missing {fn} (freeze first)"
    proj = req("POST", f"{API}/nodes/", {"data": {"type": "nodes", "attributes": {
        "title": TITLE, "category": "project", "description": DESCRIPTION[:990], "public": False}}})
    nid = proj["data"]["id"]; print("private project:", f"https://osf.io/{nid}")
    for fn in FILES:
        up = f"https://files.osf.io/v1/resources/{nid}/providers/osfstorage/?kind=file&name={urllib.parse.quote(fn)}"
        req("PUT", up, raw=open(os.path.join(HERE, fn), "rb").read(), ctype="text/markdown")
        print("uploaded:", fn)
    draft = req("POST", f"{API}/nodes/{nid}/draft_registrations/", {"data": {"type": "draft_registrations",
        "relationships": {"registration_schema": {"data": {"type": "registration-schemas", "id": SCHEMA_OPEN_ENDED}}}}})
    did = draft["data"]["id"]
    summary = open(os.path.join(HERE, "registration_summary.md")).read()[:4990]
    req("PATCH", f"{API}/draft_registrations/{did}/", {"data": {"type": "draft_registrations", "id": did,
        "attributes": {"registration_responses": {"summary": summary}}}})
    json.dump({"node_id": nid, "draft_id": did, "created": datetime.datetime.now().isoformat()},
              open(STATE, "w"), indent=1)
    print("draft registration:", did, "(NOT submitted; state saved to osf_state.json)")


def register():
    st = json.load(open(STATE)); nid, did = st["node_id"], st["draft_id"]
    if "--embargo-months" in sys.argv:
        m = int(sys.argv[sys.argv.index("--embargo-months") + 1])
        lift = (datetime.date.today() + datetime.timedelta(days=30 * m)).isoformat()
        attrs = {"draft_registration": did, "registration_choice": "embargo", "lift_embargo": lift + "T00:00:00"}
    elif "--public" in sys.argv:
        lift, attrs = None, {"draft_registration": did, "registration_choice": "immediate"}
    else:
        sys.exit("choose --embargo-months N or --public")
    reg = req("POST", f"{API}/nodes/{nid}/registrations/", {"data": {"type": "registrations", "attributes": attrs}})
    rid = reg["data"]["id"]
    st.update({"registration_id": rid, "registration_url": f"https://osf.io/{rid}", "embargo_until": lift})
    json.dump(st, open(STATE, "w"), indent=1)
    print("REGISTRATION:", f"https://osf.io/{rid}", f"(embargo until {lift})" if lift else "(public immediately)")


if __name__ == "__main__":
    if "--draft-only" in sys.argv:
        draft_only()
    elif "--register" in sys.argv:
        register()
    else:
        sys.exit(__doc__)

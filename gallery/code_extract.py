"""Extract each model's core algorithm method source (authentic, via inspect)."""
import inspect, textwrap, importlib

# pid -> (module, symbol)  symbol = class name or function name
CODE_MAP = {
 "P1":("epc.models.schelling","run_schelling"),
 "P2":("epc.models.active_brownian_particles","ActiveBrownianParticles"),
 "P3":("epc.models.gray_scott","GrayScott"),
 "P4":("epc.models.territoriality","ScentMarkingModel"),
 "P5":("epc.models.vicsek","VicsekModel"),
 "P6":("epc.models.dorsogna_spp","DOrsognaSPPModel"),
 "P7":("epc.models.lane_formation","LaneFormationModel"),
 "P8":("epc.models.nagel_schreckenberg","NagelSchreckenberg"),
 "P9":("epc.models.kuramoto","KuramotoModel"),
 "P10":("epc.models.kuramoto_nonlocal","KuramotoNonlocal"),
 "P11":("epc.models.lotka_volterra_lattice","LotkaVolterraLattice"),
 "P12":("epc.models.rps_spatial","RPSSpatialModel"),
 "P13":("epc.models.greenberg_hastings","GreenbergHastings"),
 "P14":("epc.models.btw_sandpile","run_sandpile"),
 "P15":("epc.models.game_of_life","GameOfLife"),
 "P16":("epc.models.hopfield","HopfieldNetwork"),
 "P17":("epc.models.collective_sensing","CollectiveSensingModel"),
 "P18":("epc.models.voter","VoterModel"),
 "P19":("epc.models.informed_minority","InformedMinorityModel"),
 "P20":("epc.models.quorum_sensing","AutoinducerQuorum"),
 "P21":("epc.models.hegselmann_krause","HegselmannKrauseModel"),
 "P22":("epc.models.sir_epidemic","SIREpidemicModel"),
 "P23":("epc.models.minority_game","MinorityGame"),
 "P24":("epc.models.homeostasis","ProportionalHomeostat"),
 "P25":("epc.models.canalization","CanalizedLandscape"),
 "P26":("epc.models.stochastic_resonance","BistableDoubleWell"),
 "P27":("epc.models.nowak_may","NowakMayModel"),
 "P28":("epc.models.yard_sale","YardSale"),
 "P29":("epc.models.trail_network","AntTrailModel"),
 "P30":("epc.models.autopoiesis","AutopoiesisModel"),
 "P31":("epc.models.cell_view_sorting","CellViewSorting"),
 "P32":("epc.models.division_of_labor","ResponseThresholdModel"),
}
_METHODS = ["step","_step","update","_update","_transaction","_round","simulate","run","run_to_completion"]

def extract_code(pid):
    mod_name, sym = CODE_MAP[pid]
    mod = importlib.import_module(mod_name)
    obj = getattr(mod, sym)
    try:
        if inspect.isfunction(obj):
            src = inspect.getsource(obj); where = f"{sym}()"
        else:  # class: pick the first core method present
            meth = next((m for m in _METHODS if callable(getattr(obj, m, None))
                         and m in obj.__dict__), None)
            if meth is None:
                src = inspect.getsource(obj.__init__); where = f"{sym}.__init__"
            else:
                src = inspect.getsource(getattr(obj, meth)); where = f"{sym}.{meth}()"
        src = textwrap.dedent(src)
        # cap very long methods for display
        lines = src.splitlines()
        if len(lines) > 60:
            lines = lines[:60] + ["    # … (truncated; full source in the repo)"]
        return {"source": "\n".join(lines), "where": where, "module": mod_name}
    except Exception as e:
        return {"source": f"# (source unavailable: {e})", "where": sym, "module": mod_name}

if __name__ == "__main__":
    import json
    man = json.load(open("gallery/manifest.json"))
    for e in man:
        c = extract_code(e["id"])
        e["code"] = c["source"]; e["code_where"] = c["where"]; e["code_module"] = c["module"]
        nl = len(c["source"].splitlines())
        print(f"  {e['id']:<4} {c['where']:<32} {nl} lines")
    json.dump(man, open("gallery/manifest.json","w"), indent=2)
    print("manifest updated with code")

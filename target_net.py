#!/usr/bin/env python3
"""
target_net.py -- network-based target discovery around a core gene (test case: AKT1)

Modes
  upstream    core gene -> receptors            (which inputs drive the core gene; alternative intervention points)
  downstream  core gene -> TFs -> target genes  (what the core gene does; signed, to gene expression via DoRothEA)
  sl          synthetic lethality with core-gene loss (DepMap), separate volcano view
  crosstalk   does gene X influence the pathway of causal gene C?  (shared nodes, interference edges, endpoint agreement)

Examples
  python target_net.py upstream   --gene AKT1 --demo
  python target_net.py downstream --gene AKT1 --regulon dorothea --conf ABC --tissue liver
  python target_net.py sl         --gene AKT1 --depmap ~/data/depmap/
  python target_net.py crosstalk  --causal AKT1 --goi MAPK1
Outputs go to out_<GENE>/<mode>/ (CSV tables + standalone HTML views).
"""
import argparse, json, os, sys
import pandas as pd
from targetnet import data as D
from targetnet import walk as Wk
from targetnet import config as Cf

# symbols people type that have no signaling edges; map to the functional gene and say why
SYMBOL_NOTES = {"HSP90AA2": ("HSP90AA1", "HSP90AA2 is the pseudogene HSP90AA2P (UniProt Q14568, putative); its peptides are shared with HSP90AA1"),
                "HSP90AA2P": ("HSP90AA1", "pseudogene of HSP90AA1"), "HSP90A": ("HSP90AA1", "old alias"), "HSP90B": ("HSP90AB1", "old alias"),
                "AKT": ("AKT1", "family name"),
                "PDK1_NOTE": ("PDK1", "in SIGNOR tables 'PDK1' is pyruvate dehydrogenase kinase 1, NOT PDPK1 (3-phosphoinositide-dependent kinase)"), "ERK2": ("MAPK1", "alias"), "ERK1": ("MAPK3", "alias"), "MEK1": ("MAP2K1", "alias"), "P53": ("TP53", "alias")}

def resolve_symbol(g):
    g = g.upper()
    if g in SYMBOL_NOTES:
        new, why = SYMBOL_NOTES[g]; print(f"[symbol] {g} -> {new}: {why}"); return new
    return g

TPL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "targetnet", "templates")
COMPARTMENTS = ["extracellular", "membrane", "cytosol", "lysosome", "mitochondrion", "ER_golgi", "nucleus"]


# ----------------------------------------------------------------------------- shared setup
def build_signaling(args):
    seeds = [resolve_symbol(g) for g in (getattr(args, "gene", None), getattr(args, "causal", None), getattr(args, "goi", None)) if g]
    seeds += list(getattr(args, "_extra_seeds", []) or [])          # cross-talk screen: every candidate target
    sig = D.load_signor(args.demo)
    parts = [sig]; und = pd.DataFrame(columns=["a", "b", "source"])
    if not args.no_reactome:
        rd, und = D.load_reactome_fi(args.reactome_fi, args.demo); parts.append(rd)
    omni = D.load_omnipath_signed(args.demo, signor_df=sig, seeds=seeds)
    if len(omni): parts.append(omni)
    cust = D.load_custom_edges(args.demo)
    if len(cust): parts.append(cust)
    directed = D.merge_directed(*parts)
    expressed = None if args.demo else D.load_expressed(args.tissue, args.expressed_file)
    G = Wk.build_graph(directed, expressed)
    print(f"[graph] {G.number_of_nodes()} nodes, {G.number_of_edges()} directed edges" + (f" (filtered to {args.tissue})" if args.tissue else ""))
    return G, und


def run_dir(args, tag):
    """One folder per analysis: <out_root>/<YYYYMMDD_HHMMSS>_<tag>/ ; reused within a run (api sets args.run_dir once)."""
    import datetime
    if getattr(args, "run_dir", None):
        return args.run_dir
    root = getattr(args, "out", None) or "runs"
    d = os.path.join(root, datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + tag)
    os.makedirs(d, exist_ok=True); args.run_dir = d
    with open(os.path.join(d, "run_info.json"), "w") as f:
        json.dump({k: v for k, v in vars(args).items() if not k.startswith("_") and k != "run_dir"} | {"run_dir": d, "started": datetime.datetime.now().isoformat(timespec="seconds")},
                  f, indent=1, default=str)
    print(f"[run] {d}")
    return d


def write_html(tpl, token, payload, path, d3=False):
    """Write the view payload as JSON and the standalone HTML view.
    The HTML is the same Plotly figure the notebooks show (plotly.js from CDN), so file and notebook never disagree;
    with d3=True the original D3 page is written as well (<name>_d3.html: force layout, draggable)."""
    jp = os.path.splitext(path)[0] + ".json"
    with open(jp, "w") as f:
        json.dump(payload, f)
    try:
        from targetnet import plots as P
        kind = os.path.basename(jp)
        if kind.startswith("hop1"):
            fig = P.flow_figure(jp, max_nodes=10**6, sub_max=6)
        elif kind.startswith("cell_flow"):
            fig = P.cell_flow_figure(jp)
        else:
            fn = P.flow_figure if kind.startswith(("flow", "core_impact")) else P.cell_map_figure if kind.startswith("cell_map") else \
                 P.crosstalk_figure if kind.startswith("crosstalk") else P.volcano_figure
            fig = fn(jp)
        fig.write_html(path, include_plotlyjs="cdn", full_html=True)
    except Exception as e:   # plotly missing or figure error: fall back to the D3 page so a view always exists
        print(f"[view] Plotly export failed ({e.__class__.__name__}: {e}); writing D3 page instead"); d3 = True
        open(path, "w").write(open(os.path.join(TPL, tpl)).read().replace(token, json.dumps(payload)))
    if d3:
        open(os.path.splitext(path)[0] + "_d3.html", "w").write(open(os.path.join(TPL, tpl)).read().replace(token, json.dumps(payload)))


def annotate(df, G, gene, args, und):
    """Add co-essentiality, physical contact with seed, degree, compartment to a walk table."""
    cand = set(df.node) | {gene}
    coess, _ = (D.demo_depmap(gene, cand) if args.demo else D.load_depmap(args.depmap, gene, cand)) if (args.demo or args.depmap) else ({}, None)
    phys = set()
    try:
        st = D.load_string(cand, demo=args.demo)
        st = pd.concat([st, und[und.a.isin(cand) & und.b.isin(cand)]], ignore_index=True) if len(und) else st
        for r in st.itertuples():
            if r.a == gene: phys.add(r.b)
            if r.b == gene: phys.add(r.a)
    except Exception as e:
        print(f"[STRING] skipped: {e}")
    comp = D.load_compartments(cand, args.demo)
    deg = dict(G.degree())
    df = df.copy()
    import networkx as nx
    anc = nx.ancestors(G, gene) if gene in G else set(); desc = nx.descendants(G, gene) if gene in G else set()
    df["feedback"] = df.node.map(lambda n: n in anc and n in desc)     # on a loop through the seed
    df["coessentiality_r"] = df.node.map(lambda n: round(coess[n], 3) if n in coess else "")
    df["physical_with_seed"] = df.node.isin(phys)
    df["degree"] = df.node.map(deg)
    df["compartment"] = df.node.map(lambda n: comp.get(n, "cytosol"))
    comp[gene] = comp.get(gene, "cytosol")
    return df, comp


def seed_partners(gene, G, und, df, out, receptors=None, ligands=None, quiet=False):
    """
    Physical / complex partners of the seed that the walk does NOT follow (undirected Reactome FI / STRING context edges),
    e.g. ACVRL1 - BMPR2 / ACVR2A / ENG / TGFBR2. For each partner: is it already in the walk, where does it sit, and what signed
    output does it have in the graph. Written to partners.csv as a hypothesis table (complex-mediated routes to check by hand or
    to encode in custom_edges); nothing is added to the walk automatically, because a partner's output is not the seed's output
    (a type II receptor phosphorylates the type I receptor, which then signals; borrowing TGFBR2's SMAD2/3 output for ACVRL1 would be wrong).
    """
    if not len(und): return pd.DataFrame()
    ps = set(und.b[und.a == gene]) | set(und.a[und.b == gene]); ps.discard(gene)
    if not ps: return pd.DataFrame()
    inwalk = set(df.node); receptors = receptors or set(); ligands = ligands or set()
    comp = D.load_compartments(ps | {gene})
    rows = []
    for p in sorted(ps):
        outs = [(v, d["effect"]) for _, v, d in G.out_edges(p, data=True)] if p in G else []
        ins = [(u, d["effect"]) for u, _, d in G.in_edges(p, data=True)] if p in G else []
        rows.append(dict(partner=p, in_walk=p in inwalk, role="receptor" if p in receptors else "ligand" if p in ligands else "",
                         compartment=comp.get(p, ""), direct_edge_with_seed=G.has_edge(p, gene) or G.has_edge(gene, p),
                         n_signed_out=len(outs), n_signed_in=len(ins),
                         top_outputs=", ".join(f"{v}{'+' if e == 1 else '-' if e == -1 else '?'}" for v, e in sorted(outs, key=lambda t: -G[p][t[0]]['weight'])[:6])))
    tab = pd.DataFrame(rows).sort_values(["in_walk", "n_signed_out"], ascending=[False, False])
    tab.to_csv(os.path.join(out, "partners.csv"), index=False)
    if not quiet:
        mem = tab[(tab.compartment == "membrane") | (tab.role != "")]
        print(f"[partners] {len(tab)} physical/complex partners of {gene} in the context layer ({int(tab.in_walk.sum())} already in the walk); "
              f"{len(mem)} at the membrane -> partners.csv (complex-mediated routes are NOT walked; see top_outputs)")
        if len(mem): print(mem.head(8)[["partner", "in_walk", "role", "n_signed_out", "top_outputs"]].to_string(index=False))
    return tab


def cell_map_payload(gene, df, G, comp, mode, und):
    nodes = [dict(id=gene, cls="seed", compartment=comp.get(gene, "cytosol"), score=1.0, sign=1, hops=0)]
    for r in df.itertuples():
        cls = (("ligand" if getattr(r, "endpoint_type", "") == "ligand" else "receptor") if mode == "upstream" else "tf") if r.is_endpoint else mode
        nodes.append(dict(id=r.node, cls=cls, compartment=r.compartment, score=float(r.score), sign=int(r.sign), hops=int(r.hops),
                          note=f"best path: {r.best_path}"))
    keep = {n["id"] for n in nodes}
    edges = [dict(source=a, target=b, kind="signaling", effect=d["effect"]) for a, b, d in G.edges(data=True) if a in keep and b in keep]
    for r in und.itertuples():
        if r.a in keep and r.b in keep and not G.has_edge(r.a, r.b) and not G.has_edge(r.b, r.a):
            edges.append(dict(source=r.a, target=r.b, kind="physical", effect=0))
    return dict(seed=gene, mode=mode, compartments=COMPARTMENTS, nodes=nodes, edges=edges)


# ----------------------------------------------------------------------------- mode: upstream
def _stop_if_empty(df, gene, direction, args, exclude=()):
    """A walk that reaches nothing cannot be annotated or drawn: stop with the reason and what to change."""
    if len(df): return
    msg = (f"{gene} reaches no node {direction} within max_depth {args.max_depth} at min_score {args.min_score}"
           + (f" once {', '.join(exclude)} is removed from the graph" if exclude else "") + ".\n"
           "  Likely causes: the seed has few outgoing edges in SIGNOR (check it is the intended symbol / not a complex member), "
           "or its only neighbours are hubs whose penalty pushes the first hop below min_score.\n"
           "  Try: --min-score 0.2 (or crosstalk: min_score in target_net.yaml), --hub-penalty 0.1, a larger --max-depth; "
           "run `python target_net.py downstream --gene " + gene + "` alone to see how far it gets."
           + (f"\n  If {gene} reaches TFs only through {', '.join(exclude)}, that dependency is itself the cross-talk result: "
              f"{gene}'s whole downstream output passes through {', '.join(exclude)}." if exclude else ""))
    print("[walk] " + msg); raise SystemExit(msg)


def run_upstream(args):
    gene = resolve_symbol(args.gene); out = os.path.join(run_dir(args, f"upstream_{gene}"), "upstream"); os.makedirs(out, exist_ok=True)
    G, und = build_signaling(args)
    receptors = D.load_receptors(args.demo); ligands = D.load_ligands(args.demo) - receptors
    Wk.prepare_seed(G, gene, "up", D.CFG, D.REACTOME_DIRECTED)
    df = Wk.signed_walk(G, gene, receptors | ligands, direction="up", max_depth=args.max_depth, min_score=args.min_score, hub_strength=args.hub_penalty, extend_to=ligands)
    _stop_if_empty(df, gene, "upstream", args)
    df, comp = annotate(df, G, gene, args, und)
    df["endpoint_type"] = ["ligand" if n in ligands else "receptor" if n in receptors else "" for n in df.node]
    df.to_csv(os.path.join(out, "upstream_nodes.csv"), index=False)
    rec = df[df.endpoint_type == "receptor"].sort_values("score", ascending=False); rec.to_csv(os.path.join(out, "receptors.csv"), index=False)
    lig = df[df.endpoint_type == "ligand"].sort_values("score", ascending=False); lig.to_csv(os.path.join(out, "ligands.csv"), index=False)
    seed_partners(gene, G, und, df, out, receptors, ligands)
    flow = dict(mode="upstream", seed=gene, endpoint_label="receptors", regulon="", genes=[],
                subtitle=f"{len(df)} nodes, {len(rec)} receptors and {len(lig)} ligands reached within {args.max_depth} hops (path score ≥ {args.min_score}). Drawn right-to-left: ligands and receptors on the left, {gene} on the right.",
                nodes=[dict(id=r.node, hops=int(r.hops), score=float(r.score), sign=int(r.sign), is_endpoint=bool(r.is_endpoint), endpoint_type=r.endpoint_type, best_path=r.best_path, feedback=bool(r.feedback),
                            coess=r.coessentiality_r, compartment=r.compartment, cls="node") for r in df.itertuples()] +
                      [dict(id=gene, hops=0, score=1.0, sign=1, is_endpoint=False, best_path=gene, cls="seed")],
                edges=_edges(G, set(df.node) | {gene}))
    write_html("flow.html", "/*__FLOW_JSON__*/", flow, os.path.join(out, "flow.html"), d3=getattr(args, "d3", False))
    write_html("cell_map.html", "/*__NETWORK_JSON__*/", cell_map_payload(gene, df, G, comp, "upstream", und), os.path.join(out, "cell_map.html"), d3=getattr(args, "d3", False))
    write_hop1_view(flow, out)
    print(f"\n[upstream] {len(df)} nodes, {len(rec)} receptors, {len(lig)} ligands -> {out}/")
    return df
    print(rec.head(12)[["node", "hops", "score", "influence", "sign", "sign_confidence", "best_path"]].to_string(index=False))


# ----------------------------------------------------------------------------- mode: downstream
def run_downstream(args, gene=None, quiet=False, exclude=(), G=None, und=None):
    gene = resolve_symbol(gene or args.gene)
    out = os.path.join(run_dir(args, f"downstream_{gene}"), ("downstream" if not exclude else f"downstream_{gene}_without_" + "_".join(exclude))); os.makedirs(out, exist_ok=True)
    if G is None: G, und = build_signaling(args)
    if exclude:   # cross-talk: influence of this gene that does NOT pass through the other gene
        G = G.copy(); G.remove_nodes_from([e for e in exclude if e in G])
    tfs = D.load_tfs(args.demo)
    Wk.prepare_seed(G, gene, "down", D.CFG, D.REACTOME_DIRECTED)
    through_tf = bool(D.CFG.get("walk", {}).get("through_tf", True))
    df = Wk.signed_walk(G, gene, tfs, direction="down", max_depth=args.max_depth, min_score=args.min_score, hub_strength=args.hub_penalty,
                        stop_at_endpoints=not through_tf)
    _stop_if_empty(df, gene, "downstream", args, exclude)
    regulon = D.load_regulon(args.regulon, args.conf, args.regulon_file, args.demo)
    genes = Wk.extend_regulon(df, regulon, tfs, args.max_targets)
    df, comp = annotate(df, G, gene, args, und)
    df.to_csv(os.path.join(out, "downstream_nodes.csv"), index=False)
    if not exclude: seed_partners(gene, G, und, df, out, quiet=quiet)
    tf_df = df[df.is_endpoint].sort_values("score", ascending=False); tf_df.to_csv(os.path.join(out, "transcription_factors.csv"), index=False)
    if len(genes): genes.to_csv(os.path.join(out, "target_genes.csv"), index=False)
    flow = dict(mode="downstream", seed=gene, endpoint_label="transcription factors", regulon=f"{args.regulon} {args.conf}", through_tf=through_tf,
                subtitle=f"{len(df)} nodes, {len(tf_df)} TFs reached within {args.max_depth} hops; {genes.gene.nunique() if len(genes) else 0} target genes via {args.regulon} ({args.conf}).",
                genes=json.loads(genes.to_json(orient="records")) if len(genes) else [],
                nodes=[dict(id=r.node, hops=int(r.hops), score=float(r.score), sign=int(r.sign), is_endpoint=bool(r.is_endpoint), best_path=r.best_path, feedback=bool(r.feedback),
                            coess=r.coessentiality_r, compartment=r.compartment, cls="node") for r in df.itertuples()] +
                      [dict(id=gene, hops=0, score=1.0, sign=1, is_endpoint=False, best_path=gene, cls="seed")],
                edges=_edges(G, set(df.node) | {gene}))
    write_html("flow.html", "/*__FLOW_JSON__*/", flow, os.path.join(out, "flow.html"), d3=getattr(args, "d3", False))
    write_html("cell_map.html", "/*__NETWORK_JSON__*/", cell_map_payload(gene, df, G, comp, "downstream", und), os.path.join(out, "cell_map.html"), d3=getattr(args, "d3", False))
    write_hop1_view(flow, out)
    if not quiet:
        print(f"\n[downstream] {len(df)} nodes, {len(tf_df)} TFs, {genes.gene.nunique() if len(genes) else 0} genes -> {out}/")
        print(tf_df.head(12)[["node", "hops", "score", "influence", "sign", "sign_confidence", "best_path"]].to_string(index=False))
    return G, df, genes


# ----------------------------------------------------------------------------- mode: cross-talk screen (core vs a list of targets)
def run_crosstalk_screen(args):
    from targetnet import screen as Sc
    C = resolve_symbol(args.causal)
    targets = [resolve_symbol(t) for t in Sc.parse_targets(args.goi_list)]
    out = os.path.join(run_dir(args, f"crosstalk_screen_{C}_{len(targets)}targets"), "screen"); os.makedirs(out, exist_ok=True)
    args._extra_seeds = targets
    G, und = build_signaling(args)
    xt = D.CFG.get("crosstalk", {})
    if not getattr(args, "_explicit_depth", False): args.max_depth = xt.get("max_depth", 3)
    if not getattr(args, "_explicit_score", False): args.min_score = xt.get("min_score", 0.3)
    tfs = D.load_tfs(args.demo); regulon = D.load_regulon(args.regulon, args.conf, args.regulon_file, args.demo)
    print(f"[screen] core {C} vs {len(targets)} targets; walks max_depth {args.max_depth}, min_score {args.min_score}; baseline {args.n_random} random genes")
    res = Sc.screen(G, C, targets, tfs, regulon, max_depth=args.max_depth, min_score=args.min_score, hub=args.hub_penalty,
                    min_conf=xt.get("min_sign_confidence", 0.5), max_targets=args.max_targets, n_random=args.n_random, seed=args.seed)
    tab = res["table"]; tab.attrs["core"] = C
    tab.to_csv(os.path.join(out, "screen.csv"), index=False)
    if len(res["edges"]): res["edges"].to_csv(os.path.join(out, "screen_edges.csv"), index=False)
    if len(res["baseline"]): res["baseline"].to_csv(os.path.join(out, "baseline_random.csv"), index=False)
    json.dump(dict(core=C, targets=targets, core_nodes=res["core_nodes"], skipped=res["skipped"], max_depth=args.max_depth, min_score=args.min_score),
              open(os.path.join(out, "screen_info.json"), "w"), indent=1)
    if len(tab):
        from targetnet import plots as P
        P.screen_figure(tab).write_html(os.path.join(out, "screen.html"), include_plotlyjs="cdn")
        print(f"\n[screen] ranked {len(tab)} targets -> {out}/screen.csv")
        print(tab.head(15)[["rank", "target", "interference_score", "n_interference", "tf_enrichment", "tf_p", "dependency", "net_direction", "baseline_z", "composite"]].to_string(index=False))
    return tab


# ----------------------------------------------------------------------------- mode: pathway interference (CORE network + TARGET inhibition)
def run_pathway_interference(args):
    """1) CORE's full downstream network (receptor -> TFs -> genes) as a cell cross-section;
       2) how that flow changes when TARGET is inhibited (entry edges from TARGET's side, propagated inside the CORE network)."""
    from targetnet import interference as I
    C, X = resolve_symbol(args.causal), resolve_symbol(args.goi)
    root = run_dir(args, f"pathway_{C}_vs_{X}"); out = os.path.join(root, "interference"); os.makedirs(out, exist_ok=True)
    G, und = build_signaling(args)
    # step 1: CORE network on the full graph (writes downstream/ with flow, cell_map, hop1, tables)
    G, dC, gC = run_downstream(args, gene=C, quiet=True, G=G, und=und)
    cell = json.load(open(os.path.join(root, "downstream", "flow.json")))
    for n in cell["nodes"]: n.setdefault("compartment", "cytosol")
    cell["subtitle"] = f"{len(dC)} nodes, {int(dC.is_endpoint.sum())} TFs within {args.max_depth} hops (min_score {args.min_score}); {gC.gene.nunique() if len(gC) else 0} DoRothEA genes"
    write_html("flow.html", "/*__FLOW_JSON__*/", cell, os.path.join(out, "cell_flow.html"), d3=False)
    # step 2: TARGET's reach (CORE removed, so only TARGET's own routes count) and the inhibition perturbation
    tfs = D.load_tfs(args.demo)
    H = G.copy(); H.remove_nodes_from([C] if C in H else [])
    Wk.prepare_seed(H, X, "down", D.CFG, D.REACTOME_DIRECTED)
    dX = Wk.signed_walk(H, X, tfs, direction="down", max_depth=args.max_depth, min_score=args.min_score, hub_strength=args.hub_penalty,
                        stop_at_endpoints=not bool(D.CFG.get("walk", {}).get("through_tf", True))) if X in H else Wk.EMPTY_WALK.copy()
    res = I.interference_flow(G, dC, C, X, dX, genes_df=gC, max_depth=args.max_depth, hub_strength=args.hub_penalty)
    res["nodes"].to_csv(os.path.join(out, "flow_change.csv"), index=False)
    res["edges"].to_csv(os.path.join(out, "entry_edges.csv"), index=False)
    res["genes"].to_csv(os.path.join(out, "flow_change_genes.csv"), index=False)
    json.dump(res["summary"], open(os.path.join(out, "summary.json"), "w"), indent=1)
    pert = dict(target=X, summary=res["summary"],
                nodes={r.node: dict(delta=float(r.delta), direction=r.direction, entry=r.entry, via=r.via) for r in res["nodes"].itertuples()},
                genes={r.gene: float(r.delta) for r in res["genes"].itertuples()} if len(res["genes"]) else {},
                edges=json.loads(res["edges"].to_json(orient="records")) if len(res["edges"]) else [])
    cellp = dict(cell, perturbation=pert)
    write_html("flow.html", "/*__FLOW_JSON__*/", cellp, os.path.join(out, "cell_flow_inhibition.html"), d3=False)
    s_ = res["summary"]
    print(f"\n[pathway] {C} network: {s_['core_nodes']} nodes, {s_['core_tfs']} TFs -> {root}/downstream/ and {out}/cell_flow.html")
    print(f"[interference] inhibiting {X}: {s_['entry_points']} entry points ({s_['entry_edges']} edges from {X}'s side); "
          f"{s_['nodes_changed']} nodes change ({s_['nodes_decreased']} down / {s_['nodes_increased']} up); TFs {s_['tfs_changed']}/{s_['core_tfs']} "
          f"({s_['tfs_decreased']} down / {s_['tfs_increased']} up); genes {s_['genes_down']} down / {s_['genes_up']} up; flow-change index {s_['flow_change_index']}"
          + (f"; {X} is itself on the {C} network" if s_['target_in_core_network'] else ""))
    if len(res["nodes"]):
        print(res["nodes"][res["nodes"].direction != "unchanged"].head(12)[["node", "hops_from_core", "is_tf", "delta", "direction", "entry", "via"]].to_string(index=False))
    return res


# ----------------------------------------------------------------------------- calibrate: node counts vs min_score (walk only)
def run_calibrate(args):
    """Walk the seed up and down at several min_score / hub_penalty values (no annotation, no DepMap) and print the sizes,
    so walk.min_score can be set on the actual graph (SIGNOR + OmniPath + ...) instead of guessed."""
    gene = resolve_symbol(args.gene)
    G, und = build_signaling(args)
    receptors = D.load_receptors(args.demo); ligands = D.load_ligands(args.demo) - receptors; tfs = D.load_tfs(args.demo)
    Wk.prepare_seed(G, gene, "up", D.CFG, D.REACTOME_DIRECTED); Wk.prepare_seed(G, gene, "down", D.CFG, D.REACTOME_DIRECTED)
    scores = [float(x) for x in args.scores.split(",")]; hubs = [float(x) for x in args.hubs.split(",")]
    print(f"\n[calibrate] {gene}: graph {G.number_of_nodes()} nodes / {G.number_of_edges()} edges; max_depth {args.max_depth}")
    print(f"{'hub_penalty':>11} {'min_score':>9} | {'up nodes':>8} {'receptors':>9} {'ligands':>7} {'<=2 hops':>8} | {'down nodes':>10} {'TFs':>5} {'<=2 hops':>8}")
    for h in hubs:
        for ms in scores:
            up = Wk.signed_walk(G, gene, receptors | ligands, direction="up", max_depth=args.max_depth, min_score=ms, hub_strength=h, extend_to=ligands)
            dn = Wk.signed_walk(G, gene, tfs, direction="down", max_depth=args.max_depth, min_score=ms, hub_strength=h)
            nrec = int(up.node.isin(receptors).sum()); nlig = int(up.node.isin(ligands).sum())
            print(f"{h:>11.2f} {ms:>9.2f} | {len(up):>8} {nrec:>9} {nlig:>7} {int((up.hops <= 2).sum()):>8} | {len(dn):>10} {int(dn.is_endpoint.sum()):>5} {int((dn.hops <= 2).sum()):>8}")
    print("\nPick the min_score where the receptor list is a few dozen textbook inputs (not hundreds) and set walk.min_score in target_net.yaml.")


# ----------------------------------------------------------------------------- mode: synthetic lethality
def run_sl(args):
    gene = resolve_symbol(args.gene); out = os.path.join(run_dir(args, f"sl_{gene}"), "sl"); os.makedirs(out, exist_ok=True)
    if args.demo:
        coess, full = D.demo_depmap(gene, set())
    else:
        if not args.depmap: sys.exit("sl mode needs --depmap <dir> (see fetch_depmap.py)")
        coess, full = D.load_depmap(args.depmap, gene, set())
    if not len(full): sys.exit(D.SL_NOTE["text"] or "no SL table produced")
    # which hits are already on the pathway? (run a quick downstream+upstream walk to mark them)
    on_path = set()
    try:
        G, _ = build_signaling(args)
        on_path = set(Wk.signed_walk(G, gene, D.load_tfs(args.demo), "down", args.max_depth, args.min_score, args.hub_penalty).node) | \
                  set(Wk.signed_walk(G, gene, D.load_receptors(args.demo), "up", args.max_depth, args.min_score, args.hub_penalty).node)
    except SystemExit:
        pass
    full["on_pathway"] = full.gene.isin(on_path); full["coessentiality_r"] = full.gene.map(lambda g: round(coess.get(g, float("nan")), 3))
    full.sort_values("delta").to_csv(os.path.join(out, "synthetic_lethality.csv"), index=False)
    sig = full[full.significant].sort_values("delta")
    payload = dict(seed=gene, demo=bool(args.demo), thresholds=dict(padj=0.05, delta=0.15), n_lost=int(full.n_lost.iloc[0]),
                   points=json.loads(full[["gene", "delta", "padj", "n_lost", "significant", "on_pathway", "coessentiality_r"]].to_json(orient="records")),
                   top=json.loads(sig.head(30).to_json(orient="records")))
    write_html("sl_plot.html", "/*__SL_JSON__*/", payload, os.path.join(out, "synthetic_lethality.html"), d3=getattr(args, "d3", False))
    print(f"\n[sl] {len(sig)} hits -> {out}/"); print(sig.head(10)[["gene", "delta", "padj", "n_lost", "on_pathway"]].to_string(index=False))


# ----------------------------------------------------------------------------- mode: crosstalk
def run_crosstalk(args):
    C, X = resolve_symbol(args.causal), resolve_symbol(args.goi)
    out = os.path.join(run_dir(args, f"crosstalk_{C}_x_{X}"), "crosstalk"); os.makedirs(out, exist_ok=True)
    G, und = build_signaling(args)
    xt = D.CFG.get("crosstalk", {})            # local walks for cross-talk (unless the user overrode on the command line)
    if not getattr(args, "_explicit_depth", False): args.max_depth = xt.get("max_depth", 3)
    if not getattr(args, "_explicit_score", False): args.min_score = xt.get("min_score", 0.3)
    print(f"[crosstalk] walks: max_depth {args.max_depth}, min_score {args.min_score} (crosstalk: section of the config)")
    _, dC, gC = run_downstream(args, gene=C, quiet=True, exclude=(X,), G=G, und=und)   # C's pathway, X removed
    _, dX, gX = run_downstream(args, gene=X, quiet=True, exclude=(C,), G=G, und=und)   # X's pathway, C removed
    res = Wk.crosstalk(G, dC, dX, C, X, gC, gX, tf_universe=D.load_tfs(args.demo), min_conf=xt.get("min_sign_confidence", 0.5))
    for k in ("shared", "interference", "tf_agreement", "gene_agreement"):
        res[k].to_csv(os.path.join(out, f"{k}.csv"), index=False)
    json.dump(res["summary"], open(os.path.join(out, "summary.json"), "w"), indent=1)
    # payload for the two-lane view
    sC = dC.set_index("node"); sX = dX.set_index("node"); shared = set(res["shared"].node) if len(res["shared"]) else set()
    nodes = [dict(id=C, lane="causal", cls="seed", hops_c=0, sign_c=1, score=1.0), dict(id=X, lane="goi", cls="seed", hops_x=0, sign_x=1, score=1.0)]
    agree = dict(zip(res["shared"].node, res["shared"].agreement)) if len(res["shared"]) else {}
    for n in set(dC.index.map(lambda i: dC.loc[i, "node"])) | set(dX.node):
        lane = "shared" if n in shared else "causal" if n in sC.index else "goi"
        d = dict(id=n, lane=lane, cls="node", score=float(max(sC.score.get(n, 0), sX.score.get(n, 0))),
                 is_endpoint=bool(sC.is_endpoint.get(n, False) or sX.is_endpoint.get(n, False)), agreement=agree.get(n))
        if n in sC.index: d.update(hops_c=int(sC.hops[n]), sign_c=int(sC.sign[n]))
        if n in sX.index: d.update(hops_x=int(sX.hops[n]), sign_x=int(sX.sign[n]))
        nodes.append(d)
    endpoints = []
    for r in (res["tf_agreement"].itertuples() if len(res["tf_agreement"]) else []):
        endpoints.append(dict(name=r.node, layer="tf", sign_from_causal=int(r.sign_from_causal), sign_from_goi=int(r.sign_from_goi), agreement=r.agreement))
    for r in (res["gene_agreement"].itertuples() if len(res["gene_agreement"]) else []):
        endpoints.append(dict(name=r.gene, layer="gene", tf_from_causal=r.tf_from_causal, tf_from_goi=r.tf_from_goi,
                              sign_from_causal=int(r.sign_from_causal), sign_from_goi=int(r.sign_from_goi), agreement=r.agreement))
    # physical contacts (STRING experimental/database channels + Reactome undirected) between X / X's hop-1 nodes and C's pathway:
    # the evidence layer for chaperone-client and scaffold relations that signed databases rarely curate (e.g. HSP90 - FGFR3)
    xside = {X} | set(dX.node[dX.hops == 1]); cside = {C} | set(dC.node)
    contacts = []
    try:
        st = D.load_string(xside | cside, demo=args.demo)
        st = pd.concat([st, und[und.a.isin(xside | cside) & und.b.isin(xside | cside)]], ignore_index=True) if len(und) else st
        for r in st.drop_duplicates(["a", "b"]).itertuples():
            for u, v in ((r.a, r.b), (r.b, r.a)):
                if u in xside and v in cside and not G.has_edge(u, v) and not G.has_edge(v, u):
                    contacts.append(dict(goi_side=u, core_side=v, hops_goi=0 if u == X else 1, hops_causal=0 if v == C else int(dC.set_index("node").hops[v]), source=r.source))
    except Exception as e:
        print(f"[contacts] skipped: {e}")
    contacts = pd.DataFrame(contacts).drop_duplicates(["goi_side", "core_side"]) if contacts else pd.DataFrame(columns=["goi_side", "core_side", "hops_goi", "hops_causal", "source"])
    if len(contacts): contacts = contacts.sort_values(["hops_goi", "hops_causal"]).reset_index(drop=True)
    contacts.to_csv(os.path.join(out, "physical_contacts.csv"), index=False)
    direct = contacts[(contacts.hops_goi == 0) & (contacts.hops_causal == 0)]
    print(f"[contacts] {len(contacts)} physical contacts between {X} (+hop-1) and the {C} pathway without a signed edge"
          + (f"; DIRECT {X} - {C} contact in {', '.join(direct.source)}" if len(direct) else f"; no direct {X} - {C} contact in STRING/Reactome"))
    payload = dict(causal=C, goi=X, summary=res["summary"], nodes=nodes, edges=_edges(G, {n["id"] for n in nodes}),
                   interference=json.loads(res["interference"].to_json(orient="records")) if len(res["interference"]) else [], endpoints=endpoints,
                   contacts=json.loads(contacts.to_json(orient="records")) if len(contacts) else [])
    write_html("crosstalk.html", "/*__XT_JSON__*/", payload, os.path.join(out, "crosstalk.html"), d3=getattr(args, "d3", False))
    # impact of X on C's main pathway: table + the core flow view recoloured by impact
    imp = Wk.impact_on_core(G, dC, C, X, res); imp.to_csv(os.path.join(out, "impact.csv"), index=False)
    core_flow = json.load(open(os.path.join(run_dir(args, f"crosstalk_{C}_x_{X}"), f"downstream_{C}_without_{X}", "flow.json")))
    core_flow["impact"] = {r.node: dict(impact=r.impact, net=int(r.net_sign_from_goi), agreement=r.agreement, via=r.via) for r in imp.itertuples()}
    core_flow["impact_from"] = X
    write_html("flow.html", "/*__FLOW_JSON__*/", core_flow, os.path.join(out, "core_impact.html"), d3=False)
    ni = imp.impact.value_counts().to_dict()
    print(f"\n[impact] {X} on the {C} pathway: direct {ni.get('direct', 0)}, shared {ni.get('shared', 0)}, indirect {ni.get('indirect', 0)}, untouched {ni.get('none', 0)} -> {out}/impact.csv, core_impact.html")
    s = res["summary"]
    print(f"\n[crosstalk] {X} vs {C} pathway -> {out}/   (each pathway walked with the other gene removed)")
    print(f"  cross-talk index {s['crosstalk_index']}  (shared endpoints {s['n_shared_endpoints']} / causal endpoints {s['n_endpoints_causal']})")
    print(f"  TF overlap {s['tf_overlap']} of {s['tf_causal']} causal TFs; {X} reaches {s['tf_goi']}/{s['tf_universe']} TFs in the graph "
          f"({(s['goi_coverage_of_all_tfs'] or 0)*100:.0f}%) -> expected overlap {s['tf_overlap_expected']}, enrichment {s['tf_enrichment']}, p = {s['tf_overlap_p']:.2g}")
    print(f"  shared nodes {s['n_shared_nodes']}, interference edges {s['n_interference_edges']}, endpoints syn/ant/mixed/low-conf = "
          f"{s['synergistic']}/{s['antagonistic']}/{s['mixed']}/{s['low_confidence']}")
    if len(res["interference"]):
        print("  top interference edges (ranked by proximity to both genes):")
        print(res["interference"].head(15)[["rank", "from_node", "to_node", "hops_goi_to_from", "hops_causal_to_to", "mechanism", "residue", "net_effect_of_goi_on_target", "causal_effect_on_target", "agreement"]].to_string(index=False))


def apply_config(args):
    """Load the YAML config and fill CLI/API options that were not given explicitly."""
    cfg = Cf.load_config(getattr(args, "config", None)); D.set_config(cfg)
    print(f"[config] {cfg['_file'] or 'built-in defaults (download on demand)'}")
    w = cfg["walk"]
    if args.max_depth is None: args.max_depth = w["max_depth"]
    if args.min_score is None: args.min_score = w["min_score"]
    if args.hub_penalty is None: args.hub_penalty = w["hub_penalty"]
    if getattr(args, "depmap", None) is None: args.depmap = Cf.db(cfg, "depmap", "dir")
    if getattr(args, "regulon", None) is None: args.regulon = Cf.db(cfg, "regulon", "source") or "dorothea"
    if getattr(args, "conf", None) is None: args.conf = Cf.db(cfg, "regulon", "conf") or "ABC"
    if getattr(args, "reactome_fi", None) is None: args.reactome_fi = Cf.db(cfg, "reactome_fi")
    if getattr(args, "regulon_file", None) is None: args.regulon_file = Cf.db(cfg, "regulon")
    return cfg


def check_config(cfg):
    """Print, per database, what the loaders will actually use with this config."""
    import glob as _g
    def human(n): return f"{n/1e6:,.1f} MB" if n > 1e6 else f"{n/1e3:,.0f} kB"
    def status(path):
        if not path: return ("not set", "")
        pth = os.path.expanduser(path)
        if os.path.isdir(pth): return ("folder", "")
        return ("ok", human(os.path.getsize(pth))) if os.path.exists(pth) else ("MISSING", "")
    dl = cfg.get("allow_download", True)
    print(f"config file : {cfg['_file'] or 'none (built-in defaults)'}")
    print(f"cache_dir   : {cfg['cache_dir']}   allow_download: {dl}")
    rows = [
      ("signor", "walk graph (required)", Cf.db(cfg, "signor"), f"source={Cf.db(cfg,'signor','source') or 'file'}, indirect={Cf.db(cfg,'signor','indirect') or 'drop'}", "download from SIGNOR, else OmniPath mirror"),
      ("reactome_fi", "context layer (mode=context) / walk (mode=directed)", Cf.db(cfg, "reactome_fi"), f"mode={Cf.db(cfg,'reactome_fi','mode') or 'context'}", "download 2023 release"),
      ("uniprot", "receptors + compartments (bulk, offline)", Cf.db(cfg, "uniprot"), "", "per-gene UniProt API (slow)"),
      ("receptors", "optional override of the receptor list", Cf.db(cfg, "receptors"), "", "derived from uniprot"),
      ("ligands", "optional override of the secreted-ligand list", Cf.db(cfg, "ligands"), "", "derived from uniprot"),
      ("compartments", "optional override of compartments", Cf.db(cfg, "compartments"), "", "derived from uniprot"),
      ("tf_list", "downstream endpoints (required for downstream/crosstalk)", Cf.db(cfg, "tf_list"), "", "download Lambert 2018"),
      ("regulon", "TF -> target genes", Cf.db(cfg, "regulon"), f"source={Cf.db(cfg,'regulon','source')}, conf={Cf.db(cfg,'regulon','conf')}, license={Cf.db(cfg,'regulon','license')}", "download via OmniPath"),
      ("omnipath", f"gap filler for SIGNOR (mode={Cf.db(cfg,'omnipath','mode') or 'gaps'}, licence-filtered OmniPath)", Cf.db(cfg, "omnipath"), f"datasets={Cf.db(cfg,'omnipath','datasets')}", "download via OmniPath API"),
      ("custom_edges", "optional user-curated edges merged into the walk graph", Cf.db(cfg, "custom_edges"), "", "none (optional)"),
      ("hpa_tissue", "expression filter (--tissue only)", Cf.db(cfg, "hpa_tissue"), "", "download HPA consensus"),
    ]
    print(f"\n{'database':<14}{'used for':<52}{'status':<9}{'size':<10}path / options")
    for name, use, path, opts, fallback in rows:
        st, size = status(path)
        note = "" if st in ("ok", "folder") else (f"-> {fallback}" if dl else "-> STOP (allow_download is false)")
        print(f"{name:<14}{use:<52}{st:<9}{size:<10}{path or ''}  {opts}  {note}")
    d = Cf.db(cfg, "depmap", "dir")
    print(f"\ndepmap dir    : {d or 'not set'}")
    if d and os.path.isdir(os.path.expanduser(d)):
        for kind, use in [("effect", "co-essentiality (required for genetics)"), ("mutations", "SL: seed-lost lines by LoF mutation"), ("expression", "SL fallback: lowest-quartile expression")]:
            f = D.depmap_file(os.path.expanduser(d), kind)
            print(f"  {kind:<11}{use:<46}{'ok' if f else 'MISSING':<9}{human(os.path.getsize(f)) if f else '':<10}{os.path.basename(f) if f else ''}")
    print(f"\nstring        : API only (min_score {Cf.db(cfg,'string','min_score')}) -> {'used' if dl else 'skipped (offline)'}")
    print(f"walk          : {cfg['walk']}\ncrosstalk     : {cfg.get('crosstalk')}\nsl            : {cfg['sl']}\nview          : {cfg.get('view')}")


def signor_lookup(gene, cfg):
    """Print every SIGNOR row that mentions a gene (as entity A, entity B, or inside a complex / family name),
    so the user can see why a seed has or lacks edges in the walk graph (complexes and DIRECT=NO rows are not in it)."""
    import io
    loc = D._local("signor")
    if not loc: raise SystemExit("signor-lookup needs a local SIGNOR file (databases.signor.path)")
    df = pd.read_csv(io.StringIO(open(loc, encoding="utf-8", errors="replace").read()), sep="\t", dtype=str, quoting=3, on_bad_lines="skip")
    df.columns = [c.upper().replace(" ", "") for c in df.columns]
    g = gene.upper()
    hit = df[df.ENTITYA.str.upper().str.contains(rf"\b{g}\b", na=False) | df.ENTITYB.str.upper().str.contains(rf"\b{g}\b", na=False)]
    if "TAX_ID" in hit.columns: hit = hit[hit.TAX_ID.astype(str).str.strip().isin(["9606", "9606.0", "-1"])]
    cols = [c for c in ["ENTITYA", "TYPEA", "EFFECT", "MECHANISM", "RESIDUE", "ENTITYB", "TYPEB", "DIRECT", "SCORE"] if c in hit.columns]
    print(f"[SIGNOR] {len(hit)} human rows mention {g} ({loc})")
    if not len(hit): return
    print(hit[cols].to_string(index=False, max_colwidth=28))
    isprot = lambda t: t.astype(str).str.lower().str.strip() == "protein"
    asA = hit[hit.ENTITYA.str.upper() == g]; asB = hit[hit.ENTITYB.str.upper() == g]
    direct = hit.DIRECT.astype(str).str.upper().isin(["YES", "TRUE", "T", "1"]) if "DIRECT" in hit.columns else pd.Series(True, index=hit.index)
    print(f"\nin the walk graph (protein -> protein, DIRECT=YES): outgoing {int((isprot(asA.TYPEB) & direct[asA.index]).sum())}, incoming {int((isprot(asB.TYPEA) & direct[asB.index]).sum())}")
    other = hit[(hit.ENTITYA.str.upper() != g) & (hit.ENTITYB.str.upper() != g)]
    if len(other): print(f"rows where {g} appears only inside a complex / family name (dropped by the loader): {len(other)}  e.g. {sorted(set(other.ENTITYA) | set(other.ENTITYB))[:6]}")
    nd = hit[~direct]; 
    if len(nd): print(f"indirect (DIRECT=NO) rows: {len(nd)} (kept only with signor.indirect: weight|keep)")


def write_hop1_view(flow, out, d3=False):
    """Second view of the same walk: the seed and EVERY direct neighbour (hop 1), nothing else, no node cap.
    The full flow view is best-path oriented and capped; this one answers 'what touches my gene directly?' exhaustively."""
    keep = {n["id"] for n in flow["nodes"] if n.get("hops", 0) <= 1}
    h1 = dict(flow, hop1_only=True, nodes=[n for n in flow["nodes"] if n["id"] in keep],
              edges=[e for e in flow["edges"] if e["source"] in keep and e["target"] in keep], genes=[])
    h1["subtitle"] = f"every direct neighbour in the signed graph ({len(keep) - 1}); edge sign and weight as in the full view; endpoints (receptors / TFs) in their own column"
    write_html("flow.html", "/*__FLOW_JSON__*/", h1, os.path.join(out, "hop1.html"), d3=d3)


def _edges(G, keep):
    return [dict(source=a, target=b, effect=d["effect"], weight=d["weight"], residue=d["residue"], db=d["source"])
            for a, b, d in G.edges(data=True) if a in keep and b in keep]


# ----------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    def common(p):
        p.add_argument("--d3", action="store_true", help="also write the original D3 pages (<view>_d3.html, force layout, draggable)")
        p.add_argument("--out", default="runs", help="root folder for run outputs (default: runs/<YYYYMMDD_HHMMSS>_<mode>_<genes>/)")
        p.add_argument("--config", default=None, help="YAML with database paths (default: ./target_net.yaml, ~/.target_net.yaml, $TARGETNET_CONFIG)")
        p.add_argument("--demo", action="store_true", help="offline curated AKT1/MAPK1 axis (no network)")
        p.add_argument("--max-depth", type=int, default=None); p.add_argument("--min-score", type=float, default=None)
        p.add_argument("--hub-penalty", type=float, default=None, help="0 = none, 0.5 = strong down-weighting of hub nodes")
        p.add_argument("--tissue", default=None, help="HPA tissue name to filter expressed genes (e.g. liver)")
        p.add_argument("--expressed-file", default=None, help="text file, one expressed gene symbol per line (overrides --tissue)")
        p.add_argument("--no-reactome", action="store_true"); p.add_argument("--reactome-fi", default=None)
        p.add_argument("--depmap", default=None, help="DepMap directory (adds co-essentiality annotation)")
    def regulon_opts(p):
        p.add_argument("--regulon", default=None, help="dorothea | collectri (default from config)")
        p.add_argument("--conf", default=None, help="DoRothEA confidence levels to keep (default from config)")
        p.add_argument("--regulon-file", default=None, help="local TSV/CSV with columns tf,target,mor,confidence")
        p.add_argument("--max-targets", type=int, default=None, help="cap targets per TF (default: all)")
    p = sub.add_parser("upstream"); p.add_argument("--gene", required=True); common(p)
    p = sub.add_parser("downstream"); p.add_argument("--gene", required=True); common(p); regulon_opts(p)
    p = sub.add_parser("sl"); p.add_argument("--gene", required=True); common(p)
    p = sub.add_parser("pathway-interference", help="CORE network as a cell cross-section + predicted change when TARGET is inhibited")
    p.add_argument("--causal", required=True, help="CORE gene (usually a receptor)"); p.add_argument("--goi", required=True, help="TARGET gene (to be inhibited)")
    common(p); regulon_opts(p)
    p = sub.add_parser("calibrate", help="node / receptor / TF counts of the walk at several min_score and hub_penalty values (walk only, fast)")
    p.add_argument("--gene", required=True); p.add_argument("--scores", default="0.2,0.25,0.3,0.35,0.4,0.5"); p.add_argument("--hubs", default="0.25")
    common(p)
    p = sub.add_parser("check-config", help="show which file each database loader will use with the current YAML")
    p.add_argument("--config", default=None)
    p = sub.add_parser("signor-lookup", help="show every SIGNOR row mentioning a gene (why does the seed have / lack edges?)")
    p.add_argument("--gene", required=True); p.add_argument("--config", default=None)
    p = sub.add_parser("omnipath-url", help="print (or --save FILE) the OmniPath URL for the signed-interaction gap-filler layer")
    p.add_argument("--save", default=None); p.add_argument("--config", default=None)
    p = sub.add_parser("regulon-url", help="print (or --save FILE) the OmniPath URL for the configured regulon")
    p.add_argument("--save", default=None); p.add_argument("--config", default=None)
    p.add_argument("--source", default=None); p.add_argument("--conf", default=None); p.add_argument("--license", default=None)
    p = sub.add_parser("crosstalk"); p.add_argument("--causal", required=True); p.add_argument("--goi", required=True); common(p); regulon_opts(p)
    p = sub.add_parser("crosstalk-screen", help="rank a list of candidate target genes by their cross-talk with one core gene")
    p.add_argument("--causal", required=True); p.add_argument("--goi-list", required=True, help="genes: 'A B C', 'A,B,C' or a text/CSV file (first column)")
    p.add_argument("--n-random", type=int, default=20, help="degree-matched random genes for the baseline (0 = off)"); p.add_argument("--seed", type=int, default=0)
    common(p); regulon_opts(p)
    args = ap.parse_args()
    args._explicit_depth = getattr(args, "max_depth", None) is not None
    args._explicit_score = getattr(args, "min_score", None) is not None
    if args.mode == "check-config":
        cfg = Cf.load_config(args.config); D.set_config(cfg); check_config(cfg); return
    if args.mode == "signor-lookup":
        cfg = Cf.load_config(args.config); D.set_config(cfg); signor_lookup(args.gene, cfg); return
    if args.mode == "omnipath-url":
        cfg = Cf.load_config(args.config); D.set_config(cfg)
        datasets = Cf.db(cfg, "omnipath", "datasets") or "omnipath,pathwayextra,kinaseextra,ligrecextra"
        url = (D.OMNIPATH_URL + f"?datasets={datasets}&organisms=9606&genesymbols=1&directed=1&types=post_translational"
               "&fields=sources,references,type,is_directed,is_stimulation,is_inhibition")
        print(url)
        if args.save:
            import requests; os.makedirs(os.path.dirname(os.path.expanduser(args.save)) or ".", exist_ok=True)
            r = requests.get(url, timeout=900); r.raise_for_status(); open(os.path.expanduser(args.save), "w").write(r.text)
            D.cached_text("omnipath_resources", "https://omnipathdb.org/resources?format=json")   # licence list, cached next to it
            print(f"saved {args.save}: {len(r.text.splitlines()) - 1} rows (set databases.omnipath.path to this file)")
        return
    if args.mode == "regulon-url":
        cfg = Cf.load_config(args.config)
        url = D.omnipath_regulon_url(args.source or Cf.db(cfg, "regulon", "source") or "dorothea", args.conf or Cf.db(cfg, "regulon", "conf") or "ABC",
                                     args.license or Cf.db(cfg, "regulon", "license") or "commercial")
        print(url)
        if args.save:
            import requests; os.makedirs(os.path.dirname(os.path.expanduser(args.save)) or ".", exist_ok=True)
            r = requests.get(url, timeout=600); r.raise_for_status(); open(os.path.expanduser(args.save), "w").write(r.text)
            D.set_config(cfg); print(f"saved {args.save}: {len(D.load_regulon(path=os.path.expanduser(args.save)))} edges")
        return
    apply_config(args)
    try:
        {"upstream": run_upstream, "downstream": run_downstream, "sl": run_sl, "crosstalk": run_crosstalk, "crosstalk-screen": run_crosstalk_screen,
         "calibrate": run_calibrate, "pathway-interference": run_pathway_interference}[args.mode](args)
    except D.DownloadError as e:
        sys.exit(str(e))


if __name__ == "__main__":
    main()

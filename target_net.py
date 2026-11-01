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
                "AKT": ("AKT1", "family name"), "ERK2": ("MAPK1", "alias"), "ERK1": ("MAPK3", "alias"), "MEK1": ("MAP2K1", "alias"), "P53": ("TP53", "alias")}

def resolve_symbol(g):
    g = g.upper()
    if g in SYMBOL_NOTES:
        new, why = SYMBOL_NOTES[g]; print(f"[symbol] {g} -> {new}: {why}"); return new
    return g

TPL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "targetnet", "templates")
COMPARTMENTS = ["membrane", "cytosol", "lysosome", "mitochondrion", "ER_golgi", "nucleus"]


# ----------------------------------------------------------------------------- shared setup
def build_signaling(args):
    sig = D.load_signor(args.demo)
    parts = [sig]; und = pd.DataFrame(columns=["a", "b", "source"])
    if not args.no_reactome:
        rd, und = D.load_reactome_fi(args.reactome_fi, args.demo); parts.append(rd)
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


def write_html(tpl, token, payload, path):
    """Write the standalone D3 view and, next to it, the same payload as JSON (used by targetnet.plots for Plotly figures)."""
    html = open(os.path.join(TPL, tpl)).read().replace(token, json.dumps(payload))
    open(path, "w").write(html)
    with open(os.path.splitext(path)[0] + ".json", "w") as f:
        json.dump(payload, f)


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


def cell_map_payload(gene, df, G, comp, mode, und):
    nodes = [dict(id=gene, cls="seed", compartment=comp.get(gene, "cytosol"), score=1.0, sign=1, hops=0)]
    for r in df.itertuples():
        cls = ("receptor" if mode == "upstream" else "tf") if r.is_endpoint else mode
        nodes.append(dict(id=r.node, cls=cls, compartment=r.compartment, score=float(r.score), sign=int(r.sign), hops=int(r.hops),
                          note=f"best path: {r.best_path}"))
    keep = {n["id"] for n in nodes}
    edges = [dict(source=a, target=b, kind="signaling", effect=d["effect"]) for a, b, d in G.edges(data=True) if a in keep and b in keep]
    for r in und.itertuples():
        if r.a in keep and r.b in keep and not G.has_edge(r.a, r.b) and not G.has_edge(r.b, r.a):
            edges.append(dict(source=r.a, target=r.b, kind="physical", effect=0))
    return dict(seed=gene, mode=mode, compartments=COMPARTMENTS, nodes=nodes, edges=edges)


# ----------------------------------------------------------------------------- mode: upstream
def run_upstream(args):
    gene = resolve_symbol(args.gene); out = os.path.join(run_dir(args, f"upstream_{gene}"), "upstream"); os.makedirs(out, exist_ok=True)
    G, und = build_signaling(args)
    receptors = D.load_receptors(args.demo)
    df = Wk.signed_walk(G, gene, receptors, direction="up", max_depth=args.max_depth, min_score=args.min_score, hub_strength=args.hub_penalty)
    df, comp = annotate(df, G, gene, args, und)
    df.to_csv(os.path.join(out, "upstream_nodes.csv"), index=False)
    rec = df[df.is_endpoint].sort_values("score", ascending=False)
    rec.to_csv(os.path.join(out, "receptors.csv"), index=False)
    flow = dict(mode="upstream", seed=gene, endpoint_label="receptors", regulon="", genes=[],
                subtitle=f"{len(df)} nodes, {len(rec)} receptors reached within {args.max_depth} hops (path score ≥ {args.min_score}). Drawn right-to-left: receptors on the left, {gene} on the right.",
                nodes=[dict(id=r.node, hops=int(r.hops), score=float(r.score), sign=int(r.sign), is_endpoint=bool(r.is_endpoint), best_path=r.best_path, feedback=bool(r.feedback),
                            coess=r.coessentiality_r, compartment=r.compartment, cls="node") for r in df.itertuples()] +
                      [dict(id=gene, hops=0, score=1.0, sign=1, is_endpoint=False, best_path=gene, cls="seed")],
                edges=_edges(G, set(df.node) | {gene}))
    write_html("flow.html", "/*__FLOW_JSON__*/", flow, os.path.join(out, "flow.html"))
    write_html("cell_map.html", "/*__NETWORK_JSON__*/", cell_map_payload(gene, df, G, comp, "upstream", und), os.path.join(out, "cell_map.html"))
    print(f"\n[upstream] {len(df)} nodes, {len(rec)} receptors -> {out}/")
    return df
    print(rec.head(12)[["node", "hops", "score", "sign", "coessentiality_r", "best_path"]].to_string(index=False))


# ----------------------------------------------------------------------------- mode: downstream
def run_downstream(args, gene=None, quiet=False, exclude=(), G=None, und=None):
    gene = resolve_symbol(gene or args.gene)
    out = os.path.join(run_dir(args, f"downstream_{gene}"), ("downstream" if not exclude else f"downstream_{gene}_without_" + "_".join(exclude))); os.makedirs(out, exist_ok=True)
    if G is None: G, und = build_signaling(args)
    if exclude:   # cross-talk: influence of this gene that does NOT pass through the other gene
        G = G.copy(); G.remove_nodes_from([e for e in exclude if e in G])
    tfs = D.load_tfs(args.demo)
    df = Wk.signed_walk(G, gene, tfs, direction="down", max_depth=args.max_depth, min_score=args.min_score, hub_strength=args.hub_penalty)
    regulon = D.load_regulon(args.regulon, args.conf, args.regulon_file, args.demo)
    genes = Wk.extend_regulon(df, regulon, tfs, args.max_targets)
    df, comp = annotate(df, G, gene, args, und)
    df.to_csv(os.path.join(out, "downstream_nodes.csv"), index=False)
    tf_df = df[df.is_endpoint].sort_values("score", ascending=False); tf_df.to_csv(os.path.join(out, "transcription_factors.csv"), index=False)
    if len(genes): genes.to_csv(os.path.join(out, "target_genes.csv"), index=False)
    flow = dict(mode="downstream", seed=gene, endpoint_label="transcription factors", regulon=f"{args.regulon} {args.conf}",
                subtitle=f"{len(df)} nodes, {len(tf_df)} TFs reached within {args.max_depth} hops; {genes.gene.nunique() if len(genes) else 0} target genes via {args.regulon} ({args.conf}).",
                genes=json.loads(genes.to_json(orient="records")) if len(genes) else [],
                nodes=[dict(id=r.node, hops=int(r.hops), score=float(r.score), sign=int(r.sign), is_endpoint=bool(r.is_endpoint), best_path=r.best_path, feedback=bool(r.feedback),
                            coess=r.coessentiality_r, compartment=r.compartment, cls="node") for r in df.itertuples()] +
                      [dict(id=gene, hops=0, score=1.0, sign=1, is_endpoint=False, best_path=gene, cls="seed")],
                edges=_edges(G, set(df.node) | {gene}))
    write_html("flow.html", "/*__FLOW_JSON__*/", flow, os.path.join(out, "flow.html"))
    write_html("cell_map.html", "/*__NETWORK_JSON__*/", cell_map_payload(gene, df, G, comp, "downstream", und), os.path.join(out, "cell_map.html"))
    if not quiet:
        print(f"\n[downstream] {len(df)} nodes, {len(tf_df)} TFs, {genes.gene.nunique() if len(genes) else 0} genes -> {out}/")
        print(tf_df.head(12)[["node", "hops", "score", "sign", "best_path"]].to_string(index=False))
    return G, df, genes


# ----------------------------------------------------------------------------- mode: synthetic lethality
def run_sl(args):
    gene = resolve_symbol(args.gene); out = os.path.join(run_dir(args, f"sl_{gene}"), "sl"); os.makedirs(out, exist_ok=True)
    if args.demo:
        coess, full = D.demo_depmap(gene, set())
    else:
        if not args.depmap: sys.exit("sl mode needs --depmap <dir> (see fetch_depmap.py)")
        coess, full = D.load_depmap(args.depmap, gene, set())
    if not len(full): sys.exit("no SL table produced")
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
    write_html("sl_plot.html", "/*__SL_JSON__*/", payload, os.path.join(out, "synthetic_lethality.html"))
    print(f"\n[sl] {len(sig)} hits -> {out}/"); print(sig.head(10)[["gene", "delta", "padj", "n_lost", "on_pathway"]].to_string(index=False))


# ----------------------------------------------------------------------------- mode: crosstalk
def run_crosstalk(args):
    C, X = resolve_symbol(args.causal), resolve_symbol(args.goi)
    out = os.path.join(run_dir(args, f"crosstalk_{C}_x_{X}"), "crosstalk"); os.makedirs(out, exist_ok=True)
    G, und = build_signaling(args)
    _, dC, gC = run_downstream(args, gene=C, quiet=True, exclude=(X,), G=G, und=und)   # C's pathway, X removed
    _, dX, gX = run_downstream(args, gene=X, quiet=True, exclude=(C,), G=G, und=und)   # X's pathway, C removed
    res = Wk.crosstalk(G, dC, dX, C, X, gC, gX)
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
    payload = dict(causal=C, goi=X, summary=res["summary"], nodes=nodes, edges=_edges(G, {n["id"] for n in nodes}),
                   interference=json.loads(res["interference"].to_json(orient="records")) if len(res["interference"]) else [], endpoints=endpoints)
    write_html("crosstalk.html", "/*__XT_JSON__*/", payload, os.path.join(out, "crosstalk.html"))
    s = res["summary"]
    print(f"\n[crosstalk] {X} vs {C} pathway -> {out}/   (each pathway walked with the other gene removed)")
    print(f"  cross-talk index {s['crosstalk_index']}  (shared endpoints {s['n_shared_endpoints']} / causal endpoints {s['n_endpoints_causal']})")
    print(f"  shared nodes {s['n_shared_nodes']}, interference edges {s['n_interference_edges']}, endpoints syn/ant/mixed = {s['synergistic']}/{s['antagonistic']}/{s['mixed']}")
    if len(res["interference"]):
        print(res["interference"][["from_node", "to_node", "mechanism", "residue", "net_effect_of_goi_on_target", "causal_effect_on_target", "agreement", "via_feedback_loop"]].to_string(index=False))


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


def _edges(G, keep):
    return [dict(source=a, target=b, effect=d["effect"], weight=d["weight"], residue=d["residue"], db=d["source"])
            for a, b, d in G.edges(data=True) if a in keep and b in keep]


# ----------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    def common(p):
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
    p = sub.add_parser("regulon-url", help="print (or --save FILE) the OmniPath URL for the configured regulon")
    p.add_argument("--save", default=None); p.add_argument("--config", default=None)
    p.add_argument("--source", default=None); p.add_argument("--conf", default=None); p.add_argument("--license", default=None)
    p = sub.add_parser("crosstalk"); p.add_argument("--causal", required=True); p.add_argument("--goi", required=True); common(p); regulon_opts(p)
    args = ap.parse_args()
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
        {"upstream": run_upstream, "downstream": run_downstream, "sl": run_sl, "crosstalk": run_crosstalk}[args.mode](args)
    except D.DownloadError as e:
        sys.exit(str(e))


if __name__ == "__main__":
    main()

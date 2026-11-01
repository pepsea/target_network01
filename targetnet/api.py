"""
Notebook-friendly API for target_net.

    from targetnet.api import research_gene, research_crosstalk, show

    r = research_gene("AKT1", config="target_net.yaml") # upstream + downstream (+ sl if depmap dir is configured)
    r["receptors"], r["tfs"], r["target_genes"], r["sl"] # pandas DataFrames
    figure(r["views"]["downstream_flow"]).show()        # Plotly figure inline (small; plotly.js from CDN)
    show(r["views"]["downstream_flow"])                 # or the full D3 view in an IFrame

    x = research_crosstalk(core="FGFR3", target="HSP90AA1", demo=True)
    x["summary"], x["interference"], x["endpoints"]
    show(x["views"]["crosstalk"])
"""
import os, sys, json
from argparse import Namespace
import pandas as pd

_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import target_net as cli  # noqa: E402

DEFAULTS = dict(out="runs", config=None, demo=False, max_depth=None, min_score=None, hub_penalty=None, tissue=None, expressed_file=None,
                no_reactome=False, reactome_fi=None, depmap=None, regulon=None, conf=None, regulon_file=None, max_targets=None)


def _args(**kw):
    a = dict(DEFAULTS); a.update({k: v for k, v in kw.items() if k in DEFAULTS})
    unknown = set(kw) - set(DEFAULTS)
    if unknown: raise TypeError(f"unknown options: {unknown}")
    return Namespace(**a)


def _read(path):
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()


def research_gene(gene, run_sl=None, **options):
    """
    Single-gene research: upstream (to receptors), downstream (to TFs and target genes), and synthetic lethality
    when a DepMap directory is given (or demo=True). Options = CLI options as keywords (demo, depmap, tissue, max_depth, out, ...).
    Outputs go to <out>/<YYYYMMDD_HHMMSS>_gene_<GENE>/{upstream,downstream,sl}/ (out defaults to "runs").
    Returns dict of DataFrames + run_dir + paths of HTML views.
    """
    a = _args(**options); cli.apply_config(a); gene = cli.resolve_symbol(gene)
    a.gene = gene
    out = cli.run_dir(a, f"gene_{gene}")       # one timestamped folder for all three analyses
    cli.run_upstream(a); cli.run_downstream(a)
    if run_sl is None: run_sl = bool(a.demo or a.depmap)
    sl_note = ""
    if run_sl:
        try: cli.run_sl(a)
        except SystemExit as e: sl_note = str(e); print(f"[sl] skipped: {e}")
    return dict(
        gene=gene, run_dir=out, sl_note=sl_note,
        upstream_nodes=_read(f"{out}/upstream/upstream_nodes.csv"), receptors=_read(f"{out}/upstream/receptors.csv"),
        downstream_nodes=_read(f"{out}/downstream/downstream_nodes.csv"), tfs=_read(f"{out}/downstream/transcription_factors.csv"),
        target_genes=_read(f"{out}/downstream/target_genes.csv"), sl=_read(f"{out}/sl/synthetic_lethality.csv"),
        views=dict(upstream_flow=f"{out}/upstream/flow.html", upstream_cell_map=f"{out}/upstream/cell_map.html",
                   downstream_flow=f"{out}/downstream/flow.html", downstream_cell_map=f"{out}/downstream/cell_map.html",
                   sl=f"{out}/sl/synthetic_lethality.html"))


def research_crosstalk(core, target, **options):
    """
    Cross-talk research: does `target` (gene of interest) influence the pathway of `core` (e.g. disease-causal gene)?
    Each pathway is walked with the other gene removed. Outputs go to <out>/<YYYYMMDD_HHMMSS>_crosstalk_<CORE>_x_<TARGET>/.
    Returns summary, shared nodes, interference edges, endpoint agreement, run_dir, view paths.
    """
    a = _args(**options); a._explicit_depth = options.get("max_depth") is not None; a._explicit_score = options.get("min_score") is not None
    cli.apply_config(a); a.causal, a.goi = core, target
    C, X = cli.resolve_symbol(core), cli.resolve_symbol(target)
    root = cli.run_dir(a, f"crosstalk_{C}_x_{X}")
    cli.run_crosstalk(a)
    out = f"{root}/crosstalk"
    summary = json.load(open(f"{out}/summary.json"))
    tf, ge = _read(f"{out}/tf_agreement.csv"), _read(f"{out}/gene_agreement.csv")
    endpoints = pd.concat([tf.assign(layer="tf").rename(columns={"node": "name"}), ge.assign(layer="gene").rename(columns={"gene": "name"})], ignore_index=True) \
        if len(tf) or len(ge) else pd.DataFrame()
    return dict(core=C, target=X, run_dir=root, summary=summary, shared=_read(f"{out}/shared.csv"), interference=_read(f"{out}/interference.csv"),
                endpoints=endpoints, views=dict(crosstalk=f"{out}/crosstalk.html", core_flow=f"{root}/downstream_{C}_without_{X}/flow.html",
                                                target_flow=f"{root}/downstream_{X}_without_{C}/flow.html"))


def figure(view_path, focus=None, **kw):
    """Plotly figure for a view: pass the .html or .json path from result['views'] (e.g. figure(r['views']['downstream_flow'])).
    focus: genes to highlight — list of symbols, a text/CSV file (symbol[, value] per line), or dict/Series gene -> value
    (e.g. log2FC from a proteomics experiment). Other nodes are drawn as small grey dots; focus genes keep colour, size and label,
    and are coloured by value when values are given."""
    if focus is not None: kw["focus"] = focus
    from . import plots
    jp = os.path.splitext(view_path)[0] + ".json"
    if not os.path.exists(jp):
        print(f"(no view at {jp})"); return None
    kind = os.path.basename(jp)
    fn = plots.flow_figure if kind.startswith("flow") else plots.cell_map_figure if kind.startswith("cell_map") else \
         plots.crosstalk_figure if kind.startswith("crosstalk") else plots.volcano_figure
    if fn is plots.volcano_figure: kw.pop("focus", None)
    return fn(jp, **kw)


def show(html_path, height=650):
    """Render a saved HTML view inside the notebook."""
    from IPython.display import IFrame, display
    if not os.path.exists(html_path):
        print(f"(no view at {html_path})"); return
    display(IFrame(src=os.path.relpath(html_path), width="100%", height=height))


def summarize_gene(r, n=10):
    """Compact text summary for a research_gene() result."""
    lines = [f"== {r['gene']} =="]
    if len(r["receptors"]):
        lines.append(f"Upstream: {len(r['upstream_nodes'])} nodes, {len(r['receptors'])} receptors. Top: " +
                     ", ".join(f"{x.node}({x.hops}h,{'+' if x.sign==1 else '-' if x.sign==-1 else '±'})" for x in r["receptors"].head(n).itertuples()))
    if len(r["tfs"]):
        lines.append(f"Downstream: {len(r['downstream_nodes'])} nodes, {len(r['tfs'])} TFs, {r['target_genes'].gene.nunique() if len(r['target_genes']) else 0} target genes. Top TFs: " +
                     ", ".join(f"{x.node}({'↑' if x.sign==1 else '↓' if x.sign==-1 else '±'})" for x in r["tfs"].head(n).itertuples()))
    if len(r["sl"]):
        hits = r["sl"][r["sl"].significant]
        lines.append(f"Synthetic lethality: {len(hits)} hits. Top: " + ", ".join(hits.sort_values('delta').head(n).gene))
    elif r.get("sl_note"):
        lines.append("Synthetic lethality: not run — " + r["sl_note"].replace("[DepMap] ", ""))
    return "\n".join(lines)


def summarize_crosstalk(x, n=10):
    s = x["summary"]
    lines = [f"== does {x['target']} influence the {x['core']} pathway? ==",
             f"cross-talk index {s['crosstalk_index']} (shared endpoints {s['n_shared_endpoints']}/{s['n_endpoints_causal']})",
             f"TF overlap {s['tf_overlap']}/{s['tf_causal']} vs expected {s['tf_overlap_expected']} -> enrichment {s['tf_enrichment']}, p = {s['tf_overlap_p']:.2g} "
             f"({x['target']} reaches {(s['goi_coverage_of_all_tfs'] or 0)*100:.0f}% of all TFs in the graph)",
             f"shared nodes {s['n_shared_nodes']}; interference edges {s['n_interference_edges']}; endpoints syn/ant/mixed/low-conf "
             f"{s['synergistic']}/{s['antagonistic']}/{s['mixed']}/{s['low_confidence']}"]
    if len(x["interference"]):
        lines.append("top interference edges (by proximity):")
        for e in x["interference"].head(n).itertuples():
            lines.append(f"  #{e.rank} {e.from_node} -> {e.to_node} (hops {e.hops_goi_to_from}/{e.hops_causal_to_to}) [{e.residue or e.mechanism}] {e.agreement}"
                         + (" (feedback loop)" if e.via_feedback_loop else ""))
    return "\n".join(lines)

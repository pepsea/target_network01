"""
Notebook-friendly API for target_net.

    from targetnet.api import research_gene, research_crosstalk, show

    r = research_gene("AKT1", config="target_net.yaml") # upstream + downstream (+ sl if depmap dir is configured)
    r["receptors"], r["tfs"], r["target_genes"], r["sl"] # pandas DataFrames
    show(r["views"]["downstream_flow"])                 # render an HTML view inline (IFrame)

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

DEFAULTS = dict(config=None, demo=False, max_depth=None, min_score=None, hub_penalty=None, tissue=None, expressed_file=None,
                no_reactome=False, reactome_fi=None, depmap=None, regulon=None, conf=None, regulon_file=None, max_targets=None)


def _args(**kw):
    a = dict(DEFAULTS); a.update({k: v for k, v in kw.items() if k in DEFAULTS})
    unknown = set(kw) - set(DEFAULTS)
    if unknown: raise TypeError(f"unknown options: {unknown}")
    return Namespace(**a)


def _read(path):
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()


def research_gene(gene, run_sl=None, workdir=None, **options):
    """
    Single-gene research: upstream (to receptors), downstream (to TFs and target genes), and synthetic lethality
    when a DepMap directory is given (or demo=True). Options = CLI options as keywords (demo, depmap, tissue, max_depth, ...).
    Returns dict of DataFrames + paths of HTML views.
    """
    if workdir: os.makedirs(workdir, exist_ok=True); os.chdir(workdir)
    a = _args(**options); cli.apply_config(a); gene = cli.resolve_symbol(gene)
    a.gene = gene
    cli.run_upstream(a); cli.run_downstream(a)
    if run_sl is None: run_sl = bool(a.demo or a.depmap)
    if run_sl:
        try: cli.run_sl(a)
        except SystemExit as e: print(f"[sl] skipped: {e}")
    out = f"out_{gene}"
    return dict(
        gene=gene,
        upstream_nodes=_read(f"{out}/upstream/upstream_nodes.csv"), receptors=_read(f"{out}/upstream/receptors.csv"),
        downstream_nodes=_read(f"{out}/downstream/downstream_nodes.csv"), tfs=_read(f"{out}/downstream/transcription_factors.csv"),
        target_genes=_read(f"{out}/downstream/target_genes.csv"), sl=_read(f"{out}/sl/synthetic_lethality.csv"),
        views=dict(upstream_flow=f"{out}/upstream/flow.html", upstream_cell_map=f"{out}/upstream/cell_map.html",
                   downstream_flow=f"{out}/downstream/flow.html", downstream_cell_map=f"{out}/downstream/cell_map.html",
                   sl=f"{out}/sl/synthetic_lethality.html"))


def research_crosstalk(core, target, workdir=None, **options):
    """
    Cross-talk research: does `target` (gene of interest) influence the pathway of `core` (e.g. disease-causal gene)?
    Each pathway is walked with the other gene removed. Returns summary, shared nodes, interference edges, endpoint agreement.
    """
    if workdir: os.makedirs(workdir, exist_ok=True); os.chdir(workdir)
    a = _args(**options); cli.apply_config(a); a.causal, a.goi = core, target
    cli.run_crosstalk(a)
    C, X = cli.resolve_symbol(core), cli.resolve_symbol(target)
    out = f"out_{C}/crosstalk_{X}"
    summary = json.load(open(f"{out}/summary.json"))
    tf, ge = _read(f"{out}/tf_agreement.csv"), _read(f"{out}/gene_agreement.csv")
    endpoints = pd.concat([tf.assign(layer="tf").rename(columns={"node": "name"}), ge.assign(layer="gene").rename(columns={"gene": "name"})], ignore_index=True) \
        if len(tf) or len(ge) else pd.DataFrame()
    return dict(core=C, target=X, summary=summary, shared=_read(f"{out}/shared.csv"), interference=_read(f"{out}/interference.csv"),
                endpoints=endpoints, views=dict(crosstalk=f"{out}/crosstalk.html", core_flow=f"out_{C}/downstream_without_{X}/flow.html",
                                                target_flow=f"out_{X}/downstream_without_{C}/flow.html"))


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
    return "\n".join(lines)


def summarize_crosstalk(x, n=10):
    s = x["summary"]
    lines = [f"== does {x['target']} influence the {x['core']} pathway? ==",
             f"cross-talk index {s['crosstalk_index']} (shared endpoints {s['n_shared_endpoints']}/{s['n_endpoints_causal']}); "
             f"shared nodes {s['n_shared_nodes']}; interference edges {s['n_interference_edges']}; endpoints syn/ant/mixed {s['synergistic']}/{s['antagonistic']}/{s['mixed']}"]
    if len(x["interference"]):
        for e in x["interference"].head(n).itertuples():
            lines.append(f"  {e.from_node} -> {e.to_node} [{e.residue or e.mechanism}] {e.agreement}" + (" (feedback loop)" if e.via_feedback_loop else ""))
    return "\n".join(lines)

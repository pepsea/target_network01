"""
Cross-talk screen: one disease-causal core gene C against a LIST of candidate target genes X1..Xn.

For each X the same local cross-talk analysis as research_crosstalk is run (walk of C with X removed, walk of X with
C removed, shared nodes, interference edges, TF-overlap enrichment), then every X is scored on four independent axes
and ranked. Walks only (no DepMap / STRING / compartments), so a screen of ~50 genes takes a few minutes.

Score axes (all per target X, higher = stronger cross-talk with C):
  interference_score  sum over interference edges (X-only node -> C-pathway node) of edge_weight / (1 + proximity);
                      proximity = hops from X to the edge source + hops from C to the edge target. Rewards edges close to
                      both genes; a hub that touches C's pathway only 3 hops away on both sides scores little.
  tf_enrichment       observed / expected overlap of the TFs reached by X and by C (hypergeometric p in tf_p). Corrects for
                      X reaching a large share of all TFs.
  dependency          1 - |D(C without X)| / |D(C)|: fraction of C's downstream pathway that disappears when X is removed.
                      This is the "bypass" question: if C's output needs X, inhibiting X cuts C's signal.
  net_direction       (syn - ant) / (syn + ant) over interference edges (syn_edges / ant_edges) AND confidence-gated shared
                      endpoints (syn_endpoints / ant_endpoints; these are the counts research_crosstalk prints): +1 = X pushes C's pathway the same way, -1 = X opposes it. Not part of the rank
                      (direction is a property, not a magnitude) but reported next to it.
  baseline_z          z-score of interference_score against N degree-matched random genes run against the same C
                      (the automatic baseline). Answers "is this more than any gene of similar connectivity would show?".

composite_rank: mean of the percentile ranks of interference_score, tf_enrichment (0 when p > 0.5) and dependency, with
baseline_z used as a tie-breaker. Transparent on purpose: every component is in the table, change the weights if you must.
"""
import math, random
import numpy as np
import pandas as pd
from . import walk as Wk


def _light_downstream(G, gene, tfs, regulon, max_depth, min_score, hub, exclude=(), max_targets=None):
    """Walk + regulon only (no annotation layers). Returns (nodes_df, genes_df); empty frames if the seed cannot start."""
    from . import data as D
    if gene in G: Wk.prepare_seed(G, gene, "down", D.CFG, D.REACTOME_DIRECTED)   # identical seed handling to research_crosstalk
    H = G
    if exclude:
        H = G.copy(); H.remove_nodes_from([e for e in exclude if e in H])
    if gene not in H:
        return Wk.EMPTY_WALK.copy(), pd.DataFrame()
    df = Wk.signed_walk(H, gene, tfs, direction="down", max_depth=max_depth, min_score=min_score, hub_strength=hub)
    genes = Wk.extend_regulon(df, regulon, tfs, max_targets) if len(df) and regulon is not None and len(regulon) else pd.DataFrame()
    return df, genes


def _interference_score(inter_df):
    if inter_df is None or not len(inter_df): return 0.0
    return float((inter_df.weight / (1.0 + inter_df.proximity)).sum())


def _direction(res):
    """Direction over BOTH evidence kinds: interference edges (X-side node -> C-pathway node) and shared endpoints (TFs / genes).
    Returned separately so the numbers can be matched with research_crosstalk, whose summary counts endpoints only."""
    s = res["summary"]; inter = res["interference"]
    syn_e = int((inter.agreement == "synergistic").sum()) if len(inter) else 0
    ant_e = int((inter.agreement == "antagonistic").sum()) if len(inter) else 0
    syn_p, ant_p = s.get("synergistic", 0), s.get("antagonistic", 0)
    syn, ant = syn_e + syn_p, ant_e + ant_p
    return (syn - ant) / (syn + ant) if (syn + ant) else 0.0, syn_e, ant_e, syn_p, ant_p


def degree_matched_random(G, gene, n, rng, exclude=(), tol=0.5):
    """n random genes with out-degree within [d/(1+tol), d*(1+tol)] of `gene` (widened until enough candidates exist)."""
    d = G.out_degree(gene) if gene in G else 1
    pool = []
    t = tol
    while len(pool) < max(3 * n, 10) and t < 50:
        lo, hi = d / (1 + t), d * (1 + t)
        pool = [g for g, k in G.out_degree() if lo <= k <= hi and g != gene and g not in exclude]
        t *= 2
    return rng.sample(pool, min(n, len(pool))) if pool else []


def screen(G, core, targets, tfs, regulon, max_depth=3, min_score=0.3, hub=0.25, min_conf=0.5, max_targets=None,
           n_random=20, seed=0, top_edges=5, log=print):
    """
    Returns dict(table=ranked DataFrame (one row per target), edges=top interference edges per target (long table),
                 core_nodes=|D(C)|, baseline=DataFrame of random-gene scores).
    """
    rng = random.Random(seed)
    tfs = set(tfs)
    dC_full, _ = _light_downstream(G, core, tfs, regulon, max_depth, min_score, hub)
    n_full = len(dC_full)
    if n_full == 0:
        log(f"[screen] {core} reaches no node downstream at depth {max_depth} / min_score {min_score}: nothing to screen against")
    tf_universe = tfs

    def one(X):
        dC, gC = _light_downstream(G, core, tfs, regulon, max_depth, min_score, hub, exclude=(X,), max_targets=max_targets)
        dX, gX = _light_downstream(G, X, tfs, regulon, max_depth, min_score, hub, exclude=(core,), max_targets=max_targets)
        if len(dX) == 0:
            return None, None
        res = Wk.crosstalk(G, dC, dX, core, X, gC, gX, tf_universe=tf_universe, min_conf=min_conf)
        s = res["summary"]; nd, syn_e, ant_e, syn_p, ant_p = _direction(res)
        inter = res["interference"]
        row = dict(target=X, n_nodes_target=len(dX), n_tfs_target=int(dX.is_endpoint.sum()),
                   n_core_nodes_without_target=len(dC), dependency=round(1 - len(dC) / n_full, 3) if n_full else 0.0,
                   shared_nodes=s["n_shared_nodes"], shared_endpoints=s["n_shared_endpoints"],
                   tf_overlap=s["tf_overlap"], tf_expected=s["tf_overlap_expected"], tf_enrichment=s["tf_enrichment"] or 0.0,
                   tf_p=s["tf_overlap_p"], n_interference=len(inter), interference_score=round(_interference_score(inter), 3),
                   syn_edges=syn_e, ant_edges=ant_e, syn_endpoints=syn_p, ant_endpoints=ant_p, net_direction=round(nd, 2),
                   top_edges="; ".join(f"{r.from_node}->{r.to_node}({r.agreement[:3]})" for r in inter.head(3).itertuples()) if len(inter) else "")
        edges = inter.head(top_edges).assign(target=X) if len(inter) else None
        return row, edges

    rows, edge_tabs, skipped = [], [], []
    for i, X in enumerate(targets, 1):
        if X == core: continue
        if X not in G:
            skipped.append((X, "not in graph")); continue
        row, edges = one(X)
        if row is None:
            skipped.append((X, "no downstream walk (no outgoing signed edge)")); continue
        rows.append(row)
        if edges is not None: edge_tabs.append(edges)
        log(f"[screen] {i}/{len(targets)} {X}: interference {row['interference_score']:.2f} ({row['n_interference']} edges), "
            f"TF enrichment {row['tf_enrichment']:.2f} (p={row['tf_p']:.2g}), dependency {row['dependency']:.2f}, direction {row['net_direction']:+.2f}")
    for X, why in skipped: log(f"[screen] skipped {X}: {why}")
    tab = pd.DataFrame(rows)
    if not len(tab):
        return dict(table=tab, edges=pd.DataFrame(), core_nodes=n_full, baseline=pd.DataFrame(), skipped=skipped)

    # automatic baseline: degree-matched random genes against the same core
    base_rows = []
    if n_random:
        done = set(targets) | {core}
        for X in tab.target:
            for R in degree_matched_random(G, X, max(1, n_random // max(1, len(tab)) + 1), rng, exclude=done):
                if R in done: continue
                done.add(R)
                r, _ = one(R)
                if r: base_rows.append(dict(random_gene=R, matched_to=X, interference_score=r["interference_score"],
                                            tf_enrichment=r["tf_enrichment"], dependency=r["dependency"]))
        if len(base_rows) >= 5:
            b = pd.DataFrame(base_rows); mu, sd = b.interference_score.mean(), b.interference_score.std(ddof=1) or 1.0
            tab["baseline_z"] = ((tab.interference_score - mu) / sd).round(2)
            tab["baseline_pct"] = tab.interference_score.map(lambda v: round(100 * (b.interference_score < v).mean(), 1))
            log(f"[screen] baseline: {len(b)} degree-matched random genes vs {core}: interference {mu:.2f} ± {sd:.2f}")
        else:
            tab["baseline_z"] = np.nan; tab["baseline_pct"] = np.nan
    else:
        tab["baseline_z"] = np.nan; tab["baseline_pct"] = np.nan

    # composite: mean percentile rank of the three magnitude axes; enrichment counts only when p <= 0.5
    enr = tab.tf_enrichment.where(tab.tf_p <= 0.5, 0.0)
    comp = (tab.interference_score.rank(pct=True) + enr.rank(pct=True) + tab.dependency.rank(pct=True)) / 3
    tab["composite"] = comp.round(3)
    tab = tab.sort_values(["composite", "baseline_z"], ascending=[False, False]).reset_index(drop=True)
    tab.insert(0, "rank", range(1, len(tab) + 1))
    edges = pd.concat(edge_tabs, ignore_index=True) if edge_tabs else pd.DataFrame()
    return dict(table=tab, edges=edges, core_nodes=n_full, baseline=pd.DataFrame(base_rows), skipped=skipped)


def parse_targets(targets):
    """list | plain text ('AKT1 MTOR, GSK3B' or one per line) | path to a text/CSV file (first column)."""
    import os
    if isinstance(targets, (list, tuple, set)): return [str(t).upper().strip() for t in targets if str(t).strip()]
    s = str(targets)
    if os.path.exists(os.path.expanduser(s)):
        p = os.path.expanduser(s)
        if p.lower().endswith((".csv", ".tsv", ".txt")):
            try:
                df = pd.read_csv(p, sep=None, engine="python", header=None, comment="#")
                col = df.iloc[:, 0].astype(str)
                if col.iloc[0].lower() in ("gene", "symbol", "target", "node"): col = col.iloc[1:]
                return [g.upper().strip() for g in col if g.strip()]
            except Exception:
                s = open(p).read()
        else:
            s = open(p).read()
    import re
    return [g.upper() for g in re.split(r"[\s,;]+", s) if g and not g.startswith("#")]

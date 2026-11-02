"""
Pathway interference: how does inhibiting a TARGET gene change the signal flow of a CORE gene's network?

Step 1 (baseline)   : the CORE network = downstream walk from CORE (receptor -> ... -> TFs -> DoRothEA genes) on the full graph.
Step 2 (perturbation): a "TARGET inhibited" perturbation enters the CORE network through every edge u -> v where
                      u is TARGET itself or a node TARGET reaches (walk of TARGET with CORE removed), and v is in the CORE network.
                      Entry strength at v = - sign(TARGET -> u) x effect(u -> v) x weight(u -> v) x score_TARGET(u)
                      (minus because TARGET is INHIBITED: what TARGET normally pushes up goes down, and vice versa).
                      The perturbation then propagates INSIDE the CORE network along CORE's own signed edges
                      (value' = sum_in value(u) x effect(u -> v) x weight(u -> v) x hub penalty), summed over 1..max_depth steps.
Unit perturbation   : the entry strengths are scaled to total mass 1 before propagation (entry_mass = the raw total, kept as the
                      footprint). The SAME linear propagation started with a unit at CORE gives baseline_flow per node (how much of
                      CORE's own unit signal arrives there). delta_rel = delta / baseline_flow therefore compares like with like:
                      one unit of 'TARGET inhibited' against one unit of CORE signal at that node.
                      A node counts as changed only when |delta_rel| >= interference.min_rel_change (YAML, default 0.1 = 10 %).
                      Without this relative gate every node of a large live network ends up 'changed' (absolute deltas never reach 0
                      on a hub-rich network), and flow_change_index saturates at 1.0 for every target.
Result per CORE node : delta (signed, same units as baseline_flow), baseline_flow, delta_rel, direction (decreased / increased / unchanged),
                      entry (direct / propagated), via (TARGET-side nodes that act on it). Per DoRothEA gene: delta = delta_rel(TF) x mode of regulation.
Interpretation       : delta_rel < 0 = this node's activity (or gene's expression) is predicted to FALL when TARGET is inhibited;
                      > 0 = to RISE (TARGET was a brake on it). -0.30 = the CORE-driven flow at this node drops by about 30 %.
                      TF sums (tf_delta_down / up) and the index are computed on delta_rel, so they are comparable between runs.
"""
import math
from collections import defaultdict
import numpy as np
import pandas as pd
from . import walk as Wk


def _propagate(G, core_nodes, idx, pen, x0, core, max_depth, min_abs, pin_core):
    """Linear signed propagation of x0 along the CORE network's own edges for max_depth steps; returns the summed influence."""
    N = len(core_nodes); total = x0.copy(); x = x0.copy()
    edges = [(idx[a], idx[b], d["effect"] * d["weight"] * pen.get(b, 1.0)) for a, b, d in G.edges(data=True)
             if a in idx and b in idx and d["effect"] != 0]
    for _ in range(max_depth):
        nxt = np.zeros(N)
        for ia, ib, w in edges:
            if x[ia] != 0: nxt[ib] += x[ia] * w
        if pin_core: nxt[idx[core]] = 0.0                     # the seed itself is only moved by direct entry (no feedback onto the seed)
        x = nxt; total += x
        if np.abs(x).sum() < min_abs: break
    return total


def interference_flow(G, core_df, core, target, target_df, genes_df=None, max_depth=5, hub_strength=0.25, min_abs=1e-4, min_rel=None):
    if min_rel is None:
        from . import data as Dt
        min_rel = float(Dt.CFG.get("interference", {}).get("min_rel_change", 0.1))
    core_nodes = [core] + [n for n in core_df.node if n != core]
    cset = set(core_nodes)
    sX = dict(zip(target_df.node, target_df.sign)); sX[target] = 1
    scX = dict(zip(target_df.node, target_df.score)); scX[target] = 1.0
    hX = dict(zip(target_df.node, target_df.hops)); hX[target] = 0
    xside = set(sX)
    pen = Wk.hub_penalty(G, strength=hub_strength)
    # 1) entry points: edges from the TARGET side into the CORE network
    entry = defaultdict(float); via = defaultdict(set); entry_edges = []
    for u in xside:
        if u not in G: continue
        if u in cset:
            # shared node (or TARGET itself sitting on the CORE pathway): TARGET's own effect on it is the entry; inhibition flips it
            if u == core and u != target: continue
            val = -sX[u] * scX[u]
            if val == 0: continue                      # TARGET's sign on this shared node is ambiguous: no directional entry
            entry[u] += val; via[u].add(target)
            entry_edges.append(dict(from_node=target, to_node=u, hops_target_to_from=0, sign_target_to_from=1, edge_effect=int(sX[u]),
                                    weight=round(scX[u], 3), mechanism="TARGET pathway (shared node)" if u != target else "TARGET itself",
                                    residue="", db="", effect_of_target_inhibition_on_to=int(np.sign(val)), entry_strength=round(val, 4)))
            continue
        for v in G.successors(u):
            if v not in cset: continue
            d = G[u][v]; eff = d["effect"]
            if sX[u] == 0 or eff == 0: continue
            val = -sX[u] * eff * d["weight"] * scX[u]
            entry[v] += val; via[v].add(u)
            entry_edges.append(dict(from_node=u, to_node=v, hops_target_to_from=hX[u], sign_target_to_from=int(sX[u]), edge_effect=int(eff),
                                    weight=round(d["weight"], 3), mechanism=d.get("mechanism", ""), residue=d.get("residue", ""), db=d.get("source", ""),
                                    effect_of_target_inhibition_on_to=int(np.sign(val)), entry_strength=round(val, 4)))
    # 2) propagate inside the CORE network along CORE's own edges; the same operator with a unit at CORE = baseline flow
    idx = {n: i for i, n in enumerate(core_nodes)}; N = len(core_nodes)
    x0 = np.zeros(N)
    for v, val in entry.items(): x0[idx[v]] = val
    # unit perturbation: the raw entries are scaled to total mass 1 so that the perturbation and CORE's own unit signal (baseline below)
    # are in the same units; the raw total (entry_mass) is reported separately as TARGET's mechanistic footprint
    entry_mass = float(np.abs(x0).sum())
    if entry_mass > 0: x0 = x0 / entry_mass
    total = _propagate(G, core_nodes, idx, pen, x0, core, max_depth, min_abs, pin_core=core not in entry)
    b0 = np.zeros(N); b0[idx[core]] = 1.0
    base = np.abs(_propagate(G, core_nodes, idx, pen, b0, core, max_depth, min_abs, pin_core=True))
    base[idx[core]] = 1.0
    rows = []
    hC = dict(zip(core_df.node, core_df.hops)); hC[core] = 0
    end = dict(zip(core_df.node, core_df.is_endpoint)); end[core] = False
    for n in core_nodes:
        dlt = float(total[idx[n]]); bf = float(base[idx[n]])
        rel = dlt / bf if bf >= min_abs else 0.0
        direct = n in entry
        changed = abs(rel) >= min_rel
        direction = "unchanged" if not changed else "decreased" if rel < 0 else "increased"
        rows.append(dict(node=n, hops_from_core=hC.get(n, 0), is_tf=bool(end.get(n, False)), delta=round(dlt, 4), baseline_flow=round(bf, 4),
                         delta_rel=round(rel, 4), direction=direction,
                         entry="direct" if direct else ("propagated" if changed else ""), via=";".join(sorted(via.get(n, ()))),
                         entry_strength=round(entry.get(n, 0.0), 4)))
    nodes = pd.DataFrame(rows)
    nodes["abs_delta"] = nodes.delta_rel.abs()
    nodes = nodes.sort_values(["abs_delta", "hops_from_core"], ascending=[False, True]).drop(columns="abs_delta").reset_index(drop=True)
    genes = pd.DataFrame()
    if genes_df is not None and len(genes_df):
        dmap = dict(zip(nodes.node, nodes.delta_rel.where(nodes.direction != "unchanged", 0.0)))
        g = genes_df.copy()
        g["delta_tf"] = g.tf.map(dmap).fillna(0.0)
        g["delta"] = (g.delta_tf * g.mor).round(4)
        g["direction"] = g.delta.map(lambda v: "unchanged" if v == 0 else "down" if v < 0 else "up")
        genes = g[g.delta != 0].groupby("gene", as_index=False).agg(delta=("delta", "sum"), via_tfs=("tf", lambda s: ";".join(sorted(set(s)))))
        genes["direction"] = genes.delta.map(lambda v: "down" if v < 0 else "up")
        genes = genes.reindex(genes.delta.abs().sort_values(ascending=False).index).reset_index(drop=True)
    edges = pd.DataFrame(entry_edges)
    if len(edges): edges = edges.reindex(edges.entry_strength.abs().sort_values(ascending=False).index).reset_index(drop=True)
    n_tf = int(nodes.is_tf.sum()); ch = nodes[(nodes.is_tf) & (nodes.direction != "unchanged")]
    tfd = ch.delta_rel; d_down = abs(float(tfd[tfd < 0].sum())); d_up = float(tfd[tfd > 0].sum())
    summary = dict(core=core, target=target, core_nodes=N - 1, core_tfs=n_tf, entry_points=len(entry), entry_edges=len(edges), entry_mass=round(entry_mass, 3),
                   nodes_changed=int((nodes.direction != "unchanged").sum()), nodes_decreased=int((nodes.direction == "decreased").sum()),
                   nodes_increased=int((nodes.direction == "increased").sum()), tfs_changed=len(ch),
                   tfs_decreased=int((ch.direction == "decreased").sum()), tfs_increased=int((ch.direction == "increased").sum()),
                   genes_changed=len(genes), genes_down=int((genes.direction == "down").sum()) if len(genes) else 0,
                   genes_up=int((genes.direction == "up").sum()) if len(genes) else 0,
                   flow_change_index=round(len(ch) / n_tf, 3) if n_tf else 0.0,
                   tf_delta_down=round(d_down, 3), tf_delta_up=round(d_up, 3), net_tf_delta=round(d_up - d_down, 3),
                   net_direction=round((d_up - d_down) / (d_up + d_down), 2) if (d_up + d_down) else 0.0,
                   target_in_core_network=target in cset, min_rel_change=min_rel)
    return dict(nodes=nodes, genes=genes, edges=edges, summary=summary)


def route_dependency(G, core, core_df, target, tfs, max_depth=5, min_score=0.35, hub_strength=0.3, through_tf=True, score_mode="relative", weak_ratio=None):
    """
    Necessity test (Boolean / OR logic complement to the linear amplitude model): CORE is walked AGAIN with TARGET removed from the graph,
    with the same walk parameters, and every CORE TF is compared before / after:
      lost      TARGET's removal leaves no route from CORE to this TF  -> a TARGET inhibitor blocks this output
      weakened  the TF is still reached but its best-path score falls below weak_ratio x the original (only a poorer detour remains)
      kept      an equally good route exists without TARGET -> the TF stays active (redundant / parallel input; e.g. EGFR for the FGFR trunk)
    necessity_index = (lost + weakened) / CORE TFs;  route_score_loss = 1 - sum(score after) / sum(score before) over CORE TFs.
    A TARGET can have a high amplitude score (flow_change_index) and necessity 0: it feeds the same trunk in parallel, so removing it
    lowers the total input but CORE's own route is intact. Interpretation: parallel input, not a block.
    """
    if weak_ratio is None:
        from . import data as Dt
        weak_ratio = float(Dt.CFG.get("interference", {}).get("weak_ratio", 0.5))
    before = core_df[core_df.is_endpoint.astype(bool)][["node", "score", "hops"]]
    n_tf = len(before)
    H = G.copy()
    if target in H and target != core: H.remove_node(target)
    after_df = Wk.signed_walk(H, core, tfs, direction="down", max_depth=max_depth, min_score=min_score, hub_strength=hub_strength,
                              stop_at_endpoints=not through_tf, score_mode=score_mode) if core in H else Wk.EMPTY_WALK.copy()
    after = dict(zip(after_df.node, after_df.score)) if len(after_df) else {}
    rows = []
    for r in before.itertuples():
        a = float(after.get(r.node, 0.0)); b = float(r.score) or 1e-9; ratio = a / b
        status = "lost" if r.node not in after else ("weakened" if ratio < weak_ratio else "kept")
        rows.append(dict(node=r.node, hops_from_core=int(r.hops), score_with_target=round(b, 3), score_without_target=round(a, 3), ratio=round(ratio, 3), status=status))
    t = pd.DataFrame(rows, columns=["node", "hops_from_core", "score_with_target", "score_without_target", "ratio", "status"])
    lost = int((t.status == "lost").sum()) if len(t) else 0; weak = int((t.status == "weakened").sum()) if len(t) else 0
    sb = float(before.score.sum()) if n_tf else 0.0; sa = float(sum(after.get(n, 0.0) for n in before.node)) if n_tf else 0.0
    nodes_before = set(core_df.node); nodes_after = set(after_df.node) if len(after_df) else set()
    summ = dict(core_tfs=n_tf, tfs_lost=lost, tfs_weakened=weak, tfs_kept=n_tf - lost - weak,
                necessity_index=round((lost + weak) / n_tf, 3) if n_tf else 0.0,
                route_score_loss=round(max(0.0, 1 - sa / sb), 3) if sb else 0.0,   # clipped: a hub-penalty artefact can raise scores slightly
                nodes_lost=len(nodes_before - nodes_after - {target}), weak_ratio=weak_ratio)
    if len(t): t = t.sort_values(["status", "ratio"], key=lambda c: c.map({"lost": 0, "weakened": 1, "kept": 2}) if c.name == "status" else c).reset_index(drop=True)
    return dict(table=t, summary=summ)


def classify_target(necessity_index, flow_change_index, net_direction, strong=0.2, weak=0.05):
    """Plain-language class from the two axes: necessity (route loss) x amplitude (linear flow change) and the direction."""
    if necessity_index >= strong: return "bottleneck: inhibition blocks CORE output"
    if pd.notna(net_direction) and net_direction >= 0.3 and flow_change_index >= weak: return "brake: inhibition raises CORE output"
    if flow_change_index >= strong: return "parallel input: lowers amplitude, CORE route intact (redundant)"
    if flow_change_index >= weak or necessity_index > 0: return "minor modulator"
    return "no effect on CORE output"


def screen_targets(G, core, core_df, core_genes, targets, tfs, max_depth=5, min_score=0.35, hub_strength=0.3, through_tf=True, score_mode="relative",
                   n_random=20, seed=0, log=print):
    """
    Rank candidate TARGETs by how much their inhibition changes the CORE network (the CORE network is computed once).
    Per target: the TARGET walk with CORE removed (same thresholds), then interference_flow. Score axes:
      flow_change_index  fraction of CORE's TFs whose activity changes
      tf_delta_total     sum of |delta_rel| over CORE's changed TFs (relative change of CORE-driven flow; 1.0 = one TF fully silenced)
      direct_entries     CORE nodes receiving a direct edge from TARGET's side (mechanistic footprint); entry_mass = raw summed entry strength
      tf_delta_down/up   signed parts of that total (sum of decreases / sum of increases at the TFs); net_tf_delta = up - down;
                         net_direction = (up - down)/(up + down): -1 = inhibition silences CORE output, +1 = it raises it (brake release)
      down_share         share of changed TFs that go DOWN
      necessity_index    route_dependency: fraction of CORE TFs that lose their route (or keep only a much poorer one) when TARGET is
                         removed from the graph -> can a TARGET inhibitor BLOCK CORE's output? (EGFR vs FGFR3: amplitude high, necessity 0)
      route_score_loss   1 - Σ best-path score of CORE TFs without TARGET / with TARGET
      class              plain-language reading of necessity x amplitude x direction (bottleneck / parallel input / brake / minor / none)
      baseline_z         z of tf_delta_total vs degree-matched random genes run against the same CORE
      composite          mean percentile rank of necessity_index, flow_change_index, tf_delta_total and direct_entries (baseline_z breaks ties)
    """
    import random
    from . import screen as Sc
    rng = random.Random(seed)
    rows, skipped, details = [], [], {}
    def one(X):
        H = G.copy()
        if core in H: H.remove_node(core)
        if X not in H: return None
        from . import data as Dt
        Wk.prepare_seed(H, X, "down", Dt.CFG, Dt.REACTOME_DIRECTED)      # same seed handling as pathway-interference (seed rescue)
        dX = Wk.signed_walk(H, X, tfs, direction="down", max_depth=max_depth, min_score=min_score, hub_strength=hub_strength,
                            stop_at_endpoints=not through_tf, score_mode=score_mode)
        res = interference_flow(G, core_df, core, X, dX, genes_df=core_genes, max_depth=max_depth, hub_strength=hub_strength)
        dep = route_dependency(G, core, core_df, X, tfs, max_depth, min_score, hub_strength, through_tf, score_mode)
        res["route"] = dep["table"]; res["summary"].update({k: v for k, v in dep["summary"].items() if k != "core_tfs"})
        s_ = res["summary"]; nd = res["nodes"]; tf = nd[nd.is_tf]
        ch = tf[tf.direction != "unchanged"]
        d_down = abs(float(ch.delta_rel[ch.delta_rel < 0].sum())); d_up = float(ch.delta_rel[ch.delta_rel > 0].sum())
        row = dict(target=X, target_nodes=len(dX), entry_points=s_["entry_points"], entry_edges=s_["entry_edges"],
                   direct_entries=int((nd.entry == "direct").sum()), entry_mass=s_["entry_mass"], nodes_changed=s_["nodes_changed"],
                   tfs_changed=s_["tfs_changed"], tfs_down=s_["tfs_decreased"], tfs_up=s_["tfs_increased"], flow_change_index=s_["flow_change_index"],
                   tf_delta_total=round(d_down + d_up, 3), tf_delta_down=round(d_down, 3), tf_delta_up=round(d_up, 3), net_tf_delta=round(d_up - d_down, 3),
                   net_direction=round((d_up - d_down) / (d_up + d_down), 2) if (d_up + d_down) else float("nan"),
                   down_share=round(s_["tfs_decreased"] / len(ch), 2) if len(ch) else float("nan"),
                   necessity_index=s_["necessity_index"], tfs_lost=s_["tfs_lost"], tfs_weakened=s_["tfs_weakened"], route_score_loss=s_["route_score_loss"],
                   genes_down=s_["genes_down"], genes_up=s_["genes_up"], target_on_core_network=s_["target_in_core_network"],
                   top_tfs=", ".join(f"{r.node}{'↓' if r.delta < 0 else '↑'}" for r in ch.head(6).itertuples()))
        row["class"] = classify_target(row["necessity_index"], row["flow_change_index"], row["net_direction"])
        return row, res
    for i, X in enumerate(targets, 1):
        if X == core: continue
        if X not in G: skipped.append((X, "not in graph")); continue
        r = one(X)
        if r is None or r[0]["target_nodes"] == 0 and r[0]["entry_points"] == 0:
            skipped.append((X, "no downstream walk / no entry into the CORE network")); continue
        row, res = r; rows.append(row); details[X] = res
        log(f"[iscreen] {i}/{len(targets)} {X}: TFs changed {row['tfs_changed']} ({row['tfs_down']}↓ {row['tfs_up']}↑), |delta| {row['tf_delta_total']:.2f}, "
            f"direct entries {row['direct_entries']}; route loss: {row['tfs_lost']} TFs lost / {row['tfs_weakened']} weakened (necessity {row['necessity_index']:.2f}) -> {row['class']}")
    for X, why in skipped: log(f"[iscreen] skipped {X}: {why}")
    tab = pd.DataFrame(rows)
    if not len(tab): return dict(table=tab, baseline=pd.DataFrame(), skipped=skipped, details=details)
    base = []
    if n_random:
        done = set(targets) | {core}
        per = max(1, n_random // max(1, len(tab)) + 1)
        for X in tab.target:
            for R in Sc.degree_matched_random(G, X, per, rng, exclude=done):
                if R in done: continue
                done.add(R); r = one(R)
                if r: base.append(dict(random_gene=R, matched_to=X, tf_delta_total=r[0]["tf_delta_total"], flow_change_index=r[0]["flow_change_index"]))
        if len(base) >= 5:
            b = pd.DataFrame(base); mu, sd = b.tf_delta_total.mean(), b.tf_delta_total.std(ddof=1) or 1.0
            tab["baseline_z"] = ((tab.tf_delta_total - mu) / sd).round(2)
            log(f"[iscreen] baseline: {len(b)} degree-matched random genes vs {core}: |delta| at TFs {mu:.2f} ± {sd:.2f}")
        else:
            tab["baseline_z"] = float("nan")
    else:
        tab["baseline_z"] = float("nan")
    comp = (tab.necessity_index.rank(pct=True) + tab.flow_change_index.rank(pct=True) + tab.tf_delta_total.rank(pct=True) + tab.direct_entries.rank(pct=True)) / 4
    tab["composite"] = comp.round(3)
    tab = tab.sort_values(["composite", "baseline_z"], ascending=[False, False]).reset_index(drop=True)
    tab.insert(0, "rank", range(1, len(tab) + 1))
    return dict(table=tab, baseline=pd.DataFrame(base), skipped=skipped, details=details)

"""
Pathway interference: how does inhibiting a TARGET gene change the signal flow of a CORE gene's network?

Step 1 (baseline)   : the CORE network = downstream walk from CORE (receptor -> ... -> TFs -> DoRothEA genes) on the full graph.
Step 2 (perturbation): a "TARGET inhibited" perturbation enters the CORE network through every edge u -> v where
                      u is TARGET itself or a node TARGET reaches (walk of TARGET with CORE removed), and v is in the CORE network.
                      Entry strength at v = - sign(TARGET -> u) x effect(u -> v) x weight(u -> v) x score_TARGET(u)
                      (minus because TARGET is INHIBITED: what TARGET normally pushes up goes down, and vice versa).
                      The perturbation then propagates INSIDE the CORE network along CORE's own signed edges
                      (value' = sum_in value(u) x effect(u -> v) x weight(u -> v) x hub penalty), summed over 1..max_depth steps.
Result per CORE node : delta (signed magnitude), direction (decreased / increased / mixed / unchanged), entry (direct / propagated),
                      via (TARGET-side nodes that act on it). Per DoRothEA gene: delta = delta(TF) x mode of regulation.
Interpretation       : delta < 0  = this node's activity (or gene's expression) is predicted to FALL when TARGET is inhibited;
                      delta > 0 = to RISE (TARGET was a brake on it). Magnitudes are relative within one run (no units).
"""
import math
from collections import defaultdict
import numpy as np
import pandas as pd
from . import walk as Wk


def interference_flow(G, core_df, core, target, target_df, genes_df=None, max_depth=5, hub_strength=0.25, min_abs=1e-4):
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
    # 2) propagate inside the CORE network along CORE's own edges
    idx = {n: i for i, n in enumerate(core_nodes)}; N = len(core_nodes)
    x0 = np.zeros(N)
    for v, val in entry.items(): x0[idx[v]] = val
    total = x0.copy(); x = x0.copy()
    for _ in range(max_depth):
        nxt = np.zeros(N)
        for a, b, d in G.edges(data=True):
            if a in idx and b in idx and d["effect"] != 0 and x[idx[a]] != 0:
                nxt[idx[b]] += x[idx[a]] * d["effect"] * d["weight"] * pen.get(b, 1.0)
        nxt[idx[core]] = 0.0 if core not in entry else nxt[idx[core]]       # the seed itself is only moved by direct entry
        x = nxt; total += x
        if np.abs(x).sum() < min_abs: break
    rows = []
    hC = dict(zip(core_df.node, core_df.hops)); hC[core] = 0
    end = dict(zip(core_df.node, core_df.is_endpoint)); end[core] = False
    for n in core_nodes:
        dlt = float(total[idx[n]])
        direct = n in entry
        direction = "unchanged" if abs(dlt) < min_abs else "decreased" if dlt < 0 else "increased"
        rows.append(dict(node=n, hops_from_core=hC.get(n, 0), is_tf=bool(end.get(n, False)), delta=round(dlt, 4), direction=direction,
                         entry="direct" if direct else ("propagated" if abs(dlt) >= min_abs else ""), via=";".join(sorted(via.get(n, ()))),
                         entry_strength=round(entry.get(n, 0.0), 4)))
    nodes = pd.DataFrame(rows)
    nodes["abs_delta"] = nodes.delta.abs()
    nodes = nodes.sort_values(["abs_delta", "hops_from_core"], ascending=[False, True]).drop(columns="abs_delta").reset_index(drop=True)
    genes = pd.DataFrame()
    if genes_df is not None and len(genes_df):
        dmap = dict(zip(nodes.node, nodes.delta))
        g = genes_df.copy()
        g["delta_tf"] = g.tf.map(dmap).fillna(0.0)
        g["delta"] = (g.delta_tf * g.mor).round(4)
        g["direction"] = g.delta.map(lambda v: "unchanged" if abs(v) < min_abs else "down" if v < 0 else "up")
        genes = g[g.delta != 0].groupby("gene", as_index=False).agg(delta=("delta", "sum"), via_tfs=("tf", lambda s: ";".join(sorted(set(s)))))
        genes["direction"] = genes.delta.map(lambda v: "down" if v < 0 else "up")
        genes = genes.reindex(genes.delta.abs().sort_values(ascending=False).index).reset_index(drop=True)
    edges = pd.DataFrame(entry_edges)
    if len(edges): edges = edges.reindex(edges.entry_strength.abs().sort_values(ascending=False).index).reset_index(drop=True)
    n_tf = int(nodes.is_tf.sum()); ch = nodes[(nodes.is_tf) & (nodes.direction != "unchanged")]
    summary = dict(core=core, target=target, core_nodes=N - 1, core_tfs=n_tf, entry_points=len(entry), entry_edges=len(edges),
                   nodes_changed=int((nodes.direction != "unchanged").sum()), nodes_decreased=int((nodes.direction == "decreased").sum()),
                   nodes_increased=int((nodes.direction == "increased").sum()), tfs_changed=len(ch),
                   tfs_decreased=int((ch.direction == "decreased").sum()), tfs_increased=int((ch.direction == "increased").sum()),
                   genes_changed=len(genes), genes_down=int((genes.direction == "down").sum()) if len(genes) else 0,
                   genes_up=int((genes.direction == "up").sum()) if len(genes) else 0,
                   flow_change_index=round(len(ch) / n_tf, 3) if n_tf else 0.0,
                   target_in_core_network=target in cset)
    return dict(nodes=nodes, genes=genes, edges=edges, summary=summary)

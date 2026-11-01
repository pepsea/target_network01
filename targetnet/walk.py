"""
Signed pathway walk.

A node is reached with a *path score* (product of edge weights x hub penalties) and a *sign*
(product of edge effects: +1 activation, -1 inhibition). We run a best-path search over the state
space (node, sign) so that we know, for every node, the best positive path AND the best negative
path from the seed. If both are comparably good the node is 'mixed'.

Stop rules: the walk does not expand past an *endpoint* (receptor for upstream, TF for downstream),
past max_depth, or below min_score.
"""
import heapq, math
import networkx as nx
import pandas as pd


def build_graph(directed_df, expressed=None):
    G = nx.DiGraph()
    for r in directed_df.itertuples():
        if expressed is not None and (r.src not in expressed or r.dst not in expressed):
            continue
        G.add_edge(r.src, r.dst, effect=int(r.effect), weight=float(r.weight), source=r.source,
                   mechanism=r.mechanism, residue=r.residue)
    return G


def hub_penalty(G, max_deg=None, strength=0.5):
    deg = dict(G.degree()); m = max(deg.values()) if deg else 1
    return {n: 1.0 - strength * math.log1p(d) / math.log1p(m) for n, d in deg.items()}


def signed_walk(G, seed, endpoints, direction="down", max_depth=6, min_score=0.15, hub_strength=0.25):
    """
    Returns DataFrame: node, hops, score, sign(+1/-1/0 mixed), ambiguous(bool), is_endpoint, best_path (list)
    direction 'down' follows edges src->dst; 'up' follows them backwards (src is the upstream regulator).
    """
    if seed not in G:
        raise SystemExit(f"{seed} not in the signaling graph")
    pen = hub_penalty(G, strength=hub_strength)
    nbrs = (lambda n: G.successors(n)) if direction == "down" else (lambda n: G.predecessors(n))
    edge = (lambda a, b: G[a][b]) if direction == "down" else (lambda a, b: G[b][a])

    # best[(node, sign)] = (score, hops, path, ambiguous)
    best = {(seed, 1): (1.0, 0, [seed], False)}
    heap = [(-1.0, 0, seed, 1)]
    while heap:
        negs, hops, n, sgn = heapq.heappop(heap)
        score = -negs
        if best.get((n, sgn), (0,))[0] > score: continue
        if n != seed and (n in endpoints or hops >= max_depth): continue   # do not expand past endpoints
        for m in nbrs(n):
            if m == seed: continue
            e = edge(n, m)
            ns = score * e["weight"] * pen.get(m, 1.0)
            if ns < min_score: continue
            eff = e["effect"]; amb = best[(n, sgn)][3] or eff == 0
            nsgn = sgn * (eff if eff != 0 else 1)
            key = (m, nsgn)
            if ns > best.get(key, (0,))[0]:
                best[key] = (ns, hops + 1, best[(n, sgn)][2] + [m], amb)
                heapq.heappush(heap, (-ns, hops + 1, m, nsgn))

    rows = []
    nodes = {n for n, _ in best if n != seed}
    for n in nodes:
        pos, neg = best.get((n, 1)), best.get((n, -1))
        cands = [(s, 1) for s in [pos] if s] + [(s, -1) for s in [neg] if s]
        (sc, hp, path, amb), sg = max(cands, key=lambda t: t[0][0])
        other = neg if sg == 1 else pos
        mixed = other is not None and other[0] >= 0.8 * sc
        rows.append(dict(node=n, hops=hp, score=round(sc, 4), sign=0 if mixed else sg, ambiguous=amb or mixed,
                         is_endpoint=n in endpoints, best_path=" > ".join(path if direction == "down" else path[::-1])))
    df = pd.DataFrame(rows).sort_values(["hops", "score"], ascending=[True, False]).reset_index(drop=True)
    return df


def extend_regulon(down_df, regulon, tfs, max_targets_per_tf=None):
    """Add TF -> target-gene layer. Gene sign = TF sign x mode of regulation; gene score = TF score x regulon weight."""
    tf_rows = down_df[down_df.node.isin(tfs)]
    rows = []
    for r in tf_rows.itertuples():
        sub = regulon[regulon.tf == r.node].sort_values("weight", ascending=False)
        if max_targets_per_tf: sub = sub.head(max_targets_per_tf)
        for t in sub.itertuples():
            rows.append(dict(tf=r.node, gene=t.target, mor=int(t.mor), confidence=t.confidence,
                             sign=0 if r.sign == 0 else int(r.sign * t.mor), score=round(r.score * t.weight, 4),
                             hops=r.hops + 1, ambiguous=bool(r.ambiguous)))
    genes = pd.DataFrame(rows)
    if len(genes):
        # a gene regulated by several reached TFs: keep best-scoring TF, note the count
        genes["n_tfs"] = genes.groupby("gene").gene.transform("size")
        genes = genes.sort_values("score", ascending=False)
    return genes


def subgraph_edges(G, nodes, direction="down"):
    """Edges among a node set, as a DataFrame with sign/residue/source (for drawing)."""
    rows = []
    for a, b, d in G.edges(nodes, data=True):
        if a in nodes and b in nodes:
            rows.append(dict(source=a, target=b, effect=d["effect"], residue=d["residue"], db=d["source"], weight=d["weight"]))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------- crosstalk (mode 2)
def crosstalk(G, causal_df, goi_df, causal, goi, causal_genes=None, goi_genes=None):
    """
    causal_df / goi_df: outputs of signed_walk(direction='down') for the causal gene C and the gene of interest X.
    causal_genes / goi_genes: regulon-extended gene tables (optional).
    Returns dict with shared nodes, interference edges, endpoint agreement, and summary indices.
    """
    DC = {causal} | set(causal_df.node); DX = {goi} | set(goi_df.node)
    sC = dict(zip(causal_df.node, causal_df.sign)); sX = dict(zip(goi_df.node, goi_df.sign))
    hC = dict(zip(causal_df.node, causal_df.hops)); hX = dict(zip(goi_df.node, goi_df.hops))
    sC[causal] = 1; sX[goi] = 1; hC[causal] = 0; hX[goi] = 0

    # 1) shared protein nodes
    shared = sorted((DC & DX) - {causal, goi})
    shared_rows = [dict(node=n, hops_from_causal=hC[n], hops_from_goi=hX[n], sign_from_causal=sC[n], sign_from_goi=sX[n],
                        agreement=_agree(sC[n], sX[n])) for n in shared]

    # 2) interference edges: X-side node -> C-side node, where the X-side node is NOT reachable from C
    #    (from-nodes that are also ancestors of X sit on a feedback loop through X; flagged, not removed)
    anc_X = nx.ancestors(G, goi)
    inter = []
    for u in DX - DC:
        for v in G.successors(u):
            if v in DC:
                e = G[u][v]
                # net effect of X on v via this edge = sign(X->u) x effect(u->v)
                net = 0 if sX[u] == 0 or e["effect"] == 0 else sX[u] * e["effect"]
                inter.append(dict(from_node=u, to_node=v, hops_goi_to_from=hX[u], hops_causal_to_to=hC[v], edge_effect=e["effect"],
                                  residue=e["residue"], mechanism=e["mechanism"], db=e["source"], net_effect_of_goi_on_target=net,
                                  causal_effect_on_target=sC[v], agreement=_agree(sC[v], net), via_feedback_loop=u in anc_X))
    inter_df = pd.DataFrame(inter)

    # 3) endpoint agreement (TF layer and gene layer)
    tf_rows = [r for r in shared_rows if causal_df.set_index("node").is_endpoint.get(r["node"], False)]
    gene_rows = []
    if causal_genes is not None and goi_genes is not None and len(causal_genes) and len(goi_genes):
        gc = causal_genes.drop_duplicates("gene").set_index("gene"); gx = goi_genes.drop_duplicates("gene").set_index("gene")
        for g in sorted(set(gc.index) & set(gx.index)):
            gene_rows.append(dict(gene=g, tf_from_causal=gc.loc[g, "tf"], tf_from_goi=gx.loc[g, "tf"],
                                  sign_from_causal=int(gc.loc[g, "sign"]), sign_from_goi=int(gx.loc[g, "sign"]),
                                  agreement=_agree(int(gc.loc[g, "sign"]), int(gx.loc[g, "sign"]))))
    n_endpoints_C = int(causal_df.is_endpoint.sum()) + (len(causal_genes.drop_duplicates("gene")) if causal_genes is not None and len(causal_genes) else 0)
    n_shared_end = len(tf_rows) + len(gene_rows)
    agreements = [r["agreement"] for r in tf_rows + gene_rows]
    summary = dict(causal=causal, goi=goi, n_nodes_causal=len(DC) - 1, n_nodes_goi=len(DX) - 1, n_shared_nodes=len(shared),
                   n_interference_edges=len(inter_df), n_endpoints_causal=n_endpoints_C, n_shared_endpoints=n_shared_end,
                   crosstalk_index=round(n_shared_end / n_endpoints_C, 3) if n_endpoints_C else 0.0,
                   synergistic=agreements.count("synergistic"), antagonistic=agreements.count("antagonistic"),
                   mixed=agreements.count("mixed"))
    return dict(shared=pd.DataFrame(shared_rows), interference=inter_df, tf_agreement=pd.DataFrame(tf_rows),
                gene_agreement=pd.DataFrame(gene_rows), summary=summary)


def _agree(a, b):
    if a == 0 or b == 0: return "mixed"
    return "synergistic" if a == b else "antagonistic"

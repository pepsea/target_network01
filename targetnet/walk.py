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
import numpy as np
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
    df = pd.DataFrame(rows)
    if not len(df):
        # nothing reached: report why instead of failing downstream (empty frame keeps the full column set)
        first = [(m, round(edge(seed, m)["weight"] * pen.get(m, 1.0), 3)) for m in nbrs(seed)]
        first.sort(key=lambda t: -t[1])
        print(f"[walk] {seed} ({direction}): 0 nodes reached. {len(first)} direct neighbours in this direction; "
              + (f"best first-hop score {first[0][1]} ({first[0][0]}) < min_score {min_score} -> lower min_score or hub_penalty"
                 if first else "the seed has no edges in this direction in the loaded graph (or they were removed)"))
        return pd.DataFrame(columns=["node", "hops", "score", "influence", "sign", "sign_confidence", "sign_bestpath", "ambiguous", "is_endpoint", "best_path"])
    if len(df):
        inf = signed_influence(G, seed, [seed] + list(df.node), direction, max_depth, hub_strength)
        df = df.rename(columns={"sign": "sign_bestpath"}).merge(inf, on="node", how="left")
        # consensus sign over all paths is the reported sign; best-path sign kept for reference
        df["sign"] = df.sign_consensus.fillna(df.sign_bestpath).astype(int)
        df["ambiguous"] = df.ambiguous | (df.sign == 0)
        df = df[["node", "hops", "score", "influence", "sign", "sign_confidence", "sign_bestpath", "ambiguous", "is_endpoint", "best_path"]]
    df = df.sort_values(["hops", "score"], ascending=[True, False]).reset_index(drop=True)
    return df


def signed_influence(G, seed, nodes, direction="down", max_depth=5, hub_strength=0.25):
    """
    Consensus sign over ALL walks (not just the best path): two-state influence propagation.
      p = activating influence, n = inhibiting influence reaching each node from the seed.
      step: p' = Wpos^T p + Wneg^T n ; n' = Wpos^T n + Wneg^T p   (ambiguous edges feed both states at half weight)
      summed over 1..max_depth steps; the seed is not re-entered.
    Returns DataFrame node, influence (p+n), sign_consensus (+1/-1/0), sign_confidence (|p-n|/(p+n)).
    """
    from scipy import sparse
    idx = {n: i for i, n in enumerate(nodes)}; N = len(nodes)
    if seed not in idx or N == 0: return pd.DataFrame(columns=["node", "influence", "sign_consensus", "sign_confidence"])
    pen = hub_penalty(G, strength=hub_strength)
    rows_p, cols_p, w_p, rows_n, cols_n, w_n, rows_a, cols_a, w_a = [], [], [], [], [], [], [], [], []
    for a, b, d in G.edges(data=True):
        if a not in idx or b not in idx: continue
        src, dst = (a, b) if direction == "down" else (b, a)          # walk direction
        w = d["weight"] * pen.get(dst, 1.0)
        tgt = (rows_p, cols_p, w_p) if d["effect"] == 1 else (rows_n, cols_n, w_n) if d["effect"] == -1 else (rows_a, cols_a, w_a)
        tgt[0].append(idx[src]); tgt[1].append(idx[dst]); tgt[2].append(w)
    M = lambda r, c, w: sparse.csr_matrix((w, (r, c)), shape=(N, N))
    Wp, Wn, Wa = M(rows_p, cols_p, w_p), M(rows_n, cols_n, w_n), M(rows_a, cols_a, w_a)
    p = np.zeros(N); n = np.zeros(N); p[idx[seed]] = 1.0
    P = np.zeros(N); Nn = np.zeros(N)
    for _ in range(max_depth):
        amb = 0.5 * Wa.T @ (p + n)
        p, n = Wp.T @ p + Wn.T @ n + amb, Wp.T @ n + Wn.T @ p + amb
        p[idx[seed]] = 0; n[idx[seed]] = 0
        P += p; Nn += n
        if P.sum() + Nn.sum() == 0: break
    tot = P + Nn
    with np.errstate(divide="ignore", invalid="ignore"):
        conf = np.where(tot > 0, np.abs(P - Nn) / tot, 0.0)
    sign = np.where(P > 1.25 * Nn, 1, np.where(Nn > 1.25 * P, -1, 0))
    out = pd.DataFrame({"node": nodes, "influence": np.round(tot, 4), "sign_consensus": sign, "sign_confidence": np.round(conf, 3)})
    return out[out.node != seed]


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
def crosstalk(G, causal_df, goi_df, causal, goi, causal_genes=None, goi_genes=None, tf_universe=None, min_conf=0.5):
    """
    causal_df / goi_df: outputs of signed_walk(direction='down') for the causal gene C and the gene of interest X.
    causal_genes / goi_genes: regulon-extended gene tables (optional).
    Returns dict with shared nodes, interference edges, endpoint agreement, and summary indices.
    """
    DC = {causal} | set(causal_df.node); DX = {goi} | set(goi_df.node)
    confC = dict(zip(causal_df.node, causal_df.get("sign_confidence", pd.Series(1.0, index=causal_df.index))))
    confX = dict(zip(goi_df.node, goi_df.get("sign_confidence", pd.Series(1.0, index=goi_df.index))))
    sC = dict(zip(causal_df.node, causal_df.sign)); sX = dict(zip(goi_df.node, goi_df.sign))
    hC = dict(zip(causal_df.node, causal_df.hops)); hX = dict(zip(goi_df.node, goi_df.hops))
    sC[causal] = 1; sX[goi] = 1; hC[causal] = 0; hX[goi] = 0

    # 1) shared protein nodes
    shared = sorted((DC & DX) - {causal, goi})
    def agree_conf(n):
        return _agree(sC[n], sX[n]) if (confC.get(n, 1) >= min_conf and confX.get(n, 1) >= min_conf) else "low-confidence"
    shared_rows = [dict(node=n, hops_from_causal=hC[n], hops_from_goi=hX[n], sign_from_causal=sC[n], sign_from_goi=sX[n],
                        conf_causal=round(confC.get(n, 1), 2), conf_goi=round(confX.get(n, 1), 2), agreement=agree_conf(n)) for n in shared]

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
    if len(inter_df):   # rank: proximal interference first (few hops on both sides), then edge weight
        inter_df["weight"] = [G[a][b]["weight"] for a, b in zip(inter_df.from_node, inter_df.to_node)]
        inter_df["proximity"] = inter_df.hops_goi_to_from + inter_df.hops_causal_to_to
        inter_df = inter_df.sort_values(["proximity", "weight"], ascending=[True, False]).reset_index(drop=True)
        inter_df.insert(0, "rank", range(1, len(inter_df) + 1))

    # 3) endpoint agreement (TF layer and gene layer)
    is_end = causal_df.set_index("node").is_endpoint
    tf_rows = [r for r in shared_rows if is_end.get(r["node"], False)]
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
    # TF-level enrichment: is the overlap of C's TFs and X's TFs larger than expected given how many TFs X reaches at all?
    tfC = {n for n in causal_df[causal_df.is_endpoint].node}; tfX = {n for n in goi_df[goi_df.is_endpoint].node}
    universe = set(tf_universe) & set(G.nodes) if tf_universe else (tfC | tfX)
    M, nX, nC, k = len(universe), len(tfX & universe), len(tfC & universe), len(tfC & tfX & universe)
    expected = nC * nX / M if M else 0.0
    try:
        from scipy.stats import hypergeom
        p_hyper = float(hypergeom.sf(k - 1, M, nX, nC)) if M and nC and nX else 1.0
    except Exception:
        p_hyper = float("nan")
    summary = dict(causal=causal, goi=goi, n_nodes_causal=len(DC) - 1, n_nodes_goi=len(DX) - 1, n_shared_nodes=len(shared),
                   n_interference_edges=len(inter_df), n_endpoints_causal=n_endpoints_C, n_shared_endpoints=n_shared_end,
                   crosstalk_index=round(n_shared_end / n_endpoints_C, 3) if n_endpoints_C else 0.0,
                   tf_overlap=k, tf_causal=nC, tf_goi=nX, tf_universe=M, tf_overlap_expected=round(expected, 2),
                   tf_enrichment=round(k / expected, 2) if expected else None, tf_overlap_p=p_hyper,
                   goi_coverage_of_all_tfs=round(nX / M, 3) if M else None,
                   synergistic=agreements.count("synergistic"), antagonistic=agreements.count("antagonistic"),
                   mixed=agreements.count("mixed"), low_confidence=agreements.count("low-confidence"))
    return dict(shared=pd.DataFrame(shared_rows), interference=inter_df, tf_agreement=pd.DataFrame(tf_rows),
                gene_agreement=pd.DataFrame(gene_rows), summary=summary)


def _agree(a, b):
    if a == 0 or b == 0: return "mixed"
    return "synergistic" if a == b else "antagonistic"

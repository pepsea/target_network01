"""
Plotly versions of the four views, for inline use in Jupyter.

Each function reads the JSON payload written next to the corresponding D3 HTML view
(flow.json, cell_map.json, crosstalk.json, synthetic_lethality.json) and returns a plotly Figure.

    from targetnet.plots import flow_figure, cell_map_figure, crosstalk_figure, volcano_figure
    flow_figure("runs/out_AKT1/downstream/flow.json").show()

Note on size: plotly.js itself is ~3.5 MB. Use the CDN renderer so it is NOT embedded in every notebook:
    import plotly.io as pio; pio.renderers.default = "notebook_connected"
(`use_cdn_renderer()` below does this). The figures themselves are small (a few hundred KB at most).
"""
import json, math
from collections import defaultdict
import plotly.graph_objects as go

# flat, muted palette on white (owner's convention)
COL = dict(seed="#1f3a5f", endpoint="#2d5f8f", node="#8a8a8a", shared="#6b5b95", causal="#4a7fb5", goi="#b36b3c",
           syn="#2e7d32", ant="#b23b3b", mix="#888888", edge="#777777", amb="#bbbbbb", phys="#c8c8c8")
FONT = dict(family="Helvetica, Arial, sans-serif", size=11, color="#222")


def use_cdn_renderer():
    """Keep notebooks small: load plotly.js from the CDN instead of embedding it in every output cell."""
    import plotly.io as pio
    pio.renderers.default = "notebook_connected"


def _load(p):
    return json.load(open(p))


def _base(title, w=1000, h=650):
    fig = go.Figure()
    fig.update_layout(title=dict(text=title, x=0.01, font=dict(size=14)), width=w, height=h, plot_bgcolor="white", paper_bgcolor="white",
                      font=FONT, showlegend=True, legend=dict(orientation="h", y=-0.02, x=0, font=dict(size=10)), margin=dict(l=20, r=20, t=50, b=40),
                      xaxis=dict(visible=False), yaxis=dict(visible=False), hovermode="closest")
    return fig


def _edge_traces(fig, edges, pos, rad, color_fn=None, dash_fn=None, name="edge", width=1.2, head_px=7):
    """Edges as cubic S-curves, one line trace per (colour, dash) style; arrow/bar heads as ONE rotated-marker trace.
    (No per-edge annotations: those make Plotly quadratic and take minutes for a few thousand edges.)"""
    groups = defaultdict(list)
    for e in edges:
        s, t = e["source"], e["target"]
        if s not in pos or t not in pos: continue
        col = color_fn(e) if color_fn else (COL["amb"] if e.get("effect", 1) == 0 else COL["edge"])
        dash = dash_fn(e) if dash_fn else ("dot" if e.get("effect", 1) == 0 else "solid")
        groups[(col, dash)].append(e)
    hx, hy, hang, hsym, hcol = [], [], [], [], []
    T = [i / 12 for i in range(13)]
    for (col, dash), es in groups.items():
        xs, ys = [], []
        for e in es:
            (x0, y0), (x1, y1) = pos[e["source"]], pos[e["target"]]
            dx, dy = x1 - x0, y1 - y0; dd = math.hypot(dx, dy) or 1; r = rad(e["target"]) * 0.75
            x1e, y1e = x1 - dx / dd * r, y1 - dy / dd * r
            mx = (x0 + x1e) / 2
            px = py = None
            for tt in T:
                u = 1 - tt
                bx = u**3 * x0 + 3 * u**2 * tt * mx + 3 * u * tt**2 * mx + tt**3 * x1e
                by = u**3 * y0 + 3 * u**2 * tt * y0 + 3 * u * tt**2 * y1e + tt**3 * y1e
                if tt == T[-2]: px, py = bx, by
                xs.append(bx); ys.append(by)
            xs.append(None); ys.append(None)
            # head: angle of the last curve segment (plotly marker.angle is clockwise from "up")
            ang = math.degrees(math.atan2(x1e - px, y1e - py))
            hx.append(x1e); hy.append(y1e); hang.append(ang)
            hsym.append("triangle-up" if e.get("effect", 1) >= 0 else "line-ew"); hcol.append(col)
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=col, width=width, dash=dash), hoverinfo="skip", name=name, showlegend=False))
    if hx:
        fig.add_trace(go.Scatter(x=hx, y=hy, mode="markers", hoverinfo="skip", showlegend=False,
                                 marker=dict(symbol=hsym, angle=hang, size=head_px + (2 if width > 2 else 0), color=hcol,
                                             line=dict(width=max(1.5, width), color=hcol))))


def _cap_nodes(nodes, edges, max_nodes, keep_ids=()):
    """Keep the seed, the nodes on the best paths of the top-scoring endpoints (so hop structure survives), then the
    best remaining nodes up to max_nodes; drop edges to removed nodes."""
    if len(nodes) <= max_nodes: return nodes, edges, False
    byid = {n["id"]: n for n in nodes}
    keep = {n["id"] for n in nodes if n.get("cls") == "seed" or n["id"] in keep_ids}
    ends = sorted([n for n in nodes if n.get("is_endpoint")], key=lambda n: -n.get("score", 0))
    for n in ends:                                   # whole best path of each endpoint, best endpoints first
        path = [p.strip() for p in str(n.get("best_path", "")).split(">")]
        if len(keep | set(path)) > max_nodes: break
        keep |= {p for p in path if p in byid}
    for n in sorted(nodes, key=lambda n: -n.get("score", 0)):
        if len(keep) >= max_nodes: break
        keep.add(n["id"])
    nodes = [n for n in nodes if n["id"] in keep]
    edges = [e for e in edges if e.get("source", e.get("from_node")) in keep and e.get("target", e.get("to_node")) in keep]
    return nodes, edges, True


def parse_focus(focus):
    """focus: None | list/set of gene symbols | path to a text file (one symbol per line, optional 2nd column = value)
       | dict {gene: value} | pandas Series (index = gene). Returns (set_of_genes, {gene: value} or {})."""
    if focus is None: return set(), {}
    try:
        import pandas as pd
        if isinstance(focus, pd.Series): focus = focus.to_dict()
        elif isinstance(focus, pd.DataFrame): focus = dict(zip(focus.iloc[:, 0], focus.iloc[:, 1])) if focus.shape[1] > 1 else list(focus.iloc[:, 0])
    except ImportError:
        pass
    if isinstance(focus, str):
        vals = {}
        for line in open(focus):
            parts = line.replace(",", "\t").split("\t")
            g = parts[0].strip().upper()
            if not g or g.startswith("#"): continue
            try: vals[g] = float(parts[1]) if len(parts) > 1 and parts[1].strip() else None
            except ValueError: vals[g] = None
        return set(vals), {g: v for g, v in vals.items() if v is not None}
    if isinstance(focus, dict):
        d = {str(k).upper(): v for k, v in focus.items()}
        return set(d), {g: float(v) for g, v in d.items() if v is not None}
    return {str(g).upper() for g in focus}, {}


def _value_colour(v, vmax):
    """diverging blue (negative) - light grey (0) - red (positive), symmetric around 0."""
    if vmax == 0: return COL["goi"]
    t = max(-1.0, min(1.0, v / vmax))
    if t >= 0: r, g, b = 220, int(220 - 160 * t), int(220 - 170 * t)
    else: r, g, b = int(220 + 150 * t), int(220 + 100 * t), 220
    return f"rgb({r},{g},{b})"


def _node_trace(fig, nodes, pos, rad, fill, line_col, hover, name, line_w=1.5, dash=None, show=True, focus=None, values=None):
    """focus = set of highlighted genes (None = no focus mode). In focus mode non-focus nodes are small grey dots without labels."""
    if not nodes: return
    if focus is None:
        fig.add_trace(go.Scatter(x=[pos[n["id"]][0] for n in nodes], y=[pos[n["id"]][1] for n in nodes], mode="markers+text",
                                 text=[n["id"] for n in nodes], textposition="middle right", textfont=dict(size=10),
                                 marker=dict(size=[rad(n["id"]) * 7 for n in nodes], color=[fill(n) for n in nodes],
                                             line=dict(color=[line_col(n) for n in nodes], width=line_w)),
                                 hovertext=[hover(n) for n in nodes], hoverinfo="text", name=name, showlegend=show))
        return
    values = values or {}; vmax = max([abs(v) for v in values.values()] or [0])
    hi = [n for n in nodes if n["id"] in focus or n.get("cls") == "seed"]; lo = [n for n in nodes if n not in hi]
    if lo:
        fig.add_trace(go.Scatter(x=[pos[n["id"]][0] for n in lo], y=[pos[n["id"]][1] for n in lo], mode="markers",
                                 marker=dict(size=6, color="#d5d5d5", line=dict(width=0)), hovertext=[hover(n) for n in lo], hoverinfo="text",
                                 name=f"{name} (not in focus list)", showlegend=show))
    if hi:
        col = [(_value_colour(values[n["id"]], vmax) if n["id"] in values else fill(n)) for n in hi]
        fig.add_trace(go.Scatter(x=[pos[n["id"]][0] for n in hi], y=[pos[n["id"]][1] for n in hi], mode="markers+text",
                                 text=[f"<b>{n['id']}</b>" + (f" {values[n['id']]:+.2f}" if n["id"] in values else "") for n in hi],
                                 textposition="middle right", textfont=dict(size=11),
                                 marker=dict(size=[rad(n["id"]) * 9 for n in hi], color=col, line=dict(color="#111", width=2)),
                                 hovertext=[hover(n) + (f"<br>focus value {values[n['id']]:+.3f}" if n["id"] in values else "<br>focus gene") for n in hi],
                                 hoverinfo="text", name=f"{name} — focus" + (" (colour = value)" if values else ""), showlegend=show))


def _focus_edge_style(focus):
    """Edges touching a focus gene keep their colour; others fade."""
    if focus is None: return None
    def col(e):
        a, b = e.get("source", e.get("from_node")), e.get("target", e.get("to_node"))
        base = COL["amb"] if e.get("effect", 1) == 0 else COL["edge"]
        return base if (a in focus or b in focus) else "#e3e3e3"
    return col


# --------------------------------------------------------------------------- flow (upstream / downstream)
def flow_figure(json_path, width=1000, height=650, max_nodes=150, focus=None):
    """focus: genes to highlight (list / file / dict gene->value); other nodes become small grey dots."""
    F = _load(json_path)
    fset, fvals = parse_focus(focus); fset_or_none = fset if focus is not None else None
    nodes, edges, capped = _cap_nodes(F["nodes"], F["edges"], max_nodes, keep_ids=fset)
    F["edges"] = edges; byid = {n["id"]: n for n in nodes}
    if capped: F["subtitle"] = f"showing top {max_nodes} of {len(F['nodes'])} nodes by path score (full table in CSV / D3 view). " + F.get("subtitle", "")
    F["genes"] = [g for g in F.get("genes", []) if g["tf"] in byid]
    # columns: intermediates by hop, endpoints in one final column; upstream mode is drawn right-to-left
    max_mid = max([n["hops"] for n in nodes if not n.get("is_endpoint")] or [0])
    for n in nodes: n["col"] = max_mid + 1 if n.get("is_endpoint") else n["hops"]
    ncol = max_mid + 2
    col_x = lambda c: (c + 0.5) / ncol * 100
    if F["mode"] == "upstream": col_x = lambda c, f=col_x: 100 - f(c)
    pos = {}
    for c in range(ncol):
        arr = [n for n in nodes if n["col"] == c]
        if c > 0:   # barycenter ordering on already placed neighbours
            for n in arr:
                prev = [pos[e["source"]][1] for e in F["edges"] if e["target"] == n["id"] and e["source"] in pos]
                n["_bc"] = sum(prev) / len(prev) if prev else 50
            arr.sort(key=lambda n: (n.get("_bc", 50), -n["score"]))
        step = 90 / max(len(arr), 1)
        for i, n in enumerate(arr): pos[n["id"]] = (col_x(c), 95 - step * (i + 0.5))
    rad = lambda nid: 2.2 if byid[nid]["cls"] == "seed" else 1.0 + 1.6 * byid[nid]["score"]
    fig = _base(f"{F['mode']} pathway of {F['seed']}  —  {F.get('subtitle', '')}", width, max(height, 14 * max(len(a) for a in [[n for n in nodes if n['col'] == c] for c in range(ncol)])))
    _edge_traces(fig, F["edges"], pos, rad, color_fn=_focus_edge_style(fset_or_none))
    if focus is not None: fig.update_layout(title_text=fig.layout.title.text + f"  |  focus: {len(fset & set(byid))}/{len(fset)} focus genes in view")
    hover = lambda n: (f"<b>{n['id']}</b><br>hops {n['hops']} · path score {n['score']:.2f}<br>net sign: "
                       f"{'activated' if n['sign'] == 1 else 'inhibited' if n['sign'] == -1 else 'mixed'}"
                       + ("<br>endpoint" if n.get("is_endpoint") else "") + ("<br>feedback node" if n.get("feedback") else "")
                       + (f"<br>co-essentiality r {n['coess']}" if n.get("coess") not in (None, "") else "") + f"<br>{n.get('best_path', '')}")
    outline = lambda n: COL["syn"] if n["sign"] == 1 else COL["ant"] if n["sign"] == -1 else COL["mix"]
    fk = dict(focus=fset_or_none, values=fvals)
    _node_trace(fig, [n for n in nodes if n["cls"] == "seed"], pos, rad, lambda n: COL["seed"], lambda n: "#111", hover, "seed", **fk)
    _node_trace(fig, [n for n in nodes if n.get("is_endpoint")], pos, rad, lambda n: COL["endpoint"], outline, hover,
                "endpoint (receptor / TF); outline green = activated, red = inhibited", line_w=2.2, **fk)
    _node_trace(fig, [n for n in nodes if n["cls"] != "seed" and not n.get("is_endpoint")], pos, rad, lambda n: COL["node"], outline, hover, "intermediate", line_w=2.2, **fk)
    # gene summary per TF (downstream)
    if F.get("genes"):
        per = defaultdict(list)
        for g in F["genes"]: per[g["tf"]].append(g)
        xs, ys, txt, hv = [], [], [], []
        for tf, gs in per.items():
            if tf not in pos: continue
            up = sum(g["sign"] == 1 for g in gs); dn = sum(g["sign"] == -1 for g in gs)
            xs.append(pos[tf][0] + 13); ys.append(pos[tf][1]); txt.append(f"{len(gs)} genes ↑{up} ↓{dn}")
            hv.append(f"<b>{tf} targets ({F.get('regulon', '')})</b><br>" + "<br>".join(f"{g['gene']} {'↑' if g['sign'] == 1 else '↓' if g['sign'] == -1 else '±'} ({g['confidence']})" for g in gs[:30]))
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="text", text=txt, textposition="middle right", textfont=dict(size=9, color="#444"),
                                 hovertext=hv, hoverinfo="text", name="target genes (hover)", showlegend=True))
    for c in range(ncol):
        fig.add_annotation(x=col_x(c), y=99, text=("seed" if c == 0 else ("receptors" if F["mode"] == "upstream" else "TFs") if c == ncol - 1 else f"hop {c}"),
                           showarrow=False, font=dict(size=10, color="#888"))
    fig.update_xaxes(range=[-2, 130 if F.get("genes") else 104]); fig.update_yaxes(range=[0, 102])
    return fig


# --------------------------------------------------------------------------- cell map (location view)
_ZONES = dict(membrane=(50, 50, 46, 50), cytosol=(50, 50, 22, 32), nucleus=(50, 50, 0, 13), lysosome=(20, 32, 0, 9),
              mitochondrion=(80, 32, 0, 9), ER_golgi=(76, 74, 0, 9))   # cx, cy, r_inner, r_outer


def cell_map_figure(json_path, width=900, height=700, max_nodes=150, focus=None):
    N = _load(json_path); fset, fvals = parse_focus(focus); fset_or_none = fset if focus is not None else None
    nodes, edges, capped = _cap_nodes(N["nodes"], N["edges"], max_nodes, keep_ids=fset); N["edges"] = edges
    byid = {n["id"]: n for n in nodes}
    if capped: N["mode"] = f"{N.get('mode', '')} (top {max_nodes} of {len(N['nodes'])} nodes)"
    pos = {}
    for comp in set(n["compartment"] for n in nodes):
        cx, cy, ri, ro = _ZONES.get(comp, _ZONES["cytosol"]); arr = [n for n in nodes if n["compartment"] == comp]
        for i, n in enumerate(arr):   # deterministic spiral inside the zone
            a = i * 2.399; r = ri + (ro - ri) * (0.15 + 0.8 * ((i % 7) / 7)) if ro > ri else 0
            pos[n["id"]] = (cx + r * math.cos(a) * 1.0, cy + r * math.sin(a) * 0.85)
    rad = lambda nid: 2.0 if byid[nid]["cls"] == "seed" else 0.9 + 1.4 * byid[nid]["score"]
    fig = _base(f"{N.get('mode', '')} of {N['seed']} (location view)", width, height)
    for comp, (cx, cy, ri, ro) in _ZONES.items():
        if comp == "cytosol": continue
        fig.add_shape(type="circle", x0=cx - ro, y0=cy - ro * 0.85, x1=cx + ro, y1=cy + ro * 0.85, line=dict(color="#999", width=1.2))
        fig.add_annotation(x=cx, y=cy + ro * 0.85 + 2, text=comp.replace("_", "/"), showarrow=False, font=dict(size=10, color="#888"))
    fig.add_annotation(x=50, y=84, text="cytosol", showarrow=False, font=dict(size=10, color="#888"))
    _edge_traces(fig, [e for e in N["edges"] if e["kind"] == "signaling"], pos, rad, color_fn=_focus_edge_style(fset_or_none))
    phys = [e for e in N["edges"] if e["kind"] != "signaling"]
    if phys:
        xs, ys = [], []
        for e in phys:
            if e["source"] in pos and e["target"] in pos:
                xs += [pos[e["source"]][0], pos[e["target"]][0], None]; ys += [pos[e["source"]][1], pos[e["target"]][1], None]
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=COL["phys"], width=1, dash="dot"), hoverinfo="skip", name="physical contact"))
    cmap = dict(seed=COL["seed"], receptor=COL["endpoint"], tf=COL["shared"], upstream=COL["causal"], downstream=COL["node"])
    hover = lambda n: f"<b>{n['id']}</b><br>{n['cls']} · {n['compartment']}<br>hops {n.get('hops', 0)} · score {n['score']:.2f}<br>{n.get('note', '')}"
    for cls in ["seed", "receptor", "tf", "upstream", "downstream"]:
        _node_trace(fig, [n for n in nodes if n["cls"] == cls], pos, rad, lambda n, c=cls: cmap[c],
                    lambda n: "#999" if n.get("sign") == 0 else "#111", hover, cls, line_w=1, focus=fset_or_none, values=fvals)
    fig.update_xaxes(range=[0, 100], scaleanchor="y"); fig.update_yaxes(range=[0, 100])
    return fig


# --------------------------------------------------------------------------- cross-talk (three lanes)
def crosstalk_figure(json_path, width=1050, height=700, max_nodes=200, focus=None):
    X = _load(json_path); fset, fvals = parse_focus(focus); fset_or_none = fset if focus is not None else None
    keep = {e["from_node"] for e in X["interference"]} | {e["to_node"] for e in X["interference"]} | fset
    nodes, edges, capped = _cap_nodes(X["nodes"], X["edges"], max_nodes, keep_ids=keep); X["edges"] = edges
    byid = {n["id"]: n for n in nodes}
    lane_x = dict(causal=18, shared=50, goi=82)
    max_hop = max([max(n.get("hops_c", 0) or 0, n.get("hops_x", 0) or 0) for n in nodes] or [1])
    groups = defaultdict(list)
    for n in nodes: groups[(n["lane"], n.get("hops_x", 0) if n["lane"] == "goi" else n.get("hops_c", 0))].append(n)
    pos = {}
    for (lane, hp), arr in groups.items():
        arr.sort(key=lambda n: -n["score"]); span = 26 if lane != "shared" else 24
        for i, n in enumerate(arr):
            x = lane_x[lane] if len(arr) == 1 else lane_x[lane] - span / 2 + span * (i + 0.5) / len(arr)
            pos[n["id"]] = (x, 94 - 84 * (hp or 0) / (max_hop + 1) + (1.5 if i % 2 else -1.5))
    rad = lambda nid: 2.0 if byid[nid]["cls"] == "seed" else 0.9 + 1.2 * byid[nid]["score"]
    s = X["summary"]
    fig = _base((f"(top {max_nodes} of {len(X['nodes'])} nodes) " if capped else "") + f"Cross-talk: does {X['goi']} influence the {X['causal']} pathway?  index {s['crosstalk_index']} · "
                f"interference edges {s['n_interference_edges']} · endpoints syn/ant/mixed {s['synergistic']}/{s['antagonistic']}/{s['mixed']}", width, height)
    _edge_traces(fig, X["edges"], pos, rad, color_fn=(lambda e: "#bbbbbb") if focus is None else (lambda e: "#bbbbbb" if (e["source"] in fset or e["target"] in fset) else "#eeeeee"),
                 dash_fn=lambda e: "dot" if e.get("effect", 1) == 0 else "solid", width=0.8)
    inter = [dict(source=e["from_node"], target=e["to_node"], effect=e["edge_effect"], agreement=e["agreement"], mech=e.get("mechanism", ""), res=e.get("residue", "")) for e in X["interference"]]
    _edge_traces(fig, inter, pos, rad, color_fn=lambda e: COL["syn"] if e["agreement"] == "synergistic" else COL["ant"] if e["agreement"] == "antagonistic" else COL["mix"],
                 dash_fn=lambda e: "longdash" if "chaperone" in (e.get("mech") or "") else "solid", width=2.4)
    # invisible midpoints carry the hover text of interference edges
    mx = [(pos[e["source"]][0] + pos[e["target"]][0]) / 2 for e in inter if e["source"] in pos and e["target"] in pos]
    my = [(pos[e["source"]][1] + pos[e["target"]][1]) / 2 for e in inter if e["source"] in pos and e["target"] in pos]
    ht = [f"<b>{e['source']} → {e['target']}</b><br>{e['res'] or e['mech']}<br>{e['agreement']}" for e in inter if e["source"] in pos and e["target"] in pos]
    fig.add_trace(go.Scatter(x=mx, y=my, mode="markers", marker=dict(size=9, color="rgba(0,0,0,0)"), hovertext=ht, hoverinfo="text",
                             name="interference edge (hover midpoint): green synergistic, red antagonistic, long-dash chaperone", showlegend=True))
    hover = lambda n: (f"<b>{n['id']}</b> · {n['lane']}" + (f"<br>hops from {X['causal']}: {n['hops_c']}" if n.get("hops_c") is not None else "")
                       + (f"<br>hops from {X['goi']}: {n['hops_x']}" if n.get("hops_x") is not None else "") + (f"<br>{n['agreement']}" if n.get("agreement") else ""))
    for lane, name in [("causal", f"{X['causal']} pathway"), ("shared", "shared"), ("goi", f"{X['goi']} pathway")]:
        _node_trace(fig, [n for n in nodes if n["lane"] == lane and n["cls"] != "seed"], pos, rad, lambda n, l=lane: COL[l],
                    lambda n: "#111" if n.get("is_endpoint") else "rgba(0,0,0,0)", hover, name, line_w=1.5, focus=fset_or_none, values=fvals)
        fig.add_annotation(x=lane_x[lane], y=99, text=name.upper(), showarrow=False, font=dict(size=11, color="#888"))
    _node_trace(fig, [n for n in nodes if n["cls"] == "seed"], pos, rad, lambda n: COL["seed"], lambda n: "#111", hover, "seed genes", focus=fset_or_none, values=fvals)
    fig.update_xaxes(range=[0, 100]); fig.update_yaxes(range=[0, 102])
    return fig


# --------------------------------------------------------------------------- synthetic lethality (volcano)
def volcano_figure(json_path, width=850, height=600, label_top=20):
    S = _load(json_path); pts = [p for p in S["points"] if p["padj"] and p["padj"] > 0]
    y = lambda p: -math.log10(max(p["padj"], 1e-300))
    fig = _base(f"Synthetic lethality with {S['seed']} loss" + ("  (DEMO: simulated values)" if S.get("demo") else ""), width, height)
    fig.update_layout(showlegend=True, xaxis=dict(visible=True, title=f"Δ gene effect ({S['seed']}-lost − intact)", gridcolor="#eee", zerolinecolor="#ddd"),
                      yaxis=dict(visible=True, title="−log10 FDR", gridcolor="#eee"))
    for key, name, col, size in [("other", "other genes", "#cccccc", 5), ("on", "on the pathway network", COL["causal"], 8), ("sig", "synthetic-lethal hit", COL["goi"], 10)]:
        sub = [p for p in pts if (key == "sig" and p["significant"]) or (key == "on" and p["on_pathway"] and not p["significant"]) or (key == "other" and not p["significant"] and not p["on_pathway"])]
        fig.add_trace(go.Scatter(x=[p["delta"] for p in sub], y=[y(p) for p in sub], mode="markers", name=name,
                                 marker=dict(size=size, color=col, opacity=0.6 if key == "other" else 0.95, line=dict(width=0.5, color="#333") if key == "sig" else None),
                                 hovertext=[f"<b>{p['gene']}</b><br>Δ {p['delta']:.3f} · FDR {p['padj']:.1e}<br>co-ess r {p.get('coessentiality_r')}" for p in sub], hoverinfo="text"))
    lab = sorted([p for p in pts if p["significant"]], key=lambda p: p["delta"])[:label_top]
    fig.add_trace(go.Scatter(x=[p["delta"] for p in lab], y=[y(p) for p in lab], mode="text", text=[p["gene"] for p in lab], textposition="top right",
                             textfont=dict(size=9), hoverinfo="skip", showlegend=False))
    th = S["thresholds"]
    fig.add_hline(y=-math.log10(th["padj"]), line=dict(color="#999", dash="dash", width=1)); fig.add_vline(x=-th["delta"], line=dict(color="#999", dash="dash", width=1))
    return fig

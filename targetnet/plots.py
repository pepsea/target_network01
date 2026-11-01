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
import json, math, os
import pandas as pd
from collections import defaultdict
import plotly.graph_objects as go

# flat, muted palette on white (owner's convention)
COL = dict(seed="#1a1a1a", endpoint="#2d5f8f", ligand="#3f8f8f", node="#8a8a8a", shared="#6b5b95", causal="#4a7fb5", goi="#b36b3c",
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


CAPTIONS = {
 "impact": [
  "IMPACT VIEW  the core gene's own downstream pathway (walked with the target gene removed), recoloured by what the target gene does to each node.",
  "NODES  green = target activates this node directly (interference edge into it), red = inhibits directly, dark grey = mixed; purple = shared (target's own pathway also reaches it);",
  "       light green / light red = indirect (inherits the effect along the core pathway's own edges); pale grey = not influenced at this depth.",
  "OUTLINE green = target pushes the node the same way as the core gene (synergistic), red = opposite way (antagonistic), grey = mixed. Hover: which target-side node acts on it.",
  "LINES  the core pathway's own signaling edges (▶ activation, ┃ inhibition). Columns = hops from the core gene. Full table: impact.csv (node, impact, net_sign_from_goi, agreement, via).",
  "READ   red OUTLINE = the target opposes the core gene there (inhibiting the target releases that node); green outline = same direction. Direct nodes at hop 1-2 are the mechanistically closest points; the TF column shows where the influence reaches transcription."],
 "flow": [
  "LINES  dark grey curve = signaling edge from the signed graph (SIGNOR direct relations, plus custom_edges if configured; Reactome FI only if mode=directed or as seed rescue — hover an edge for its source); it runs from regulator to target.",
  "       ▶ arrowhead = activation / up-regulation · ┃ bar head = inhibition / down-regulation · dotted grey = relation of unknown sign.",
  "       Line thickness ∝ edge confidence (SIGNOR score). Only edges among drawn nodes are shown; the walk itself used the whole graph.",
  "NODES  black diamond = the analysed gene (seed) · blue = endpoint (receptor in upstream view, transcription factor in downstream view) · grey = intermediate.",
  "       Fill saturation ∝ path score. Outline: green = seed activates this node (consensus over all paths), red = inhibits, grey dashed = mixed; long dashes = feedback node (also upstream of the seed).",
  "COLUMNS hop N = shortest number of edges from the seed; when the walk exceeds the node cap, each hop shows its highest-scoring genes (FOCUS genes are highlighted when they are among them; all genes are in the tables); endpoints are collected in the last column(s): upstream = receptors, then their secreted ligands (teal; UniProt Secreted + Cytokine/Growth factor/Hormone); downstream = TFs. Right-hand boxes: DoRothEA target genes per TF (↑ induced, ↓ repressed by the seed); hover to list them.",
  "FOCUS  with a focus list: focus genes are orange with bold dark-orange labels (blue→red fill instead when a value such as log2FC is supplied); all other nodes small grey; edges not touching a focus gene are lighter. Focus genes outside the drawn high-score set are listed in the title, not forced into the figure."],
 "cell_map": [
  "LINES  dark grey curve with ▶ (activation) or ┃ (inhibition) = signaling edge; dotted light grey straight line = physical contact only (STRING / Reactome undirected), no direction or sign.",
  "NODES  placed by UniProt subcellular location (secreted ligands outside the membrane ring, membrane ring, cytosol ring, nucleus, lysosome, mitochondrion, ER/Golgi); position inside a zone is arbitrary.",
  "       black diamond = analysed gene · blue = receptor · purple = TF · light blue = upstream node · grey = downstream node; size ∝ path score; grey outline = mixed sign."],
 "crosstalk": [
  "LANES  left = nodes reached only from the core (causal) gene · centre = reached from both · right = reached only from the gene of interest (GOI). Row ≈ hops from the respective seed.",
  "THIN GREY LINES = ordinary signaling edges of the two pathways (▶ activation, ┃ inhibition, dotted = unknown sign).",
  "DOTTED BLUE straight line = physical contact between the GOI side and the core pathway with NO signed edge in any database (STRING experimental/database channels, Reactome undirected) — e.g. chaperone-client; listed in physical_contacts.csv.",
  "COLOURED LINES = interference edges: a GOI-side node acting on a core-pathway node. green = pushes that node the same way as the core gene (synergistic) · red = opposite way (antagonistic) · grey = undetermined.",
  "       long-dashed = chaperone / protein-stability edge rather than phosphorylation. Only the top-ranked edges (by proximity to both genes) are drawn; all are in interference.csv. Hover the midpoint for residue and mechanism.",
  "NODES  black diamonds = the two analysed genes (core, target) · purple = shared · light blue = core pathway · orange = GOI pathway; black outline = TF endpoint; size ∝ path score.",
  "TITLE  TF-overlap enrichment = observed / expected overlap of the two TF sets given how many TFs the GOI reaches at all (hypergeometric p)."],
 "volcano": [
  "x = Δ CRISPR gene effect (seed-lost − seed-intact cell lines; more negative = more essential when the seed is lost) · y = −log10 FDR.",
  "orange = passes both thresholds (synthetic-lethal candidate) · blue = gene already on the seed's signaling network (context, not SL class) · grey = other genes. Dashed lines = thresholds."],
}


def add_caption(fig, kind, extra=None):
    """Write the key to the figure's lines and nodes under the plot (travels with PNG exports)."""
    lines = CAPTIONS.get(kind, []) + (extra or [])
    txt = "<br>".join(lines)
    fig.update_layout(margin=dict(b=40 + 16 * len(lines)))
    fig.add_annotation(xref="paper", yref="paper", x=0, y=-0.02 - 0.0 , xanchor="left", yanchor="top", align="left", showarrow=False,
                       text=f"<span style='font-size:10px;color:#555'>{txt}</span>", bgcolor="rgba(255,255,255,0)")
    fig.update_layout(height=fig.layout.height + 16 * len(lines) + 30, legend=dict(orientation="h", x=0, y=-0.02 - (16 * len(lines) + 20) / fig.layout.height))
    return fig


def legend_text(kind):
    """Plain-text key for a view kind ('flow' | 'cell_map' | 'crosstalk' | 'volcano'), e.g. for notebooks or reports."""
    return "\n".join(CAPTIONS.get(kind, []))


def _edge_legend(fig, kinds=("act", "inh", "amb")):
    """Legend-only traces explaining the line styles (drawn off-canvas)."""
    items = {"act": ("signaling edge, activation (▶)", dict(color=COL["edge"], width=1.6, dash="solid")),
             "inh": ("signaling edge, inhibition (┃)", dict(color=COL["edge"], width=1.6, dash="solid")),
             "amb": ("signaling edge, unknown sign", dict(color=COL["amb"], width=1.6, dash="dot")),
             "phys": ("physical contact only (no direction)", dict(color=COL["phys"], width=1.2, dash="dot")),
             "syn": ("interference: synergistic", dict(color=COL["syn"], width=1.6)),
             "ant": ("interference: antagonistic", dict(color=COL["ant"], width=1.6)),
             "chap": ("interference via chaperone / stability", dict(color=COL["mix"], width=1.6, dash="longdash")),
             "faded": ("edge not touching a focus gene", dict(color="#c9c9c9", width=1.6))}
    for k in kinds:
        name, line = items[k]
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines", line=line, name=name, showlegend=True, hoverinfo="skip"))


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


def _cap_nodes(nodes, edges, max_nodes, keep_ids=(), per_hop=False):
    """Cap the drawn nodes at max_nodes so the view stays readable; the full tables keep everything.
    per_hop=True (gene research flow / cell map): the budget is shared across hop columns and filled with the HIGHEST-SCORING genes
    of each hop; FOCUS genes are highlighted only if they are among them (keep_ids is not forced).
    per_hop=False (legacy): seeds and keep_ids first, then hop-1 by score, then best paths of the top endpoints, then best remaining."""
    if len(nodes) <= max_nodes: return nodes, edges, False
    byid = {n["id"]: n for n in nodes}
    keep = {n["id"] for n in nodes if n.get("cls") == "seed"}
    if per_hop:
        keep |= {n["id"] for n in nodes if n["id"] in keep_ids}          # analysis results that the view exists to show (e.g. interference-edge ends); never FOCUS
        cols = defaultdict(list)
        for n in nodes:
            if n["id"] in keep: continue
            lane = n.get("lane", "")
            hp = n.get("hops_x", 0) if lane == "goi" else n.get("hops_c", n.get("hops", 0))
            cols[(lane, "end" if (n.get("is_endpoint") and not lane) else hp)].append(n)
        for arr in cols.values(): arr.sort(key=lambda n: -n.get("score", 0))
        budget = max_nodes - len(keep)
        quota = {k: 0 for k in cols}
        remaining = budget
        while remaining > 0 and any(quota[k] < len(cols[k]) for k in cols):      # round-robin over columns, best score first
            for k in cols:
                if remaining <= 0: break
                if quota[k] < len(cols[k]): quota[k] += 1; remaining -= 1
        for k, arr in cols.items(): keep |= {n["id"] for n in arr[:quota[k]]}
    else:
        keep |= {n["id"] for n in nodes if n["id"] in keep_ids}
        hop1 = lambda n: 1 in (n.get("hops"), n.get("hops_c"), n.get("hops_x"))
        for n in sorted([n for n in nodes if hop1(n)], key=lambda n: -n.get("score", 0)):
            if len(keep) >= max_nodes: break
            keep.add(n["id"])
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
    """focus: None | plain text ("GSK3B, FOXO3 MTOR" or one gene per line, optional value after a gene: "GSK3B -1.8")
       | path to such a text file | list/set of symbols | dict {gene: value} | pandas Series (index = gene).
       Returns (set_of_genes, {gene: value} or {})."""
    if focus is None: return set(), {}
    try:
        import pandas as pd
        if isinstance(focus, pd.Series): focus = focus.to_dict()
        elif isinstance(focus, pd.DataFrame): focus = dict(zip(focus.iloc[:, 0], focus.iloc[:, 1])) if focus.shape[1] > 1 else list(focus.iloc[:, 0])
    except ImportError:
        pass
    if isinstance(focus, str):
        import os, re
        text = open(os.path.expanduser(focus)).read() if os.path.exists(os.path.expanduser(focus)) else focus
        # plain text: genes separated by newlines / commas / spaces / semicolons; an optional number after a gene is its value
        vals, last = {}, None
        for tok in re.split(r"[\s,;]+", text):
            if not tok or tok.startswith("#"): continue
            try:
                v = float(tok)
                if last is not None: vals[last] = v
                continue
            except ValueError:
                pass
            last = tok.strip().upper(); vals.setdefault(last, None)
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


FOCUS_COL = "#E8862E"   # focus-gene fill (muted orange, distinct from the blue/grey class colours)
FOCUS_TXT = "#9A4B00"   # focus-gene label and outline


_is_seed = lambda n: n.get("cls") == "seed"
_sym = lambda nodes: ["diamond" if _is_seed(n) else "circle" for n in nodes]          # analysed gene(s): black diamond, bold label
_lbl = lambda n, extra="": (f"<b>{n['id']}</b>" if _is_seed(n) else n["id"]) + extra


def _node_trace(fig, nodes, pos, rad, fill, line_col, hover, name, line_w=1.5, dash=None, show=True, focus=None, values=None, label_size=10):
    """focus = set of highlighted genes (None = no focus mode). In focus mode non-focus nodes are small grey dots without labels.
    The analysed gene(s) (cls == 'seed') are always drawn as a black diamond with a bold label, in every view."""
    if not nodes: return
    if focus is None:
        fig.add_trace(go.Scatter(x=[pos[n["id"]][0] for n in nodes], y=[pos[n["id"]][1] for n in nodes], mode="markers+text",
                                 text=[_lbl(n) for n in nodes], textposition="middle right",
                                 textfont=dict(size=[label_size + 2 if _is_seed(n) else label_size for n in nodes]),
                                 marker=dict(size=[rad(n["id"]) * 7 for n in nodes], color=[fill(n) for n in nodes], symbol=_sym(nodes),
                                             line=dict(color=[line_col(n) for n in nodes], width=line_w)),
                                 hovertext=[hover(n) for n in nodes], hoverinfo="text", name=name, showlegend=show))
        return
    values = values or {}; vmax = max([abs(v) for v in values.values()] or [0])
    hi = [n for n in nodes if n["id"] in focus or n.get("cls") == "seed"]; lo = [n for n in nodes if n not in hi]
    if lo:   # context nodes: small grey circles, still labelled in small grey text so the pathway stays readable
        fig.add_trace(go.Scatter(x=[pos[n["id"]][0] for n in lo], y=[pos[n["id"]][1] for n in lo], mode="markers+text",
                                 text=[n["id"] for n in lo], textposition="middle right", textfont=dict(size=8, color="#8a8a8a"),
                                 marker=dict(size=[max(7, rad(n["id"]) * 4) for n in lo], color="#c4c4c4", line=dict(color="#9a9a9a", width=0.8)),
                                 hovertext=[hover(n) for n in lo], hoverinfo="text", name=f"{name} (not in focus list)", showlegend=show))
    if hi:
        # focus genes: dedicated highlight colour (FOCUS_COL) so they stand out from every class colour; with a supplied value the
        # fill becomes the blue->red value scale instead; the seed keeps its own colour
        def fcol(n):
            if n["id"] in values: return _value_colour(values[n["id"]], vmax)
            return fill(n) if (n.get("cls") == "seed" and n["id"] not in focus) else FOCUS_COL
        col = [fcol(n) for n in hi]
        fig.add_trace(go.Scatter(x=[pos[n["id"]][0] for n in hi], y=[pos[n["id"]][1] for n in hi], mode="markers+text",
                                 text=[f"<b>{n['id']}</b>" + (f" {values[n['id']]:+.2f}" if n["id"] in values else "") for n in hi],
                                 textposition="middle right", textfont=dict(size=11, color=[FOCUS_TXT if n["id"] in focus else "#111" for n in hi]),
                                 marker=dict(size=[rad(n["id"]) * 6.5 for n in hi], color=col, symbol=_sym(hi),
                                             line=dict(color=[FOCUS_TXT if n["id"] in focus else "#111" for n in hi], width=1.8)),
                                 hovertext=[hover(n) + (f"<br>focus value {values[n['id']]:+.3f}" if n["id"] in values else "<br>focus gene") for n in hi],
                                 hoverinfo="text", name=f"{name} — focus" + (" (colour = value)" if values else ""), showlegend=show))


def _focus_edge_style(focus):
    """Edges touching a focus gene keep their colour; others fade."""
    if focus is None: return None
    def col(e):
        a, b = e.get("source", e.get("from_node")), e.get("target", e.get("to_node"))
        base = COL["amb"] if e.get("effect", 1) == 0 else COL["edge"]
        return base if (a in focus or b in focus) else "#c9c9c9"
    return col


# --------------------------------------------------------------------------- flow (upstream / downstream)
def _view_cfg(key, default):
    """Figure budgets from the YAML `view:` section (max_nodes, max_nodes_crosstalk, max_edges_crosstalk, rows_per_column)."""
    try:
        from . import data as D
        v = (D.CFG or {}).get("view") or {}
        return int(v.get(key, default))
    except Exception:
        return default


def flow_figure(json_path, width=1000, height=650, max_nodes=None, focus=None, caption=True, sub_max=3):
    """focus: genes to highlight (list / file / dict gene->value); other nodes become small grey dots."""
    F = _load(json_path)
    if max_nodes is None: max_nodes = _view_cfg("max_nodes", 150)
    fset, fvals = parse_focus(focus); fset_or_none = fset if focus is not None else None
    nodes, edges, capped = _cap_nodes(F["nodes"], F["edges"], max_nodes, per_hop=True)
    F["edges"] = edges; byid = {n["id"]: n for n in nodes}
    if capped: F["subtitle"] = f"showing the {max_nodes} highest-scoring of {len(F['nodes'])} nodes, per hop (full table in CSV). " + F.get("subtitle", "")
    F["genes"] = [g for g in F.get("genes", []) if g["tf"] in byid]
    # columns: intermediates by hop, endpoints in one final column; upstream mode is drawn right-to-left
    max_mid = max([n["hops"] for n in nodes if not n.get("is_endpoint")] or [0])
    is_lig = lambda n: n.get("endpoint_type") == "ligand"
    has_lig = any(is_lig(n) for n in nodes)
    for n in nodes: n["col"] = max_mid + 2 if is_lig(n) else max_mid + 1 if n.get("is_endpoint") else n["hops"]
    ncol = max_mid + (3 if has_lig else 2)
    col_x = lambda c: (c + 0.5) / ncol * 100
    if F["mode"] == "upstream": col_x = lambda c, f=col_x: 100 - f(c)
    # a hop with many genes is wrapped into side-by-side sub-columns (at most ROWS genes each) instead of one very tall column
    ROWS, SUBMAX = _view_cfg("rows_per_column", 28), sub_max
    slot = 100 / ncol
    pos = {}; nsub = {}
    for c in range(ncol):
        arr = [n for n in nodes if n["col"] == c]
        if c > 0:   # barycenter ordering on already placed neighbours
            for n in arr:
                prev = [pos[e["source"]][1] for e in F["edges"] if e["target"] == n["id"] and e["source"] in pos]
                n["_bc"] = sum(prev) / len(prev) if prev else 50
            arr.sort(key=lambda n: (n.get("_bc", 50), -n["score"]))
        k = max(1, min(SUBMAX, math.ceil(len(arr) / ROWS))); nsub[c] = k
        per = math.ceil(len(arr) / k)
        for j in range(k):
            part = arr[j * per:(j + 1) * per]
            x = col_x(c) + (0 if k == 1 else (j - (k - 1) / 2) * slot * 0.75 / (k - 1 if k > 1 else 1) * (1 if F["mode"] != "upstream" else -1))
            step = 90 / max(len(part), 1)
            for i, n in enumerate(part): pos[n["id"]] = (x, 95 - step * (i + 0.5))
    tallest = max(min(ROWS, math.ceil(len([n for n in nodes if n["col"] == c]) / nsub[c])) for c in range(ncol))
    wrapped = any(k > 1 for k in nsub.values())
    label_size = 9 if wrapped else 10
    rad = lambda nid: 2.2 if byid[nid]["cls"] == "seed" else 1.0 + 1.6 * byid[nid]["score"]
    title = (f"direct {'regulators' if F['mode'] == 'upstream' else 'targets'} of {F['seed']} (hop 1, all {len(nodes) - 1} genes)  —  {F.get('subtitle', '')}"
             if F.get("hop1_only") else f"{F['mode']} pathway of {F['seed']}  —  {F.get('subtitle', '')}")
    fig = _base(title, width if not wrapped else int(width * (1.0 + 0.15 * min(nsub_max := max(nsub.values()), 4))), max(height, 16 * tallest + 120))
    _edge_traces(fig, F["edges"], pos, rad, color_fn=_focus_edge_style(fset_or_none))
    if focus is not None:
        missing = sorted(fset - set(byid))
        fig.update_layout(title_text=fig.layout.title.text + f"  |  focus: {len(fset & set(byid))}/{len(fset)} focus genes among the drawn high-score nodes"
                          + (f" (not drawn / not reached: {', '.join(missing[:8])}{'…' if len(missing) > 8 else ''})" if missing else ""))
    hover = lambda n: (f"<b>{n['id']}</b><br>hops {n['hops']} · path score {n['score']:.2f}<br>net sign: "
                       f"{'activated' if n['sign'] == 1 else 'inhibited' if n['sign'] == -1 else 'mixed'}"
                       + ("<br>endpoint" if n.get("is_endpoint") else "") + ("<br>feedback node" if n.get("feedback") else "")
                       + (f"<br>co-essentiality r {n['coess']}" if n.get("coess") not in (None, "") else "") + f"<br>{n.get('best_path', '')}")
    outline = lambda n: COL["syn"] if n["sign"] == 1 else COL["ant"] if n["sign"] == -1 else COL["mix"]
    fk = dict(focus=fset_or_none, values=fvals)
    _node_trace(fig, [n for n in nodes if n["cls"] == "seed"], pos, rad, lambda n: COL["seed"], lambda n: "#111", hover, "seed", **fk)
    imp = F.get("impact")
    if imp and focus is None:
        # IMPACT MODE (crosstalk core_impact view): the core pathway recoloured by what the target gene does to each node
        X = F.get("impact_from", "target")
        IMP_COL = {("direct", 1): "#2e7d32", ("direct", -1): "#b23b3b", ("direct", 0): "#6d6d6d",
                   ("shared", 1): "#6b5b95", ("shared", -1): "#6b5b95", ("shared", 0): "#6b5b95",
                   ("indirect", 1): "#9fcf9f", ("indirect", -1): "#e3a3a3", ("indirect", 0): "#c9c9c9"}
        def icol(n):
            d = imp.get(n["id"]); 
            return "#e9e9e9" if not d or d["impact"] == "none" else IMP_COL.get((d["impact"], int(d["net"])), "#c9c9c9")
        def ihover(n):
            d = imp.get(n["id"]) or dict(impact="none", net=0, via="", agreement="")
            tag = {"direct": f"DIRECT: {X}-side node(s) {d['via']} act on it", "shared": f"SHARED: {X} reaches it independently",
                   "indirect": f"INDIRECT: inherits {X}'s effect along the core pathway", "none": f"not reached by {X}"}[d["impact"]]
            eff = {1: "activates", -1: "inhibits", 0: "mixed"}[int(d["net"])] if d["impact"] != "none" else ""
            return hover(n) + f"<br><b>{tag}</b>" + (f"<br>{X} {eff} this node -> {d['agreement']} with the core gene" if eff else "")
        AGR = {"synergistic": COL["syn"], "antagonistic": COL["ant"], "mixed": COL["mix"]}
        def iline(n):
            d = imp.get(n["id"]); return AGR.get(d["agreement"], "#999") if d and d["impact"] != "none" else "#bbb"
        groups = [("direct", f"direct target of {X}: fill green = {X} activates it, red = inhibits, grey = mixed; OUTLINE green = same direction as {F['seed']} (synergistic), red = opposite (antagonistic)", 3.0),
                  ("shared", f"shared: reached by {X}'s own pathway too (purple)", 2.4),
                  ("indirect", f"indirect: downstream of an influenced node inside the {F['seed']} pathway (light tint)", 1.6),
                  ("none", f"not influenced by {X}", 1.0)]
        for kind, name, lw in groups:
            sel = [n for n in nodes if n["cls"] != "seed" and (imp.get(n["id"], {"impact": "none"})["impact"] == kind)]
            _node_trace(fig, sel, pos, rad, icol, iline, ihover, name, line_w=lw)
        fig.update_layout(title_text=f"impact of {X} on the {F['seed']} pathway  —  " + F.get("subtitle", ""))
    else:
        _node_trace(fig, [n for n in nodes if n.get("is_endpoint") and not is_lig(n)], pos, rad, lambda n: COL["endpoint"], outline, hover,
                    "endpoint (receptor / TF); outline green = activated, red = inhibited", line_w=2.2, label_size=label_size, **fk)
        _node_trace(fig, [n for n in nodes if is_lig(n)], pos, rad, lambda n: COL["ligand"], outline, hover,
                    "ligand (secreted cytokine / growth factor / hormone)", line_w=2.2, label_size=label_size, **fk)
        _node_trace(fig, [n for n in nodes if n["cls"] != "seed" and not n.get("is_endpoint")], pos, rad, lambda n: COL["node"], outline, hover, "intermediate", line_w=2.2, label_size=label_size, **fk)
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
        label = "seed" if c == 0 else ("ligands" if (has_lig and c == ncol - 1) else ("receptors" if F["mode"] == "upstream" else "TFs") if c == ncol - (2 if has_lig else 1) else f"hop {c}")
        fig.add_annotation(x=col_x(c), y=99, text=label,
                           showarrow=False, font=dict(size=10, color="#888"))
    fig.update_xaxes(range=[-2, 130 if F.get("genes") else 104]); fig.update_yaxes(range=[0, 102])
    _edge_legend(fig, ("act", "inh", "amb") + (("faded",) if focus is not None else ()))
    if caption: add_caption(fig, "impact" if (imp and focus is None) else "flow")
    return fig


# --------------------------------------------------------------------------- cell map (location view)
_ZONES = dict(extracellular=(50, 50, 52, 56), membrane=(50, 50, 46, 50), cytosol=(50, 50, 22, 32), nucleus=(50, 50, 0, 13), lysosome=(20, 32, 0, 9),
              mitochondrion=(80, 32, 0, 9), ER_golgi=(76, 74, 0, 9))   # cx, cy, r_inner, r_outer


def cell_map_figure(json_path, width=900, height=700, max_nodes=None, focus=None, caption=True):
    N = _load(json_path); fset, fvals = parse_focus(focus); fset_or_none = fset if focus is not None else None
    if max_nodes is None: max_nodes = _view_cfg("max_nodes", 150)
    nodes, edges, capped = _cap_nodes(N["nodes"], N["edges"], max_nodes, per_hop=True); N["edges"] = edges
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
        if comp in ("cytosol", "extracellular"): continue      # extracellular = outside the membrane ring (no circle of its own)
        fig.add_shape(type="circle", x0=cx - ro, y0=cy - ro * 0.85, x1=cx + ro, y1=cy + ro * 0.85, line=dict(color="#999", width=1.2))
        fig.add_annotation(x=cx, y=cy + ro * 0.85 + 2, text=comp.replace("_", "/"), showarrow=False, font=dict(size=10, color="#888"))
    if any(n["compartment"] == "extracellular" for n in nodes):
        fig.add_annotation(x=50, y=1.5, text="extracellular (secreted ligands)", showarrow=False, font=dict(size=10, color="#888"))
    fig.add_annotation(x=50, y=84, text="cytosol", showarrow=False, font=dict(size=10, color="#888"))
    _edge_traces(fig, [e for e in N["edges"] if e["kind"] == "signaling"], pos, rad, color_fn=_focus_edge_style(fset_or_none))
    phys = [e for e in N["edges"] if e["kind"] != "signaling"]
    if phys:
        xs, ys = [], []
        for e in phys:
            if e["source"] in pos and e["target"] in pos:
                xs += [pos[e["source"]][0], pos[e["target"]][0], None]; ys += [pos[e["source"]][1], pos[e["target"]][1], None]
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=COL["phys"], width=1, dash="dot"), hoverinfo="skip", name="physical contact"))
    cmap = dict(seed=COL["seed"], receptor=COL["endpoint"], ligand=COL["ligand"], tf=COL["shared"], upstream=COL["causal"], downstream=COL["node"])
    hover = lambda n: f"<b>{n['id']}</b><br>{n['cls']} · {n['compartment']}<br>hops {n.get('hops', 0)} · score {n['score']:.2f}<br>{n.get('note', '')}"
    for cls in ["seed", "receptor", "ligand", "tf", "upstream", "downstream"]:
        _node_trace(fig, [n for n in nodes if n["cls"] == cls], pos, rad, lambda n, c=cls: cmap[c],
                    lambda n: "#999" if n.get("sign") == 0 else "#111", hover, cls, line_w=1, focus=fset_or_none, values=fvals)
    fig.update_xaxes(range=[-8, 108], scaleanchor="y"); fig.update_yaxes(range=[-2, 102])
    _edge_legend(fig, ("act", "inh", "amb", "phys"))
    if caption: add_caption(fig, "cell_map")
    return fig


# --------------------------------------------------------------------------- cross-talk screen (ranked targets)
def screen_figure(table_path, width=900, height=600, caption=True, label_top=25):
    """
    One point per candidate target gene screened against the core gene.
    x = interference score (edges from the target's pathway into the core pathway, proximity-weighted),
    y = TF-overlap enrichment (observed / expected shared TFs), size = dependency (share of the core pathway lost when the
    target is removed), fill = net direction (teal = synergistic, brown = antagonistic, grey = balanced / none).
    Top-right, large, strongly coloured = strong, directional cross-talk.
    """
    tab = pd.read_csv(table_path) if isinstance(table_path, str) else table_path
    if not len(tab): return _base("cross-talk screen: no targets scored", width, height)
    core = tab.attrs.get("core", "") if hasattr(tab, "attrs") else ""
    fig = _base(f"cross-talk screen{(' against ' + core) if core else ''}: {len(tab)} candidate targets", width, height)
    def col(d):
        if d > 0.2: return "#3f8f8f"
        if d < -0.2: return "#b36b3c"
        return "#9a9a9a"
    size = 8 + 30 * tab.dependency.clip(0, 1)
    txt = [r.target if r["rank"] <= label_top else "" for _, r in tab.iterrows()]
    hv = [f"<b>{r.target}</b>  rank {r['rank']}<br>interference {r.interference_score:.2f} ({r.n_interference} edges)"
          f"<br>TF enrichment {r.tf_enrichment:.2f}, p = {r.tf_p:.2g}<br>dependency {r.dependency:.2f}"
          f"<br>interference edges syn/ant {r.syn_edges}/{r.ant_edges}; shared endpoints syn/ant {r.syn_endpoints}/{r.ant_endpoints} -> direction {r.net_direction:+.2f}"
          + (f"<br>baseline z {r.baseline_z:+.2f} (pct {r.baseline_pct})" if pd.notna(r.get("baseline_z", float("nan"))) else "")
          + f"<br>{r.top_edges}" for _, r in tab.iterrows()]
    fig.add_trace(go.Scatter(x=tab.interference_score, y=tab.tf_enrichment, mode="markers+text", text=txt, textposition="top center",
                             textfont=dict(size=10), marker=dict(size=size, color=[col(d) for d in tab.net_direction], line=dict(color="#333", width=1)),
                             hovertext=hv, hoverinfo="text", name="candidate target", showlegend=False))
    for name, c in [("synergistic with core", "#3f8f8f"), ("antagonistic to core", "#b36b3c"), ("balanced / none", "#9a9a9a")]:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", marker=dict(size=10, color=c, line=dict(color="#333", width=1)), name=name))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", marker=dict(size=22, color="#fff", line=dict(color="#333", width=1)), name="size = dependency of core on target"))
    xmax = max(0.5, float(tab.interference_score.max()) * 1.15); ymax = max(1.5, float(tab.tf_enrichment.max()) * 1.2)
    fig.update_xaxes(visible=True, showticklabels=True, title="interference score (proximity-weighted edges into the core pathway)", showgrid=True, gridcolor="#eee", zeroline=False, range=[-0.03 * xmax, xmax], linecolor="#999")
    fig.update_yaxes(visible=True, showticklabels=True, title="TF-overlap enrichment (observed / expected shared TFs)", showgrid=True, gridcolor="#eee", zeroline=False, range=[-0.05 * ymax, ymax], linecolor="#999")
    fig.add_hline(y=1.0, line=dict(color="#bbb", width=1, dash="dot"), annotation_text="enrichment = 1 (chance)", annotation_font=dict(size=9, color="#888"), annotation_position="bottom right")
    fig.update_layout(showlegend=True, legend=dict(orientation="h", y=-0.18))
    if caption:
        fig.add_annotation(xref="paper", yref="paper", x=0, y=-0.3, showarrow=False, align="left", font=dict(size=10, color="#555"),
                           text="Each point = one candidate target screened against the core gene (local walks). Ranking table: screen.csv (composite = mean percentile of "
                                "interference, enrichment, dependency; baseline_z vs degree-matched random genes). Direction is a property, not a magnitude.")
        fig.update_layout(margin=dict(b=150))
    return fig


# --------------------------------------------------------------------------- cross-talk (three lanes)
def crosstalk_figure(json_path, width=1100, height=820, max_nodes=None, max_edges=None, focus=None, caption=True):
    """Draws the top `max_edges` interference edges (ranked by proximity to both genes) and the nodes around them."""
    X = _load(json_path); fset, fvals = parse_focus(focus); fset_or_none = fset if focus is not None else None
    if max_nodes is None: max_nodes = _view_cfg("max_nodes_crosstalk", 120)
    if max_edges is None: max_edges = _view_cfg("max_edges_crosstalk", 40)   # (crosstalk.max_edges_plot in older yamls is ignored)
    n_all_edges = len(X["interference"])
    X["interference"] = sorted(X["interference"], key=lambda e: (e.get("rank", 0)))[:max_edges]
    keep = {e["from_node"] for e in X["interference"]} | {e["to_node"] for e in X["interference"]}      # the drawn interference edges need their ends
    nodes, edges, capped = _cap_nodes(X["nodes"], X["edges"], max_nodes, keep_ids=keep, per_hop=True); X["edges"] = edges
    byid = {n["id"]: n for n in nodes}
    lane_x = dict(causal=18, shared=50, goi=82)
    max_hop = max([max(n.get("hops_c", 0) or 0, n.get("hops_x", 0) or 0) for n in nodes] or [1])
    groups = defaultdict(list)
    for n in nodes: groups[(n["lane"], n.get("hops_x", 0) if n["lane"] == "goi" else n.get("hops_c", 0))].append(n)
    pos = {}
    NCOLS = 5                                                  # names per grid row inside a lane (label collisions start above this)
    rows_at = defaultdict(int)                                 # rows needed at each hop = the tallest lane at that hop
    for (lane, hp), arr in groups.items():
        ncols = 1 if len(arr) == 1 else min(NCOLS, len(arr))
        rows_at[hp or 0] = max(rows_at[hp or 0], math.ceil(len(arr) / ncols))
    total_rows = sum(rows_at[h] for h in range(max_hop + 1)) + max_hop      # one blank row between hops
    row_h = min(4.2, 88 / max(total_rows, 1))
    y_top = {}; y = 94
    for h in range(max_hop + 1):
        y_top[h] = y; y -= (rows_at[h] + 1) * row_h
    for (lane, hp), arr in groups.items():
        arr.sort(key=lambda n: -n["score"])
        span = 26 if lane != "shared" else 24
        ncols = 1 if len(arr) == 1 else min(NCOLS, len(arr))
        for i, n in enumerate(arr):
            r, c = divmod(i, ncols)
            x = lane_x[lane] if ncols == 1 else lane_x[lane] - span / 2 + span * (c + 0.5) / ncols
            pos[n["id"]] = (x, y_top[hp or 0] - r * row_h)
    y_bottom = min([v[1] for v in pos.values()] or [0]) - 4   # crop the empty space under the last row
    height = max(560, int((102 - y_bottom) * 7.5) + 260)       # ~28 px per grid row + title / caption; grows with the rows
    rad = lambda nid: 2.0 if byid[nid]["cls"] == "seed" else 0.9 + 1.2 * byid[nid]["score"]
    s = X["summary"]
    enr = f"TF-overlap enrichment {s.get('tf_enrichment')} (p={s.get('tf_overlap_p', float('nan')):.1g})" if s.get("tf_enrichment") is not None else f"index {s['crosstalk_index']}"
    fig = _base(f"Cross-talk: does {X['goi']} influence the {X['causal']} pathway?  {enr} · showing top {len(X['interference'])} of {n_all_edges} interference edges"
                + (f" · {max_nodes} of {len(X['nodes'])} nodes" if capped else "") + f" · endpoints syn/ant/mixed {s['synergistic']}/{s['antagonistic']}/{s['mixed']}", width, height)
    contacts = [c for c in X.get("contacts", []) if c["goi_side"] in pos and c["core_side"] in pos]
    if contacts:   # physical contact without a signed edge (STRING / Reactome): dotted straight line, no direction
        xs, ys = [], []
        for c in contacts:
            xs += [pos[c["goi_side"]][0], pos[c["core_side"]][0], None]; ys += [pos[c["goi_side"]][1], pos[c["core_side"]][1], None]
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color="#7a9cc6", width=1.6, dash="dot"), hoverinfo="skip",
                                 name=f"physical contact {X['goi']}-side – {X['causal']} pathway (no signed edge; STRING / Reactome)"))
        mx = [(pos[c["goi_side"]][0] + pos[c["core_side"]][0]) / 2 for c in contacts]; my = [(pos[c["goi_side"]][1] + pos[c["core_side"]][1]) / 2 for c in contacts]
        fig.add_trace(go.Scatter(x=mx, y=my, mode="markers", marker=dict(size=5, color="#7a9cc6"), hoverinfo="text", showlegend=False,
                                 hovertext=[f"physical contact {c['goi_side']} – {c['core_side']} ({c['source']}); no direction / sign in the databases" for c in contacts]))
    _edge_traces(fig, X["edges"], pos, rad, color_fn=(lambda e: "#bbbbbb") if focus is None else (lambda e: "#999999" if (e["source"] in fset or e["target"] in fset) else "#d4d4d4"),
                 dash_fn=lambda e: "dot" if e.get("effect", 1) == 0 else "solid", width=0.8)
    inter = [dict(source=e["from_node"], target=e["to_node"], effect=e["edge_effect"], agreement=e["agreement"], mech=e.get("mechanism", ""), res=e.get("residue", "")) for e in X["interference"]]
    _edge_traces(fig, inter, pos, rad, color_fn=lambda e: COL["syn"] if e["agreement"] == "synergistic" else COL["ant"] if e["agreement"] == "antagonistic" else COL["mix"],
                 dash_fn=lambda e: "longdash" if "chaperone" in (e.get("mech") or "") else "solid", width=1.4)
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
                    lambda n: "#111" if n.get("is_endpoint") else "rgba(0,0,0,0)", hover, name, line_w=1.5, focus=fset_or_none, values=fvals, label_size=9)
        fig.add_annotation(x=lane_x[lane], y=99, text=name.upper(), showarrow=False, font=dict(size=11, color="#888"))
    _node_trace(fig, [n for n in nodes if n["cls"] == "seed"], pos, rad, lambda n: COL["seed"], lambda n: "#111", hover, "seed genes", focus=fset_or_none, values=fvals)
    if focus is not None:
        missing = sorted(fset - set(byid))
        fig.update_layout(title_text=fig.layout.title.text + f"  |  focus: {len(fset & set(byid))}/{len(fset)} among drawn nodes"
                          + (f" (not drawn / not reached: {', '.join(missing[:8])}{'…' if len(missing) > 8 else ''})" if missing else ""))
    fig.update_xaxes(range=[0, 100]); fig.update_yaxes(range=[y_bottom, 102])
    _edge_legend(fig, ("act", "inh", "syn", "ant", "chap"))
    if caption: add_caption(fig, "crosstalk")
    return fig


# --------------------------------------------------------------------------- synthetic lethality (volcano)
def volcano_figure(json_path, width=850, height=600, label_top=20, caption=True):
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
    if caption: add_caption(fig, "volcano")
    return fig

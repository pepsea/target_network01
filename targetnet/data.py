"""
Data loaders for target_net. Every loader has a live route (public APIs / downloads, cached on disk)
and a demo route (curated AKT1 axis, offline). Licenses noted per source.

  SIGNOR        CC-BY 4.0   directed, signed, residue-level signaling edges
  Reactome FI   CC0         directed reaction-derived edges (+ undirected context)
  STRING        CC-BY 4.0   physical edges (experimental + database channels only)
  UniProt       CC-BY 4.0   receptors (KW-0675 + cell membrane), subcellular location
  Lambert 2018  CC-BY       human transcription factor list
  DoRothEA      via OmniPath (regulon A-C recommended; check source licenses before commercial redistribution)
  DepMap        CC-BY 4.0   CRISPR gene effect, mutations, expression
  HPA           CC BY-SA 3.0 tissue expression (share-alike: used as a filter only, not redistributed)
"""
import io, json, os, time, zipfile
import numpy as np
import pandas as pd
import requests

from . import config as _cfgmod
CFG = _cfgmod.load_config()          # defaults; CLI/API call set_config() to override
CACHE = CFG["cache_dir"]

def set_config(cfg):
    global CFG, CACHE
    CFG = cfg; CACHE = cfg["cache_dir"]

def _local(name, key="path"):
    """Local file for a database if configured and present; otherwise None (download if allowed, else error)."""
    p = _cfgmod.db(CFG, name, key)
    if p and os.path.exists(p): return p
    if p:   # a path was configured but is wrong: stop, never silently download something else
        raise SystemExit(f"[config] {name}: configured path not found: {p}\n        (config file: {CFG.get('_file')}; fix databases.{name}.{key})")
    if not CFG.get("allow_download", True):
        raise SystemExit(f"[config] {name}: no local file and allow_download is false (set databases.{name}.{key} in the YAML)")
    print(f"[config] {name}: no local path configured -> downloading into {CACHE} (config file: {CFG.get('_file') or 'none found'})")
    return None
HEADERS = {"User-Agent": "target-net/0.2"}
REACTOME_FI_URL = "https://reactome.org/download/tools/ReatomeFIs/FIsInGene_070323_with_annotations.txt.zip"
LAMBERT_URL = "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"
OMNIPATH_URL = "https://omnipathdb.org/interactions"
HPA_URL = "https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip"

EDGE_WEIGHT = {"SIGNOR": 0.9, "Reactome": 0.7, "SIGNOR+Reactome": 0.95}


class DownloadError(RuntimeError):
    """A network fetch failed; required layers stop the run, optional layers (STRING) are skipped."""


# ------------------------------------------------------------------ cache helpers
def _cache_path(key, ext):
    os.makedirs(CACHE, exist_ok=True)
    return os.path.join(CACHE, key.replace("/", "_") + ext)


def cached_text(key, url, params=None, sleep=0.25, binary=False):
    fn = _cache_path(key, ".bin" if binary else ".txt")
    if os.path.exists(fn):
        return open(fn, "rb" if binary else "r").read()
    try:
        r = requests.get(url, params=params, headers=HEADERS, timeout=600); r.raise_for_status()
    except requests.RequestException as e:
        raise DownloadError(f"[download] could not fetch {url.split('?')[0]} ({e.__class__.__name__}: {e}). "
                            f"Set the local path in target_net.yaml (config file used: {CFG.get('_file') or 'none found'}) or check the network.")
    time.sleep(sleep)
    data = r.content if binary else r.text
    open(fn, "wb" if binary else "w").write(data)
    return data


def cached_json(key, url, params=None, sleep=0.25):
    fn = _cache_path(key, ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    r = requests.get(url, params=params, headers=HEADERS, timeout=120)
    r.raise_for_status(); time.sleep(sleep)
    js = r.json(); json.dump(js, open(fn, "w")); return js


# ------------------------------------------------------------------ signaling edges
def load_signor(demo=False):
    """
    DataFrame: src, dst, effect(+1/-1/0), mechanism, residue, source, score.
    Accepts either the per-organism API dump (getData.php?organism=9606) or SIGNOR's full download
    'all_data_DD_MM_YY.tsv' (all organisms; filtered here to human, TAX_ID 9606).
    Only protein -> protein rows are kept; complexes / protein families / chemicals are dropped (see note in README).
    """
    if demo:
        return _demo_signor()
    loc = _local("signor")
    txt = open(loc, encoding="utf-8", errors="replace").read() if loc else cached_text("signor_9606", _cfgmod.db(CFG, "signor", "url"))
    print(f"[SIGNOR] source: {loc or 'download'}")
    df = pd.read_csv(io.StringIO(txt), sep="\t", dtype=str, quoting=3, on_bad_lines="skip")
    # normalise header variants: 'ENTITYA' / 'ENTITY A' / 'entitya' / 'TAX_ID' / 'TAXID'
    norm = {c: c.upper().replace(" ", "").replace("_", "") for c in df.columns}
    df = df.rename(columns=norm)
    need = ["ENTITYA", "TYPEA", "ENTITYB", "TYPEB", "EFFECT"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise SystemExit(f"[SIGNOR] unexpected header, missing {missing}; got {list(df.columns)[:12]}")
    n_all = len(df)
    if "TAXID" in df.columns:                       # all_data file: keep human rows only
        df = df[df.TAXID.astype(str).str.strip() == "9606"]
    types = df.TYPEA.str.lower().str.strip(), df.TYPEB.str.lower().str.strip()
    n_human = len(df)
    df = df[(types[0] == "protein") & (types[1] == "protein")]
    eff = df.EFFECT.str.lower().map(lambda e: 1 if "up-regulates" in e else -1 if "down-regulates" in e else 0)
    score = pd.to_numeric(df["SCORE"], errors="coerce") if "SCORE" in df.columns else pd.Series(float("nan"), index=df.index)
    out = pd.DataFrame({"src": df.ENTITYA.str.upper().str.strip(), "dst": df.ENTITYB.str.upper().str.strip(), "effect": eff,
                        "mechanism": df.get("MECHANISM", pd.Series("", index=df.index)).fillna(""),
                        "residue": df.get("RESIDUE", pd.Series("", index=df.index)).fillna(""),
                        "source": "SIGNOR", "score": score})
    out = out[(out.src != out.dst) & (out.src != "") & (out.dst != "")]
    # same pair reported several times: keep the row with a sign (non-zero) and the highest SIGNOR score
    out = out.sort_values(["effect", "score"], key=lambda s: (s == 0) if s.name == "effect" else -s.fillna(0)).drop_duplicates(["src", "dst"])
    print(f"[SIGNOR] rows {n_all} -> human {n_human} -> protein-protein unique edges {len(out)}"
          + (f" (SIGNOR score median {out.score.median():.2f})" if out.score.notna().any() else ""))
    return out


def load_reactome_fi(path=None, demo=False):
    """Returns (directed DataFrame like load_signor, undirected DataFrame a,b)."""
    if demo:
        return _demo_reactome()
    path = path or _local("reactome_fi")
    if path is None:
        path = _cache_path("reactome_fi", ".zip")
        if not os.path.exists(path):
            print("[Reactome] downloading FI network ..."); open(path, "wb").write(
                requests.get(_cfgmod.db(CFG, "reactome_fi", "url"), headers=HEADERS, timeout=600).content)
    print(f"[Reactome] source: {path}")
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            name = [n for n in z.namelist() if n.endswith(".txt")][0]
            df = pd.read_csv(z.open(name), sep="\t", dtype=str)
    else:
        df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    mode = _cfgmod.db(CFG, "reactome_fi", "mode") or "directed"
    min_score = float(_cfgmod.db(CFG, "reactome_fi", "min_score") or 1.0)
    n0 = len(df)
    ann = df["Annotation"].fillna("").str.lower() if "Annotation" in df.columns else pd.Series("", index=df.index)
    sc = pd.to_numeric(df["Score"], errors="coerce").fillna(0) if "Score" in df.columns else pd.Series(1.0, index=df.index)
    df = df[(sc >= min_score) & ~ann.str.contains("predicted")]          # drop NBC-predicted FIs
    ann = ann.loc[df.index]
    # direction is trusted only when the annotation states a regulatory mechanism (not mere co-membership in a reaction/complex)
    mech = ann.str.contains("activat|inhibit|catalyz|catalys|phosphorylat|ubiquitinat|expression regulat|regulat", regex=True)
    rows, und = [], []
    for a, b, d, m in zip(df.Gene1.str.upper(), df.Gene2.str.upper(), df.Direction.fillna("-").str.strip(), mech):
        if a == b: continue
        if mode == "off": break
        if mode == "directed" and m and d != "-":
            if "->" in d or d == "<->": rows.append((a, b, 1))
            if "<-" in d or d == "<->": rows.append((b, a, 1))
            if "-|" in d: rows.append((a, b, -1))
            if "|-" in d: rows.append((b, a, -1))
        else:
            und.append((a, b))
    dirdf = pd.DataFrame(rows, columns=["src", "dst", "effect"]).assign(mechanism="reactome_FI", residue="", source="Reactome").drop_duplicates(["src", "dst"])
    unddf = pd.DataFrame(und, columns=["a", "b"]).assign(source="Reactome")
    print(f"[Reactome] FI rows {n0} -> after score/predicted filter {len(df)} -> directed (mechanistic) {len(dirdf)}, context {len(unddf)}  [mode={mode}]")
    return dirdf, unddf


def merge_directed(*dfs):
    """Merge directed edge tables; same pair from two DBs -> source 'A+B', weight from EDGE_WEIGHT."""
    allr = pd.concat(dfs, ignore_index=True)
    def agg(g):
        srcs = "+".join(sorted(set(g.source)))
        eff = g.effect[g.effect != 0]
        effect = int(eff.iloc[0]) if len(eff) else 0
        if len(set(eff)) > 1: effect = 0           # databases disagree -> ambiguous
        res = ";".join(sorted({r for r in g.residue if r}))
        w = EDGE_WEIGHT.get(srcs, 0.7)
        if srcs == "Reactome": w = float(_cfgmod.db(CFG, "reactome_fi", "weight") or 0.5)
        if "score" in g and g.score.notna().any():        # SIGNOR confidence score (0-1): scale within [0.6, 1.0] x base weight
            w = w * (0.6 + 0.4 * float(g.score.max()))
        return pd.Series(dict(effect=effect, mechanism=";".join(sorted(set(g.mechanism))), residue=res, source=srcs, weight=round(w, 3)))
    out = allr.groupby(["src", "dst"], sort=False).apply(agg, include_groups=False).reset_index()
    return out


def load_string(genes, min_score=700, demo=False, batch=500):
    """Physical edges among `genes` (experimental + database channels), via the STRING API (POST, batched).
    Optional layer: any failure returns an empty table with a warning instead of stopping the run."""
    if demo:
        return _demo_string()
    empty = pd.DataFrame(columns=["a", "b", "source"])
    min_score = _cfgmod.db(CFG, "string", "min_score") or min_score
    if not CFG.get("allow_download", True):
        print("[STRING] offline mode: physical layer skipped"); return empty
    genes = sorted(set(genes))
    key = "string_" + str(abs(hash((tuple(genes), min_score))))
    fn = _cache_path(key, ".tsv")
    if os.path.exists(fn):
        df = pd.read_csv(fn, sep="\t")
    else:
        frames = []
        try:
            for i in range(0, len(genes), batch):
                chunk = genes[i:i + batch]
                r = requests.post("https://string-db.org/api/tsv/network",
                                  data={"identifiers": "\r".join(chunk), "species": 9606, "required_score": min_score,
                                        "network_type": "physical", "caller_identity": HEADERS["User-Agent"]},
                                  headers=HEADERS, timeout=300)
                r.raise_for_status(); time.sleep(0.5)
                frames.append(pd.read_csv(io.StringIO(r.text), sep="\t"))
        except Exception as e:
            print(f"[STRING] skipped (physical-contact column will be empty): {e.__class__.__name__}: {str(e)[:120]}"); return empty
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        if len(df): df.to_csv(fn, sep="\t", index=False)
    if not len(df) or "escore" not in df.columns:
        print("[STRING] no edges returned"); return empty
    df = df[(df.escore > 0) | (df.dscore > 0)]
    out = pd.DataFrame({"a": df.preferredName_A.str.upper(), "b": df.preferredName_B.str.upper(), "source": "STRING"}).drop_duplicates()
    print(f"[STRING] {len(out)} physical edges among {len(genes)} genes")
    return out


# ------------------------------------------------------------------ bulk UniProt (receptors + compartments offline)
_UNIPROT = {}

def load_uniprot_bulk():
    """Human reviewed UniProt tsv with columns Entry, Gene Names (primary), Keywords, Subcellular location [CC]. Returns DataFrame or None."""
    if "df" in _UNIPROT: return _UNIPROT["df"]
    loc = _cfgmod.db(CFG, "uniprot")
    if loc and os.path.exists(loc):
        df = pd.read_csv(loc, sep="\t", dtype=str).fillna("")
        cols = {c.lower(): c for c in df.columns}
        gene = next((cols[c] for c in cols if c.startswith("gene names")), None)
        kw = next((cols[c] for c in cols if c.startswith("keywords")), None)
        sl = next((cols[c] for c in cols if c.startswith("subcellular location")), None)
        if not (gene and sl):
            print(f"[UniProt] bulk file lacks gene / subcellular-location columns: {list(df.columns)}"); _UNIPROT["df"] = None; return None
        out = pd.DataFrame({"gene": df[gene].str.split(";").str[0].str.strip().str.upper(), "keywords": df[kw] if kw else "", "location": df[sl]})
        out = out[out.gene != ""].drop_duplicates("gene")
        print(f"[UniProt] bulk table: {len(out)} genes from {loc}")
        _UNIPROT["df"] = out; return out
    _UNIPROT["df"] = None; return None


# ------------------------------------------------------------------ endpoint sets
def load_receptors(demo=False):
    """Gene symbols with UniProt keyword Receptor (KW-0675) and cell-membrane location."""
    if demo:
        return _demo_receptors()
    loc = _cfgmod.db(CFG, "receptors")
    if loc and os.path.exists(loc):
        genes = {l.strip().upper() for l in open(loc) if l.strip() and not l.startswith("#")}
        print(f"[receptors] {len(genes)} from {loc}"); return genes
    up = load_uniprot_bulk()
    if up is not None:
        m = up.keywords.str.contains(r"\bReceptor\b", case=False) & up.location.str.contains("Cell membrane|Plasma membrane", case=False)
        genes = set(up[m].gene); print(f"[receptors] {len(genes)} from bulk UniProt (keyword Receptor + cell membrane)"); return genes
    if not CFG.get("allow_download", True):
        raise SystemExit("[config] receptors: set databases.uniprot.path (bulk tsv) or databases.receptors.path, or allow_download")
    genes = set(); cursor_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": "keyword:KW-0675 AND organism_id:9606 AND reviewed:true AND cc_scl_term:SL-0039",
              "fields": "gene_primary", "format": "tsv", "size": 500}
    r = requests.get(cursor_url, params=params, headers=HEADERS, timeout=120)
    while True:
        r.raise_for_status()
        for line in r.text.splitlines()[1:]:
            for g in line.split(";"):
                if g.strip(): genes.add(g.strip().upper())
        nxt = r.links.get("next", {}).get("url")
        if not nxt: break
        r = requests.get(nxt, headers=HEADERS, timeout=120)
    print(f"[UniProt] {len(genes)} membrane receptors")
    return genes


def load_tfs(demo=False):
    if demo:
        return _demo_tfs()
    loc = _local("tf_list")
    txt = open(loc, encoding="utf-8", errors="replace").read() if loc else cached_text("lambert_tfs", _cfgmod.db(CFG, "tf_list", "url"))
    df = pd.read_csv(io.StringIO(txt))
    col_tf = [c for c in df.columns if c.strip().lower().startswith("is tf")][0]
    col_sym = [c for c in df.columns if "HGNC" in c][0]
    tfs = set(df[df[col_tf].astype(str).str.strip() == "Yes"][col_sym].str.upper())
    print(f"[Lambert] {len(tfs)} TFs")
    return tfs


def omnipath_regulon_url(source="dorothea", conf="ABC", license="commercial"):
    """URL that returns the regulon as tsv (same thing the automatic route downloads)."""
    q = {"datasets": source, "organisms": 9606, "genesymbols": 1, "fields": "sources,references,dorothea_level", "license": license}
    if source == "dorothea": q["dorothea_levels"] = ",".join(conf)
    return OMNIPATH_URL + "?" + "&".join(f"{k}={v}" for k, v in q.items())


def _regulon_from_table(df):
    """Accept OmniPath export (source_genesymbol, target_genesymbol, is_stimulation, is_inhibition, dorothea_level)
    or a simple table (tf, target, mor, confidence)."""
    df.columns = [c.lower() for c in df.columns]
    if "source_genesymbol" in df.columns:
        out = pd.DataFrame({"tf": df.source_genesymbol, "target": df.target_genesymbol,
                            "mor": np.where(pd.to_numeric(df.get("is_inhibition", 0), errors="coerce").fillna(0) == 1, -1, 1),
                            "confidence": df.get("dorothea_level", pd.Series("A", index=df.index)).astype(str).str[0]})
    elif {"tf", "target"} <= set(df.columns):
        out = df.rename(columns={"source": "tf"})[["tf", "target"]].copy()
        out["mor"] = pd.to_numeric(df.get("mor", 1), errors="coerce").fillna(1).astype(int)
        out["confidence"] = df.get("confidence", pd.Series("A", index=df.index)).astype(str).str[0]
    elif {"source", "target"} <= set(df.columns):   # decoupler / CollecTRI style
        out = pd.DataFrame({"tf": df.source, "target": df.target, "mor": np.sign(pd.to_numeric(df.get("weight", 1), errors="coerce").fillna(1)).astype(int),
                            "confidence": df.get("confidence", pd.Series("A", index=df.index)).astype(str).str[0]})
    else:
        raise SystemExit(f"[regulon] unrecognised columns: {list(df.columns)[:10]}")
    return out


def load_regulon(source=None, conf=None, path=None, demo=False):
    """DataFrame: tf, target, mor(+1/-1), confidence, weight."""
    if demo:
        return _demo_regulon()
    source = source or _cfgmod.db(CFG, "regulon", "source") or "dorothea"
    conf = conf or _cfgmod.db(CFG, "regulon", "conf") or "ABC"
    lic = _cfgmod.db(CFG, "regulon", "license") or "commercial"
    path = path or _local("regulon")
    if path:
        print(f"[regulon] source: {path}")
        df = pd.read_csv(path, sep=None, engine="python")
    else:
        txt = cached_text(f"omnipath_{source}_{conf}_{lic}", omnipath_regulon_url(source, conf, lic))
        print(f"[regulon] source: OmniPath ({source}, levels {conf}, license={lic})")
        df = pd.read_csv(io.StringIO(txt), sep="\t")
    df = _regulon_from_table(df)
    df["tf"] = df.tf.astype(str).str.upper(); df["target"] = df.target.astype(str).str.upper()
    df = df[df.confidence.isin(list(conf))].drop_duplicates(["tf", "target"])
    df["weight"] = df.confidence.map({"A": 1.0, "B": 0.8, "C": 0.6, "D": 0.4, "E": 0.2}).fillna(0.5)
    print(f"[{source}] {len(df)} TF-target edges (conf {conf}), {df.tf.nunique()} TFs")
    return df


def load_expressed(tissue=None, expressed_file=None, min_ntpm=1.0):
    """Set of expressed genes, or None (no filter)."""
    if expressed_file:
        return {l.strip().upper() for l in open(expressed_file) if l.strip()}
    if not tissue:
        return None
    loc = _local("hpa_tissue")
    if loc and not loc.endswith(".zip"):
        df = pd.read_csv(loc, sep="\t")
    else:
        raw = open(loc, "rb").read() if loc else cached_text("hpa_tissue_consensus", _cfgmod.db(CFG, "hpa_tissue", "url"), binary=True)
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            df = pd.read_csv(z.open(z.namelist()[0]), sep="\t")
    df.columns = [c.lower() for c in df.columns]
    sub = df[df.tissue.str.lower() == tissue.lower()]
    if sub.empty:
        print(f"[HPA] tissue '{tissue}' not found; available e.g.: {sorted(df.tissue.unique())[:12]} ..."); return None
    exp = set(sub[sub.ntpm >= min_ntpm]["gene name"].str.upper())
    print(f"[HPA] {len(exp)} genes expressed in {tissue} (nTPM >= {min_ntpm})")
    return exp


def load_compartments(genes, demo=False):
    if demo:
        return _demo_compartments()
    loc = _cfgmod.db(CFG, "compartments")
    if loc and os.path.exists(loc):
        tab = pd.read_csv(loc, sep="\t", header=None, names=["gene", "compartment"], comment="#")
        m = dict(zip(tab.gene.str.upper(), tab.compartment)); return {g: m.get(g, "cytosol") for g in genes}
    kw = [("plasma membrane", "membrane"), ("cell membrane", "membrane"), ("lysosome", "lysosome"),
          ("late endosome", "lysosome"), ("mitochondri", "mitochondrion"), ("endoplasmic reticulum", "ER_golgi"),
          ("golgi", "ER_golgi"), ("nucleus", "nucleus"), ("nucleolus", "nucleus"), ("cytoplasm", "cytosol"), ("cytosol", "cytosol")]
    def first_mentioned(txt):   # UniProt lists the primary location first -> take the earliest keyword hit
        hits = [(txt.find(k), c) for k, c in kw if k in txt]
        return min(hits)[1] if hits else "cytosol"
    up = load_uniprot_bulk()
    if up is not None:
        locmap = dict(zip(up.gene, up.location.str.lower()))
        return {g: first_mentioned(locmap.get(g, "")) for g in genes}
    if not CFG.get("allow_download", True):
        print("[UniProt] offline and no bulk file: all compartments set to cytosol"); return {g: "cytosol" for g in genes}
    out = {}
    for g in sorted(genes):
        try:
            js = cached_json("uniprot_loc_" + g, "https://rest.uniprot.org/uniprotkb/search",
                             {"query": f"gene_exact:{g} AND organism_id:9606 AND reviewed:true",
                              "fields": "cc_subcellular_location", "format": "json", "size": 1})
            txt = json.dumps(js).lower()
            out[g] = first_mentioned(txt)
        except Exception:
            out[g] = "cytosol"
    return out


# ------------------------------------------------------------------ DepMap
_DEPMAP = {}

def load_depmap(depmap_dir, gene, candidates, min_lost=None, padj_th=None, delta_th=None):
    """Returns coess {gene: r}, sl_full DataFrame(gene, delta, p, padj, n_lost, significant)."""
    empty = ({}, pd.DataFrame())
    depmap_dir = depmap_dir or _cfgmod.db(CFG, "depmap", "dir")
    slc = CFG.get("sl", {}); min_lost = min_lost or slc.get("min_lost", 20); padj_th = padj_th or slc.get("padj", 0.05); delta_th = delta_th or slc.get("delta", 0.15)
    if not depmap_dir: return empty
    ge_path = os.path.join(depmap_dir, "CRISPRGeneEffect.csv")
    if not os.path.exists(ge_path):
        print("[DepMap] CRISPRGeneEffect.csv not found"); return empty
    if _DEPMAP.get("path") == ge_path:
        ge = _DEPMAP["ge"]
    else:
        print("[DepMap] loading gene effect (once per session) ...")
        ge = pd.read_csv(ge_path, index_col=0); ge.columns = [c.split(" (")[0].upper() for c in ge.columns]
        _DEPMAP.update(path=ge_path, ge=ge)
    if gene not in ge.columns:
        print(f"[DepMap] {gene} not in matrix"); return empty
    a = ge[gene]
    cols = [g for g in candidates if g in ge.columns]
    corr = ge[cols].corrwith(a) if cols else pd.Series(dtype=float)
    coess = {g: float(r) for g, r in corr.items() if not np.isnan(r)}
    lost = set()
    mut_path = os.path.join(depmap_dir, "OmicsSomaticMutations.csv")
    if os.path.exists(mut_path):
        mut = pd.read_csv(mut_path, usecols=lambda c: c in ("ModelID", "HugoSymbol", "LikelyLoF"), low_memory=False)
        m = mut[mut.HugoSymbol.str.upper() == gene]
        if "LikelyLoF" in m: m = m[m.LikelyLoF.astype(str).str.lower().isin(["true", "1"])]
        lost = set(m.ModelID) & set(ge.index); print(f"[DepMap] {len(lost)} likely-LoF lines for {gene}")
    if len(lost) < min_lost:
        ex_path = os.path.join(depmap_dir, "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
        if os.path.exists(ex_path):
            ex = pd.read_csv(ex_path, index_col=0); ex.columns = [c.split(" (")[0].upper() for c in ex.columns]
            common = ex.index.intersection(ge.index); q = ex.loc[common, gene].quantile(0.25)
            lost = set(ex.loc[common][ex.loc[common, gene] <= q].index)
            print(f"[DepMap] fallback: lowest-quartile {gene} expression = 'lost' ({len(lost)} lines)")
    if len(lost) < min_lost:
        print("[DepMap] too few lost lines -> no SL test"); return coess, pd.DataFrame()
    from scipy import stats
    li = ge.index.isin(lost)
    X, Y = ge.loc[li], ge.loc[~li]                       # vectorised Welch t-test over all genes at once
    nx_, ny_ = X.notna().sum(), Y.notna().sum()
    ok = (nx_ >= 10) & (ny_ >= 10)
    mx, my = X.mean(), Y.mean(); vx, vy = X.var(ddof=1), Y.var(ddof=1)
    se = np.sqrt(vx / nx_ + vy / ny_)
    t = (mx - my) / se
    dfree = (vx / nx_ + vy / ny_) ** 2 / ((vx / nx_) ** 2 / (nx_ - 1) + (vy / ny_) ** 2 / (ny_ - 1))
    pvals = 2 * stats.t.sf(np.abs(t), dfree)
    full = pd.DataFrame({"gene": ge.columns, "delta": (mx - my).values, "p": pvals, "n_lost": nx_.values})[ok.values].reset_index(drop=True)
    full["padj"] = bh_fdr(full.p.values)
    full["significant"] = (full.padj < padj_th) & (full.delta < -delta_th)
    print(f"[DepMap] {int(full.significant.sum())} SL candidates")
    return coess, full


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    r = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(r, 0, 1); return out


# ================================================================== DEMO DATA (curated AKT1 axis; offline smoke test)
def _demo_signor():
    e = [  # upstream of AKT1 (RTK / GPCR / integrin inputs)
        ("EGFR", "PIK3CA", 1), ("ERBB2", "PIK3CA", 1), ("PDGFRB", "PIK3CA", 1), ("FGFR3", "FRS2", 1), ("FRS2", "GAB1", 1),
        ("GAB1", "PIK3CA", 1), ("INSR", "IRS1", 1), ("IGF1R", "IRS1", 1), ("IRS1", "PIK3CA", 1), ("CXCR4", "GNB1", 1),
        ("GNB1", "PIK3CG", 1), ("ITGB1", "PTK2", 1), ("PTK2", "PIK3CA", 1), ("PIK3CA", "PDPK1", 1), ("PIK3CG", "PDPK1", 1),
        ("PTEN", "PDPK1", -1), ("PDPK1", "AKT1", 1), ("MTOR", "AKT1", 1), ("PHLPP1", "AKT1", -1), ("PPP2CA", "AKT1", -1),
        ("TRIB3", "AKT1", -1), ("KRAS", "PIK3CA", 1), ("EGFR", "KRAS", 1),
        # downstream of AKT1
        ("AKT1", "FOXO1", -1), ("AKT1", "FOXO3", -1), ("AKT1", "GSK3B", -1), ("AKT1", "TSC2", -1), ("AKT1", "AKT1S1", -1),
        ("AKT1", "MDM2", 1), ("AKT1", "BAD", -1), ("AKT1", "CHUK", 1), ("AKT1", "CREB1", 1), ("AKT1", "NOS3", 1),
        ("AKT1", "CDKN1B", -1), ("AKT1", "YAP1", -1),
        ("GSK3B", "CTNNB1", -1), ("GSK3B", "MYC", -1), ("GSK3B", "SNAI1", -1), ("GSK3B", "NFE2L2", -1),
        ("CTNNB1", "TCF7L2", 1), ("MDM2", "TP53", -1), ("CHUK", "NFKBIA", -1), ("NFKBIA", "RELA", -1),
        ("TSC2", "RHEB", -1), ("RHEB", "MTOR", 1), ("AKT1S1", "MTOR", -1), ("MTOR", "RPS6KB1", 1), ("MTOR", "EIF4EBP1", -1),
        ("MTOR", "SREBF1", 1), ("MTOR", "TFEB", -1), ("MTOR", "ULK1", -1), ("RPS6KB1", "RPS6", 1), ("RPS6KB1", "IRS1", -1),
        ("EIF4EBP1", "EIF4E", -1), ("BAD", "BCL2L1", -1), ("YAP1", "TEAD1", 1),
        # MAPK1 (gene of interest in crosstalk demo) and its downstream, with cross-edges into the AKT1 axis
        ("MAP2K1", "MAPK1", 1), ("BRAF", "MAP2K1", 1), ("KRAS", "BRAF", 1), ("MAPK1", "RPS6KA1", 1), ("MAPK1", "ELK1", 1),
        ("MAPK1", "ETV4", 1), ("MAPK1", "FOS", 1), ("MAPK1", "MYC", 1), ("MAPK1", "TSC2", -1), ("MAPK1", "FOXO3", -1),
        ("MAPK1", "MKNK1", 1), ("MKNK1", "EIF4E", 1), ("RPS6KA1", "GSK3B", -1), ("RPS6KA1", "BAD", -1), ("RPS6KA1", "CREB1", 1),
        ("RPS6KA1", "RPS6", 1), ("MAPK1", "SOS1", -1), ("SOS1", "KRAS", 1), ("MAPK1", "FRS2", -1), ("AKT1", "RAF1", -1),
        ("RAF1", "MAP2K1", 1),
        # FGFR3 axis (demo): FRS2/GRB2/SOS-RAS-RAF-MEK-ERK, PLCG1-PKC, PI3K, STATs
        ("FGFR3", "PLCG1", 1), ("PLCG1", "PRKCA", 1), ("PRKCA", "RAF1", 1), ("FRS2", "GRB2", 1), ("GRB2", "SOS1", 1),
        ("KRAS", "RAF1", 1), ("FGFR3", "STAT1", 1), ("FGFR3", "STAT3", 1), ("MAPK1", "STAT3", 1), ("SPRY2", "GRB2", -1),
        ("MAPK1", "SPRY2", 1), ("FGFR3", "SHC1", 1), ("SHC1", "GRB2", 1),
    ]
    # HSP90 chaperone-client edges (demo; stabilization = positive regulation of client protein level)
    chap = [("HSP90AA1", c) for c in ["FGFR3", "RAF1", "BRAF", "AKT1", "CDK4", "ERBB2", "EGFR", "HIF1A", "SRC", "IKBKB", "MET", "PDPK1"]]
    chap += [("CDC37", c) for c in ["RAF1", "BRAF", "AKT1", "CDK4", "SRC", "FGFR3"]] + [("CDC37", "HSP90AA1")]
    chap += [("STUB1", "FGFR3"), ("STUB1", "RAF1"), ("STUB1", "AKT1")]      # CHIP E3: degrades un-chaperoned clients
    chap = [(a, b, -1 if a == "STUB1" else 1) for a, b in chap]
    res = {("AKT1", "GSK3B"): "Ser9", ("AKT1", "FOXO1"): "Thr24;Ser256", ("AKT1", "FOXO3"): "Thr32;Ser253", ("AKT1", "TSC2"): "Ser939;Thr1462",
           ("AKT1", "BAD"): "Ser136", ("AKT1", "MDM2"): "Ser166;Ser186", ("RPS6KA1", "GSK3B"): "Ser9", ("RPS6KA1", "BAD"): "Ser112",
           ("MAPK1", "TSC2"): "Ser664", ("MAPK1", "FOXO3"): "Ser294;Ser344;Ser425", ("AKT1", "RAF1"): "Ser259", ("MAPK1", "SOS1"): "Ser1132",
           ("MAPK1", "FRS2"): "Thr-sites", ("PDPK1", "AKT1"): "Thr308", ("MTOR", "AKT1"): "Ser473", ("RPS6KA1", "CREB1"): "Ser133",
           ("AKT1", "CREB1"): "Ser133", ("MAPK1", "MYC"): "Ser62", ("GSK3B", "MYC"): "Thr58"}
    df = pd.DataFrame(e, columns=["src", "dst", "effect"])
    df["mechanism"] = "phosphorylation"; df["residue"] = [res.get((a, b), "") for a, b in zip(df.src, df.dst)]
    cdf = pd.DataFrame(chap, columns=["src", "dst", "effect"]); cdf["mechanism"] = "chaperone/stability"; cdf["residue"] = ""
    df = pd.concat([df, cdf], ignore_index=True); df["source"] = "SIGNOR"; return df

def _demo_reactome():
    d = [("INSR", "IRS1", 1), ("IRS1", "PIK3CA", 1), ("PIK3CA", "PDPK1", 1), ("PDPK1", "AKT1", 1), ("AKT1", "FOXO1", -1),
         ("AKT1", "GSK3B", -1), ("RHEB", "MTOR", 1), ("MTOR", "RPS6KB1", 1), ("MAPK1", "ELK1", 1), ("MAP2K1", "MAPK1", 1),
         ("AKT1", "TSC2", -1), ("GSK3B", "CTNNB1", -1), ("MDM2", "TP53", -1), ("ERBB3", "PIK3CA", 1), ("MET", "GAB1", 1)]
    u = [("AKT1", "PDPK1"), ("MTOR", "RICTOR"), ("TSC1", "TSC2"), ("PIK3CA", "PIK3R1")]
    return (pd.DataFrame(d, columns=["src", "dst", "effect"]).assign(mechanism="reactome_FI", residue="", source="Reactome"),
            pd.DataFrame(u, columns=["a", "b"]).assign(source="Reactome"))

def _demo_string():
    e = [("AKT1", "PDPK1"), ("AKT1", "MTOR"), ("AKT1", "GSK3B"), ("AKT1", "FOXO3"), ("AKT1", "TSC2"), ("PIK3CA", "PIK3R1"),
         ("MTOR", "RPTOR"), ("MTOR", "RICTOR"), ("TSC1", "TSC2"), ("AKT1", "PHLPP1"), ("IRS1", "PIK3CA"), ("MAPK1", "MAP2K1"), ("MAPK1", "RPS6KA1")]
    return pd.DataFrame(e, columns=["a", "b"]).assign(source="STRING")

def _demo_receptors():
    return {"EGFR", "ERBB2", "ERBB3", "PDGFRB", "FGFR3", "INSR", "IGF1R", "CXCR4", "ITGB1", "MET", "KDR", "FLT3", "KIT"}

def _demo_tfs():
    return {"FOXO1", "FOXO3", "MYC", "TP53", "RELA", "CREB1", "SREBF1", "TFEB", "SNAI1", "NFE2L2", "TCF7L2", "TEAD1",
            "ELK1", "ETV4", "FOS", "JUN", "EGR1", "STAT3", "STAT1", "NFKB1", "HIF1A"}

def _demo_regulon():
    reg = {"FOXO1": [("CDKN1B", 1), ("SOD2", 1), ("G6PC1", 1), ("PCK1", 1), ("BCL2L11", 1), ("GADD45A", 1)],
           "FOXO3": [("CDKN1B", 1), ("BCL2L11", 1), ("SOD2", 1), ("GADD45A", 1), ("PUMA", 1)],
           "MYC": [("CDK4", 1), ("LDHA", 1), ("ODC1", 1), ("CDKN1A", -1), ("CCND2", 1), ("NCL", 1)],
           "RELA": [("NFKBIA", 1), ("TNF", 1), ("IL6", 1), ("BCL2L1", 1), ("CXCL8", 1)],
           "TP53": [("CDKN1A", 1), ("MDM2", 1), ("BAX", 1), ("BBC3", 1), ("GADD45A", 1)],
           "CREB1": [("FOS", 1), ("NR4A1", 1), ("PPARGC1A", 1), ("BDNF", 1)],
           "SREBF1": [("FASN", 1), ("SCD", 1), ("HMGCR", 1), ("LDLR", 1)],
           "TFEB": [("CTSD", 1), ("LAMP1", 1), ("ATP6V1H", 1), ("SQSTM1", 1)],
           "SNAI1": [("CDH1", -1), ("VIM", 1)], "NFE2L2": [("NQO1", 1), ("HMOX1", 1), ("GCLC", 1)],
           "TCF7L2": [("MYC", 1), ("CCND1", 1), ("AXIN2", 1)], "TEAD1": [("CTGF", 1), ("CYR61", 1)],
           "ELK1": [("FOS", 1), ("EGR1", 1)], "ETV4": [("MMP1", 1), ("MMP9", 1)], "FOS": [("CCND1", 1), ("MMP1", 1)],
           "STAT3": [("BCL2L1", 1), ("MYC", 1), ("CCND1", 1), ("SOCS3", 1), ("MCL1", 1)], "STAT1": [("IRF1", 1), ("CDKN1A", 1), ("CXCL10", 1)],
           "HIF1A": [("VEGFA", 1), ("SLC2A1", 1), ("LDHA", 1), ("CA9", 1)]}
    rows = [(tf, t, m, "A" if i < 3 else "B") for tf, ts in reg.items() for i, (t, m) in enumerate(ts)]
    df = pd.DataFrame(rows, columns=["tf", "target", "mor", "confidence"])
    df["weight"] = df.confidence.map({"A": 1.0, "B": 0.8}); return df

def _demo_compartments():
    mem = ["EGFR", "ERBB2", "ERBB3", "PDGFRB", "FGFR3", "INSR", "IGF1R", "CXCR4", "ITGB1", "MET", "IRS1", "FRS2", "GAB1", "PIK3CA",
           "PIK3CG", "PIK3R1", "PTEN", "PDPK1", "GNB1", "PTK2", "KRAS", "SOS1", "NOS3"]
    lys = ["TSC1", "TSC2", "RHEB", "MTOR", "RPTOR", "RICTOR", "AKT1S1"]
    nuc = ["FOXO1", "FOXO3", "MYC", "TP53", "RELA", "CREB1", "TFEB", "SNAI1", "NFE2L2", "TCF7L2", "TEAD1", "ELK1", "ETV4", "FOS",
           "JUN", "EGR1", "MDM2", "CTNNB1", "YAP1", "CDKN1B"]
    er = ["SREBF1"]; mito = ["BAD", "BCL2L1"]
    mem += ["PLCG1", "GRB2", "SHC1", "SRC"]; nuc += ["STAT1", "STAT3", "HIF1A", "CDK4"]
    c = {g: "membrane" for g in mem}; c.update({g: "lysosome" for g in lys}); c.update({g: "nucleus" for g in nuc})
    c.update({g: "ER_golgi" for g in er}); c.update({g: "mitochondrion" for g in mito}); return c

def demo_depmap(gene, cand):
    """SIMULATED genetic layer for offline smoke tests (not real DepMap values)."""
    rng = np.random.default_rng(1); n = 1500
    genes = [f"GENE{i:04d}" for i in range(n)]
    delta = rng.normal(0, 0.06, n); p = np.exp(-np.abs(delta) / 0.06 * 2.5) * rng.uniform(0.3, 1, n)
    hits = {"AKT2": -0.45, "PIK3CB": -0.33, "MCL1": -0.30, "SGK1": -0.26, "CDK4": -0.22, "BCL2L1": -0.20}
    for i, (g, d) in enumerate(hits.items()):
        genes[i] = g; delta[i] = d; p[i] = 10 ** rng.uniform(-9, -5)
    full = pd.DataFrame(dict(gene=genes, delta=delta, p=p, n_lost=38)); full["padj"] = bh_fdr(full.p.values)
    full["significant"] = (full.padj < 0.05) & (full.delta < -0.15)
    coess = {g: round(float(rng.uniform(0.15, 0.6)), 3) for g in cand if g != gene}
    coess.update({g: round(float(rng.uniform(-0.1, 0.1)), 3) for g in hits})
    return coess, full

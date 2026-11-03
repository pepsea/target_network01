"""
YAML configuration: where each database lives and whether downloads are allowed.

Search order for the config file:
  1. explicit path (CLI --config / api config=)
  2. $TARGETNET_CONFIG
  3. ./target_net.yaml, or target_net.yaml in any parent folder (notebooks live in notebooks/)
  4. target_net.yaml in the project root (next to the targetnet package), wherever the notebook runs
  5. ~/.target_net.yaml
  6. built-in defaults (everything downloaded on demand into cache_dir)

See config.example.yaml for all keys. Every 'path' may use ~ and environment variables.
"""
import os, copy
import yaml

DEFAULTS = {
    "cache_dir": "cache",
    "allow_download": True,
    "databases": {
        "signor":      {"path": None, "url": "https://signor.uniroma2.it/getData.php?organism=9606",
                        "indirect": "drop",      # SIGNOR DIRECT column: drop (default) | weight (keep at low weight) | keep
                        "source": "file"},       # file (path or SIGNOR API) | omnipath (SIGNOR mirror via OmniPath; no SCORE/DIRECT columns)
        "reactome_fi": {"path": None, "url": "https://reactome.org/download/tools/ReatomeFIs/FIsInGene_070323_with_annotations.txt.zip",
                        "mode": "context",       # context (default) = undirected context layer only; directed = mechanistic FIs join the walk; off
                        "min_score": 1.0,        # FI score: 1.0 = pathway-annotated; <1 = NBC-predicted (excluded by default)
                        "weight": 0.5},          # edge weight for Reactome-only directed edges (SIGNOR = 0.9)
        "tf_list":     {"path": None, "url": "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"},
        "uniprot":     {"path": None, "url": "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&fields=accession,gene_primary,keyword,cc_subcellular_location&format=tsv"},
                                                             # bulk human reviewed UniProt tsv -> receptors + compartments offline
        "receptors":   {"path": None},                       # txt, one gene symbol per line (overrides uniprot bulk; else UniProt API)
        "ligands":     {"path": None},                       # txt, one gene symbol per line (overrides the UniProt-derived secreted-ligand list)
        "regulon":     {"source": "dorothea", "conf": "ABC", "path": None, "license": "commercial"},
                                                             # local tsv (OmniPath export or tf,target,mor,confidence); else OmniPath query
                                                             # license: commercial | academic | ignore  (OmniPath per-record licence filter)
        "string":      {"enabled": True, "min_score": 700},  # STRING API only; enabled: false = ignore STRING completely
        "compartments":{"path": None},                       # tsv: gene<TAB>compartment (else UniProt API)
        "depmap":      {"dir": None},                        # folder with CRISPRGeneEffect.csv etc.
        "custom_edges": {"path": None},
        "omnipath": {"mode": "full", "datasets": "omnipath,pathwayextra,kinaseextra,ligrecextra", "license": "commercial", "path": None, "signed_only": True, "min_references": 1},
        "hpa_tissue":  {"path": None, "url": "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"},
    },
    "walk": {"max_depth": 5, "min_score": 0.35, "hub_penalty": 0.3, "seed_rescue": True, "seed_rescue_weight": 0.3,
             "through_tf": True, "score_mode": "relative"},
    "view": {"max_nodes": 150, "max_nodes_crosstalk": 120, "max_edges_crosstalk": 40, "rows_per_column": 28},   # figure budgets only; tables are never capped
    "crosstalk": {"max_depth": 3, "min_score": 0.4, "min_sign_confidence": 0.5},   # local walks: cross-talk is proximal interference (plot budget -> view.max_edges_crosstalk)
    "interference": {"min_rel_change": 0.1, "weak_ratio": 0.5},    # a CORE node counts as changed when |delta / baseline_flow| >= this (10 %); gates flow_change_index
    "edge_score": {"mode": "percentile", "combine": "noisy_or"},   # edge weight = percentile-normalised score within each database; legacy = old fixed base weights
    "sl":   {"min_lost": 20, "padj": 0.05, "delta": 0.15},
}


def _expand(v):
    return os.path.expanduser(os.path.expandvars(v)) if isinstance(v, str) else v


def _merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def _upward(name, start=None):
    """target_net.yaml in the current folder or any parent (so notebooks in notebooks/ find the project-root file)."""
    d = os.path.abspath(start or os.getcwd())
    while True:
        c = os.path.join(d, name)
        if os.path.exists(c): return c
        nd = os.path.dirname(d)
        if nd == d: return None
        d = nd


def load_config(path=None):
    root_yaml = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "target_net.yaml")   # next to the package: found even when the notebook kernel runs outside the project
    cands = [path, os.environ.get("TARGETNET_CONFIG"), _upward("target_net.yaml"), root_yaml, os.path.expanduser("~/.target_net.yaml")]
    used = next((c for c in cands if c and os.path.exists(os.path.expanduser(c))), None)
    cfg = _merge(DEFAULTS, yaml.safe_load(open(os.path.expanduser(used))) or {}) if used else copy.deepcopy(DEFAULTS)
    for db in cfg["databases"].values():
        for k in ("path", "dir"):
            if k in db: db[k] = _expand(db[k])
    cfg["cache_dir"] = _expand(cfg["cache_dir"]); cfg["_file"] = used
    if path and not used:
        raise SystemExit(f"config file not found: {path}")
    return cfg


def db(cfg, name, key="path"):
    return cfg["databases"].get(name, {}).get(key)

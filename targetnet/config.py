"""
YAML configuration: where each database lives and whether downloads are allowed.

Search order for the config file:
  1. explicit path (CLI --config / api config=)
  2. $TARGETNET_CONFIG
  3. ./target_net.yaml
  4. ~/.target_net.yaml
  5. built-in defaults (everything downloaded on demand into cache_dir)

See config.example.yaml for all keys. Every 'path' may use ~ and environment variables.
"""
import os, copy
import yaml

DEFAULTS = {
    "cache_dir": "cache",
    "allow_download": True,
    "databases": {
        "signor":      {"path": None, "url": "https://signor.uniroma2.it/getData.php?organism=9606"},
        "reactome_fi": {"path": None, "url": "https://reactome.org/download/tools/ReatomeFIs/FIsInGene_070323_with_annotations.txt.zip"},
        "tf_list":     {"path": None, "url": "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"},
        "receptors":   {"path": None},                       # txt, one gene symbol per line (else UniProt API)
        "regulon":     {"source": "dorothea", "conf": "ABC", "path": None},   # tsv: tf,target,mor,confidence (else OmniPath)
        "string":      {"min_score": 700},                   # STRING API only
        "compartments":{"path": None},                       # tsv: gene<TAB>compartment (else UniProt API)
        "depmap":      {"dir": None},                        # folder with CRISPRGeneEffect.csv etc.
        "hpa_tissue":  {"path": None, "url": "https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip"},
    },
    "walk": {"max_depth": 6, "min_score": 0.15, "hub_penalty": 0.25},
    "sl":   {"min_lost": 20, "padj": 0.05, "delta": 0.15},
}


def _expand(v):
    return os.path.expanduser(os.path.expandvars(v)) if isinstance(v, str) else v


def _merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def load_config(path=None):
    cands = [path, os.environ.get("TARGETNET_CONFIG"), "target_net.yaml", os.path.expanduser("~/.target_net.yaml")]
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

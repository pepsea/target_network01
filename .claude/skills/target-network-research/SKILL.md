---
name: target-network-research
description: Gene-network target research in a Jupyter notebook with the target_net package. Single candidate gene -> upstream receptors, downstream TFs/target genes, synthetic lethality. Core (disease-causal) gene + target gene -> cross-talk analysis. Database locations come from a YAML config.
---

# Target network research (Jupyter notebook, Claude Code)

Use when the user asks to research a gene as a drug target, to find upstream alternatives / downstream effects /
synthetic-lethal partners of a gene, or whether one gene influences the pathway of another (cross-talk).
The analysis is done by the `target_net` package (this folder's parent contains `target_net.py`, `targetnet/`,
`notebooks/`, `config.example.yaml`). Never re-implement the analysis; drive the package from a notebook.

## 1. Decide the entry point from the inputs

| user gives | analysis | template | API |
|---|---|---|---|
| **one gene** (candidate / target) | upstream → receptors; downstream → TFs → target genes; synthetic lethality (DepMap) | `notebooks/gene_research.ipynb` | `research_gene(GENE, config=...)` |
| **two genes**: core (disease-causal) + target (gene of interest) | cross-talk: does TARGET influence the CORE pathway? | `notebooks/crosstalk_research.ipynb` | `research_crosstalk(core=CORE, target=TARGET, config=...)` |

If two genes are given and it is not clear which is the core gene, ask once. The analysis is asymmetric; do not guess.

## 2. Database configuration (YAML) — do this before any run

All database paths are read from a YAML file, searched in this order: `--config`/`config=` argument →
`$TARGETNET_CONFIG` → `./target_net.yaml` → `~/.target_net.yaml` → built-in defaults (download on demand).

0. Database files: `docs/DATA_SETUP.md` says where each file comes from; `bash scripts/setup_data.sh <root>` fetches the automatic ones.
1. Look for `target_net.yaml` in the project root or `~/.target_net.yaml`. If none exists, copy `config.example.yaml`
   to `target_net.yaml` and fill in the paths the user already has locally. Ask the user for locations of local copies
   only if the listing / `~/data` search finds nothing; do not download gigabytes without saying so.
2. Keys (all optional; `null` = download into `cache_dir` when `allow_download: true`):
   ```yaml
   cache_dir: ~/data/target_net_cache
   allow_download: true            # false = strict offline; every needed database must have a path
   databases:
     signor:       {path: ~/data/signor/signor_human.tsv}            # CC-BY 4.0
     reactome_fi:  {path: ~/data/reactome/FIsInGene_*.txt.zip}       # CC0
     tf_list:      {path: ~/data/tf/DatabaseExtract_v_1.01.csv}      # Lambert 2018, CC-BY
     receptors:    {path: null}                                       # txt, one symbol per line; null = UniProt API
     regulon:      {source: dorothea, conf: ABC, path: ~/data/dorothea/dorothea_hs_ABC.tsv}   # tsv tf,target,mor,confidence
     string:       {min_score: 700}                                   # API only (skipped when allow_download is false)
     compartments: {path: null}                                       # tsv gene<TAB>compartment; null = UniProt API
     depmap:       {dir: ~/data/depmap}                               # CRISPRGeneEffect.csv [+ mutations, expression]
     hpa_tissue:   {path: ~/data/hpa/rna_tissue_consensus.tsv.zip}    # CC BY-SA, filter only
   walk: {max_depth: 6, min_score: 0.15, hub_penalty: 0.25}
   sl:   {min_lost: 20, padj: 0.05, delta: 0.15}
   ```
3. Validate: `python target_net.py upstream --gene AKT1 --config target_net.yaml` prints `[config] <file>` and one
   `source:` line per database. A message like `signor: no local file and allow_download is false` names the exact
   key to fix. DepMap missing → co-essentiality and the `sl` mode are skipped; offer `python fetch_depmap.py --out <dir>` (~1.5 GB).
4. Record the config path you used in the notebook's parameters cell (`CONFIG = "../target_net.yaml"`).

## 3. Create the notebook

Copy the template to the user's working folder as `research_<GENE>.ipynb` or `crosstalk_<CORE>_<TARGET>.ipynb`.
Edit **only the parameters cell**: `GENE` (or `CORE`, `TARGET`), `CONFIG`, `DEMO = False`, `TISSUE` (HPA tissue name for
the disease context, else `None`), `TARGET_NET_DIR` (relative path to the folder containing `target_net.py`).
Leave `MAX_DEPTH / MIN_SCORE / HUB_PENALTY = None` so values come from the YAML `walk:` section unless the user asks to tune them.

Conventions (user's standing preferences):
- every code cell has a markdown cell above it explaining what it does and why; keep this when adding cells
- first version minimal, run it, look, then extend
- Python only; `DEMO = True` only for a smoke test and say so (demo DepMap values are simulated)

## 4. Execute and check

`jupyter nbconvert --to notebook --execute <nb>.ipynb --inplace --ExecutePreprocessor.timeout=1800`
(first live run with downloads can take minutes). Read the printed `summarize_gene` / `summarize_crosstalk` output and the tables.

Symbols: the package resolves `HSP90AA2 → HSP90AA1` (pseudogene), `ERK2 → MAPK1`, `MEK1 → MAP2K1`, `P53 → TP53`.
"not in the signaling graph" → check the HGNC symbol, pseudogene status, or whether the protein is non-signaling, before anything else.

Sanity checks before reporting novel hits (state the outcome in one line):
- textbook recovery: AKT1 upstream → RTKs via PI3K/PDPK1; AKT1 downstream → FOXO ↓, GSK3B ↓, mTORC1 ↑;
  MAPK1 × AKT1 cross-talk → RSK→GSK3B S9, ERK⊣SOS1 feedback. If missing, lower `min_score` or `hub_penalty` and re-run.
- too many TFs (hundreds) → set `TISSUE` or `expressed_file`
- cross-talk index needs a baseline: run the same TARGET against 3–5 unrelated core genes; hub-like targets (HSP90, broad kinases) score high against everything

## 5. Report

Answer first, then evidence. Single gene: top receptors (hops, sign, co-essentiality), TF programmes switched on/off,
synthetic-lethal hits with `on_pathway = False`. Cross-talk: index with baseline, interference edges (from → to, residue or
mechanism, synergistic/antagonistic), antagonistic edges called out (combination-therapy relevance). Always state the databases
and config used, whether DepMap / tissue filter were on, that database edges are context-free hypotheses, and that
interference edges with a residue are testable by phosphoproteomics (core inhibitor ± target inhibitor).
Point to the run folder (`r["run_dir"]` / `x["run_dir"]`, printed as `[run] …`) and its views (`flow.html`, `cell_map.html`, `crosstalk.html`, `synthetic_lethality.html`).

## Output layout
```
runs/<YYYYMMDD_HHMMSS>_gene_<GENE>/            run_info.json + upstream/ downstream/ sl/   (notebook / api.research_gene)
runs/<YYYYMMDD_HHMMSS>_<mode>_<GENE>/          CLI single mode
runs/<YYYYMMDD_HHMMSS>_crosstalk_<CORE>_x_<TARGET>/  crosstalk/ (summary.json shared.csv interference.csv tf_agreement.csv gene_agreement.csv crosstalk.html)
                                               + downstream_<CORE>_without_<TARGET>/ + downstream_<TARGET>_without_<CORE>/
each view (.html) has a .json twin for targetnet.plots; `--out DIR` / `out=` changes the root
```

## Licences (say so when results go to a client)
SIGNOR, STRING, UniProt, DepMap CC-BY 4.0 · Reactome CC0 · Lambert TF list CC-BY · HPA CC BY-SA (filter only, never redistribute) ·
DoRothEA aggregates sources — confirm the A–C subset before commercial redistribution · KEGG and PhosphoSitePlus are not used (commercial licences).

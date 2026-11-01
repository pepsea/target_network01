# target_net — network-based target discovery around a core gene (test case: AKT1)

Four modes, one signed signaling graph (SIGNOR ∪ Reactome FI), endpoints defined biologically
(receptors upstream, transcription factors downstream), DoRothEA for the TF → gene layer, DepMap for genetics.

```
python target_net.py upstream   --gene AKT1               # 1-1  core gene -> receptors
python target_net.py downstream --gene AKT1               # 1-2  core gene -> TFs -> target genes (DoRothEA A-C)
python target_net.py sl         --gene AKT1 --depmap DIR  # 1-3  synthetic lethality with core-gene loss
python target_net.py crosstalk  --causal AKT1 --goi MAPK1 # 2    does gene X influence the pathway of causal gene C?
```
Add `--demo` to any command for an offline smoke test on a curated AKT1 / MAPK1 axis (no internet needed;
the DepMap part of the demo is **simulated** and labelled as such).

## Install on a Linux server
`docs/INSTALL_LINUX.md` — venv or Docker (`Dockerfile`, `docker-compose.yml`), data folder, config, JupyterLab, systemd, troubleshooting.

## Data files
`docs/DATA_SETUP.md` lists every database (source page, command, expected file, licence, update cadence);
`bash scripts/setup_data.sh ~/Documents/database` downloads everything that can be fetched automatically.

## Configuration (YAML)
Copy `config.example.yaml` to `target_net.yaml` (project root) or `~/.target_net.yaml` and set the local paths of each
database; anything left `null` is downloaded into `cache_dir` when `allow_download: true`. `--config FILE` / `config=` overrides.
`allow_download: false` makes runs strictly offline and reports the exact YAML key that is missing.

## Claude Code skill
`.claude/skills/target-network-research/SKILL.md` — Claude Code picks it up automatically when this folder is the project
(or copy the folder to `~/.claude/skills/` for global use). It drives the notebooks from a single gene, or a core + target gene pair.

## Jupyter notebooks (recommended entry point)
```
notebooks/gene_research.ipynb       input: one candidate gene  -> upstream / downstream / synthetic lethality
notebooks/crosstalk_research.ipynb  input: CORE (disease-causal) + TARGET gene -> cross-talk
```
Both use `targetnet.api` (`research_gene`, `research_crosstalk`, `show`) and write results to `runs/<YYYYMMDD_HHMMSS>_gene_<GENE>/` or `runs/<YYYYMMDD_HHMMSS>_crosstalk_<CORE>_x_<TARGET>/` (`run_info.json` records the parameters; `out=` changes the root).

## Setup
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install requests pandas numpy networkx scipy pyyaml jupyter
python fetch_depmap.py --out ~/data/depmap/        # optional, ~1.5 GB, CC-BY 4.0 (needed for sl and co-essentiality)
```
Live runs download and cache (folder `cache/`): SIGNOR dump, Reactome FI network (~20 MB), Lambert TF list,
DoRothEA regulons via OmniPath, UniProt receptor list / locations, STRING physical edges, HPA tissue table if `--tissue` is given.

## Common options
| option | meaning | default |
|---|---|---|
| `--max-depth N` | maximum hops from the core gene | 6 |
| `--min-score S` | drop nodes whose best path score (product of edge weights × hub penalty) is below S | 0.15 |
| `--hub-penalty H` | 0 = none, 0.5 = strong down-weighting of highly connected nodes | 0.25 |
| `--tissue NAME` | keep only genes expressed in an HPA tissue (nTPM ≥ 1), e.g. `liver` | off |
| `--expressed-file F` | your own list of expressed genes (one symbol per line) | off |
| `--depmap DIR` | adds a co-essentiality column (Pearson r of CRISPR profiles with the core gene) | off |
| `--no-reactome`, `--reactome-fi FILE` | skip Reactome / use a local FI file | |
| downstream & crosstalk: `--regulon dorothea\|collectri --conf ABC --regulon-file F --max-targets N` | regulon source and confidence | dorothea ABC |

## How each mode works
**Signed walk.** Best-path search from the core gene over states (node, sign). Edge weight: SIGNOR 0.9, Reactome 0.7,
both 0.95; each node multiplies by a hub penalty. Sign = product of edge effects (+1 activation, −1 inhibition); a node
whose best positive and best negative paths are within 20 % of each other is *mixed*. The walk never expands past an
endpoint, past `--max-depth`, or below `--min-score`.

| mode | walk | endpoints | extra | outputs (inside `runs/<YYYYMMDD_HHMMSS>_<mode>_<genes>/`) |
|---|---|---|---|---|
| upstream | backwards | receptors (UniProt KW-0675 ∩ cell membrane) | — | `upstream_nodes.csv`, `receptors.csv`, `flow.html`, `cell_map.html` |
| downstream | forwards | TFs (Lambert 2018) | TF → genes via DoRothEA; gene sign = TF sign × mode of regulation | `downstream_nodes.csv`, `transcription_factors.csv`, `target_genes.csv`, `flow.html`, `cell_map.html` |
| sl | — | — | Welch t-test of CRISPR gene effect, core-lost vs intact lines, BH-FDR; hits flagged if already on the pathway | `synthetic_lethality.csv/.html` (volcano) |
| crosstalk | forwards ×2 | TFs | D(C) walked **with X removed**, D(X) walked **with C removed** (only routes independent of the other gene count) | `shared.csv`, `interference.csv`, `tf_agreement.csv`, `gene_agreement.csv`, `summary.json`, `crosstalk.html` |

**Cross-talk definitions**
- *shared nodes*: reached by both walks; sign agreement reported per node.
- *interference edge*: a node reachable only from X acting directly on a node of C's pathway (or on C itself), with the
  SIGNOR residue when known; `net_effect_of_goi_on_target` = sign(X→u) × effect(u→v), compared with C's own effect on v
  → synergistic / antagonistic / mixed. `via_feedback_loop` marks from-nodes that are also upstream of X.
- *cross-talk index*: shared endpoints (TFs + target genes) / all endpoints of C.

Every node table carries: hops, path score, sign, ambiguous, best_path, feedback (node is also upstream of the seed),
coessentiality_r, physical_with_seed (STRING / Reactome undirected), degree, compartment (UniProt).

## Views
In notebooks use the Plotly versions (`targetnet.plots`, via `api.figure(view_path)`), rendered inline with plotly.js from the CDN so the .ipynb stays small; the D3 HTML files are the full standalone views (`api.show`).
- `flow.html` — layered flow: seed → hops → endpoints column (+ gene summary boxes per TF in downstream mode; click to list genes).
  Outline colour = net sign from the seed (green activated, red inhibited, dashed mixed); long dashes = feedback node.
- `cell_map.html` — same nodes placed by subcellular compartment (location view).
- `crosstalk.html` — three lanes (C pathway / shared / X pathway); thick green/red edges = interference; tables of edges and endpoint agreement.
- `synthetic_lethality.html` — volcano plot + ranked hits (orange), pathway members in blue.

## Offline UniProt (recommended)
One bulk file replaces hundreds of UniProt API calls (receptors + compartments):
```bash
mkdir -p ~/Documents/database/uniprot && curl -o ~/Documents/database/uniprot/uniprot_human_reviewed.tsv \
  "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&fields=accession,gene_primary,keyword,cc_subcellular_location&format=tsv"
```
and set `databases.uniprot.path` to it.

## DoRothEA / CollecTRI regulons
Downloaded from OmniPath with a per-record licence filter (`regulon.license: commercial` by default — only records whose
sources allow commercial use). To keep a local copy in your database folder:
```bash
python target_net.py regulon-url --save ~/Documents/database/dorothea/dorothea_ABC_commercial.tsv
```
then set `databases.regulon.path` to that file. (`regulon-url` alone just prints the URL for curl.)
The loader also accepts the simple format `tf,target,mor,confidence` and decoupler-style `source,target,weight`.

## Reactome FI input (important for result quality)
Reactome "functional interactions" collapse whole pathways into gene–gene edges (e.g. hundreds of GPCRs `->` AKT1 in one hop).
By default the loader therefore (a) drops NBC-*predicted* FIs (`min_score: 1.0`), (b) uses an FI's direction only when its
annotation names a mechanism (activate / inhibit / catalyse / phosphorylate / expression-regulate) — other FIs become
undirected context edges — and (c) gives Reactome-only directed edges weight 0.5 (SIGNOR 0.9), so curated SIGNOR routes
win. Default is now `mode: context` (walk on SIGNOR direct relations only; Reactome as context layer): on the first live run,
`mode: directed` still contributed 109k directed FIs vs 15k SIGNOR edges and the AKT1 upstream walk reached 2,667 nodes.

## SIGNOR input
If the SIGNOR site is down and no local file exists, the loader falls back to the **OmniPath mirror of SIGNOR** automatically
(`databases.signor.source: omnipath` forces it). The mirror has direction, sign, references and residues (enzyme-substrate
table) but not SIGNOR's SCORE / DIRECT columns, so edge weights are flat (0.9) and the indirect-relation filter cannot apply.
Either the API dump or SIGNOR's full download `all_data_DD_MM_YY.tsv` works (`databases.signor.path`). The full file is
filtered to human (TAX_ID 9606) protein→protein rows; the SIGNOR confidence `SCORE` scales the edge weight (×0.6–1.0).
Relations with `DIRECT = NO` (pathway-level statements such as "PDGFRA up-regulates AKT1") are dropped by default
(`signor.indirect: drop | weight | keep`) because they short-circuit the mechanistic route.
Rows whose entity is a **complex** or **protein family** (e.g. mTORC1, "AKT family") are currently dropped — expanding them to
member proteins (SIGNOR complexes/families files) is an open item.

## Licenses
SIGNOR CC-BY 4.0 · Reactome CC0 · STRING CC-BY 4.0 · UniProt CC-BY 4.0 · Lambert TF list CC-BY · DepMap CC-BY 4.0 ·
HPA CC BY-SA 3.0 (used only as an expression filter, not redistributed) · DoRothEA: regulons aggregate several sources;
confirm the A–C subset is clear for commercial redistribution (its successor CollecTRI is available via `--regulon collectri`).
KEGG and PhosphoSitePlus are deliberately not used (commercial licenses required).

## Known limits / next steps
1. Database direction is context-free; `--tissue` helps but does not replace measured data. Phosphoproteomics after
   perturbation maps directly onto the interference edges (which sites follow C, which follow X).
2. Benchmark before trusting weights: `upstream --gene MAPK1` should return RTKs via GRB2/SOS/RAS/RAF; `crosstalk --causal AKT1 --goi MAPK1`
   should recover RSK→GSK3B S9, RSK→BAD S112, ERK⊣SOS1 (negative feedback) — it does on the demo graph.
3. Logic simulation (Boolean model of the signed graph) is the natural next layer for "what if X is inhibited".

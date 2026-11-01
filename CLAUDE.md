# target_net — project notes for Claude Code

Network-based drug-target research around a core gene. Entry points: `target_net.py` (CLI), `targetnet/api.py` (notebook API),
`notebooks/*.ipynb` (templates), `.claude/skills/target-network-research/SKILL.md` (the workflow skill — follow it for analysis requests).

## Layout
- `targetnet/config.py`  YAML config: database paths, download policy, walk/SL thresholds (`config.example.yaml` -> `target_net.yaml`)
- `targetnet/data.py`    loaders (SIGNOR, Reactome FI, STRING, UniProt, Lambert TFs, DoRothEA/OmniPath, DepMap, HPA) + curated demo data
- `targetnet/walk.py`    signed best-path walk (node set, score) + signed influence propagation (consensus sign, influence); regulon extension; cross-talk logic
- `targetnet/plots.py`   the four views as Plotly figures (notebook inline AND the .html files in run folders; built from the .json next to each view)
- `targetnet/templates/` original D3 pages, written only with --d3 as <view>_d3.html
- `fetch_depmap.py`      DepMap downloader (portal API, Figshare fallback)
- `docs/INSTALL_LINUX.md` server install (venv / Docker); `Dockerfile`, `docker-compose.yml`
- `docs/DATA_SETUP.md`   manual for every database file; `scripts/setup_data.sh` downloads the automatic ones

## Conventions (owner's preferences)
- Python; minimal first version, then extend. Notebook/py code: a comment or markdown cell above each block explaining it.
- Commercial-use data only (no KEGG, no PhosphoSitePlus); state licences when adding a source.
- Demo mode (`--demo`) must stay runnable offline; demo DepMap values are simulated and must be labelled as such.
- Plots: flat, muted colours, white background, English labels.

## Test
`for m in "upstream --gene AKT1" "downstream --gene AKT1" "sl --gene AKT1" "crosstalk --causal FGFR3 --goi HSP90AA2"; do python target_net.py $m --demo; done`
Notebooks: `jupyter nbconvert --to notebook --execute notebooks/gene_research.ipynb --output /tmp/t.ipynb`
Expected textbook recovery: AKT1 upstream -> RTKs via PI3K/PDPK1; MAPK1 x AKT1 -> RSK->GSK3B S9, ERK-|SOS1 feedback.

## Open items
- Live-data calibration (AKT1, 2026-11-01): SIGNOR 15,428 direct+indirect pp edges; Reactome FI 'mechanistic' still 109k -> default mode=context; SIGNOR DIRECT=NO dropped; walk defaults max_depth 5 / min_score 0.2. Second run (SIGNOR direct only, Reactome context): upstream 903 nodes / 79 receptors; receptor ranking sensible (RET, EGFR, FGFR1, IGF1R, PDGFRs). Best-path sign misled for EGFR (via PTPN1) -> consensus sign added. min_score calibration on live AKT1 upstream: 0.2 -> 903 nodes/79 receptors, 0.3 -> 222/13 -> default set to 0.25. Note: SIGNOR entity 'PDK1' = pyruvate dehydrogenase kinase, not PDPK1
- SIGNOR complexes / protein families: expand to member proteins (currently dropped by load_signor; `signor-lookup --gene X` shows what is lost). Live case ACVRL1 (2026-11-01): 0 outgoing protein edges -> downstream empty; research_gene now continues (upstream + SL) and records downstream_note; mitigations added: walk.seed_rescue (seed's own mechanistic Reactome FI edges when it has none) and databases.custom_edges (user tsv); membership available offline-free via OmniPath `complexes?resources=SIGNOR`
- Ligand layer (2026-11-01, from GDF2 -> ACVRL1): upstream endpoints = receptors + secreted ligands (UniProt Secreted + Cytokine/Growth factor/Hormone); receptors may step one hop into ligands (not score-gated); ligands.csv, teal 'ligands' column, extracellular zone in cell map, role=ligand in gene_table
- Seed partners (2026-11-01): partners.csv = physical/complex partners of the seed from the context layer with their own signed outputs; NOT walked (type II -> type I receptor logic makes output-borrowing wrong). Complex-mediated routes go through custom_edges
- Chaperone/client edge layer (Taipale 2012, Picard list) beyond SIGNOR
- Cross-talk: live AKT1 x ACVR1B gave index 1.0 / 2,265 nodes (hub saturation) -> local walks (depth 3, min_score 0.3), hypergeometric TF-overlap enrichment, confidence-gated agreement, proximity-ranked edges. Still to do: automatic baseline (same target vs N random core genes)
- Boolean/logic simulation of the signed graph

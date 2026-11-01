# target_net — project notes for Claude Code

Network-based drug-target research around a core gene. Entry points: `target_net.py` (CLI), `targetnet/api.py` (notebook API),
`notebooks/*.ipynb` (templates), `.claude/skills/target-network-research/SKILL.md` (the workflow skill — follow it for analysis requests).

## Layout
- `targetnet/config.py`  YAML config: database paths, download policy, walk/SL thresholds (`config.example.yaml` -> `target_net.yaml`)
- `targetnet/data.py`    loaders (SIGNOR, Reactome FI, STRING, UniProt, Lambert TFs, DoRothEA/OmniPath, DepMap, HPA) + curated demo data
- `targetnet/walk.py`    signed best-path walk over (node, sign) states; regulon extension; cross-talk logic
- `targetnet/templates/` D3 views: flow.html, cell_map.html, crosstalk.html, sl_plot.html (each view also written as .json)
- `targetnet/plots.py`   Plotly versions of the four views for inline notebook use (read the .json next to each view)
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
- Live-data calibration (AKT1, 2026-11-01): SIGNOR 15,428 direct+indirect pp edges; Reactome FI 'mechanistic' still 109k -> default mode=context; SIGNOR DIRECT=NO dropped; walk defaults max_depth 5 / min_score 0.2. Re-check: upstream/downstream node counts should land ~100-400; receptors for AKT1 via PIK3CA/PDPK1
- SIGNOR complexes / protein families: expand to member proteins (currently dropped by load_signor); membership available offline-free via OmniPath `complexes?resources=SIGNOR`
- Chaperone/client edge layer (Taipale 2012, Picard list) beyond SIGNOR
- Cross-talk baseline: same target vs several unrelated core genes
- Boolean/logic simulation of the signed graph

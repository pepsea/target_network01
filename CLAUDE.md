# target_net — project notes for Claude Code

Network-based drug-target research around a core gene. Entry points: `target_net.py` (CLI), `targetnet/api.py` (notebook API),
`notebooks/*.ipynb` (templates), `.claude/skills/target-network-research/SKILL.md` (the workflow skill — follow it for analysis requests).

## Layout
- `targetnet/config.py`  YAML config: database paths, download policy, walk/SL thresholds (`config.example.yaml` -> `target_net.yaml`)
- `targetnet/data.py`    loaders (SIGNOR, Reactome FI, STRING, UniProt, Lambert TFs, DoRothEA/OmniPath, DepMap, HPA) + curated demo data
- `targetnet/walk.py`    signed best-path walk (node set, score) + signed influence propagation (consensus sign, influence); regulon extension; cross-talk logic
- `targetnet/screen.py`  cross-talk screen: one core gene vs a list of targets -> interference score, TF enrichment, dependency, direction, degree-matched random baseline, composite rank (notebooks/crosstalk_screen.ipynb, CLI crosstalk-screen)
- `targetnet/interference.py` pathway interference: CORE network + 'TARGET inhibited' perturbation (entry edges from TARGET's side, propagated along CORE's signed edges) -> flow_change per node / gene (notebooks/pathway_interference.ipynb, CLI pathway-interference, plots.cell_flow_figure = cell cross-section layout)
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
`for m in "upstream --gene AKT1" "downstream --gene AKT1" "sl --gene AKT1" "crosstalk --causal FGFR3 --goi HSP90AA2" "crosstalk-screen --causal FGFR3 --goi-list AKT1,MAPK1,HSP90AA2 --n-random 5" "pathway-interference --causal FGFR3 --goi MAPK1"; do python target_net.py $m --demo; done`
Notebooks: `jupyter nbconvert --to notebook --execute notebooks/gene_research.ipynb --output /tmp/t.ipynb`
Expected textbook recovery: AKT1 upstream -> RTKs via PI3K/PDPK1; MAPK1 x AKT1 -> RSK->GSK3B S9, ERK-|SOS1 feedback.

## Open items
- Live-data calibration (AKT1, 2026-11-01): SIGNOR 15,428 direct+indirect pp edges; Reactome FI 'mechanistic' still 109k -> default mode=context; SIGNOR DIRECT=NO dropped; walk defaults max_depth 5 / min_score 0.2. Second run (SIGNOR direct only, Reactome context): upstream 903 nodes / 79 receptors; receptor ranking sensible (RET, EGFR, FGFR1, IGF1R, PDGFRs). Best-path sign misled for EGFR (via PTPN1) -> consensus sign added. min_score calibration on live AKT1 upstream: 0.2 -> 903 nodes/79 receptors, 0.3 -> 222/13 -> default set to 0.25. Note: SIGNOR entity 'PDK1' = pyruvate dehydrogenase kinase, not PDPK1
- SIGNOR complexes / protein families: expand to member proteins (currently dropped by load_signor; `signor-lookup --gene X` shows what is lost). Live case ACVRL1 (2026-11-01): 0 outgoing protein edges -> downstream empty; research_gene now continues (upstream + SL) and records downstream_note; mitigations added: walk.seed_rescue (seed's own mechanistic Reactome FI edges when it has none) and databases.custom_edges (user tsv); membership available offline-free via OmniPath `complexes?resources=SIGNOR`
- Ligand layer (2026-11-01, from GDF2 -> ACVRL1): upstream endpoints = receptors + secreted ligands (UniProt Secreted + Cytokine/Growth factor/Hormone); receptors may step one hop into ligands (not score-gated); ligands.csv, teal 'ligands' column, extracellular zone in cell map, role=ligand in gene_table
- OmniPath layer (2026-11-01): databases.omnipath, licence-filtered (purpose=commercial, local filter via /resources), signed_only. Default switched gaps -> full on the owner's principle 'more interactions; do not restrict the analysis for the sake of the figure' (figures wrap crowded hops, node cap is view-only). First live run under full (2026-11-01): graph 8,492 nodes / 81,296 edges (SIGNOR 13k + OmniPath 78,851); AKT1 upstream at min_score 0.25 -> 3,846 nodes / 392 receptors / 262 ligands = flooded -> `python target_net.py calibrate --gene AKT1` sweeps min_score/hub_penalty (walk only); calibrate table (AKT1, max_depth 6): hub 0.25 -> receptors 375/307/203/110/42/17 at min_score 0.20/0.25/0.30/0.35/0.40/0.50; hub 0.40 -> 220/82/38/26/16/9. Defaults set to walk min_score 0.35 / hub_penalty 0.3, confirmed live: AKT1 upstream 932 nodes / 49 receptors / 128 ligands (max_depth 6); crosstalk min_score 0.4; downstream still ~2,000 nodes / 300 TFs at these values (AKT1 is a hub; acceptable per owner: more interactions, view handles it). Signs: OmniPath consensus_* fields resolve stim+inh conflicts (ACVRL1->SMAD1 was dropped as 'unsigned' before); conflicting-without-consensus kept as ambiguous. COMPLEX:A_B entities expanded to member proteins (ACVR2A_ACVRL1 -> SMAD5 etc.). Seed rescue weight lowered to 0.3 after ACVRL1 -> PPP1CA (Reactome FI artefact) dragged 1,921 nodes in. After the fixes ACVRL1 downstream = SMAD1 (0.60), SMAD9, SMAD5 at hop 1 (all Lambert TFs -> walk ends there by design) + 27 DoRothEA genes: correct and compact
- Seed partners (2026-11-01): partners.csv = physical/complex partners of the seed from the context layer with their own signed outputs; NOT walked (type II -> type I receptor logic makes output-borrowing wrong). Complex-mediated routes go through custom_edges
- walk.through_tf (2026-11-01, default true): downstream walk continues past TFs along the graph's transcriptional-regulation edges (SMAD1 -> ID1 ...); TFs stay is_endpoint for tables/cross-talk; flow view then places TFs at their hop (no TF column) with DoRothEA counts under the node. false = old behaviour (stop at first TF)
- View node cap (2026-11-01, owner's rule for ALL views): per hop/lane, the highest-scoring genes fill the budget (round-robin over columns); FOCUS genes are highlighted only if they are among them, never forced in (missing ones listed in the title). Cross-talk additionally keeps the ends of the drawn interference edges
- Hop-1 views (2026-11-01): upstream/hop1.html and downstream/hop1.html = seed + every direct neighbour, no node cap (views upstream_hop1 / downstream_hop1); the flow views stay best-path oriented and capped
- Impact view (2026-11-01): crosstalk writes impact.csv + core_impact.html = C's flow recoloured by X's effect (direct / shared / indirect / none; fill = X's sign, outline = agreement); walk.impact_on_core
- Chaperone/client edge layer (Taipale 2012, Picard list) beyond SIGNOR. Live FGFR3 x HSP90AA1 (2026-11-01): no HSP90AA1->FGFR3 edge in SIGNOR/OmniPath-gaps although FGFR3 is a known HSP90 client -> added (a) OmniPath edges incident to the analysed genes always taken in mode=gaps, (b) omnipath.signed_only (unsigned directed records flooded FGFR3's hop 1), (c) crosstalk physical_contacts.csv + dotted blue lines (STRING/Reactome contact without signed edge)
- Cross-talk: live AKT1 x ACVR1B gave index 1.0 / 2,265 nodes (hub saturation) -> local walks (depth 3, min_score 0.3), hypergeometric TF-overlap enrichment, confidence-gated agreement, proximity-ranked edges. Automatic baseline now exists in the screen (degree-matched random targets vs the same core, baseline_z); the single-pair notebook still has none
- Pathway interference (2026-11-01): third analysis, see interference.py; demo FGFR3 vs MAPK1 inhibited -> STAT3 down, SOS1/GRB2/GAB1/FRS2/RAF1 up (feedback release). Propagation is linear signed influence restricted to the CORE network (no saturation / logic) -> the Boolean simulation item would refine it
- Boolean/logic simulation of the signed graph

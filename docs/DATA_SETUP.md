# Data setup manual — where every file comes from

All databases live under one folder (default `~/Documents/database/`) and are declared in `target_net.yaml`.
Anything left `null` in the YAML is downloaded once into `cache_dir` when `allow_download: true`.
`bash scripts/setup_data.sh` performs every automatic download below in one go.

Legend: **required** = the walk cannot run without it · *optional* = adds a layer or a filter.

| # | database | needed for | YAML key | expected file | size | licence |
|---|---|---|---|---|---|---|
| 1 | SIGNOR | signaling edges (**required**) | `signor.path` | `all_data_DD_MM_YY.tsv` | ~30 MB | CC-BY 4.0 |
| 2 | Reactome FI | signaling edges (**required**, with 1) | `reactome_fi.path` | `FIsInGene_MMDDYY_with_annotations.txt[.zip]` | ~20 MB | CC0 |
| 3 | UniProt (bulk) | receptors + compartments | `uniprot.path` | `uniprot_human_reviewed.tsv` | ~30 MB | CC-BY 4.0 |
| 4 | Lambert TF list | TF endpoints (**required** for downstream) | `tf_list.path` | `DatabaseExtract_v_1.01.csv` | 2 MB | CC-BY |
| 5 | DoRothEA (OmniPath) | TF → target genes | `regulon.path` | `dorothea_ABC_commercial.tsv` | ~5 MB | per-record (filtered to commercial) |
| 6 | DepMap | co-essentiality, synthetic lethality | `depmap.dir` | `CRISPRGeneEffect.csv`, `OmicsSomaticMutations.csv`, (`OmicsExpressionProteinCodingGenesTPMLogp1.csv` (newer releases: `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`)) | 1 GB + 0.5 GB (+ 0.5 GB) | CC-BY 4.0 |
| 7 | HPA tissue | `--tissue` expression filter | `hpa_tissue.path` | `rna_tissue_consensus.tsv[.zip]` | ~50 MB | CC BY-SA 3.0 |
| 8 | STRING | physical-contact annotation | — (API) | — | — | CC-BY 4.0 |

---

## 1. SIGNOR — signed, directed, residue-level signaling (required)
- Page: https://signor.uniroma2.it/downloads.php → **"Download all data"** → `all_data_DD_MM_YY.tsv` (date of the release).
- Or per organism via API: `curl -o signor_human.tsv "https://signor.uniroma2.it/getData.php?organism=9606"`.
- The loader keeps human (`TAX_ID 9606`) protein→protein rows and uses the `SCORE` column as edge confidence.
- Optional, for the open item "expand complexes/families": from the same page download **Complexes** (`SIGNOR_complexes.csv`) and **Protein families** (`SIGNOR_PF.csv`).
- Updated: monthly. Re-download a few times a year; keep the date in the file name.
- Site down? The pipeline falls back to the OmniPath mirror of SIGNOR (`databases.signor.source: omnipath`); it lacks SCORE/DIRECT, so prefer the native file when available.

## 2. Reactome Functional Interaction network (required)
- Page: https://reactome.org/download-data → section **"Functional Interactions"** (Reactome FI) → `FIsInGene_<MMDDYY>_with_annotations.txt.zip`.
- Columns used: `Gene1 Gene2 Annotation Direction Score`. Direction `->`/`<-` activation, `-|`/`|-` inhibition, `-` undirected (goes to the physical layer).
- The `.txt` and the `.zip` are both accepted. Updated: about yearly.

## 3. UniProt bulk table (recommended; replaces hundreds of API calls)
```bash
mkdir -p ~/Documents/database/uniprot
curl -o ~/Documents/database/uniprot/uniprot_human_reviewed.tsv \
  "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&fields=accession,gene_primary,keyword,cc_subcellular_location&format=tsv"
```
- Gives ~20,500 reviewed human proteins. Receptors = keyword *Receptor* + location *Cell membrane*; compartments = first-mentioned subcellular location.
- Without it the tool calls the UniProt API per gene (slow, needs internet). Updated: every 8 weeks; once a year is fine.

## 4. Lambert 2018 human transcription factors (required for downstream / crosstalk)
- Page: http://humantfs.ccbr.utoronto.ca/download.php → **"Human TFs – Database Extract (v1.01)"** → `DatabaseExtract_v_1.01.csv`.
- Direct: `curl -o DatabaseExtract_v_1.01.csv http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv`
- Rows with `Is TF? == Yes` (1,639 TFs) are the downstream endpoints. Static (2018).

## 5. DoRothEA regulons via OmniPath (TF → target gene)
```bash
python target_net.py regulon-url --save ~/Documents/database/dorothea/dorothea_ABC_commercial.tsv
```
- Builds the URL from the YAML (`source: dorothea`, `conf: ABC`, `license: commercial`), downloads, and verifies it loads.
- `license: commercial` keeps only records whose source databases allow commercial use; `academic` returns everything.
- Alternative: `source: collectri` (DoRothEA's 2023 successor, same command). Also accepted: a tsv with `tf,target,mor,confidence`.

## 6. DepMap (CRISPR gene effect, mutations, expression)
```bash
python fetch_depmap.py --out ~/Documents/database/depmap --list   # shows release and URLs
python fetch_depmap.py --out ~/Documents/database/depmap          # downloads the three files
```
- Manual: https://depmap.org/portal/data_page → **All Downloads** → latest `DepMap Public YYQN` → `CRISPRGeneEffect.csv`, `OmicsSomaticMutations.csv`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (newer releases: `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`). Mirrored on Figshare.
- `CRISPRGeneEffect.csv` alone → co-essentiality. Synthetic lethality needs `OmicsSomaticMutations.csv` (mutation route) or the expression file (fallback: lowest-quartile expression = "lost").
- Updated: quarterly (YYQ1…Q4). Keep the release name in a `README` inside the folder.

## 7. Human Protein Atlas tissue expression (for `--tissue`)
- Page: https://www.proteinatlas.org/about/download → **"RNA consensus tissue gene data"** → `rna_tissue_consensus.tsv.zip`.
- Columns used: `Gene name`, `Tissue`, `nTPM` (genes with nTPM ≥ 1 in the chosen tissue are kept). Both `.tsv` and `.zip` accepted.
- Licence is CC BY-SA: used only as a filter, never redistributed in outputs. Updated: yearly.

## 8. Custom edges (optional, your own curation or in-house data)
A tsv/csv at `databases.custom_edges.path` with columns `src, dst, effect` (+1 / −1 / 0, or activate / inhibit) and optional
`mechanism, residue, score (0–1), source`. Lines starting with `#` are ignored. Edges are merged into the walk graph as source
`custom` (weight 0.8 × score; a pair also in SIGNOR becomes `SIGNOR+custom`, 0.95). Use it for relations SIGNOR lacks
(e.g. `ACVRL1 SMAD1 activate phosphorylation`) or for AirID / phosphoproteomics results you want the walk to see. Every edge you add is
your assertion — keep the file under version control next to the project.

## 9. STRING (physical contact with the seed)
- API only (`https://string-db.org/api/...`, experimental + database channels, score ≥ 700). Needs internet; with `allow_download: false` this layer is skipped.
- A local STRING file is not supported yet (the full `protein.physical.links.detailed.v12.0.txt.gz` is 1.5 GB; open item if offline use matters).

---

## Folder layout and YAML that match this manual
```
~/Documents/database/
  signor/all_data_12_07_26.tsv
  reactome/FIsInGene_061424_with_annotations.txt
  uniprot/uniprot_human_reviewed.tsv
  tf/DatabaseExtract_v_1.01.csv
  dorothea/dorothea_ABC_commercial.tsv
  depmap/CRISPRGeneEffect.csv  OmicsSomaticMutations.csv  [OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv]
  hpa/rna_tissue_consensus_tissues.tsv
  target_net_cache/            (anything the tool downloads itself)
```
`config.example.yaml` in the project root is the complete template for this layout (every key with its default and a
comment); `cp config.example.yaml target_net.yaml` and edit only the paths that differ.

## Verify
```bash
python target_net.py check-config                     # per database: path, ok / MISSING / not set, size, options, fallback
python target_net.py upstream --gene AKT1 --config target_net.yaml
```
Expected log lines: `[config] target_net.yaml` · `[SIGNOR] source: … rows N -> human M -> protein-protein unique edges K` (K ≫ 10,000) ·
`[Reactome] source: …` · `[UniProt] bulk table: ~20,000 genes` · `[receptors] ~1,000+ from bulk UniProt` · `[Lambert] 1639 TFs` (downstream) ·
`[dorothea] … TF-target edges` (downstream) · `[DepMap] loading gene effect` (if configured). Then check that `receptors.csv` for AKT1 lists RTKs
(EGFR, INSR, IGF1R, PDGFR…) within 3–4 hops.

## Strict offline mode
Set `allow_download: false` once all files are local. Any missing file then stops the run with the exact YAML key to fix
(STRING is skipped silently).

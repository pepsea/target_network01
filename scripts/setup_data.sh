#!/usr/bin/env bash
# Download every automatically downloadable database for target_net into one folder.
# Usage:  bash scripts/setup_data.sh [~/Documents/database]
# Skips files that already exist. SIGNOR and Reactome FI need a click on their download pages (see docs/DATA_SETUP.md);
# the script prints what is still missing at the end.
set -euo pipefail
ROOT="${1:-$HOME/Documents/database}"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$ROOT"/{signor,reactome,uniprot,tf,dorothea,depmap,hpa,target_net_cache}
echo "database root: $ROOT"

get() {  # get <dest> <url>
  if [ -s "$1" ]; then echo "  [skip] $(basename "$1") exists"; else echo "  [get ] $(basename "$1")"; curl -L --fail --retry 3 -o "$1" "$2"; fi
}

echo "[3] UniProt human reviewed (receptors + compartments)"
get "$ROOT/uniprot/uniprot_human_reviewed.tsv" \
  "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&fields=accession,gene_primary,keyword,cc_subcellular_location&format=tsv"

echo "[4] Lambert 2018 TF list"
mkdir -p "$ROOT/omnipath"; python3 "$(dirname "$0")/../target_net.py" omnipath-url --save "$ROOT/omnipath/omnipath_signed.tsv" || echo "OmniPath signed layer: download later with: python target_net.py omnipath-url --save $ROOT/omnipath/omnipath_signed.tsv"
get "$ROOT/tf/DatabaseExtract_v_1.01.csv" "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"

echo "[5] DoRothEA regulons via OmniPath (licence filter from target_net.yaml / default commercial)"
if [ -s "$ROOT/dorothea/dorothea_ABC_commercial.tsv" ]; then echo "  [skip] dorothea_ABC_commercial.tsv exists";
else (cd "$HERE" && python target_net.py regulon-url --save "$ROOT/dorothea/dorothea_ABC_commercial.tsv"); fi

echo "[7] HPA tissue consensus"
get "$ROOT/hpa/rna_tissue_consensus.tsv.zip" "https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip"

echo "[6] DepMap (large; portal API with Figshare fallback)"
if [ -s "$ROOT/depmap/CRISPRGeneEffect.csv" ]; then echo "  [skip] CRISPRGeneEffect.csv exists";
else (cd "$HERE" && python fetch_depmap.py --out "$ROOT/depmap"); fi

echo "[1] SIGNOR"
if ls "$ROOT"/signor/all_data_*.tsv >/dev/null 2>&1 || [ -s "$ROOT/signor/signor_human.tsv" ]; then echo "  [skip] SIGNOR file exists";
else echo "  [get ] per-organism dump (or download 'all data' manually from https://signor.uniroma2.it/downloads.php)";
     curl -L --fail --retry 3 -o "$ROOT/signor/signor_human.tsv" "https://signor.uniroma2.it/getData.php?organism=9606" \
       || echo "  [warn] SIGNOR site unreachable - the pipeline falls back to the OmniPath mirror automatically (set databases.signor.source: omnipath to force it)"; fi

echo "[2] Reactome FI"
if ls "$ROOT"/reactome/FIsInGene_*_with_annotations.txt* >/dev/null 2>&1; then echo "  [skip] Reactome FI file exists";
else echo "  [get ] trying the 2023 release URL (if it fails, download from https://reactome.org/download-data -> Functional Interactions)";
     curl -L --fail --retry 3 -o "$ROOT/reactome/FIsInGene_070323_with_annotations.txt.zip" \
       "https://reactome.org/download/tools/ReatomeFIs/FIsInGene_070323_with_annotations.txt.zip" || echo "  [warn] Reactome download failed - get it manually"; fi

echo
echo "Summary:"
for f in signor/all_data_*.tsv signor/signor_human.tsv reactome/FIsInGene_*_with_annotations.txt* uniprot/uniprot_human_reviewed.tsv \
         tf/DatabaseExtract_v_1.01.csv dorothea/dorothea_ABC_commercial.tsv depmap/CRISPRGeneEffect.csv depmap/OmicsSomaticMutations.csv hpa/rna_tissue_consensus*.tsv*; do
  if ls "$ROOT"/$f >/dev/null 2>&1; then echo "  ok      $f"; else echo "  MISSING $f"; fi
done
echo
echo "Next: write target_net.yaml with these paths (template: docs/DATA_SETUP.md), then"
echo "      python target_net.py upstream --gene AKT1 --config target_net.yaml"

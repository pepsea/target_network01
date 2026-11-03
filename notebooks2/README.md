# notebooks2 — 標的遺伝子の流れの図（最終版）

遺伝子名を 1 つ入れると、その遺伝子の**上流・下流の遺伝子の選択**と**図**ができるノートブックです。
`notebooks/network_development.ipynb`（試行錯誤）で最終的に採用した部分だけをまとめています。

| ファイル | 内容 |
|---|---|
| `flow_view.ipynb` | 本体。`config.yaml` の一番上の `gene` に標的遺伝子名を入れて、上から順に実行する |
| `config.yaml` | パラメータとデータベースの場所（標的遺伝子名以外はすべてここ） |
| `focus_demo.txt` | フォーカス遺伝子のデモのリスト（`config.yaml` の `focus.path`） |
| `results/<日付>/<時刻>_<遺伝子名>/` | 実行ごとの出力: 図（`.png`・`.html`）、遺伝子の表（`_genes.csv`）、使った設定（`config_used.yaml`）、フォーカスリストの写し。**GitHub には上げない**（`.gitignore`） |

実行: `jupyter nbconvert --to notebook --execute flow_view.ipynb --output flow_view.ipynb`（約 1 分）、またはノートブックを開いて上から実行する。

## データのダウンロード（すべて）

保存先は、既定では `~/Documents/database/`。場所を変えるときは `config.yaml` の `databases:` を書き換える。
ダウンロードしたファイルは、サイズ・出どころ・ライセンスを確かめてから使う。

| # | データ | 使い道 | `config.yaml` のキー | ライセンス |
|---|---|---|---|---|
| 1 | SIGNOR | 符号つきのシグナル伝達の辺・低分子・複合体 | `databases.signor.path` | CC BY 4.0 |
| 2 | OmniPath（相互作用） | 符号つきの辺（記録ごとにライセンス判定する） | `databases.omnipath.path` | リソースごとに異なる（ノートブックの 3 で判定） |
| 3 | OmniPath（リソース一覧） | 上のライセンス判定 | `databases.omnipath.resources` | — |
| 4 | UniProt（キーワード・局在） | 受容体・リガンド・機能・局在 | `databases.uniprot.path` | CC BY 4.0 |
| 5 | UniProt（遺伝子名） | 表とホバーの遺伝子名 | `databases.uniprot.names` | CC BY 4.0 |
| 6 | Lambert の転写因子一覧 | 転写因子の定義 | `databases.tf_list.path` | ファイルに記載なし（不明として使う） |
| 7 | DepMap | 点の大きさ（入力遺伝子との相関） | `databases.depmap.dir` | CC BY 4.0 |
| 8 | Human-GEM | 代謝の辺（酵素と代謝物） | `databases.human_gem.path` / `.genes` | CC BY 4.0（公開ページの記載。ローカルのファイルには記載なし） |

### 1. SIGNOR（必須）
- ページ: https://signor.uniroma2.it/downloads.php → **Download all data** → `all_data_DD_MM_YY.tsv`
- 保存先: `~/Documents/database/signor/`。使う列: `ENTITYA`・`TYPEA`・`ENTITYB`・`TYPEB`・`EFFECT`・`MECHANISM`・`RESIDUE`・`TAX_ID`・`DIRECT`・`SCORE`
- 毎月更新される。ファイル名に日付を残す。

### 2・3. OmniPath（必須）
```bash
mkdir -p ~/Documents/database/target_net_cache && cd ~/Documents/database/target_net_cache
curl -o omnipath_signed_omnipath_pathwayextra_kinaseextra_ligrecextra.txt \
  "https://omnipathdb.org/interactions?datasets=omnipath,pathwayextra,kinaseextra,ligrecextra&organisms=9606&genesymbols=1&directed=1&types=post_translational&fields=sources,references,type,is_directed,is_stimulation,is_inhibition,consensus_stimulation,consensus_inhibition"
curl -o omnipath_resources.txt "https://omnipathdb.org/resources?format=json"
```
- ノートブックの 3 で、記録が載せているリソースのライセンスを判定し、**要ライセンスまたは曖昧なリソースを 1 つでも含む記録は捨てる**（`config.yaml` の `license_strict`）。
- `notebooks`・`targetnet` の本体（`python target_net.py ...`）を一度動かすと、同じファイルが `target_net_cache/` に自動で保存される。

### 4. UniProt（キーワード・局在、必須）
```bash
mkdir -p ~/Documents/database/uniprot
curl -o ~/Documents/database/uniprot/uniprot_human_reviewed.tsv \
  "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&fields=accession,gene_primary,keyword,cc_subcellular_location&format=tsv"
```
受容体 = キーワード Receptor + 細胞膜、リガンド = 分泌 + Cytokine / Growth factor / Hormone。

### 5. UniProt（遺伝子名）
```bash
curl -o ~/Documents/database/uniprot/uniprot_human_names.tsv \
  "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&fields=accession,gene_primary,protein_name&format=tsv"
```
列は `Entry`・`Gene Names (primary)`・`Protein names`（約 20,000 行、約 2.4 MB）。表示には、最初の名前だけを使う。サーバーが一時的に失敗することがあるので、`wc -l` で行数を確かめて、失敗したら再実行する。

### 6. Lambert の転写因子一覧（必須）
- ページ: http://humantfs.ccbr.utoronto.ca/download.php → **Human TFs – Database Extract (v1.01)** → `DatabaseExtract_v_1.01.csv`
- 直接: `curl -L -o ~/Documents/database/tf/DatabaseExtract_v_1.01.csv http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv`
- `Is TF? == Yes` の行を転写因子とする。2018 年版で更新されない。

### 7. DepMap（点の大きさ用。無くても動く）
```bash
cd /path/to/target_network01
python fetch_depmap.py --out ~/Documents/database/depmap --list   # リリースと URL を確認
python fetch_depmap.py --out ~/Documents/database/depmap          # ダウンロード（約 1.5 GB）
```
- このノートブックで使うのは `CRISPRGeneEffect.csv` だけ。手動: https://depmap.org/portal/data_page → **All Downloads** → 最新の `DepMap Public YYQN`。
- 入力遺伝子が DepMap にないときは、点の大きさが全部最小になる。四半期ごとに更新される。

### 8. Human-GEM（代謝の辺。`metabolism.enabled: false` で使わない）
```bash
mkdir -p ~/Documents/database/Human-GEM && cd ~/Documents/database/Human-GEM
for f in Human-GEM.yml genes.tsv; do
  curl -O "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/$f"
done
```
- ページ: https://github.com/SysBioChalmers/Human-GEM（`model/` フォルダ）。`Human-GEM.yml`（反応・代謝物・遺伝子の対応、約 7.5 MB）と `genes.tsv`（Ensembl ID と遺伝子記号の対応）を使う。

## フォルダの配置

```
~/Documents/database/
  signor/all_data_12_07_26.tsv
  target_net_cache/omnipath_signed_omnipath_pathwayextra_kinaseextra_ligrecextra.txt
  target_net_cache/omnipath_resources.txt
  uniprot/uniprot_human_reviewed.tsv
  uniprot/uniprot_human_names.tsv
  tf/DatabaseExtract_v_1.01.csv
  depmap/CRISPRGeneEffect.csv
  Human-GEM/Human-GEM.yml   Human-GEM/genes.tsv
```

## 使い方の要点
- 標的遺伝子名: `config.yaml` の一番上の `gene`。
- フォーカス遺伝子: `config.yaml` の `focus.path` のテキスト（1 行に 1 遺伝子、`#` はコメント）。図で丸が青、名前が太字になる。使わないときは `null`。
- 出力: `results/<日付>/<時刻>_<遺伝子名>/`。`results/` は `.gitignore` で GitHub に上げない設定。
- パラメータ（`geo_min`・`max_hops`・`n_mid`・`n_end` など）の意味は、ノートブックの冒頭の説明にある。

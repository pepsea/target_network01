# Installing target_net on a Linux server

Two ways: **A. plain Python venv** (simplest) or **B. Docker** (reproducible, JupyterLab served on a port).
Both share the same data folder layout and `target_net.yaml`.

## 0. Requirements
- Linux x86_64 or arm64, Python 3.10+ (A) or Docker 24+ with compose (B), git, curl
- Disk: ~3 GB for databases (DepMap is 1.5–2 GB), plus a few hundred MB per 100 runs
- RAM: 8 GB is enough for everything except DepMap loading (`CRISPRGeneEffect.csv` → 3–4 GB in pandas); 16 GB recommended if `sl` mode is used
- Network: outbound HTTPS to GitHub and the database sites (only for setup; afterwards the server can run offline with `allow_download: false`)

## 1. Code
```bash
sudo apt-get install -y git curl python3 python3-venv          # Debian/Ubuntu
cd /opt && sudo mkdir -p target_net && sudo chown $USER target_net
git clone https://github.com/pepsea/target_network01.git /opt/target_net
cd /opt/target_net
```

## 2. Databases
Choose one:

**2a. Copy the folder from the Mac** (fastest; keeps exactly the versions you validated):
```bash
# run on the Mac
rsync -avh --progress ~/Documents/database/ user@server:/data/database/
```
**2b. Download on the server** (automatic files) + fetch SIGNOR / Reactome by hand (see `docs/DATA_SETUP.md`):
```bash
mkdir -p /data/database
bash scripts/setup_data.sh /data/database
```
Either way the layout is:
```
/data/database/
  signor/all_data_12_07_26.tsv
  reactome/FIsInGene_061424_with_annotations.txt
  uniprot/uniprot_human_reviewed.tsv
  tf/DatabaseExtract_v_1.01.csv
  dorothea/dorothea_ABC_commercial.tsv
  depmap/CRISPRGeneEffect.csv  OmicsSomaticMutations.csv
  hpa/rna_tissue_consensus_tissues.tsv
  target_net_cache/
```

## 3. Configuration
```bash
cp config.example.yaml target_net.yaml
```
Edit `target_net.yaml` so every path points under `/data/database` (the block in `docs/DATA_SETUP.md`, with `~/Documents/database` → `/data/database`).
Once all files are present set `allow_download: false` for a strictly offline server.

---

## A. Python venv
```bash
cd /opt/target_net
python3 -m venv .venv && source .venv/bin/activate
pip install --upgrade pip && pip install -r requirements.txt
# smoke test (offline demo graph)
python target_net.py upstream --gene AKT1 --demo
# live test against the databases
python target_net.py upstream --gene AKT1 --config target_net.yaml
```
Run JupyterLab for remote use (bind to localhost and tunnel; do not expose it to the network without a password/TLS):
```bash
jupyter lab --no-browser --ip 127.0.0.1 --port 8888 --notebook-dir /opt/target_net
# on your Mac:  ssh -N -L 8888:127.0.0.1:8888 user@server   → open http://localhost:8888
```
Optional systemd unit (`/etc/systemd/system/target_net-jupyter.service`):
```ini
[Unit]
Description=target_net JupyterLab
After=network.target
[Service]
User=youruser
WorkingDirectory=/opt/target_net
Environment=TARGETNET_CONFIG=/opt/target_net/target_net.yaml
ExecStart=/opt/target_net/.venv/bin/jupyter lab --no-browser --ip 127.0.0.1 --port 8888 --notebook-dir /opt/target_net
Restart=on-failure
[Install]
WantedBy=multi-user.target
```
`sudo systemctl enable --now target_net-jupyter`.

## B. Docker
`Dockerfile` and `docker-compose.yml` are in the repo root. The image contains the code and Python dependencies only;
databases, config and run outputs are mounted, so updating code never touches data.
```bash
cd /opt/target_net
# paths: edit docker-compose.yml if your database folder is not /data/database
docker compose build
docker compose run --rm cli upstream --gene AKT1 --demo                     # smoke test
docker compose run --rm cli upstream --gene AKT1                            # live (uses /config/target_net.yaml)
docker compose up -d jupyter                                                # JupyterLab on http://server:8888 (token printed by: docker compose logs jupyter)
```
Inside the containers: code at `/app`, databases read-only at `/data/database`, config at `/config/target_net.yaml`
(`TARGETNET_CONFIG` points to it), outputs at `/app/runs` (bind-mounted to `./runs`).

## 4. Claude Code on the server
```bash
npm install -g @anthropic-ai/claude-code    # needs Node 18+
cd /opt/target_net && claude
```
`CLAUDE.md` and `.claude/skills/target-network-research/` are in the repo, so the project context and the analysis skill load automatically.

## 5. Updating
```bash
cd /opt/target_net && git pull            # code
bash scripts/setup_data.sh /data/database # re-run after deleting a file you want refreshed (it skips existing files)
docker compose build                      # Docker only
```
Record database versions (SIGNOR date, Reactome date, DepMap release) in `/data/database/VERSIONS.txt`; `run_info.json` in each run folder records the config used.

## 6. Troubleshooting
| symptom | cause / fix |
|---|---|
| `[config] signor: no local file and allow_download is false` | path wrong in YAML, or file missing — the key named is the one to fix |
| `X not in the signaling graph` | symbol not in SIGNOR/Reactome (pseudogene? alias?); try the HGNC symbol |
| `[DepMap] too few lost lines -> no SL test` | add `OmicsSomaticMutations.csv` or the expression file; or the core gene is essential (expected, e.g. MTOR) |
| MemoryError loading DepMap | 16 GB RAM, or run `sl` on a machine with more memory; co-essentiality only needs the columns of the candidate genes (optimisation open item) |
| Plotly figures blank in JupyterLab | first cell calls `use_cdn_renderer()`; the browser needs internet to fetch plotly.js from the CDN, or switch to `pio.renderers.default = "notebook"` (embeds 3.5 MB per notebook) |
| STRING step slow/failing | API only; set `allow_download: false` to skip it, or accept the wait |

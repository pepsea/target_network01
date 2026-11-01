#!/usr/bin/env python3
"""
fetch_depmap.py -- download the DepMap files build_network.py needs, without clicking through the portal.

  python fetch_depmap.py --out ~/data/depmap/              # latest public release
  python fetch_depmap.py --out ~/data/depmap/ --release 25Q2
  python fetch_depmap.py --out ~/data/depmap/ --list       # just show what would be downloaded

Route 1: DepMap portal download API (no login needed for public releases).
Route 2: Figshare (DepMap mirrors every public release there) if route 1 fails.
Files (CC-BY 4.0): CRISPRGeneEffect.csv, OmicsSomaticMutations.csv, OmicsExpressionProteinCodingGenesTPMLogp1.csv
"""
import argparse, os, re, sys, json
import requests

WANTED = ["CRISPRGeneEffect.csv", "OmicsSomaticMutations.csv", "OmicsExpressionProteinCodingGenesTPMLogp1.csv", "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"]
HEADERS = {"User-Agent": "target-neighborhood-pipeline/0.1"}
PORTAL_API = "https://depmap.org/portal/download/api/downloads"
FIGSHARE_SEARCH = "https://api.figshare.com/v2/articles/search"


def stream_download(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 1_000_000:
        print(f"  [skip] {os.path.basename(dest)} already present"); return
    tmp = dest + ".part"
    with requests.get(url, headers=HEADERS, stream=True, timeout=600, allow_redirects=True) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0)); done = 0
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk); done += len(chunk)
                if total: print(f"\r  {os.path.basename(dest)}: {done/1e6:,.0f} / {total/1e6:,.0f} MB", end="")
    print(); os.replace(tmp, dest)


# ---------------------------------------------------------------- route 1: DepMap portal
def portal_files(release=None):
    js = requests.get(PORTAL_API, headers=HEADERS, timeout=120).json()
    # the payload is a dict with a "table" list; each row has fileName, releaseName, downloadUrl (format has changed before)
    rows = js.get("table") or js.get("downloads") or js
    if not isinstance(rows, list):
        raise RuntimeError(f"unexpected portal payload keys: {list(js)[:10]}")
    def get(r, *keys):
        for k in keys:
            if k in r and r[k]: return r[k]
        return None
    out = []
    for r in rows:
        fn = get(r, "fileName", "filename", "name"); rel = get(r, "releaseName", "release", "releaseLabel")
        url = get(r, "downloadUrl", "url", "download_url")
        if fn and rel and url: out.append((rel, fn, url))
    rels = sorted({r for r, _, _ in out if re.search(r"DepMap Public \d\dQ\d", r)},
                  key=lambda s: re.search(r"(\d\d)Q(\d)", s).groups(), reverse=True)
    if not rels: raise RuntimeError("no 'DepMap Public YYQN' releases found in portal listing")
    chosen = next((r for r in rels if release and release in r), rels[0])
    files = {fn: url for r, fn, url in out if r == chosen and fn in WANTED}
    return chosen, files


# ---------------------------------------------------------------- route 2: Figshare
def figshare_files(release=None):
    q = f"DepMap Public {release}" if release else "DepMap Public"
    arts = requests.post(FIGSHARE_SEARCH, json={"search_for": q, "page_size": 50, "order": "published_date",
                                                "order_direction": "desc"}, headers=HEADERS, timeout=120).json()
    arts = [a for a in arts if re.search(r"DepMap Public \d\dQ\d", a.get("title", ""))]
    if not arts: raise RuntimeError("no DepMap Public articles found on Figshare")
    arts.sort(key=lambda a: re.search(r"(\d\d)Q(\d)", a["title"]).groups(), reverse=True)
    art = arts[0]
    files = requests.get(f"https://api.figshare.com/v2/articles/{art['id']}/files", headers=HEADERS, timeout=120).json()
    return art["title"], {f["name"]: f["download_url"] for f in files if f["name"] in WANTED}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--release", default=None)
    ap.add_argument("--list", action="store_true"); ap.add_argument("--only", nargs="*", default=None)
    a = ap.parse_args()
    out = os.path.expanduser(a.out); os.makedirs(out, exist_ok=True)
    wanted = a.only or WANTED

    for name, fn in (("DepMap portal", portal_files), ("Figshare", figshare_files)):
        try:
            rel, files = fn(a.release)
            missing = [w for w in wanted if w not in files]
            if missing: print(f"[{name}] release {rel}: missing {missing}")
            if not files: continue
            print(f"[{name}] release: {rel}")
            for w in wanted:
                if w in files: print(f"  {w}  <-  {files[w]}")
            if a.list: return
            for w in wanted:
                if w in files: stream_download(files[w], os.path.join(out, w))
            print(f"\n[done] files in {out}\n  next: python build_network.py --gene MTOR --depmap {out}")
            return
        except Exception as e:
            print(f"[{name}] failed: {e}")
    sys.exit("Both routes failed. Manual fallback: https://depmap.org/portal/data_page -> All Downloads -> "
             + ", ".join(WANTED))


if __name__ == "__main__":
    main()

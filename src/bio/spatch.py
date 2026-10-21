"""Download SPATCH benchmark files (Nat Commun 2025; spatch.pku-genomics.org).

The site serves files through public Aliyun Drive share links. This replicates its download button:
list datasets -> get a share token -> get a signed download URL -> stream to disk.
Usage: uv run python -m bio.spatch list
       uv run python -m bio.spatch get <dataset_id> <field> <outdir>   (field e.g. transcriptome, supplement)
"""
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

LISTING = "https://spatch.pku-genomics.org/spatch/download/spatial/pageInfo"
API = "https://bj21400.api.aliyunfile.com/v2"


def _post(url, data, headers=None, form=False):
    body = urllib.parse.urlencode(data).encode() if form else json.dumps(data).encode()
    h = {"Content-Type": "application/x-www-form-urlencoded" if form else "application/json", **(headers or {})}
    with urllib.request.urlopen(urllib.request.Request(url, body, h), timeout=60) as r:
        return json.load(r)


def list_datasets():
    return _post(LISTING, {"page": 1, "limit": 100}, form=True)["data"]["records"]


def share_files(share_id):
    token = _post(f"{API}/share_link/get_share_token", {"share_id": share_id})["share_token"]
    items = _post(f"{API}/file/list", {"limit": 100, "share_id": share_id, "parent_file_id": "root",
                                       "url_expire_sec": 7200}, {"x-share-token": token})["items"]
    return items


def download(share_id, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for it in share_files(share_id):
        dest = outdir / it["name"]
        if dest.exists() and dest.stat().st_size == it.get("size"):
            print("have", dest)
            continue
        print(f"get {it['name']} ({it.get('size', 0) / 1e6:.0f} MB)")
        req = urllib.request.Request(it["download_url"], headers={"Referer": "https://www.aliyundrive.com/"})
        with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
            while chunk := r.read(1 << 20):
                f.write(chunk)


if __name__ == "__main__":
    if sys.argv[1] == "list":
        for x in list_datasets():
            print(x["id"], x["platform"], x["cancer"], x["cellOrSpotNumber"], x.get("transcriptomeSize"))
    elif sys.argv[1] == "get":
        rec = {str(x["id"]): x for x in list_datasets()}[sys.argv[2]]
        download(rec[sys.argv[3]], sys.argv[4])

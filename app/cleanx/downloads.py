import datetime as _dt
import os
import time
import urllib.parse
import urllib.request
from pathlib import Path

from . import config
from .labels import load_json, save_json

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
ENZYME_DAT_URL = "https://ftp.expasy.org/databases/enzyme/enzyme.dat"
PFAM_CLANS_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"
PFAM_RELNOTES_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt"

NEW_ENTRY_FIELDS = ["accession", "ec", "date_created", "length", "fragment", "keywordid",
                    "organism_name", "protein_existence", "xref_pfam", "sequence"]
ALL_EC_FIELDS = ["accession", "ec", "date_created", "xref_pfam", "length"]


def _fetch(url, dest, timeout=1800, retries=3):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cleanx/0.1 (research)"})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as fh:
                headers = dict(r.headers.items())
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            os.replace(tmp, dest)
            return headers
        except Exception as e:  # network errors are retried, then re-raised
            last = e
            print(f"  download attempt {attempt}/{retries} failed: {e}", flush=True)
            time.sleep(10 * attempt)
    raise RuntimeError(f"download failed: {url}") from last


def _manifest(dest, url, headers, extra=None):
    info = {"url": url, "downloaded_utc": _dt.datetime.utcnow().isoformat(timespec="seconds"),
            "uniprot_release": headers.get("X-UniProt-Release") or headers.get("x-uniprot-release"),
            "uniprot_release_date": headers.get("X-UniProt-Release-Date") or headers.get("x-uniprot-release-date"),
            "bytes": Path(dest).stat().st_size}
    info.update(extra or {})
    save_json(info, str(dest) + ".json")
    return info


def uniprot_table(query, fields, dest, max_records=None, force=False):
    """Download a UniProtKB TSV (gzip) and write <dest>.json with the release it came from.

    max_records uses the paginated search endpoint (first page only) for smoke tests."""
    dest = Path(dest)
    if dest.is_file() and not force:
        return load_json(str(dest) + ".json") if Path(str(dest) + ".json").is_file() else {}
    params = {"query": query, "fields": ",".join(fields), "format": "tsv", "compressed": "true"}
    if max_records:
        params["size"] = str(int(max_records))
        url = UNIPROT_SEARCH + "?" + urllib.parse.urlencode(params)
    else:
        url = UNIPROT_STREAM + "?" + urllib.parse.urlencode(params)
    print(f"  downloading UniProt: {query}  ->  {dest.name}", flush=True)
    t0 = time.time()
    headers = _fetch(url, dest)
    info = _manifest(dest, url, headers, {"query": query, "fields": fields, "max_records": max_records})
    print(f"  done in {time.time() - t0:.0f}s (release {info['uniprot_release']})", flush=True)
    return info


def plain_file(url, dest, force=False):
    dest = Path(dest)
    if dest.is_file() and not force:
        return load_json(str(dest) + ".json") if Path(str(dest) + ".json").is_file() else {}
    print(f"  downloading {url}", flush=True)
    headers = _fetch(url, dest)
    return _manifest(dest, url, headers)


def read_uniprot_tsv(path, fields):
    import pandas as pd
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, compression="infer")
    if len(df.columns) != len(fields):
        raise ValueError(f"{path}: expected {len(fields)} columns {fields}, got {list(df.columns)}")
    df.columns = fields
    return df

"""Download and cache UniProt EC evidence for E7.

Use --offline to parse cached responses without network access.

See docs/code_reference.md#cx-12-uniprot-evidence."""
import argparse
import datetime as _dt
import gzip
import hashlib
import json
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from cleanx import config

config.enter_app_dir()

import pandas as pd

from cleanx import ecbench as ECB
from cleanx import labels as L
from cleanx import provenance as PV
from cleanx.downloads import UNIPROT_STREAM
from cleanx.pipeline import write_rows

JSON_DIR = config.CX_DOWNLOADS / "uniprot_json"
PROV_DIR = config.CX_DATA / "provenance"
MANIFEST = PROV_DIR / "manifest.json"
SINGLE = "https://rest.uniprot.org/uniprotkb/{acc}?format=json"
SETS = ["new392", "sp_temporal_cal", "sp_temporal_test", "ecbench_test"]
BATCH = 100
fails = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""), flush=True)
    if not ok:
        fails.append(name)
    return ok


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def accessions(name):
    if name == "ecbench_test":
        path = ECB.DEST / ECB.FILES[("test", "100")]
        if not path.is_file():
            ECB.fetch()
        ECB.verify({("test", "100"): path})
        ids, _, _ = ECB.read_results(path)
        return ids
    ids, _, _ = L.load_named_dataset(name)
    return ids


def http_get(url, timeout=300, retries=3):
    """(status, headers, body bytes); HTTP 4xx is returned, not raised; network errors are retried."""
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cleanx/0.1 (research)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.status, dict(r.headers.items()), r.read()
        except urllib.error.HTTPError as e:
            if 400 <= e.code < 500:
                return e.code, dict(e.headers.items()) if e.headers else {}, b""
            last = e
        except Exception as e:  # noqa: BLE001 - network errors are retried, then re-raised
            last = e
        print(f"  attempt {attempt}/{retries} failed: {last}", flush=True)
        time.sleep(10 * attempt)
    raise RuntimeError(f"download failed: {url}") from last


def release_of(headers):
    low = {str(k).lower(): v for k, v in headers.items()}
    return low.get("x-uniprot-release") or ""


def batches(ids, d):
    """[(accessions, cache path)] in order; the file name carries the SHA-256 of the accession list."""
    out = []
    for b in range(0, len(ids), BATCH):
        chunk = ids[b:b + BATCH]
        h = hashlib.sha256("\n".join(chunk).encode("utf-8")).hexdigest()[:16]
        out.append((chunk, d / f"batch_{b // BATCH:03d}_{h}.json.gz"))
    return out


def fetch_single(acc, d, record):
    """Cache one entry from the single-entry endpoint; returns True if an entry was stored."""
    path, nf = d / f"single_{acc}.json.gz", d / f"single_{acc}.notfound"
    if path.is_file():
        return True
    if nf.is_file():
        return False
    url = SINGLE.format(acc=urllib.parse.quote(acc))
    status, headers, body = http_get(url)
    time.sleep(0.5)
    if status == 200 and body:
        try:
            json.loads(body.decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            raise SystemExit(f"single entry {acc}: response is not JSON ({e})")
        path.write_bytes(gzip.compress(body, mtime=0))
        record[path.relative_to(config.CX_DATA).as_posix()] = {"url": url, "release": release_of(headers),
                                                               "downloaded_utc": now()}
        return True
    nf.write_text(f"HTTP {status}\n", encoding="utf-8")
    record[nf.relative_to(config.CX_DATA).as_posix()] = {"url": url, "http_status": status, "downloaded_utc": now()}
    return False


def download_set(ids, d, record):
    """Fill the cache for one set: the batches, then single entries for queries the batches do not resolve, then
    the target of every inactive entry that was merged into exactly one other entry."""
    d.mkdir(parents=True, exist_ok=True)
    for b, (chunk, path) in enumerate(batches(ids, d)):
        if path.is_file():
            continue
        q = " OR ".join(f"accession:{a}" for a in chunk)
        url = UNIPROT_STREAM + "?" + urllib.parse.urlencode({"query": q, "format": "json", "compressed": "true"})
        status, headers, body = http_get(url)
        if status != 200:
            raise SystemExit(f"{d.name}: UniProt returned HTTP {status} for batch {b}")
        try:
            json.loads(gzip.decompress(body).decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            raise SystemExit(f"{d.name}: batch {b} is not gzip JSON ({e})")
        path.write_bytes(body)
        record[path.relative_to(config.CX_DATA).as_posix()] = {"url": url, "release": release_of(headers),
                                                               "n_queries": len(chunk), "downloaded_utc": now()}
        time.sleep(1)
    batch_idx = PV.index_entries([e for _, p in batches(ids, d) for e in PV.read_json_results(p)])
    for a in ids:
        if PV.resolve(a, batch_idx)[0] is None and fetch_single(a, d, record):
            for e in PV.read_json_results(d / f"single_{a}.json.gz"):
                typ, targets = PV.inactive_reason(e) if PV.is_inactive(e) else ("", [])
                if typ == "MERGED" and len(targets) == 1:
                    fetch_single(targets[0], d, record)


def now():
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_cache(ids, d):
    """(entries, via) from the cache of one set: the expected batch files and every single-entry file."""
    missing = [p.name for _, p in batches(ids, d) if not p.is_file()]
    if missing:
        raise SystemExit(f"{d}: {len(missing)} batch responses are not cached (e.g. {missing[0]}); run without "
                         f"--offline on a machine with internet access")
    batch_entries = [e for _, p in batches(ids, d) for e in PV.read_json_results(p)]
    batch_idx = PV.index_entries(batch_entries)
    singles = sorted(d.glob("single_*.json.gz"))
    entries = batch_entries + [e for p in singles for e in PV.read_json_results(p)]
    via = {}
    for a in ids:
        if PV.resolve(a, batch_idx)[0] is not None:
            via[a] = "batch"
        elif (d / f"single_{a}.json.gz").is_file():
            via[a] = "single"
        elif (d / f"single_{a}.notfound").is_file():
            via[a] = "not_found"
        else:
            via[a] = "not_downloaded"
    return entries, via


def parse_set(ids, d, tag):
    entries, via = read_cache(ids, d)
    erows, crows = PV.tables(ids, entries, via)
    ed = write_rows(erows, PROV_DIR / f"{tag}_entries.csv")
    write_rows(crows, PROV_DIR / f"{tag}_ec_evidence.csv")
    return ed, len(crows)


def smoke():
    """Known answers for the parser and the class rules on synthetic entries (no network)."""
    E, S, O = ("ECO:0000269", "PubMed", "1"), ("ECO:0000250", "UniProtKB", "P11111"), ("ECO:0000305", "", "")
    R = ("ECO:0000255", "HAMAP-Rule", "MF_00001")
    entries = [
        PV.synthetic_entry("A00001", {"1.1.1.1": [E]}),                                   # EXP in the name
        PV.synthetic_entry("A00002", {"2.2.2.2": [S]}, secondary=["B00002"]),              # SIM, found by secondary
        PV.synthetic_entry("A00003", {}, in_reaction={"3.3.3.3": [S, E]}),                  # reaction only; EXP wins
        PV.synthetic_entry("A00004", {"4.4.4.4": [O]}, contains={"4.4.4.5": [S]}),          # OTHER + SIM in a domain
        PV.synthetic_entry("A00005", {"5.5.5.5": [E], "5.5.5.6": [S]}),                     # multi-EC, one SIM
        PV.synthetic_entry("A00006", {}, inactive="DELETED"),                               # deleted entry
        PV.synthetic_entry("M00008", {}, inactive="MERGED", merged_into=["A00001"]),        # merged into A00001
        PV.synthetic_entry("D00009", {}, inactive="DEMERGED", merged_into=["A00001", "A00005"]),
        PV.synthetic_entry("A00010", {"2.1.3.2": [R], "2.1.3.3": [R, O]}),                  # family rule (RULE)
        PV.synthetic_entry("A00011", {"1.1.1.1": [E]}, secondary=["S00012"]),               # S00012 is claimed by
        PV.synthetic_entry("A00013", {"1.1.1.1": [E]}, secondary=["S00012"]),               # two entries: ambiguous
    ]
    variant = PV.synthetic_entry("A00007", {})  # serializer variants: bare-string EC, object-valued reaction EC
    variant["proteinDescription"]["alternativeNames"] = [{"fullName": {"value": "alt"}, "ecNumbers": ["7.7.7.7"]}]
    variant["comments"].append({"commentType": "CATALYTIC ACTIVITY", "reaction": {
        "name": "C = D", "ecNumber": {"value": "7.7.7.8"}, "evidences": [
            {"evidenceCode": "ECO:0000250", "source": "UniProtKB", "id": "Q22222"}]}})
    entries.append(variant)
    queries = ["A00001", "B00002", "A00003", "A00004", "A00005", "A00006", "Z99999"]
    with tempfile.TemporaryDirectory() as td:  # round trip through a gzip JSON response file
        p = Path(td) / "batch_000.json.gz"
        p.write_bytes(gzip.compress(json.dumps({"results": entries}).encode("utf-8"), mtime=0))
        entries = PV.read_json_results(p)
    erows, crows = PV.tables(queries, entries)
    found = {r["query_accession"]: r["found"] for r in erows}
    check("entries found by primary and secondary accession; inactive and unknown ones are not found",
          found == {"A00001": True, "B00002": True, "A00003": True, "A00004": True, "A00005": True,
                    "A00006": False, "Z99999": False}, str(found))
    codes = {(r["query_accession"], r["ec"]): set(r["codes"].split(";")) - {""} for r in crows}
    check("EC evidence is collected from names, reactions and domain names",
          codes == {("A00001", "1.1.1.1"): {"ECO:0000269"}, ("B00002", "2.2.2.2"): {"ECO:0000250"},
                    ("A00003", "3.3.3.3"): {"ECO:0000250", "ECO:0000269"}, ("A00004", "4.4.4.4"): {"ECO:0000305"},
                    ("A00004", "4.4.4.5"): {"ECO:0000250"}, ("A00005", "5.5.5.5"): {"ECO:0000269"},
                    ("A00005", "5.5.5.6"): {"ECO:0000250"}}, str(codes))
    src = {(r["query_accession"], r["ec"]): r["sim_sources"] for r in crows}
    check("similarity sources are the cited UniProtKB accessions", src[("B00002", "2.2.2.2")] == "P11111" and
          src[("A00001", "1.1.1.1")] == "")
    where = {(r["query_accession"], r["ec"]): (r["in_name"], r["in_reaction"]) for r in crows}
    check("where each EC was found is recorded", where[("A00003", "3.3.3.3")] == (False, True) and
          where[("A00001", "1.1.1.1")] == (True, False))
    per = {}
    for r in crows:
        per.setdefault(r["query_accession"], {})[r["ec"]] = r["codes"].split(";")
    cls = {q: PV.protein_class(t, per.get(q, {}), found[q]) for q, t in
           [("A00001", ["1.1.1.1"]), ("B00002", ["2.2.2.2"]), ("A00003", ["3.3.3.3"]), ("A00004", ["4.4.4.4"]),
            ("A00005", ["5.5.5.5", "5.5.5.6"]), ("A00006", ["6.6.6.6"]), ("Z99999", ["1.1.1.1"])]}
    check("protein classes follow the pre-registered rules",
          cls == {"A00001": "EXP", "B00002": "SIM", "A00003": "EXP", "A00004": "OTHER", "A00005": "SIM",
                  "A00006": "LABEL_ABSENT", "Z99999": "LABEL_ABSENT"}, str(cls))
    check("a benchmark label missing from the entry makes the protein LABEL_ABSENT",
          PV.protein_class(["1.1.1.1", "9.9.9.9"], per["A00001"]) == "LABEL_ABSENT")
    meta = {r["query_accession"]: r for r in erows}
    check("lineage and superkingdom are parsed", meta["A00001"]["superkingdom"] == "Bacteria")
    _, vrows = PV.tables(["A00007"], entries)
    vcodes = {r["ec"]: (r["codes"], r["sim_sources"]) for r in vrows}
    check("serializer variants (bare-string EC, object-valued reaction EC) are parsed",
          vcodes == {"7.7.7.7": ("", ""), "7.7.7.8": ("ECO:0000250", "Q22222")}, str(vcodes))
    q2 = ["M00008", "D00009", "S00012", "A00006", "A00010", "Z99999", "B00002", "A00001"]
    er2, cr2 = PV.tables(q2, entries)
    res = {r["query_accession"]: (r["found"], r["resolution"], r.get("primary_accession", "")) for r in er2}
    check("accession resolution: merged followed, demerged/deleted/ambiguous/unknown not found",
          res == {"M00008": (True, "merged", "A00001"), "D00009": (False, "inactive:DEMERGED", "D00009"),
                  "S00012": (False, "ambiguous_secondary", ""), "A00006": (False, "inactive:DELETED", "A00006"),
                  "A00010": (True, "primary", "A00010"), "Z99999": (False, "absent", ""),
                  "B00002": (True, "secondary", "A00002"), "A00001": (True, "primary", "A00001")}, str(res))
    per2 = {}
    for r in cr2:
        per2.setdefault(r["query_accession"], {})[r["ec"]] = r["codes"].split(";")
    msrc = {r["ec"]: r["model_sources"] for r in cr2 if r["query_accession"] == "A00010"}
    check("sequence-model (family rule) labels are RULE; experimental beats every other code",
          PV.protein_class(["2.1.3.2", "2.1.3.3"], per2["A00010"]) == "RULE" and
          PV.ec_class(["ECO:0000255", "ECO:0000269"]) == "EXP" and PV.ec_class(["ECO:0000255", "ECO:0000250"]) ==
          "SIM" and PV.protein_class(["1.1.1.1"], per2["M00008"]) == "EXP" and
          msrc == {"2.1.3.2": "HAMAP-Rule:MF_00001", "2.1.3.3": "HAMAP-Rule:MF_00001"}, str(msrc))
    with tempfile.TemporaryDirectory() as td:  # cache names depend on the accession list (no stale reuse)
        a = batches(["P1", "P2", "P3"], Path(td))
        b = batches(["P1", "P2"], Path(td))
        c = batches(["P1", "P2", "P3"], Path(td))
        check("batch cache files are keyed by their accession list",
              a[0][1] != b[0][1] and a[0][1] == c[0][1] and a[0][0] == ["P1", "P2", "P3"])
    check("release header is read case-insensitively",
          release_of({"x-uniprot-release": "2026_04"}) == "2026_04" and
          release_of({"X-UniProt-Release": "2026_04"}) == "2026_04" and release_of({}) == "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", nargs="+", default=SETS)
    ap.add_argument("--offline", action="store_true", help="parse cached responses only (no network)")
    ap.add_argument("--limit", type=int, default=0, help="first N accessions per set (live check)")
    ap.add_argument("--smoke", action="store_true", help="parser known-answer checks on synthetic entries")
    args = ap.parse_args()
    config.ensure_dirs()
    t0 = time.time()
    if args.smoke:
        smoke()
        print(f"\nCX12 SMOKE {'FAILED' if fails else 'DONE'} in {time.time() - t0:.0f}s")
        sys.exit(1 if fails else 0)

    PROV_DIR.mkdir(parents=True, exist_ok=True)
    manifest = L.load_json(MANIFEST) if MANIFEST.is_file() else {"files": {}, "sets": {}}
    for name in args.sets:
        ids = accessions(name)
        if args.limit:
            ids = ids[:args.limit]
        tag = name if not args.limit else f"{name}_limit{args.limit}"
        d = JSON_DIR / tag
        print(f"{name}: {len(ids)} accessions", flush=True)
        if not args.offline:
            download_set(ids, d, manifest["files"])
        ed, n_ec = parse_set(ids, d, tag)
        found = int(ed.found.sum())
        prefix = d.relative_to(config.CX_DATA).as_posix() + "/"
        releases = sorted({v.get("release", "") for k, v in manifest["files"].items()
                           if k.startswith(prefix) and v.get("release")})
        manifest["sets"][tag] = {"n_queries": len(ids), "n_found": found, "n_not_found": len(ids) - found,
                                 "n_ec_rows": n_ec, "uniprot_releases": releases,
                                 "via": ed.via.value_counts().to_dict(),
                                 "resolution": ed.resolution.value_counts().to_dict(), "parsed_utc": now()}
        print(f"  found {found}/{len(ids)} entries ({len(ids) - found} not found), {n_ec} EC rows; "
              f"resolution {manifest['sets'][tag]['resolution']}; "
              f"UniProt release(s) {releases or 'unknown (offline cache)'}", flush=True)
    for rel, info in manifest["files"].items():
        p = config.CX_DATA / rel
        if p.is_file():
            info["sha256"] = sha256(p)
    L.save_json(manifest, MANIFEST)
    rels = sorted({r for s in manifest["sets"].values() for r in s.get("uniprot_releases", [])})
    if len(rels) > 1:
        print(f"  WARNING: responses come from more than one UniProt release: {rels}; delete "
              f"{JSON_DIR} and download again for a single release", flush=True)
    print(f"\nCX12 DONE in {time.time() - t0:.0f}s -> {PROV_DIR}")


if __name__ == "__main__":
    main()

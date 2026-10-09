"""Label provenance (E7): how was each benchmark EC label on a UniProtKB entry asserted?

UniProtKB attaches Evidence and Conclusion Ontology (ECO) codes to every annotation. For an EC number the evidence
is the union of the codes on that EC in the protein description (recommended, alternative and submission names,
including 'includes'/'contains' domains) and on CATALYTIC ACTIVITY reactions that carry that EC:
  ECO:0000269  experimental evidence used in manual assertion
  ECO:0000250  sequence similarity evidence used in manual assertion ('by similarity'); its source names the
               protein the annotation was propagated from (source 'UniProtKB', id = accession)
  ECO:0000255 / ECO:0000256 / ECO:0000259  match to a sequence model (e.g. a HAMAP rule) used in manual /
               automatic assertion
Everything else (e.g. ECO:0000305 curator inference, ECO:0000303 author statement, imported annotations, no code)
is OTHER.

Classes (fixed in the E7 pre-registration, cx_13_provenance.py):
  EC level      EXP if any ECO:0000269; else SIM if any ECO:0000250; else RULE if any sequence-model code;
                else OTHER
  protein level LABEL_ABSENT if the entry was not found or any benchmark label is not an EC on the current entry;
                else SIM if any label is SIM; else EXP if every label is EXP; else RULE if any label is RULE;
                else OTHER. SIM and RULE together are the labels propagated by homology (PROPAGATED).

Accession resolution: an active entry whose primary accession is the query; else the single active entry that
lists the query as a secondary accession (more than one such entry: 'ambiguous_secondary', not found); else an
inactive entry for the query that was merged into exactly one active entry present in the responses ('merged');
otherwise not found ('inactive:<reason>' or 'absent').

This module parses the UniProt REST JSON (uniprotkb entries) into flat tables; it never downloads anything."""
import gzip
import json
from pathlib import Path

EXPERIMENTAL = frozenset({"ECO:0000269"})
SIMILARITY = frozenset({"ECO:0000250"})
SEQUENCE_MODEL = frozenset({"ECO:0000255", "ECO:0000256", "ECO:0000259"})
EC_CLASSES = ("EXP", "SIM", "RULE", "OTHER")
PROTEIN_CLASSES = ("EXP", "SIM", "RULE", "OTHER", "LABEL_ABSENT")
PROPAGATED = frozenset({"SIM", "RULE"})


def _evidences(obj):
    out = []
    for ev in obj.get("evidences") or []:
        code = str(ev.get("evidenceCode") or "").strip()
        if code:
            out.append((code, str(ev.get("source") or ""), str(ev.get("id") or "")))
    return out


def _walk_names(node, found):
    """Collect (ec, evidences) from every 'ecNumbers' list anywhere under the protein description."""
    if isinstance(node, dict):
        for k, v in node.items():
            if k == "ecNumbers" and isinstance(v, list):
                for e in v:  # uniprot-core serialises an EC as {"value", "evidences"}, a bare ECNumber as a string
                    if isinstance(e, dict) and e.get("value"):
                        found.append((str(e["value"]).strip(), _evidences(e)))
                    elif isinstance(e, str) and e.strip():
                        found.append((e.strip(), []))
            else:
                _walk_names(v, found)
    elif isinstance(node, list):
        for v in node:
            _walk_names(v, found)


def ec_evidence(entry):
    """entry (UniProt JSON) -> {ec: {'codes', 'sim_sources', 'model_sources': sets; 'in_name', 'in_reaction'}}."""
    out = {}

    def add(ec, evs, where):
        rec = out.setdefault(ec, {"codes": set(), "sim_sources": set(), "model_sources": set(),
                                  "in_name": False, "in_reaction": False})
        rec[where] = True
        for code, source, ident in evs:
            rec["codes"].add(code)
            if code in SIMILARITY and source == "UniProtKB" and ident:
                rec["sim_sources"].add(ident)
            if code in SEQUENCE_MODEL and source:
                rec["model_sources"].add(f"{source}:{ident}" if ident else source)

    names = []
    _walk_names(entry.get("proteinDescription") or {}, names)
    for ec, evs in names:
        add(ec, evs, "in_name")
    for c in entry.get("comments") or []:
        if c.get("commentType") != "CATALYTIC ACTIVITY":
            continue
        r = c.get("reaction") or {}
        ec = r.get("ecNumber") or ""
        if isinstance(ec, dict):  # defensive: an evidenced-value object instead of the usual plain string
            ec = ec.get("value") or ""
        ec = str(ec).strip()
        if ec:
            add(ec, _evidences(r), "in_reaction")
    return out


def is_inactive(entry):
    return str(entry.get("entryType") or "").lower().startswith("inactive")


def inactive_reason(entry):
    r = entry.get("inactiveReason") or {}
    return (str(r.get("inactiveReasonType") or "UNKNOWN").upper(),
            [str(t) for t in r.get("mergeDemergeTo") or []])


def entry_meta(entry):
    org = entry.get("organism") or {}
    lineage = [str(x) for x in (org.get("lineage") or [])]
    audit = entry.get("entryAudit") or {}
    return {"primary_accession": str(entry.get("primaryAccession") or ""),
            "secondary_accessions": ";".join(str(x) for x in entry.get("secondaryAccessions") or []),
            "entry_type": str(entry.get("entryType") or ""),
            "superkingdom": lineage[0] if lineage else "",
            "lineage": ";".join(lineage),
            "protein_existence": str(entry.get("proteinExistence") or ""),
            "first_public_date": str(audit.get("firstPublicDate") or ""),
            "entry_version": audit.get("entryVersion", ""),
            "annotation_score": entry.get("annotationScore", ""),
            "inactive_reason": json.dumps(entry.get("inactiveReason")) if entry.get("inactiveReason") else ""}


def index_entries(entries):
    """{'primary': acc -> active entry, 'secondary': acc -> [active entries listing it], 'inactive': acc -> entry}.
    The same entry read twice (e.g. from two responses) is kept once."""
    prim, sec, inact = {}, {}, {}
    for e in entries:
        acc = str(e.get("primaryAccession") or "")
        if not acc:
            continue
        if is_inactive(e):
            inact.setdefault(acc, e)
            continue
        prim[acc] = e
    for acc, e in prim.items():
        for s in e.get("secondaryAccessions") or []:
            lst = sec.setdefault(str(s), [])
            if acc not in [str(x.get("primaryAccession")) for x in lst]:
                lst.append(e)
    return {"primary": prim, "secondary": sec, "inactive": inact}


def resolve(q, index):
    """(entry or None, how) for one query accession; see the module docstring for the order of the rules."""
    if q in index["primary"]:
        return index["primary"][q], "primary"
    cands = index["secondary"].get(q, [])
    if len(cands) == 1:
        return cands[0], "secondary"
    if len(cands) > 1:
        return None, "ambiguous_secondary"
    e = index["inactive"].get(q)
    if e is not None:
        typ, targets = inactive_reason(e)
        if typ == "MERGED" and len(targets) == 1 and targets[0] in index["primary"]:
            return index["primary"][targets[0]], "merged"
        return None, f"inactive:{typ}"
    return None, "absent"


def tables(query_ids, entries, via=None):
    """Flat tables for a list of query accessions.

    Returns (entry_rows, ec_rows): one entry row per query (found or not) and one row per (query, EC on the
    resolved entry). via: optional query -> how the response was obtained (e.g. 'batch', 'single')."""
    idx = index_entries(entries)
    erows, crows = [], []
    for q in query_ids:
        e, how = resolve(q, idx)
        base = {"query_accession": q, "found": e is not None, "resolution": how, "via": (via or {}).get(q, "")}
        if e is None:
            gone = idx["inactive"].get(q)
            erows.append({**base, **(entry_meta(gone) if gone is not None else {})})
            continue
        erows.append({**base, **entry_meta(e)})
        for ec, rec in sorted(ec_evidence(e).items()):
            crows.append({"query_accession": q, "primary_accession": str(e.get("primaryAccession") or ""),
                          "ec": ec, "codes": ";".join(sorted(rec["codes"])),
                          "sim_sources": ";".join(sorted(rec["sim_sources"])),
                          "model_sources": ";".join(sorted(rec["model_sources"])),
                          "in_name": rec["in_name"], "in_reaction": rec["in_reaction"]})
    return erows, crows


def ec_class(codes):
    codes = set(codes)
    if codes & EXPERIMENTAL:
        return "EXP"
    if codes & SIMILARITY:
        return "SIM"
    if codes & SEQUENCE_MODEL:
        return "RULE"
    return "OTHER"


def protein_class(true_ecs, ec_codes, found=True):
    """true_ecs: the benchmark labels of one protein; ec_codes: {ec: iterable of codes} on its current entry."""
    if not found or not true_ecs or any(e not in ec_codes for e in true_ecs):
        return "LABEL_ABSENT"
    cls = [ec_class(ec_codes[e]) for e in true_ecs]
    if "SIM" in cls:
        return "SIM"
    if all(c == "EXP" for c in cls):
        return "EXP"
    if "RULE" in cls:
        return "RULE"
    return "OTHER"


def source_accession(ident):
    """Canonical accession of a cited source ('P12345-2' isoform -> 'P12345')."""
    return str(ident).strip().split("-")[0]


def read_json_results(path):
    """A saved UniProt JSON response ({'results': [...]} or a single entry), gzip or plain."""
    raw = Path(path).read_bytes()
    if str(path).endswith(".gz"):
        raw = gzip.decompress(raw)
    obj = json.loads(raw.decode("utf-8"))
    if isinstance(obj, dict) and "results" in obj:
        return list(obj["results"])
    return [obj]


def synthetic_entry(acc, ecs, secondary=(), lineage=("Bacteria",), in_reaction=None, contains=None,
                    inactive=None, merged_into=()):
    """UniProt-JSON-shaped entry for tests. ecs: {ec: [(code, source, id), ...]} placed on the recommended name;
    in_reaction: {ec: [(code, source, id)]} placed on CATALYTIC ACTIVITY comments; contains: {ec: [...]} placed
    on a 'contains' domain name; inactive: None or a reason type ('DELETED', 'MERGED', 'DEMERGED') with
    merged_into the target accessions."""
    def evs(lst):
        return [{"evidenceCode": c, "source": s, "id": i} for c, s, i in lst]

    if inactive:
        return {"entryType": "Inactive", "primaryAccession": acc,
                "inactiveReason": {"inactiveReasonType": inactive, "mergeDemergeTo": list(merged_into)}}
    entry = {"entryType": "UniProtKB reviewed (Swiss-Prot)",
             "primaryAccession": acc, "secondaryAccessions": list(secondary),
             "organism": {"scientificName": "Synthetica", "lineage": list(lineage)},
             "proteinExistence": "1: Evidence at protein level",
             "entryAudit": {"firstPublicDate": "2023-01-01", "entryVersion": 3}, "annotationScore": 3.0,
             "proteinDescription": {"recommendedName": {
                 "fullName": {"value": "synthetic"},
                 "ecNumbers": [{"value": ec, "evidences": evs(ev)} for ec, ev in ecs.items()]}},
             "comments": [{"commentType": "FUNCTION", "texts": [{"value": "x"}]}]}
    if contains:
        entry["proteinDescription"]["contains"] = [{"recommendedName": {
            "fullName": {"value": "domain"},
            "ecNumbers": [{"value": ec, "evidences": evs(ev)} for ec, ev in contains.items()]}}]
    for ec, ev in (in_reaction or {}).items():
        entry["comments"].append({"commentType": "CATALYTIC ACTIVITY",
                                  "reaction": {"name": "A = B", "ecNumber": ec, "evidences": evs(ev)}})
    return entry

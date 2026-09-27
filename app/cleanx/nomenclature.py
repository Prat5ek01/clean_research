"""EC nomenclature drift: map EC numbers that IUBMB transferred or deleted to their current successors
(ENZYME database, enzyme.dat)."""
import re

EC_RE = re.compile(r"\d+\.\d+\.\d+\.n?\d+")


def parse_enzyme_dat(path):
    release, active, transfers, deleted = None, set(), {}, set()
    cur, de = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("CC   Release of"):
                release = line[5:].strip()
            elif line.startswith("ID   "):
                cur, de = line[5:].strip(), []
            elif line.startswith("DE   ") and cur is not None:
                de.append(line[5:].strip())
            elif line.startswith("//"):
                if cur is not None:
                    text = " ".join(de)
                    if text.startswith("Transferred entry"):
                        transfers[cur] = EC_RE.findall(text)
                    elif text.startswith("Deleted entry"):
                        deleted.add(cur)
                    else:
                        active.add(cur)
                cur, de = None, []
    return {"release": release, "active": active, "transfers": transfers, "deleted": deleted}


def forward_map(nom):
    """old EC -> sorted list of current ECs (empty list = deleted). Unchanged ECs are absent."""
    transfers, deleted = nom["transfers"], nom["deleted"]

    def resolve(ec, seen):
        if ec in seen:
            return set()
        seen = seen | {ec}
        if ec in transfers:
            out = set()
            for t in transfers[ec]:
                out |= resolve(t, seen)
            return out
        if ec in deleted:
            return set()
        return {ec}

    return {ec: sorted(resolve(ec, frozenset())) for ec in list(transfers) + sorted(deleted)}


def forward_counts(ec_counts, fmap):
    """Training EC counts re-expressed in current nomenclature (a transferred EC credits all successors)."""
    out = {}
    for ec, n in ec_counts.items():
        for cur in fmap.get(ec, [ec]):
            out[cur] = out.get(cur, 0) + n
    return out

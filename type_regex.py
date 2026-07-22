import re

SV_MIN_BP = 50
CNV_MIN_BP = 1000

def clean_number(s):
    return s.replace(".", "").replace(",", "")

CHROM_RE = re.compile(r"chr([\dXYM]+)", re.IGNORECASE)

def get_chrom(s):
    m = CHROM_RE.search(s)
    return "chr" + m.group(1) if m else None

def get_positions(s):
    t = re.sub(r"chr[\dXYM]+", "", s, flags=re.IGNORECASE)
    t = re.sub(r"[ACGT]+\s*>\s*[ACGT]+", "", t, flags=re.IGNORECASE)
    t = re.sub(r"[\d.]+\s*kb", "", t, flags=re.IGNORECASE)
    nums = re.findall(r"\d[\d.,]{3,}\d|\d{4,}", t)
    return [clean_number(n) for n in nums]

def estimate_size(s, positions, indel_match):
    kb = re.search(r"([\d.]+)\s*kb", s, re.IGNORECASE)
    if kb:
        try:
            return float(kb.group(1)) * 1000
        except ValueError:
            pass
    if len(positions) >= 2:
        return abs(int(positions[1]) - int(positions[0]))
    if indel_match:
        parts = re.split(r"\s*>\s*", indel_match.group(0))
        if len(parts) == 2:
            return abs(len(parts[0]) - len(parts[1]))
    return None

SNV_RE = re.compile(r"([ACGT])\s*>\s*([ACGT])(?![ACGT])", re.IGNORECASE)

def match_snv(s):
    m = SNV_RE.search(s)
    if not m:
        return None
    return {"vtype": "SNV", "chrom": get_chrom(s),
            "pos": (get_positions(s) or [None])[0],
            "ref": m.group(1).upper(), "alt": m.group(2).upper()}

SV_KEYWORD_RE = re.compile(r"\b(del|dup|inv)\b|del$|dup$", re.IGNORECASE)
SV_INDEL_RE = re.compile(r"([ACGT]{2,}\s*>\s*[ACGT]+|[ACGT]+\s*>\s*[ACGT]{2,})", re.IGNORECASE)

def match_sv(s):
    positions = get_positions(s)
    kw = SV_KEYWORD_RE.search(s)
    indel = SV_INDEL_RE.search(s)
    if not kw and not indel:
        return None
    size = estimate_size(s, positions, indel)
    if size is not None and size >= CNV_MIN_BP:
        return None
    vtype = "SV"
    if kw:
        vtype = "SV/" + kw.group(0).lower().replace("$", "")
    return {"vtype": vtype, "chrom": get_chrom(s),
            "start": positions[0] if positions else None,
            "end": positions[1] if len(positions) > 1 else None, "size": size}

def match_cnv(s):
    positions = get_positions(s)
    size = estimate_size(s, positions, None)
    has_delsdup = re.search(r"\b(del|dup)\b|del$|dup$", s, re.IGNORECASE)
    is_range = len(positions) >= 2
    if size is not None and size >= CNV_MIN_BP and has_delsdup:
        return {"vtype": "CNV", "chrom": get_chrom(s),
                "start": positions[0] if positions else None,
                "end": positions[1] if len(positions) > 1 else None, "size": size}
    return None

def classify_variant(raw):
    s = raw.strip()
    if not s:
        return None
    for matcher in (match_cnv, match_sv, match_snv):
        result = matcher(s)
        if result:
            result["raw"] = raw
            result["status"] = "ok"
            return result
    return {"raw": raw, "vtype": None, "status": "needs_review"}

if __name__ == "__main__":
    with open("variants.tsv", encoding="utf-8") as f:
        for line in f:
            r = classify_variant(line)
            if r is None:
                continue
            print("[", r["status"], "] tip=", r.get("vtype"),
                  "chrom=", r.get("chrom"),
                  "pos=", r.get("pos") or r.get("start"),
                  "end=", r.get("end"),
                  "ref=", r.get("ref"), "alt=", r.get("alt"),
                  "size=", r.get("size"))
            print("     <-", r["raw"].strip())


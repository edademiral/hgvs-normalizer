import re

# Tek pozisyon: 123, -70, *100, 4072-1234 ya da belirsizse ?
POS = r"(?:[-*]?\d+(?:[+-]\d+)?|\?)"

# Sınır: tek pozisyon ya da parantezli belirsiz aralık
BOUND = rf"(?:\({POS}_{POS}\)|{POS})"

# Referans + gen ayırıcı: NM_004006.2(DMD): kısmını yakalar
REF_RE = re.compile(r"^([A-Za-z]{2,3}_?\d+(?:\.\d+)?)(?:\(([^)]+)\))?:(.+)$")

# Sonucun geçerli HGVS formatına uyup uymadığını kontrol eden kalıp
VALIDATE_RE = re.compile(
    rf"^[cgmnrp]\.{BOUND}(?:_{BOUND})?"
    rf"(?:del|dup|inv|delins[ACGTUN]+|ins(?:[ACGTUN]+|\d+_\d+))$"
)


def normalize_hgvs_sv(raw):
    # Tum bosluklari sil
    s = re.sub(r"\s+", "", raw.strip())

    # Referans (NM_...) ve gen (DMD) kismini varyanttan ayir
    ref = gene = None
    m = REF_RE.match(s)
    if m:
        ref, gene, s = m.group(1), m.group(2), m.group(3)

    # del/dup/inv/ins anahtar kelimelerini kucuk harfe cevir
    s = re.sub(r"(?i)(delins|del|dup|inv|ins)", lambda x: x.group(1).lower(), s)

    # del/dup sonrasi baz dizisini sil (delACT -> del)
    s = re.sub(r"(del|dup)([ACGTUN]+)", r"\1", s)

    # ins/delins'te eklenen diziyi KORU, buyuk harfe sabitle (insacgt -> insACGT)
    s = re.sub(r"(delins|ins)([acgtun]+)", lambda x: x.group(1) + x.group(2).upper(), s)

    # Parcalari standart forma birlestir
    canonical = (f"{ref}({gene}):{s}" if ref and gene
                 else f"{ref}:{s}" if ref else s)

    # Sonucu ve gecerlilik bilgisini dondur
    return {"raw": raw, "ref": ref, "gene": gene, "variant": s,
            "canonical": canonical, "valid": bool(VALIDATE_RE.match(s))}


# Dosya calistirilinca ornekleri dener
if __name__ == "__main__":
    tests = [
        "NM_004006.2(DMD): c.4072-1234_5155+1076DUP",
        "c.76_78delACT",
        "c.(?_-70)_(*100_?)del",
        "c.123_124insacgt",
        "c.123_456DelinsGGG",
    ]
    for t in tests:
        r = normalize_hgvs_sv(t)
        print(r["raw"], "->", r["canonical"], "valid=", r["valid"])
import re
from type_regex import classify_variant
from sema import kayit_olustur

GIRDI = "variants.tsv"


def _int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def cnv_tipi_belirle(raw):
    s = raw.lower()
    if re.search(r"\b(del|deletion|loss)\b|del$", s):
        return "copy_number_loss"
    if re.search(r"\b(dup|duplication|gain|amp)\b|dup$", s):
        return "copy_number_gain"
    return None


def sv_tipi_belirle(raw):
    s = raw.lower()
    for kelime, kod in [("del", "DEL"), ("dup", "DUP"),
                        ("inv", "INV"), ("ins", "INS")]:
        if re.search(r"\b" + kelime, s):
            return kod
    return None


def kayda_donustur(r):
    raw = r["raw"]
    chrom = r.get("chrom")
    start = _int(r.get("pos") or r.get("start"))
    end = _int(r.get("end"))
    ref = r.get("ref")
    alt = r.get("alt")
    vtype = r.get("vtype") or ""

    if r.get("status") == "needs_review" or not vtype:
        return kayit_olustur("BILINMEYEN", raw=raw, chrom=chrom, start=start)

    if vtype == "SNV":
        return kayit_olustur("SNV", raw=raw, chrom=chrom, start=start,
                             ref=ref, alt=alt)

    if vtype == "CNV":
        return kayit_olustur("CNV", raw=raw, chrom=chrom, start=start,
                             end=end, cnv_type=cnv_tipi_belirle(raw),
                             size=_int(r.get("size")))

    if not (ref and alt):
        m = re.search(r"([ACGT]+)\s*>\s*([ACGT]+)", raw, re.IGNORECASE)
        if m:
            ref, alt = m.group(1).upper(), m.group(2).upper()

    if ref and alt and (len(ref) > 1 or len(alt) > 1):
        return kayit_olustur("INDEL", raw=raw, chrom=chrom, start=start,
                             ref=ref, alt=alt)

    return kayit_olustur("SV", raw=raw, chrom=chrom, start=start,
                         end=end, svtype=sv_tipi_belirle(raw),
                         size=_int(r.get("size")))


if __name__ == "__main__":
    kovalar = {}

    with open(GIRDI, encoding="utf-8") as f:
        for satir in f:
            if not satir.strip():
                continue
            ham = classify_variant(satir)
            if ham is None:
                continue
            kayit = kayda_donustur(ham)
            grup = kayit.vtype if kayit.status == "ok" else "REVIEW"
            kovalar.setdefault(grup, []).append(kayit)

    for grup, kayitlar in kovalar.items():
        dosya = grup.lower() + "_db.tsv"
        with open(dosya, "w", encoding="utf-8") as out:
            out.write("\t".join(kayitlar[0].sutunlar()) + "\n")
            for k in kayitlar:
                out.write("\t".join(k.satir()) + "\n")

    print()
    print("#" * 58)
    print("#  ISLEM TAMAMLANDI")
    print("#" * 58)
    for grup in sorted(kovalar):
        kayitlar = kovalar[grup]
        print()
        print("  " + grup + "  ->  " + grup.lower() + "_db.tsv   (" +
              str(len(kayitlar)) + " kayit)")
        print("  " + "-" * 52)
        for k in kayitlar:
            not_metni = "; ".join(k.notes) if k.notes else ""
            print("    " + k.raw.strip()[:38].ljust(40) + " " + not_metni)
    print()
    print("#" * 58)

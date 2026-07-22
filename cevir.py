import requests
from type_regex import classify_variant

GIRDI = "variants.tsv"
GENOM = "GRCh38"
API = "https://rest.variantvalidator.org/VariantValidator/variantvalidator"


def hgvs_cevir(chrom, pos, ref, alt):
    if not (chrom and pos and ref and alt):
        return None
    varyant = chrom + "-" + str(pos) + "-" + ref + "-" + alt
    # DUZELTME 1: /all yerine /mane_select
    url = API + "/" + GENOM + "/" + varyant + "/mane_select?content-type=application/json"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            return None
        data = r.json()
        # DUZELTME 2: HGVS cevabin ANAHTARINDA (NM_...:c... seklinde)
        for anahtar in data.keys():
            if ":c." in anahtar or ":g." in anahtar or ":n." in anahtar:
                return anahtar
        return None
    except Exception:
        return None


kovalar = {"SNV": [], "SV": [], "CNV": [], "review": []}

with open(GIRDI, encoding="utf-8") as f:
    for line in f:
        r = classify_variant(line)
        if r is None:
            continue
        tip = r.get("vtype")
        if r["status"] == "needs_review" or tip is None:
            r["hgvs"] = None
            kovalar["review"].append(r)
            continue
        chrom = r.get("chrom")
        pos = r.get("pos") or r.get("start")
        ref = r.get("ref")
        alt = r.get("alt")
        print("Cevriliyor:", r["raw"].strip()[:40], "...")
        hgvs = hgvs_cevir(chrom, pos, ref, alt)
        r["hgvs"] = hgvs
        if hgvs is None and ref and alt:
            kovalar["review"].append(r)
            continue
        if tip == "SNV":
            kovalar["SNV"].append(r)
        elif tip == "CNV":
            kovalar["CNV"].append(r)
        else:
            kovalar["SV"].append(r)

dosya_adlari = {"SNV": "snv.tsv", "SV": "sv.tsv", "CNV": "cnv.tsv", "review": "review.tsv"}
for tip, kova in kovalar.items():
    with open(dosya_adlari[tip], "w", encoding="utf-8") as out:
        out.write("chrom\tstart\tend\tref\talt\ttype\tsize\thgvs\traw\n")
        for r in kova:
            out.write("\t".join([
                str(r.get("chrom") or ""),
                str(r.get("pos") or r.get("start") or ""),
                str(r.get("end") or ""),
                str(r.get("ref") or ""),
                str(r.get("alt") or ""),
                str(r.get("vtype") or ""),
                str(r.get("size") or ""),
                str(r.get("hgvs") or ""),
                r["raw"].strip(),
            ]) + "\n")

print()
print("#" * 52)
print("#  ISLEM TAMAMLANDI")
print("#" * 52)
for tip in ["SNV", "SV", "CNV", "review"]:
    kova = kovalar[tip]
    print()
    print("  " + tip + "  ->  " + dosya_adlari[tip] + "   (" + str(len(kova)) + " varyant)")
    print("  " + "-" * 46)
    for r in kova:
        h = r.get("hgvs") or "(cevrilemedi)"
        print("    " + r["raw"].strip()[:30].ljust(32) + " " + str(h))
print()
print("#" * 52)
print("  TOPLAM: " + str(sum(len(k) for k in kovalar.values())) + " varyant")
print("#" * 52)

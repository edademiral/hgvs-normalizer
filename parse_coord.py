import re

# Koordinat Parser - daginik gercek veriye gore kalibre edildi
# Isi: chr + pozisyon + bazlar/tip cikarmak, temizlemek.
# HGVS'ye CEVIRMEZ (o is seqrepo/VariantValidator). Genom: GRCh38.


def clean_number(s):
    # Binlik ayiricilari (nokta VE virgul) sil: 41.021.576 -> 41021576
    return s.replace(".", "").replace(",", "")


def parse_variant(raw):
    s = raw.strip()
    if not s:
        return None

    result = {"raw": raw, "gene": None, "chrom": None,
              "start": None, "end": None, "ref": None,
              "alt": None, "vtype": None, "size": None,
              "status": "ok"}

    # 1) Gen adi (bassa, harflerle) varsa ayir
    gene_m = re.match(r"^([A-Z][A-Z0-9]+)\s*(?:DEL|del)?\s*[:\s]", s)
    if gene_m and "chr" in s.lower():
        result["gene"] = gene_m.group(1)

    # 2) Kromozom
    chrom_m = re.search(r"chr([\dXYM]+)", s, re.IGNORECASE)
    if chrom_m:
        result["chrom"] = "chr" + chrom_m.group(1)

    # 3) SNV / indel: baz>baz
    sub_m = re.search(r"([ACGT]+)\s*>\s*([ACGT]+)", s, re.IGNORECASE)
    if sub_m:
        result["ref"] = sub_m.group(1).upper()
        result["alt"] = sub_m.group(2).upper()
        if len(result["ref"]) == 1 and len(result["alt"]) == 1:
            result["vtype"] = "SNV"
        elif len(result["ref"]) > len(result["alt"]):
            result["vtype"] = "deletion"
        else:
            result["vtype"] = "insertion"

    # 4) Pozisyon(lar)
    # Once KROMOZOMU cikar (chr12'deki 12 pozisyon sanilmasin!)
    s_nopos = re.sub(r"chr[\dXYM]+", "", s, flags=re.IGNORECASE)
    # Bazlari cikar
    s_nopos = re.sub(r"[ACGT]+\s*>\s*[ACGT]+", "", s_nopos, flags=re.IGNORECASE)
    # Boyut bilgisini cikar (18.6kb)
    s_nopos = re.sub(r"[\d.]+\s*kb", "", s_nopos, flags=re.IGNORECASE)
    # Sayilari yakala (en az 4 basamak = pozisyon)
    nums = re.findall(r"\d[\d.,]{3,}\d|\d{4,}", s_nopos)
    nums = [clean_number(n) for n in nums]
    if len(nums) == 1:
        result["start"] = nums[0]
    elif len(nums) >= 2:
        result["start"] = nums[0]
        result["end"] = nums[1]

    # 5) Tip: metinde del/dup/inv gecti mi
    if re.search(r"\bdel\b", s, re.IGNORECASE) or re.search(r"del$", s, re.IGNORECASE):
        result["vtype"] = "deletion"
    elif re.search(r"\bdup\b", s, re.IGNORECASE):
        result["vtype"] = "duplication"
    elif re.search(r"\binv\b", s, re.IGNORECASE):
        result["vtype"] = "inversion"

    # 6) Boyut (18.6kb gibi)
    size_m = re.search(r"([\d.]+)\s*kb", s, re.IGNORECASE)
    if size_m:
        result["size"] = size_m.group(1) + "kb"

    # 7) Tip hala bilinmiyorsa isaretle (elle kontrol)
    if result["vtype"] is None:
        if result["start"] and result["end"]:
            result["vtype"] = "range (tip belirsiz!)"
            result["status"] = "needs_review"
        else:
            result["status"] = "needs_review"

    return result


# Test: gercek dosyayla
if __name__ == "__main__":
    with open("variants.tsv", encoding="utf-8") as f:
        for line in f:
            r = parse_variant(line)
            if r is None:
                continue
            print("[", r["status"], "]",
                  "gene=", r["gene"], "chrom=", r["chrom"],
                  "start=", r["start"], "end=", r["end"],
                  "ref=", r["ref"], "alt=", r["alt"],
                  "type=", r["vtype"], "size=", r["size"])
            print("     <-", r["raw"].strip())
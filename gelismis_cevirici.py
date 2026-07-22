import re

GIRDI = "variants.tsv"

REFSEQ_MAP = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12",
    "5": "NC_000005.10", "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11",
    "9": "NC_000009.12", "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10", "16": "NC_000016.10",
    "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
    "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
    "M": "NC_012920.1", "MT": "NC_012920.1"
}

CHROM_RE = re.compile(r"(?:chr|chromosome|chr\.)\s*([\dXYM]+)|(NC_0000\d{2})", re.IGNORECASE)
SV_KEYWORD_RE = re.compile(r"\b(del|deletion|dup|duplication|inv|inversion|ins|insertion|loss|gain|amp|amplification|delins|indel)\b", re.IGNORECASE)
CHANGE_RE = re.compile(r"([ACGT]+)\s*(?:>|->|→|/|to)\s*([ACGT]+)", re.IGNORECASE)
STICKY_SNV_RE = re.compile(r"\b(\d+)\s*([ACGT])\s*(?:>|->|→|/|to|)\s*([ACGT])\b", re.IGNORECASE)

def pre_clean_text(s):
    # Girdi dosyasını satır satır okur
    s = s.strip()
    # Ok işaretlerini (->, →) standart '>' işaretine dönüştürür
    return s.replace("→", ">").replace("->", ">").replace("=>", ">")

def get_chrom(s):
    # Metinden kromozom adını (chr12, chrM vb.) çeker
    m = CHROM_RE.search(s)
    if not m: return None
    chrom_num = m.group(1) or m.group(2)
    if chrom_num.upper().startswith("NC_"):
        try:
            num = int(chrom_num[-2:])
            return f"chr{num}"
        except ValueError: return None
    return "chr" + chrom_num.upper()

def get_positions(s):
    # Metindeki gen adlarını ve 'kb' birimlerini siler
    s_clean = re.sub(r"[\d.]+\s*kb", "", s, flags=re.IGNORECASE)
    s_clean = re.sub(r"\b[A-Z]{2,}\d*\b", "", s_clean)
    s_clean = re.sub(r"chr(?:|变|omosome|.)\s*[\dXYM]+", "", s_clean, flags=re.IGNORECASE)
    s_clean = re.sub(r"[ACGT]+\s*>\s*[ACGT]+", "", s_clean, flags=re.IGNORECASE)
    
    # Sayılardaki nokta ve virgülleri temizleyip saf koordinatları ayıklar
    nums = re.findall(r"\b\d[\d.,]*\b", s_clean)
    cleaned_nums = []
    for n in nums:
        cleaned = n.replace(".", "").replace(",", "")
        if cleaned.isdigit() and len(cleaned) >= 4:
            cleaned_nums.append(cleaned)
    return cleaned_nums

def classify_variant(raw):
    s = pre_clean_text(raw)
    if not s: return None
    
    chrom = get_chrom(s)
    positions = get_positions(s)
    kw_match = SV_KEYWORD_RE.search(s)
    change_match = CHANGE_RE.search(s)
    sticky_match = STICKY_SNV_RE.search(s)
    
    # Metinde del, dup, inv, ins veya > işareti yoksa belirsiz (review) ilan eder
    if not kw_match and not change_match and not sticky_match:
        return {"raw": raw, "vtype": None, "status": "needs_review"}
        
    # Tek harf değişimi (A>G, T>C) varsa varyantı SNV olarak etiketler
    if sticky_match:
        return {
            "raw": raw, "status": "ok", "vtype": "SNV", "chrom": chrom,
            "pos": sticky_match.group(1), "ref": sticky_match.group(2).upper(), "alt": sticky_match.group(3).upper()
        }
    if change_match and len(change_match.group(1)) == 1 and len(change_match.group(2)) == 1:
        return {
            "raw": raw, "status": "ok", "vtype": "SNV", "chrom": chrom,
            "pos": positions[0] if positions else None, "ref": change_match.group(1).upper(), "alt": change_match.group(2).upper()
        }
        
    kw = kw_match.group(0).lower() if kw_match else "delins"
    if "del" in kw: kw = "del"
    elif "dup" in kw or "amp" in kw or "gain" in kw: kw = "dup"
    elif "inv" in kw: kw = "inv"
    elif "ins" in kw: kw = "ins"
    
    is_cnv = False
    if "kb" in s.lower() or "gain" in s.lower() or "loss" in s.lower() or "amplification" in s.lower():
        # Metinde del veya dup var ve aralık >= 1000 bp ise CNV olarak etiketler
        is_cnv = True
    elif len(positions) >= 2:
        size = abs(int(positions[1]) - int(positions[0]))
        if size >= 1000: is_cnv = True
        
    # Metinde del, dup, inv, ins var ve aralık < 1000 bp ise SV olarak etiketler
    # TGAA>T gibi harf sayısı değişen mutasyonları SV (Indel) olarak etiketler
    vtype = "CNV" if is_cnv else f"SV/{kw}"
    
    return {
        "raw": raw, "status": "ok", "vtype": vtype, "chrom": chrom,
        "start": positions[0] if positions else None,
        "end": positions[1] if len(positions) > 1 else positions[0] if positions else None,
        "kw": kw,
        "ref": change_match.group(1).upper() if change_match else None,
        "alt": change_match.group(2).upper() if change_match else None
    }

def yerel_hgvs_olustur(r):
    chrom = r.get("chrom")
    if not chrom: return None
    chr_clean = chrom.replace("chr", "").upper()
    # Kromozom adını resmi NCBI RefSeq numarasıyla (NC_...) eşleştirir
    # chr12 bilgisini NC_000012.12 koduna dönüştürür
    # chrM bilgisini NC_012920.1 koduna dönüştürür
    refseq_id = REFSEQ_MAP.get(chr_clean)
    if not refseq_id: return None
    
    vtype = r.get("vtype")
    if vtype == "SNV":
        pos = r.get("pos")
        if pos and r.get("ref") and r.get("alt"):
            # SNV için RefSeq no, pozisyon ve harf değişimini birleştirir
            return f"{refseq_id}:g.{pos}{r['ref']}>{r['alt']}"
    else:
        start = r.get("start")
        end = r.get("end")
        kw = r.get("kw", "del")
        if start:
            # Harf sayısı değişen mutasyonlar için pozisyonun sonuna 'delins' takısını ekler
            if r.get("ref") and r.get("alt"):
                return f"{refseq_id}:g.{start}delins{r['alt']}"
            if end and start != end:
                # del, dup, inv için başlangıç ve bitiş koordinatlarını alt çizgiyle bağlar
                return f"{refseq_id}:g.{start}_{end}{kw}"
            return f"{refseq_id}:g.{start}{kw}"
    return None

if __name__ == "__main__":
    kovalar = {"SNV": [], "SV": [], "CNV": [], "review": []}
    print("Kusursuz yerel motor calisiyor...")
    
    with open(GIRDI, encoding="utf-8") as f:
        for line in f:
            if not line.strip(): continue
            r = classify_variant(line)
            
            if r["status"] == "needs_review":
                r["hgvs"] = "(belirsiz tip - insan bakmali)"
                kovalar["review"].append(r)
                continue
                
            hgvs = yerel_hgvs_olustur(r)
            r["hgvs"] = hgvs
            
            # del, dup, inv, ins içermeyen veya eksik bilgi barındıran satırları review.tsv dosyasına ayırır
            if not hgvs:
                r["hgvs"] = "(eksik bilgi - cevrilemedi)"
                kovalar["review"].append(r)
                continue
                
            # Yapısı netleşen varyantları sütunlar halinde snv.tsv, sv.tsv veya cnv.tsv dosyalarına kaydeder
            if r["vtype"] == "SNV": kovalar["SNV"].append(r)
            elif r["vtype"] == "CNV": kovalar["CNV"].append(r)
            else: kovalar["SV"].append(r)

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
                    "",
                    str(r.get("hgvs") or ""),
                    r["raw"].strip()
                ]) + "\n")

    print("\n" + "#"*52 + "\n#  ISLEM TAMAMLANDI (LOKAL MOTOR)\n" + "#"*52)
    for tip in ["SNV", "SV", "CNV", "review"]:
        kova = kovalar[tip]
        print(f"\n  {tip}  ->  {dosya_adlari[tip]}   ({len(kova)} varyant)")
        print("  " + "-"*46)
        for r in kova:
            print("    " + r["raw"].strip()[:30].ljust(32) + " " + str(r["hgvs"]))
    print("\n" + "#"*52 + f"\n  TOPLAM: {sum(len(k) for k in kovalar.values())} varyant\n" + "#"*52)
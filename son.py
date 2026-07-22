import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.validator
from type_regex import classify_variant
from kopru import kayda_donustur, GIRDI

# UTA Veritabanı araçlarını global olarak bir kere başlatıyoruz (Performans için)
print("UTA Veritabanı bağlantısı kuruluyor...")
HP = hgvs.parser.Parser()
try:
    HDP = hgvs.dataproviders.uta.connect()
    AM = hgvs.assemblymapper.AssemblyMapper(
        HDP, assembly_name='GRCh38', alt_aln_method='splign', replace_reference=True
    )
    VR = hgvs.validator.Validator(HDP)
    UTA_AKTIF = True
    print("UTA Bağlantısı Başarılı!\n")
except Exception as e:
    UTA_AKTIF = False
    print(f"UYARI: UTA bağlantısı kurulamadı ({e}). Sadece genomik taslaklar yazılacak.\n")

def taslak_uret(k):
    acc = k.ref_accession
    if not acc or k.start is None:
        return None
    t = k.vtype
    if t == "SNV":
        return acc + ":g." + str(k.start) + k.ref + ">" + k.alt
    if t == "INDEL":
        return acc + ":g." + str(k.start) + "delins" + k.alt
    if t == "CNV" and k.end:
        yon = "dup" if k.cnv_type == "copy_number_gain" else "del"
        return acc + ":g." + str(k.start) + "_" + str(k.end) + yon
    if t == "SV" and k.end:
        yon = {"DEL": "del", "DUP": "dup", "INV": "inv"}.get(k.svtype, "del")
        return acc + ":g." + str(k.start) + "_" + str(k.end) + yon
    if t == "REPEAT":
        return (acc + ":g." + str(k.start) + k.repeat_unit +
                "[" + str(k.copy_number) + "]")
    return None

def hgvs_dogrula_ve_donustur(g_taslak):
    """
    Genomik taslağı (g.) alır, UTA üzerinden doğrular ve 
    mümkünse transkript düzeyine (c.) çevirir.
    """
    if not UTA_AKTIF or not g_taslak:
        return g_taslak, "dogrulanmadi"
    
    try:
        var_g = HP.parse_hgvs_variant(g_taslak)
        VR.validate(var_g)
        var_c = AM.g_to_c(var_g)
        return str(var_c), "dogrulandi"
    except Exception:
        return g_taslak, "dogrulanmadi"


EK = ["hgvs_taslak", "validation_status"]

kovalar = {}
with open(GIRDI, encoding="utf-8") as f:
    for satir in f:
        if not satir.strip():
            continue
        ham = classify_variant(satir)
        if ham is None:
            continue
        k = kayda_donustur(ham)
        
        # Taslak üretimi ve doğrulama adımları
        ham_taslak = taslak_uret(k) if k.status == "ok" else None
        
        if ham_taslak:
            nihai_taslak, durum = hgvs_dogrula_ve_donustur(ham_taslak)
        else:
            nihai_taslak, durum = "", "uretilemedi"
            
        k._ek = [nihai_taslak, durum]
        grup = k.vtype if k.status == "ok" else "REVIEW"
        kovalar.setdefault(grup, []).append(k)

for grup, kayitlar in kovalar.items():
    with open(grup.lower() + "_db.tsv", "w", encoding="utf-8") as out:
        out.write("\t".join(kayitlar[0].sutunlar() + EK) + "\n")
        for k in kayitlar:
            out.write("\t".join(k.satir() + k._ek) + "\n")

print()
print("=" * 58)
for grup in sorted(kovalar):
    print()
    print("  " + grup + "  (" + str(len(kovalar[grup])) + ")  ->  " +
          grup.lower() + "_db.tsv")
    for k in kovalar[grup]:
        print("    " + k.raw.strip()[:32].ljust(34) +
              "[" + k._ek[1] + "] " + k._ek[0])
print()
print("=" * 58)
print("UTA aktif mi:", UTA_AKTIF)

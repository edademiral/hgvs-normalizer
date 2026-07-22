import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.validator

def varyant_dogrula_ve_cevir(taslak_hgvs):
    hp = hgvs.parser.Parser()
    
    try:
        print(f"\n--- {taslak_hgvs} İşleniyor ---")
        var_g = hp.parse_hgvs_variant(taslak_hgvs)
        print("1. [BAŞARILI] Sözdizimi geçerli. Obje oluşturuldu.")
    except Exception as e:
        print(f"1. [HATA] Ayrıştırma başarısız: {e}")
        return "[uretilemedi]"

    hdp = hgvs.dataproviders.uta.connect()
    vr = hgvs.validator.IntrinsicValidator()
    
    try:
        vr.validate(var_g)
        print("2. [BAŞARILI] Varyant kurallara uygun (Intrinsic Validation).")
    except Exception as e:
        print(f"2. [HATA] Doğrulama başarısız: {e}")
        return "[dogrulanmadi]"

    am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh38')
    
    try:
        ilgili_transkriptler = am.relevant_transcripts(var_g)
        if ilgili_transkriptler:
            var_c = am.g_to_c(var_g, ilgili_transkriptler[0])
            print(f"3. [BAŞARILI] Transkripte dönüştürüldü: {var_c}")
            return str(var_c)
        else:
            print("3. [UYARI] Bu lokasyon için transkript bulunamadı.")
            return str(var_g)
            
    except Exception as e:
        print(f"3. [HATA] Haritalama başarısız: {e}")
        return "[dogrulanmadi]"

test_varyantlari = [
    "NC_000011.10:g.45917701T>C",  
    "NC_000099.99:g.12345A>T"      
]

for v in test_varyantlari:
    sonuc = varyant_dogrula_ve_cevir(v)
    print(f"👉 NİHAİ ÇIKTI: {sonuc}\n")

"""
sema.py - Varyant veri modeli (sema) ve dogrulama kurallari.
Her varyantin ORTAK cekirdegi var; turune gore SPESIFIK alanlar ekleniyor.
Zorunlu alan eksikse kayit otomatik needs_review olur - uydurma yapilmaz.
"""

from dataclasses import dataclass, field, asdict
from typing import Optional, List

REFSEQ_GRCH38 = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
    "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
    "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
    "M": "NC_012920.1", "MT": "NC_012920.1",
}

BAZLAR = set("ACGTN")


def refseq_bul(chrom, assembly="GRCh38"):
    if not chrom or assembly != "GRCh38":
        return None
    anahtar = str(chrom).replace("chr", "").replace("Chr", "").upper()
    return REFSEQ_GRCH38.get(anahtar)


@dataclass
class VariantRecord:
    raw: str = ""
    assembly: str = "GRCh38"
    chrom: Optional[str] = None
    start: Optional[int] = None
    ref_accession: Optional[str] = None
    hgvs: Optional[str] = None
    vtype: Optional[str] = None
    status: str = "ok"
    notes: List[str] = field(default_factory=list)

    ORTAK_SUTUNLAR = ["assembly", "chrom", "start", "ref_accession",
                      "vtype", "hgvs", "status", "notes", "raw"]
    OZEL_SUTUNLAR: List[str] = field(default_factory=list, repr=False)

    def dogrula(self):
        if not self.chrom:
            self.notes.append("kromozom yok")
        if self.start is None:
            self.notes.append("baslangic pozisyonu yok")
        if self.chrom and not self.ref_accession:
            self.ref_accession = refseq_bul(self.chrom, self.assembly)
            if not self.ref_accession:
                self.notes.append("RefSeq bulunamadi: " + str(self.chrom))
        self._dogrula_ozel()
        if self.notes:
            self.status = "needs_review"
        return self.status == "ok"

    def _dogrula_ozel(self):
        pass

    def sutunlar(self):
        return self.ORTAK_SUTUNLAR[:-1] + list(self.OZEL_SUTUNLAR) + ["raw"]

    def satir(self):
        d = asdict(self)
        d["notes"] = "; ".join(self.notes)
        return [str(d.get(s) if d.get(s) is not None else "")
                for s in self.sutunlar()]


@dataclass
class SNVRecord(VariantRecord):
    ref: Optional[str] = None
    alt: Optional[str] = None

    def __post_init__(self):
        self.vtype = "SNV"
        self.OZEL_SUTUNLAR = ["ref", "alt"]

    def _dogrula_ozel(self):
        if not self.ref or not self.alt:
            self.notes.append("SNV icin ref/alt zorunlu")
            return
        if len(self.ref) != 1 or len(self.alt) != 1:
            self.notes.append("SNV tek baz olmali (indel olabilir)")
        if set(self.ref.upper()) - BAZLAR or set(self.alt.upper()) - BAZLAR:
            self.notes.append("gecersiz baz karakteri")
        if self.ref.upper() == self.alt.upper():
            self.notes.append("ref ve alt ayni")


@dataclass
class IndelRecord(VariantRecord):
    ref: Optional[str] = None
    alt: Optional[str] = None
    size: Optional[int] = None

    def __post_init__(self):
        self.vtype = "INDEL"
        self.OZEL_SUTUNLAR = ["ref", "alt", "size"]

    def _dogrula_ozel(self):
        if not self.ref or not self.alt:
            self.notes.append("INDEL icin ref/alt zorunlu")
            return
        if set(self.ref.upper()) - BAZLAR or set(self.alt.upper()) - BAZLAR:
            self.notes.append("gecersiz baz karakteri")
        if len(self.ref) == len(self.alt) == 1:
            self.notes.append("tek baz -> SNV olmali")
        if self.size is None:
            self.size = abs(len(self.ref) - len(self.alt))


GECERLI_SVTYPE = {"DEL", "DUP", "INV", "INS", "BND", "DELINS"}


@dataclass
class SVRecord(VariantRecord):
    end: Optional[int] = None
    svtype: Optional[str] = None
    size: Optional[int] = None
    ci_start: Optional[str] = None
    ci_end: Optional[str] = None

    def __post_init__(self):
        self.vtype = "SV"
        self.OZEL_SUTUNLAR = ["end", "svtype", "size", "ci_start", "ci_end"]

    def _dogrula_ozel(self):
        if not self.svtype:
            self.notes.append("SV icin svtype zorunlu (DEL/DUP/INV/INS)")
        elif self.svtype.upper() not in GECERLI_SVTYPE:
            self.notes.append("bilinmeyen svtype: " + str(self.svtype))
        else:
            self.svtype = self.svtype.upper()
        if self.svtype and self.svtype.upper() != "INS" and self.end is None:
            self.notes.append("SV icin bitis koordinati zorunlu")
        if self.start is not None and self.end is not None:
            if self.end < self.start:
                self.notes.append("bitis baslangictan kucuk")
            elif self.size is None:
                self.size = self.end - self.start


@dataclass
class CNVRecord(VariantRecord):
    end: Optional[int] = None
    cnv_type: Optional[str] = None
    copy_number: Optional[int] = None
    size: Optional[int] = None

    def __post_init__(self):
        self.vtype = "CNV"
        self.OZEL_SUTUNLAR = ["end", "cnv_type", "copy_number", "size"]

    def _dogrula_ozel(self):
        if self.end is None:
            self.notes.append("CNV icin bitis koordinati zorunlu")
        if not self.cnv_type:
            self.notes.append("CNV icin tip zorunlu (loss/gain)")
        elif self.cnv_type not in ("copy_number_loss", "copy_number_gain"):
            self.notes.append("gecersiz cnv_type: " + str(self.cnv_type))
        if self.start is not None and self.end is not None:
            if self.end < self.start:
                self.notes.append("bitis baslangictan kucuk")
            elif self.size is None:
                self.size = self.end - self.start


@dataclass
class RepeatRecord(VariantRecord):
    repeat_unit: Optional[str] = None
    copy_number: Optional[int] = None
    ref_copy_number: Optional[int] = None

    def __post_init__(self):
        self.vtype = "REPEAT"
        self.OZEL_SUTUNLAR = ["repeat_unit", "copy_number", "ref_copy_number"]

    def _dogrula_ozel(self):
        if not self.repeat_unit:
            self.notes.append("REPEAT icin tekrar motifi zorunlu")
        elif set(self.repeat_unit.upper()) - BAZLAR:
            self.notes.append("gecersiz tekrar motifi")
        if self.copy_number is None:
            self.notes.append("REPEAT icin kopya sayisi zorunlu")


SINIFLAR = {
    "SNV": SNVRecord,
    "INDEL": IndelRecord,
    "SV": SVRecord,
    "CNV": CNVRecord,
    "REPEAT": RepeatRecord,
}


def kayit_olustur(vtype, **alanlar):
    sinif = SINIFLAR.get(str(vtype).upper())
    if sinif is None:
        r = VariantRecord(**{k: v for k, v in alanlar.items()
                             if k in VariantRecord.__annotations__})
        r.status = "needs_review"
        r.notes.append("bilinmeyen tur: " + str(vtype))
        return r
    gecerli = {k: v for k, v in alanlar.items()
               if k in sinif.__annotations__ or k in VariantRecord.__annotations__}
    r = sinif(**gecerli)
    r.dogrula()
    return r


if __name__ == "__main__":
    ornekler = [
        ("SNV", dict(raw="chr11-45917701 T>C", chrom="chr11",
                     start=45917701, ref="T", alt="C")),
        ("INDEL", dict(raw="chr6-43624673 TGAA>T", chrom="chr6",
                       start=43624673, ref="TGAA", alt="T")),
        ("CNV", dict(raw="chr12 del 18.6kb", chrom="chr12",
                     start=41021576, end=41040185,
                     cnv_type="copy_number_loss")),
        ("CNV", dict(raw="tip yok", chrom="chr5",
                     start=126489209, end=126759255)),
        ("REPEAT", dict(raw="HTT CAG(42)", chrom="chr4",
                        start=3074877, repeat_unit="CAG", copy_number=42)),
    ]
    for vtype, alanlar in ornekler:
        r = kayit_olustur(vtype, **alanlar)
        print("=" * 55)
        print(r.vtype, "-", r.status)
        print("  RefSeq  :", r.ref_accession)
        print("  Sutunlar:", r.sutunlar())
        if r.notes:
            print("  UYARI   :", "; ".join(r.notes))

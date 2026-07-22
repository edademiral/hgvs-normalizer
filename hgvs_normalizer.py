#!/usr/bin/env python3
"""
hgvs_normalizer.py - Normalize free-text variant descriptions to HGVS.

Pipeline stages
---------------
1. Parse      : messy literature/clinical text -> typed record
2. Schema     : type-specific required fields; incomplete records are flagged
3. Draft      : build an *unvalidated* genomic HGVS candidate (GRCh38)
4. Validate   : OPTIONAL - biocommons hgvs + UTA; reference-base checked

Design principle
----------------
An unvalidated string is NEVER written to the `hgvs` column. Drafts live in
`hgvs_draft` and every record carries an explicit `validation_status`.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import logging
import os
import re
import sys
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Sequence

__version__ = "0.5.0" # Sürüm güncellendi (Bug fix içerir)

ASSEMBLY = "GRCh38"
SV_MIN_BP = 50          # pending mentor confirmation
CNV_MIN_BP = 1_000      # pending mentor confirmation
VALIDATION_MAX_BP = 100_000   # larger regions are skipped (slow/unsupported)

logger = logging.getLogger("hgvs_normalizer")

# GRCh38 RefSeq accessions. NCBI versions change - re-verify periodically.
REFSEQ_GRCH38: Dict[str, str] = {
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
NUCLEOTIDES = frozenset("ACGTN")


def lookup_accession(chrom: Optional[str]) -> Optional[str]:
    if not chrom:
        return None
    return REFSEQ_GRCH38.get(str(chrom).replace("chr", "").upper())


# ===================================================================
# Stage 2: schema - common core plus type-specific fields
# ===================================================================
@dataclass
class VariantRecord:
    """Base record. Fields shared by every variant type."""

    raw: str = ""
    assembly: str = ASSEMBLY
    chrom: Optional[str] = None
    start: Optional[int] = None
    accession: Optional[str] = None
    hgvs: Optional[str] = None             # validated only
    hgvs_draft: Optional[str] = None       # unvalidated candidate
    vtype: Optional[str] = None
    status: str = "ok"                     # ok | needs_review
    validation_status: str = "not_validated"
    notes: List[str] = field(default_factory=list)

    COMMON_COLUMNS: Sequence[str] = (
        "assembly", "chrom", "start", "accession", "vtype",
        "hgvs", "hgvs_draft", "validation_status", "status", "notes",
    )
    TYPE_COLUMNS: List[str] = field(default_factory=list, repr=False)

    def check(self) -> bool:
        if not self.chrom:
            self.notes.append("missing chromosome")
        if self.start is None:
            self.notes.append("missing start position")
        if self.chrom and not self.accession:
            self.accession = lookup_accession(self.chrom)
            if not self.accession:
                self.notes.append(f"no RefSeq accession for {self.chrom}")
        self._check_type_specific()
        if self.notes:
            self.status = "needs_review"
        return self.status == "ok"

    def _check_type_specific(self) -> None:
        """Overridden by subclasses."""

    def columns(self) -> List[str]:
        return list(self.COMMON_COLUMNS) + list(self.TYPE_COLUMNS) + ["raw"]

    def as_row(self) -> List[str]:
        data = asdict(self)
        data["notes"] = "; ".join(self.notes)
        return ["" if data.get(c) is None else str(data.get(c))
                for c in self.columns()]


@dataclass
class SNVRecord(VariantRecord):
    ref: Optional[str] = None
    alt: Optional[str] = None

    def __post_init__(self):
        self.vtype = "SNV"
        self.TYPE_COLUMNS = ["ref", "alt"]

    def _check_type_specific(self):
        if not self.ref or not self.alt:
            self.notes.append("SNV requires ref and alt")
            return
        if len(self.ref) != 1 or len(self.alt) != 1:
            self.notes.append("SNV must be a single base on both sides")
        if set(self.ref.upper()) - NUCLEOTIDES or set(self.alt.upper()) - NUCLEOTIDES:
            self.notes.append("invalid nucleotide character")
        if self.ref.upper() == self.alt.upper():
            self.notes.append("ref equals alt")


@dataclass
class IndelRecord(VariantRecord):
    ref: Optional[str] = None
    alt: Optional[str] = None
    size: Optional[int] = None

    def __post_init__(self):
        self.vtype = "INDEL"
        self.TYPE_COLUMNS = ["ref", "alt", "size"]

    def _check_type_specific(self):
        if not self.ref or not self.alt:
            self.notes.append("indel requires ref and alt")
            return
        if set(self.ref.upper()) - NUCLEOTIDES or set(self.alt.upper()) - NUCLEOTIDES:
            self.notes.append("invalid nucleotide character")
        if self.size is None:
            self.size = abs(len(self.ref) - len(self.alt))


VALID_SVTYPES = frozenset({"DEL", "DUP", "INV", "INS"})


@dataclass
class SVRecord(VariantRecord):
    end: Optional[int] = None
    svtype: Optional[str] = None
    size: Optional[int] = None

    def __post_init__(self):
        self.vtype = "SV"
        self.TYPE_COLUMNS = ["end", "svtype", "size"]

    def _check_type_specific(self):
        if not self.svtype:
            self.notes.append("SV requires svtype")
        elif self.svtype.upper() not in VALID_SVTYPES:
            self.notes.append(f"unknown svtype: {self.svtype}")
        else:
            self.svtype = self.svtype.upper()
        if self.svtype and self.svtype != "INS" and self.end is None:
            self.notes.append("SV requires an end coordinate")
        if self.start is not None and self.end is not None:
            if self.end < self.start:
                self.notes.append("end precedes start")
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
        self.TYPE_COLUMNS = ["end", "cnv_type", "copy_number", "size"]

    def _check_type_specific(self):
        if self.end is None:
            self.notes.append("CNV requires an end coordinate")
        if not self.cnv_type:
            self.notes.append("CNV requires a direction (loss/gain)")
        elif self.cnv_type not in ("copy_number_loss", "copy_number_gain"):
            self.notes.append(f"invalid cnv_type: {self.cnv_type}")
        if self.start is not None and self.end is not None:
            if self.end < self.start:
                self.notes.append("end precedes start")
            elif self.size is None:
                self.size = self.end - self.start


@dataclass
class RepeatRecord(VariantRecord):
    repeat_unit: Optional[str] = None
    copy_number: Optional[int] = None

    def __post_init__(self):
        self.vtype = "REPEAT"
        self.TYPE_COLUMNS = ["repeat_unit", "copy_number"]

    def _check_type_specific(self):
        if not self.repeat_unit:
            self.notes.append("repeat requires a motif")
        elif set(self.repeat_unit.upper()) - NUCLEOTIDES:
            self.notes.append("invalid repeat motif")
        if self.copy_number is None:
            self.notes.append("repeat requires a copy number")


# ===================================================================
# Stage 1: parser (Bug fixes applied: thousands vs decimal separators)
# ===================================================================
class VariantParser:
    """Extracts structured fields from free-text variant descriptions."""

    CHROM = re.compile(r"chr[\s.]*([\dXYM]+|MT)", re.I)
    SUBSTITUTION = re.compile(r"\b([ACGT]+)\s*(?:>|->|\u2192)\s*([ACGT]+)\b", re.I)
    REPEAT = re.compile(r"\b([ACGT]{2,6})\s*[\(\[]\s*(\d+)\s*[\)\]]", re.I)
    KILOBASE = re.compile(r"([\d.,]+)\s*kb\b", re.I)
    LOSS = re.compile(r"\b(del|deletion|loss)\b|del\b", re.I)
    GAIN = re.compile(r"\b(dup|duplication|gain|amp)\b|dup\b", re.I)
    INVERSION = re.compile(r"\b(inv|inversion)\b", re.I)
    INSERTION = re.compile(r"\b(?:ins|insertion)\s*([ACGT]+)\b", re.I)
    DELINS = re.compile(r"\b(?:del|deletion)\s*([ACGT]*)\s*(?:ins|insertion)\s*([ACGT]+)\b", re.I)
    
    # Yeni kural: Yalnızca arkasında tam 3 rakam varsa (veya bu şekilde tekrarlanıyorsa) binlik ayırıcıdır.
    THOUSANDS = re.compile(r"^\d{1,3}(?:[.,]\d{3})+$")

    @classmethod
    def _parse_number(cls, token: str) -> Optional[float]:
        """Convert a string number with potential separators into a float."""
        if cls.THOUSANDS.match(token):
            # Binlik ayırıcı: tüm nokta ve virgülleri sil.
            cleaned = token.replace(".", "").replace(",", "")
            return float(cleaned)
        
        # Eğer virgül varsa ve binlik ayırıcı kuralına uymuyorsa, büyük ihtimalle ondalık ayırıcıdır.
        # Nokta zaten varsayılan ondalık ayırıcı. Virgülü noktaya çevir.
        cleaned = token.replace(",", ".")
        try:
            return float(cleaned)
        except ValueError:
            return None

    def _positions(self, text: str) -> List[int]:
        """Numeric coordinates, with chromosome/bases/sizes removed first."""
        cleaned = self.CHROM.sub(" ", text)
        cleaned = self.SUBSTITUTION.sub(" ", cleaned)
        cleaned = self.DELINS.sub(" ", cleaned)
        cleaned = self.INSERTION.sub(" ", cleaned)
        cleaned = self.KILOBASE.sub(" ", cleaned)
        cleaned = self.REPEAT.sub(" ", cleaned)
        out: List[int] = []
        for token in re.findall(r"\d[\d.,]{2,}\d|\d{3,}", cleaned):
            val = self._parse_number(token)
            if val is not None and val.is_integer(): # Koordinatlar tam sayı olmalıdır.
                out.append(int(val))
        return out

    def _size_and_mismatch(self, text: str, positions: Sequence[int]) -> tuple[Optional[int], Optional[str]]:
        """
        Calculates size primarily from coordinates.
        Uses text 'kb' description for cross-validation.
        Returns: (size_in_bp, mismatch_warning_if_any)
        """
        coord_size = None
        text_size = None

        if len(positions) >= 2:
            coord_size = abs(positions[1] - positions[0])

        match = self.KILOBASE.search(text)
        if match:
             val = self._parse_number(match.group(1))
             if val is not None:
                 text_size = int(val * 1000)

        # 1. Koordinat varsa Otorite Koordinattır.
        if coord_size is not None:
            if text_size is not None:
                # %10'dan fazla sapma varsa uyar (186000 vs 18609 gibi durumlar)
                if abs(coord_size - text_size) > (coord_size * 0.1):
                     warning = f"size mismatch: coords implies {coord_size}bp, text implies {text_size}bp. Using coords."
                     return coord_size, warning
            return coord_size, None
        
        # 2. Sadece metin varsa ona güven
        if text_size is not None:
             return text_size, None

        return None, None


    def parse(self, line: str) -> Optional[VariantRecord]:
        text = line.strip()
        if not text or text.startswith("#"):
            return None

        chrom_match = self.CHROM.search(text)
        chrom = f"chr{chrom_match.group(1).upper()}" if chrom_match else None
        positions = self._positions(text)
        start = positions[0] if positions else None
        end = positions[1] if len(positions) > 1 else None
        
        size, mismatch_warning = self._size_and_mismatch(text, positions)
        common = dict(raw=line, chrom=chrom, start=start)

        record = None

        repeat = self.REPEAT.search(text)
        if repeat:
            record = RepeatRecord(repeat_unit=repeat.group(1).upper(),
                                  copy_number=int(repeat.group(2)), **common)

        elif substitution := self.SUBSTITUTION.search(text):
            ref = substitution.group(1).upper()
            alt = substitution.group(2).upper()
            record = (SNVRecord(ref=ref, alt=alt, **common)
                      if len(ref) == 1 and len(alt) == 1
                      else IndelRecord(ref=ref, alt=alt, **common))
        elif delins := self.DELINS.search(text):
            deleted_seq = delins.group(1).upper()
            inserted_seq = delins.group(2).upper()
            record = IndelRecord(ref=deleted_seq or "N", alt=inserted_seq, **common)

        elif ins_match := self.INSERTION.search(text):
            inserted_seq = ins_match.group(1).upper()
            record = IndelRecord(ref="N", alt=inserted_seq, **common)

        else:
            if self.LOSS.search(text):
                direction, svtype = "copy_number_loss", "DEL"
            elif self.GAIN.search(text):
                direction, svtype = "copy_number_gain", "DUP"
            elif self.INVERSION.search(text):
                direction, svtype = None, "INV"
            else:
                record = VariantRecord(**common)
                record.status = "needs_review"
                record.notes.append("variant type could not be determined")
                record.accession = lookup_accession(chrom)

            if not record:
                if svtype in ("DEL", "DUP") and size is not None and size >= CNV_MIN_BP:
                    record = CNVRecord(end=end, cnv_type=direction, size=size, **common)
                else:
                    record = SVRecord(end=end, svtype=svtype, size=size, **common)

        if record:
             if mismatch_warning:
                 record.notes.append(mismatch_warning)
             record.check()
             
        return record


# ===================================================================
# Stage 3: draft builder
# ===================================================================
class HgvsDraftBuilder:
    """Builds an UNVALIDATED genomic HGVS candidate from a record."""

    def build(self, record: VariantRecord) -> Optional[str]:
        acc, pos = record.accession, record.start
        if not acc or pos is None:
            return None
        vtype = record.vtype
        if vtype == "SNV":
            return f"{acc}:g.{pos}{record.ref}>{record.alt}"
        if vtype == "INDEL":
            return f"{acc}:g.{pos}delins{record.alt}"
        if vtype == "REPEAT":
            return f"{acc}:g.{pos}{record.repeat_unit}[{record.copy_number}]"
        if vtype == "CNV" and record.end:
            op = "dup" if record.cnv_type == "copy_number_gain" else "del"
            return f"{acc}:g.{pos}_{record.end}{op}"
        if vtype == "SV" and record.end:
            op = {"DEL": "del", "DUP": "dup", "INV": "inv"}.get(record.svtype)
            if op:
                return f"{acc}:g.{pos}_{record.end}{op}"
        return None


# ===================================================================
# Stage 4: validators (interchangeable)
# ===================================================================
class NullValidator:
    name = "disabled"
    data_version = "n/a"

    def validate(self, record: VariantRecord) -> None:
        record.validation_status = "not_validated"


class UtaValidator:
    name = "uta"

    def __init__(self) -> None:
        import hgvs.parser
        import hgvs.dataproviders.uta
        import hgvs.assemblymapper
        import hgvs.validator

        self._parser = hgvs.parser.Parser()
        self._hdp = hgvs.dataproviders.uta.connect()
        self.data_version = self._hdp.data_version()
        self._mapper = hgvs.assemblymapper.AssemblyMapper(
            self._hdp, assembly_name=ASSEMBLY,
            alt_aln_method="splign", replace_reference=True,
        )
        self._validator = hgvs.validator.Validator(self._hdp)

    def validate(self, record: VariantRecord) -> None:
        draft = record.hgvs_draft
        if not draft:
            record.validation_status = "no_draft"
            return
        size = getattr(record, "size", None)
        if size and size > VALIDATION_MAX_BP:
            record.validation_status = "skipped_too_large"
            return
        try:
            parsed = self._parser.parse_hgvs_variant(draft)
        except Exception as exc:
            record.validation_status = "unparseable"
            record.notes.append(type(exc).__name__)
            return
        try:
            self._validator.validate(parsed)
        except Exception as exc:
            record.validation_status = "failed_validation"
            record.notes.append(str(exc)[:100])
            return
        try:
            record.hgvs = str(self._mapper.g_to_c(parsed))
            record.validation_status = "validated_c"
        except Exception:
            record.hgvs = str(parsed)
            record.validation_status = "validated_g"


# ===================================================================
# Orchestration
# ===================================================================
class NormalizationPipeline:
    def __init__(self, parser=None, builder=None, validator=None):
        self.parser = parser or VariantParser()
        self.builder = builder or HgvsDraftBuilder()
        self.validator = validator or NullValidator()
        self.failed_records = []

    def run(self, input_path: str) -> List[VariantRecord]:
        records: List[VariantRecord] = []
        self.failed_records = []
        with open(input_path, encoding="utf-8") as handle:
            for line in handle:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                
                try:
                    record = self.parser.parse(line)
                    if record is None:
                        continue
                    if record.status == "ok":
                        record.hgvs_draft = self.builder.build(record)
                    records.append(record)
                except Exception as e:
                    self.failed_records.append({"raw_line": text, "error": str(e)})

        logger.info("parsed %d records", len(records))

        for index, record in enumerate(records, 1):
            self.validator.validate(record)
            logger.debug("[%d/%d] %s -> %s", index, len(records),
                         record.raw.strip()[:40], record.validation_status)
        return records

    def write(self, records, output_dir: str) -> Dict[str, List[VariantRecord]]:
        os.makedirs(output_dir, exist_ok=True)
        buckets: Dict[str, List[VariantRecord]] = {}
        for record in records:
            key = record.vtype if (record.status == "ok" and record.vtype) else "REVIEW"
            buckets.setdefault(key, []).append(record)

        for key, group in buckets.items():
            path = os.path.join(output_dir, f"{key.lower()}.tsv")
            with open(path, "w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(group[0].columns())
                for record in group:
                    writer.writerow(record.as_row())

        if self.failed_records:
            failed_path = os.path.join(output_dir, "failed_records.tsv")
            with open(failed_path, "w", encoding="utf-8") as handle:
                handle.write("raw_line\terror_message\n")
                for err_rec in self.failed_records:
                    handle.write(f"{err_rec['raw_line']}\t{err_rec['error']}\n")

        manifest = os.path.join(output_dir, "run_manifest.txt")
        with open(manifest, "w", encoding="utf-8") as handle:
            handle.write(f"tool_version: {__version__}\n")
            handle.write("timestamp: %s\n"
                         % datetime.datetime.now().isoformat(timespec="seconds"))
            handle.write(f"assembly: {ASSEMBLY}\n")
            handle.write(f"sv_min_bp: {SV_MIN_BP}\n")
            handle.write(f"cnv_min_bp: {CNV_MIN_BP}\n")
            handle.write(f"validator: {self.validator.name}\n")
            handle.write(f"validator_data_version: {self.validator.data_version}\n")
            for key, group in sorted(buckets.items()):
                handle.write(f"{key}: {len(group)}\n")
            handle.write(f"FAILED: {len(self.failed_records)}\n")
            
        return buckets


# ===================================================================
# Self-test
# ===================================================================
TEST_CASES = [
    ("chr11-45917701 T>C", "SNV", "ok"),
    ("chrM-8612 T>C", "SNV", "ok"),
    ("chr6-43624673 TGAA>T", "INDEL", "ok"),
    ("CNTN1:chr12:41.021.576-41.040.185 del 18.6kb", "CNV", "ok"),
    ("CASK DEL:chrX:41772635-41895004", "CNV", "ok"),
    ("chr5-126,489,209-126,759,255", None, "needs_review"),
    ("HTT chr4:3074877 CAG(42)", "REPEAT", "ok"),
    ("chr3:1000-5000 inversion", "SV", "ok"),
    ("chr9:100.000-200.000 del 1,5kb", "CNV", "needs_review"), # Yeni test eklendi (mismatch yakalamalı)
]


def run_self_test() -> bool:
    parser = VariantParser()
    passed = failed = 0
    for text, expected_type, expected_status in TEST_CASES:
        record = parser.parse(text)
        actual_type = record.vtype if record else None
        actual_status = record.status if record else None
        if actual_type == expected_type and actual_status == expected_status:
            passed += 1
            print(f"  PASS  {text[:46]:<46} {actual_type}")
        else:
            failed += 1
            print(f"  FAIL  {text[:46]:<46} "
                  f"expected={expected_type}/{expected_status} "
                  f"got={actual_type}/{actual_status}")
    print(f"\n  {passed} passed, {failed} failed")
    return failed == 0


# ===================================================================
# CLI
# ===================================================================
def main() -> None:
    argp = argparse.ArgumentParser(
        description="Normalize free-text variant descriptions to HGVS")
    argp.add_argument("--input", default="variants.tsv")
    argp.add_argument("--output-dir", default="output")
    argp.add_argument("--validate", action="store_true",
                      help="validate against UTA (needs network, slow)")
    argp.add_argument("--self-test", action="store_true")
    argp.add_argument("--verbose", action="store_true")
    args = argp.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s")

    if args.self_test:
        print("\nKNOWN-ANSWER TESTS\n" + "-" * 62)
        sys.exit(0 if run_self_test() else 1)

    if not os.path.exists(args.input):
        logger.error("input not found: %s", args.input)
        sys.exit(1)

    validator = NullValidator()
    if args.validate:
        logger.info("connecting to UTA (first connection can be slow)...")
        try:
            validator = UtaValidator()
            logger.info("UTA ready: %s", validator.data_version)
        except Exception as exc:
            logger.warning("validation unavailable (%s); drafts stay unvalidated",
                           type(exc).__name__)

    pipeline = NormalizationPipeline(validator=validator)
    records = pipeline.run(args.input)
    buckets = pipeline.write(records, args.output_dir)

    print("\n" + "=" * 66)
    print(f"  DONE  (tool v{__version__}, {ASSEMBLY}, "
          f"validator: {validator.name} / {validator.data_version})")
    print("=" * 66)
    for key in sorted(buckets):
        print(f"\n  {key}  ({len(buckets[key])})  ->  "
              f"{args.output_dir}/{key.lower()}.tsv")
        for record in buckets[key]:
            print(f"    {record.raw.strip()[:32]:<32} "
                  f"[{record.validation_status}] "
                  f"{record.hgvs or record.hgvs_draft or ''}")
            
    if pipeline.failed_records:
        print(f"\n  FAILED  ({len(pipeline.failed_records)})  ->  {args.output_dir}/failed_records.tsv")
        
    print("\n" + "=" * 66)
    print(f"  Provenance: {args.output_dir}/run_manifest.txt")


if __name__ == "__main__":
    main()
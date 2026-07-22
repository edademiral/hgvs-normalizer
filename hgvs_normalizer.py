#!/usr/bin/env python3
"""
hgvs_normalizer.py - Normalize free-text variant descriptions to HGVS.

Pipeline stages
---------------
0. Passthrough : input that is ALREADY valid HGVS is detected and forwarded
                 to the validator untouched (never re-parsed)
1. Parse       : messy literature/clinical text -> typed record
2. Schema      : type-specific required fields; incomplete records are flagged
3. Draft       : build an *unvalidated* genomic HGVS candidate (GRCh38),
                 with VCF-style allele trimming and optional 3' shifting
4. Validate    : OPTIONAL - local FASTA (pyfaidx) and/or biocommons hgvs + UTA

Design principle
----------------
An unvalidated string is NEVER written to the `hgvs` column. Drafts live in
`hgvs_draft` and every record carries an explicit `validation_status`.
Nothing is dropped silently: every discarded token, assumption or ambiguity
is recorded in `notes`.

Changes in 0.6.0
----------------
* FIX  indel/insertion drafts were biologically wrong (TGAA>T produced
       `delinsT`). Alleles are now trimmed VCF-style before drafting.
* FIX  `del\\b` / `dup\\b` matched inside words ("Mendelian", "model", "CAMP")
       and "loss of function" was classified as a CNV loss.
* FIX  coordinates below 100 bp were silently dropped (`\\d{3,}`), which lost
       every mitochondrial position such as m.73 and m.152.
* FIX  accession digits (NM_001173464.2) were harvested as coordinates.
* FIX  REVIEW bucket wrote a header from the first record only, so mixed-type
       rows were column-shifted.
* FIX  `AssemblyMapper.g_to_c()` was called without a transcript accession and
       always raised, silently falling back to `validated_g`.
* FIX  `replace_reference=True` silently rewrote wrong reference bases.
* NEW  chromosome bound checks, Mb/bp size units, bare `12:123-456` loci,
       copy-number extraction, line numbers, run provenance hashes.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import logging
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from typing import ClassVar, Dict, List, Optional, Sequence, Tuple

__version__ = "0.6.0"

ASSEMBLY = "GRCh38"
SV_MIN_BP = 50              # below this an event is an indel, not an SV
CNV_MIN_BP = 1_000          # DEL/DUP at or above this are treated as CNV
VALIDATION_MAX_BP = 100_000  # larger regions are skipped (slow/unsupported)
SHIFT_MAX_ITER = 10_000     # guard against runaway 3' shifting

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
    "M": "NC_012920.1",
}

# GRCh38 primary assembly lengths - used to reject impossible coordinates.
CHROM_LENGTHS_GRCH38: Dict[str, int] = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
    "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
    "M": 16569,
}

NUCLEOTIDES = frozenset("ACGTN")


def norm_chrom(chrom: Optional[str]) -> Optional[str]:
    """'chr12' / 'Chr 12' / 'MT' -> canonical key ('12', 'X', 'M')."""
    if not chrom:
        return None
    key = re.sub(r"^chr[\s._]*", "", str(chrom).strip(), flags=re.I).upper()
    if key == "MT":
        key = "M"
    return key or None


def lookup_accession(chrom: Optional[str]) -> Optional[str]:
    return REFSEQ_GRCH38.get(norm_chrom(chrom) or "")


def chrom_length(chrom: Optional[str]) -> Optional[int]:
    return CHROM_LENGTHS_GRCH38.get(norm_chrom(chrom) or "")


# ===================================================================
# Allele handling - the core biological correctness fix
# ===================================================================
def trim_alleles(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """
    Reduce VCF-style left-anchored alleles to their minimal representation.

    HGVS describes only the changed bases, VCF carries a shared anchor base.
    'TGAA' > 'T' at 43624673 is a three-base deletion starting at 43624674,
    not a delins. Suffix is trimmed first, then prefix (bcftools order).
    """
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def shift_deletion_3prime(fetch, chrom: str, pos: int, ref: str) -> Tuple[int, str]:
    """Move a deletion to the most 3' position, as HGVS requires."""
    length = len(ref)
    for _ in range(SHIFT_MAX_ITER):
        nxt = fetch(chrom, pos + length, pos + length)
        if not nxt or nxt.upper() != ref[0].upper():
            break
        ref = ref[1:] + nxt.upper()
        pos += 1
    return pos, ref


def shift_insertion_3prime(fetch, chrom: str, pos: int, alt: str) -> Tuple[int, str]:
    """Move an insertion (placed after `pos`) to the most 3' position."""
    for _ in range(SHIFT_MAX_ITER):
        nxt = fetch(chrom, pos + 1, pos + 1)
        if not nxt or nxt.upper() != alt[0].upper():
            break
        alt = alt[1:] + nxt.upper()
        pos += 1
    return pos, alt


# ===================================================================
# Stage 2: schema - common core plus type-specific fields
# ===================================================================
@dataclass
class VariantRecord:
    """Base record. Fields shared by every variant type."""

    raw: str = ""
    line_no: Optional[int] = None
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

    # ClassVar so these never leak into asdict()/as_row()
    COMMON_COLUMNS: ClassVar[Sequence[str]] = (
        "line_no", "assembly", "chrom", "start", "accession", "vtype",
        "hgvs", "hgvs_draft", "validation_status", "status", "notes",
    )
    TYPE_COLUMNS: ClassVar[Sequence[str]] = ()
    REQUIRES_LOCUS: ClassVar[bool] = True

    def check(self) -> bool:
        if self.REQUIRES_LOCUS:
            if not self.chrom:
                self.notes.append("missing chromosome")
            if self.start is None:
                self.notes.append("missing start position")
            if self.chrom and not self.accession:
                self.accession = lookup_accession(self.chrom)
                if not self.accession:
                    self.notes.append(f"no RefSeq accession for {self.chrom}")
            self._check_bounds()
        self._check_type_specific()
        if self.notes:
            self.status = "needs_review"
        return self.status == "ok"

    def _check_bounds(self) -> None:
        limit = chrom_length(self.chrom)
        for label in ("start", "end"):
            value = getattr(self, label, None)
            if value is None:
                continue
            if value < 1:
                self.notes.append(f"{label} must be >= 1 (got {value})")
            elif limit and value > limit:
                self.notes.append(
                    f"{label} {value} exceeds {self.chrom} length {limit} on {ASSEMBLY}")

    def _check_type_specific(self) -> None:
        """Overridden by subclasses."""

    def columns(self) -> List[str]:
        return list(self.COMMON_COLUMNS) + list(self.TYPE_COLUMNS) + ["raw"]

    def as_dict(self) -> Dict[str, str]:
        data = asdict(self)
        data["notes"] = "; ".join(self.notes)
        return {k: ("" if v is None else str(v)) for k, v in data.items()}


@dataclass
class HgvsInputRecord(VariantRecord):
    """Input that already carries a valid-looking HGVS string."""

    REQUIRES_LOCUS: ClassVar[bool] = False

    def __post_init__(self):
        self.vtype = "HGVS_INPUT"

    def _check_type_specific(self):
        if not self.hgvs_draft:
            self.notes.append("no HGVS string extracted")
        else:
            self.notes.append("input already in HGVS notation; not re-parsed")
            # a note alone should not force review for this type
            self.status = "ok"

    def check(self) -> bool:
        had = bool(self.hgvs_draft)
        super().check()
        self.status = "ok" if had else "needs_review"
        return self.status == "ok"


@dataclass
class SNVRecord(VariantRecord):
    ref: Optional[str] = None
    alt: Optional[str] = None

    TYPE_COLUMNS: ClassVar[Sequence[str]] = ("ref", "alt")

    def __post_init__(self):
        self.vtype = "SNV"

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

    TYPE_COLUMNS: ClassVar[Sequence[str]] = ("ref", "alt", "size")

    def __post_init__(self):
        self.vtype = "INDEL"

    def _check_type_specific(self):
        if not self.ref and not self.alt:
            self.notes.append("indel requires at least one allele")
            return
        for label, allele in (("ref", self.ref), ("alt", self.alt)):
            if allele and set(allele.upper()) - NUCLEOTIDES:
                self.notes.append(f"invalid nucleotide character in {label}")
        if self.ref and self.alt and self.ref.upper() == self.alt.upper():
            self.notes.append("ref equals alt")
        if self.size is None:
            self.size = abs(len(self.ref or "") - len(self.alt or ""))


VALID_SVTYPES = frozenset({"DEL", "DUP", "INV", "INS"})


@dataclass
class SVRecord(VariantRecord):
    end: Optional[int] = None
    svtype: Optional[str] = None
    size: Optional[int] = None

    TYPE_COLUMNS: ClassVar[Sequence[str]] = ("end", "svtype", "size")

    def __post_init__(self):
        self.vtype = "SV"

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
        if self.size is not None and self.size < SV_MIN_BP:
            self.notes.append(
                f"size {self.size}bp is below SV_MIN_BP={SV_MIN_BP}; "
                "allele sequences are needed to describe it as an indel")


@dataclass
class CNVRecord(VariantRecord):
    end: Optional[int] = None
    cnv_type: Optional[str] = None
    copy_number: Optional[int] = None
    size: Optional[int] = None

    TYPE_COLUMNS: ClassVar[Sequence[str]] = ("end", "cnv_type", "copy_number", "size")

    def __post_init__(self):
        self.vtype = "CNV"

    def _check_type_specific(self):
        if self.end is None:
            self.notes.append("CNV requires an end coordinate")
        if not self.cnv_type:
            self.notes.append("CNV requires a direction (loss/gain)")
        elif self.cnv_type not in ("copy_number_loss", "copy_number_gain"):
            self.notes.append(f"invalid cnv_type: {self.cnv_type}")
        if self.copy_number is not None and self.copy_number < 0:
            self.notes.append("copy number must be >= 0")
        if self.start is not None and self.end is not None:
            if self.end < self.start:
                self.notes.append("end precedes start")
            elif self.size is None:
                self.size = self.end - self.start


@dataclass
class RepeatRecord(VariantRecord):
    repeat_unit: Optional[str] = None
    copy_number: Optional[int] = None

    TYPE_COLUMNS: ClassVar[Sequence[str]] = ("repeat_unit", "copy_number")

    def __post_init__(self):
        self.vtype = "REPEAT"

    def _check_type_specific(self):
        if not self.repeat_unit:
            self.notes.append("repeat requires a motif")
        elif set(self.repeat_unit.upper()) - NUCLEOTIDES:
            self.notes.append("invalid repeat motif")
        if self.copy_number is None:
            self.notes.append("repeat requires a copy number")


# ===================================================================
# Stage 1: parser
# ===================================================================
CHROM_TOKEN = r"(?:2[0-2]|1\d|[1-9]|X|Y|MT|M)"
NUM = r"\d{1,3}(?:[.,]\d{3})+|\d+"


class VariantParser:
    """Extracts structured fields from free-text variant descriptions."""

    # --- input that is already HGVS ---------------------------------
    HGVS_FULL = re.compile(
        r"\b((?:N[CMRGTPW]_\d+(?:\.\d+)?|ENS[TGP]\d+(?:\.\d+)?|LRG_\d+(?:t\d+)?)"
        r"\s*:\s*[cgmnopr]\.\S+)", re.I)

    # --- locus: coordinates must be anchored to a chromosome --------
    LOCUS_CHR = re.compile(
        rf"\bchr[\s._]*(?P<chrom>{CHROM_TOKEN})\b\s*[:\-_\s]\s*(?P<start>{NUM})"
        rf"(?:\s*[-\u2013\u2014_]\s*(?P<end>{NUM}))?", re.I)
    LOCUS_BARE = re.compile(
        rf"(?:^|[\s|;,(])(?P<chrom>{CHROM_TOKEN})\s*:\s*(?P<start>{NUM})"
        rf"(?:\s*[-\u2013\u2014_]\s*(?P<end>{NUM}))?", re.I)
    CHROM_ONLY = re.compile(rf"\bchr[\s._]*({CHROM_TOKEN})\b", re.I)

    # --- separators --------------------------------------------------
    THOUSANDS = re.compile(r"^\d{1,3}(?:[.,]\d{3})+$")
    AMBIGUOUS_SEP = re.compile(r"^\d{1,3}[.,]\d{3}$")

    # --- alleles -----------------------------------------------------
    SUBSTITUTION = re.compile(
        r"(?<![A-Za-z])([ACGTN]+)\s*(?:>|->|\u2192)\s*([ACGTN]+)(?![A-Za-z])", re.I)
    # slash form excludes N so that "N/A" cannot be read as a substitution
    SUBSTITUTION_SLASH = re.compile(
        r"(?<![A-Za-z])([ACGT]+)\s*/\s*([ACGT]+)(?![A-Za-z])", re.I)
    DELINS = re.compile(
        r"(?<![A-Za-z])del(?:etion)?\s*([ACGTN]*)\s*ins(?:ertion)?\s*([ACGTN]+)"
        r"(?![A-Za-z])", re.I)
    INSERTION = re.compile(
        r"(?<![A-Za-z])ins(?:ertion)?\s*([ACGTN]+)(?![A-Za-z])", re.I)

    # --- repeats -----------------------------------------------------
    REPEAT_PAREN = re.compile(
        r"\(\s*([ACGT]{2,10})\s*\)\s*(?:n\b|[x*\u00d7]?\s*(\d+))", re.I)
    REPEAT_BRACKET = re.compile(
        r"(?<![A-Za-z])([ACGT]{2,10})\s*[(\[]\s*(\d+)\s*[)\]]", re.I)
    REPEAT_TIMES = re.compile(
        r"(?<![A-Za-z])([ACGT]{2,10})\s*[x*\u00d7]\s*(\d+)\b", re.I)

    # --- dosage direction --------------------------------------------
    # stripped before classification so that "loss of function" is not a CNV
    LOF_GOF = re.compile(
        r"\b(?:loss|gain)[\s\-]?of[\s\-]?function\b|\b[LG]OF\b", re.I)
    LOSS = re.compile(
        r"(?<![A-Za-z])(?:del(?:etion|eted)?|loss)(?![A-Za-z])", re.I)
    GAIN = re.compile(
        r"(?<![A-Za-z])(?:dup(?:lication|licated)?|gain|amplification)(?![A-Za-z])",
        re.I)
    INVERSION = re.compile(
        r"(?<![A-Za-z])inv(?:ersion|erted)?(?![A-Za-z])", re.I)

    # --- sizes and copy number ---------------------------------------
    SIZE_UNIT = re.compile(r"([\d.,]+)\s*(kbp|kb|mbp|mb|bp)\b", re.I)
    MULTIPLIER = {"bp": 1, "kb": 1_000, "kbp": 1_000,
                  "mb": 1_000_000, "mbp": 1_000_000}
    COPY_NUMBER = re.compile(
        r"\b(?:cn|copy[\s\-]*number)\s*[:=]?\s*(\d+)\b"
        r"|(?<![A-Za-z\d])x\s*(\d+)\b", re.I)

    # ---------------------------------------------------------------
    @classmethod
    def parse_coordinate(cls, token: str) -> Tuple[Optional[int], Optional[str]]:
        """Genomic coordinates are integers; separators are thousands marks."""
        if cls.THOUSANDS.match(token):
            value = int(token.replace(".", "").replace(",", ""))
            if cls.AMBIGUOUS_SEP.match(token):
                return value, (f"ambiguous separator in '{token}', "
                               f"assumed thousands -> {value}")
            return value, None
        if token.isdigit():
            return int(token), None
        return None, f"unparseable coordinate token: '{token}'"

    @classmethod
    def parse_magnitude(cls, token: str) -> Optional[float]:
        """Size numbers may legitimately be decimal (18.6kb, 1,5kb)."""
        if cls.THOUSANDS.match(token):
            return float(token.replace(".", "").replace(",", ""))
        try:
            return float(token.replace(",", "."))
        except ValueError:
            return None

    def _text_size(self, text: str) -> Optional[int]:
        match = self.SIZE_UNIT.search(text)
        if not match:
            return None
        value = self.parse_magnitude(match.group(1))
        if value is None:
            return None
        return int(round(value * self.MULTIPLIER[match.group(2).lower()]))

    def _size_and_mismatch(self, text: str, start: Optional[int], end: Optional[int]
                           ) -> Tuple[Optional[int], Optional[str]]:
        """Coordinates are authoritative; the prose size is cross-validation."""
        coord_size = abs(end - start) if (start is not None and end is not None) else None
        text_size = self._text_size(text)

        if coord_size is not None:
            if text_size is not None and abs(coord_size - text_size) > coord_size * 0.1:
                return coord_size, (f"size mismatch: coordinates imply {coord_size}bp, "
                                    f"text implies {text_size}bp; using coordinates")
            return coord_size, None
        if text_size is not None:
            return text_size, "size taken from prose only; no end coordinate"
        return None, None

    def _copy_number(self, text: str) -> Optional[int]:
        match = self.COPY_NUMBER.search(text)
        if not match:
            return None
        raw = match.group(1) or match.group(2)
        return int(raw) if raw else None

    # ---------------------------------------------------------------
    def parse(self, line: str, line_no: Optional[int] = None) -> Optional[VariantRecord]:
        text = line.strip()
        if not text or text.startswith("#"):
            return None

        notes: List[str] = []

        # Stage 0 - already HGVS? do not re-parse, forward to the validator.
        hgvs_match = self.HGVS_FULL.search(text)
        if hgvs_match:
            record = HgvsInputRecord(raw=line, line_no=line_no)
            record.hgvs_draft = re.sub(r"\s+", "", hgvs_match.group(1)).rstrip(".,;")
            record.check()
            return record

        # Locus - a coordinate only counts if it is anchored to a chromosome.
        chrom = start = end = None
        locus = self.LOCUS_CHR.search(text) or self.LOCUS_BARE.search(text)
        if locus:
            chrom = "chr" + (norm_chrom(locus.group("chrom")) or "")
            start, note = self.parse_coordinate(locus.group("start"))
            if note:
                notes.append(note)
            if locus.group("end"):
                end, note = self.parse_coordinate(locus.group("end"))
                if note:
                    notes.append(note)
        else:
            chrom_only = self.CHROM_ONLY.search(text)
            if chrom_only:
                chrom = "chr" + (norm_chrom(chrom_only.group(1)) or "")
            notes.append("no chromosome-anchored coordinate found")

        size, mismatch = self._size_and_mismatch(text, start, end)
        if mismatch:
            notes.append(mismatch)

        common = dict(raw=line, line_no=line_no, chrom=chrom, start=start)
        record: Optional[VariantRecord] = None

        # Repeats first: CAG(42) must not be read as an allele.
        repeat = (self.REPEAT_BRACKET.search(text)
                  or self.REPEAT_TIMES.search(text)
                  or self.REPEAT_PAREN.search(text))
        if repeat:
            groups = repeat.groups()
            copies = groups[1] if len(groups) > 1 else None
            record = RepeatRecord(
                repeat_unit=groups[0].upper(),
                copy_number=int(copies) if copies else None, **common)
            if copies is None:
                notes.append("repeat count given as 'n'; not a specific allele")

        elif substitution := (self.SUBSTITUTION.search(text)
                              or self.SUBSTITUTION_SLASH.search(text)):
            ref = substitution.group(1).upper()
            alt = substitution.group(2).upper()
            record = (SNVRecord(ref=ref, alt=alt, **common)
                      if len(ref) == 1 and len(alt) == 1
                      else IndelRecord(ref=ref, alt=alt, **common))

        elif delins := self.DELINS.search(text):
            deleted = delins.group(1).upper()
            inserted = delins.group(2).upper()
            if not deleted and end is None:
                # length of the deleted stretch is unknown - do not guess
                record = VariantRecord(**common)
                record.status = "needs_review"
                notes.append("delins with unknown deleted allele and no end "
                             "coordinate; cannot build an unambiguous draft")
                record.accession = lookup_accession(chrom)
            else:
                if not deleted and end is not None:
                    deleted = "N" * (end - (start or end) + 1)
                    notes.append("deleted allele inferred from the coordinate range")
                record = IndelRecord(ref=deleted, alt=inserted, **common)

        elif insertion := self.INSERTION.search(text):
            record = IndelRecord(ref="", alt=insertion.group(1).upper(), **common)

        else:
            # Strip loss/gain-of-function before asking about dosage.
            dosage_text = self.LOF_GOF.sub(" ", text)
            if self.LOF_GOF.search(text):
                notes.append("loss/gain-of-function wording ignored for dosage typing")

            if self.LOSS.search(dosage_text):
                direction, svtype = "copy_number_loss", "DEL"
            elif self.GAIN.search(dosage_text):
                direction, svtype = "copy_number_gain", "DUP"
            elif self.INVERSION.search(dosage_text):
                direction, svtype = None, "INV"
            else:
                direction = svtype = None
                record = VariantRecord(**common)
                record.status = "needs_review"
                notes.append("variant type could not be determined")
                record.accession = lookup_accession(chrom)

            if record is None:
                if svtype in ("DEL", "DUP") and size is not None and size >= CNV_MIN_BP:
                    record = CNVRecord(end=end, cnv_type=direction, size=size,
                                       copy_number=self._copy_number(text), **common)
                else:
                    record = SVRecord(end=end, svtype=svtype, size=size, **common)

        if record is not None:
            record.notes.extend(notes)
            record.check()
        return record


# ===================================================================
# Stage 3: draft builder
# ===================================================================
class HgvsDraftBuilder:
    """Builds an UNVALIDATED genomic HGVS candidate from a record."""

    def __init__(self, fetch=None):
        # fetch(chrom, start, end) -> reference sequence, 1-based inclusive
        self.fetch = fetch

    def build(self, record: VariantRecord) -> Optional[str]:
        if isinstance(record, HgvsInputRecord):
            return record.hgvs_draft
        acc, pos = record.accession, record.start
        if not acc or pos is None:
            return None

        vtype = record.vtype
        if vtype in ("SNV", "INDEL"):
            return self._build_allele(record, acc, pos)
        if vtype == "REPEAT" and record.repeat_unit and record.copy_number is not None:
            end = pos + len(record.repeat_unit) - 1
            return f"{acc}:g.{pos}_{end}{record.repeat_unit}[{record.copy_number}]"
        if vtype == "CNV" and record.end:
            op = "dup" if record.cnv_type == "copy_number_gain" else "del"
            return f"{acc}:g.{pos}_{record.end}{op}"
        if vtype == "SV" and record.end:
            op = {"DEL": "del", "DUP": "dup", "INV": "inv"}.get(record.svtype)
            if op:
                return f"{acc}:g.{pos}_{record.end}{op}"
        return None

    def _build_allele(self, record, acc: str, pos: int) -> Optional[str]:
        ref, alt = (record.ref or ""), (record.alt or "")
        pos, ref, alt = trim_alleles(pos, ref, alt)

        if not ref and not alt:
            record.notes.append("ref equals alt after trimming; no change described")
            record.status = "needs_review"
            return None

        # Pure deletion
        if ref and not alt:
            if self.fetch:
                pos, ref = self._shift(record, shift_deletion_3prime, pos, ref)
            end = pos + len(ref) - 1
            return f"{acc}:g.{pos}del" if len(ref) == 1 else f"{acc}:g.{pos}_{end}del"

        # Pure insertion, placed between pos-1 and pos after trimming
        if alt and not ref:
            anchor = pos - 1
            if self.fetch:
                anchor, alt = self._shift(record, shift_insertion_3prime, anchor, alt)
            return f"{acc}:g.{anchor}_{anchor + 1}ins{alt}"

        if len(ref) == 1 and len(alt) == 1:
            return f"{acc}:g.{pos}{ref}>{alt}"
        if len(ref) == 1:
            return f"{acc}:g.{pos}delins{alt}"
        end = pos + len(ref) - 1
        return f"{acc}:g.{pos}_{end}delins{alt}"

    def _shift(self, record, shifter, pos: int, allele: str):
        try:
            new_pos, new_allele = shifter(self.fetch, record.chrom, pos, allele)
        except Exception as exc:                       # noqa: BLE001
            record.notes.append(f"3' shift skipped ({type(exc).__name__})")
            return pos, allele
        if new_pos != pos:
            record.notes.append(f"3'-shifted {new_pos - pos}bp per HGVS rules")
        return new_pos, new_allele


# ===================================================================
# Reference sequence access (optional, offline)
# ===================================================================
class FastaReference:
    """Thin pyfaidx wrapper: reference-base checking and 3' shifting."""

    def __init__(self, path: str):
        from pyfaidx import Fasta
        self._fasta = Fasta(path, as_raw=True, sequence_always_upper=True)
        self._keys = set(self._fasta.keys())
        self.name = os.path.basename(path)

    def _contig(self, chrom: str) -> Optional[str]:
        key = norm_chrom(chrom)
        for candidate in (f"chr{key}", key, f"chrM" if key == "M" else None,
                          "MT" if key == "M" else None):
            if candidate and candidate in self._keys:
                return candidate
        return None

    def fetch(self, chrom: str, start: int, end: int) -> Optional[str]:
        contig = self._contig(chrom)
        if not contig or start < 1:
            return None
        return self._fasta[contig][start - 1:end] or None

    def check_reference(self, record: VariantRecord) -> None:
        ref = getattr(record, "ref", None)
        if not ref or record.start is None or set(ref) - NUCLEOTIDES:
            return
        observed = self.fetch(record.chrom, record.start, record.start + len(ref) - 1)
        if observed is None:
            record.notes.append("reference sequence unavailable for this locus")
        elif observed.upper() != ref.upper():
            record.notes.append(
                f"REFERENCE MISMATCH: {record.chrom}:{record.start} is "
                f"'{observed}' but the input claims '{ref}'")
            record.status = "needs_review"


# ===================================================================
# Stage 4: validators (interchangeable)
# ===================================================================
class NullValidator:
    name = "disabled"
    data_version = "n/a"

    def validate(self, record: VariantRecord) -> None:
        record.validation_status = "not_validated"


class FastaValidator:
    """Offline validator: reference bases only, no transcript mapping."""

    name = "fasta"

    def __init__(self, reference: FastaReference):
        self.reference = reference
        self.data_version = reference.name

    def validate(self, record: VariantRecord) -> None:
        if not record.hgvs_draft:
            record.validation_status = "no_draft"
            return
        before = len(record.notes)
        self.reference.check_reference(record)
        mismatched = any("REFERENCE MISMATCH" in n for n in record.notes[before:])
        record.validation_status = "failed_reference_check" if mismatched \
            else "reference_checked"


class UtaValidator:
    """Full validator: biocommons hgvs + UTA, with 3' normalization."""

    name = "uta"

    def __init__(self, mane: Optional[Dict[str, str]] = None) -> None:
        import hgvs.parser
        import hgvs.dataproviders.uta
        import hgvs.assemblymapper
        import hgvs.validator
        import hgvs.normalizer

        self._parser = hgvs.parser.Parser()
        self._hdp = hgvs.dataproviders.uta.connect()
        self.data_version = self._hdp.data_version()
        self._mapper = hgvs.assemblymapper.AssemblyMapper(
            self._hdp, assembly_name=ASSEMBLY,
            alt_aln_method="splign",
            # never silently rewrite a wrong reference base
            replace_reference=False,
        )
        self._normalizer = hgvs.normalizer.Normalizer(
            self._hdp, shuffle_direction=3, cross_boundaries=False)
        self._validator = hgvs.validator.Validator(self._hdp)
        self._mane = mane or {}

    def _pick_transcript(self, transcripts: Sequence[str]) -> Optional[str]:
        """Prefer MANE Select when a mapping is supplied, else curated RefSeq."""
        if not transcripts:
            return None
        for tx in transcripts:
            if tx in self._mane.values():
                return tx
        curated = sorted((t for t in transcripts if t.startswith("NM_")), reverse=True)
        return curated[0] if curated else sorted(transcripts)[0]

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
        except Exception as exc:                       # noqa: BLE001
            record.validation_status = "unparseable"
            record.notes.append(f"parse error: {type(exc).__name__}: {str(exc)[:80]}")
            return

        try:
            parsed = self._normalizer.normalize(parsed)
        except Exception as exc:                       # noqa: BLE001
            record.notes.append(f"3' normalization skipped: {type(exc).__name__}")

        try:
            self._validator.validate(parsed)
        except Exception as exc:                       # noqa: BLE001
            record.validation_status = "failed_validation"
            record.notes.append(f"validation error: {str(exc)[:120]}")
            return

        record.hgvs = str(parsed)
        record.validation_status = "validated_g"

        # transcript projection - the transcript actually used must be recorded
        try:
            transcripts = self._mapper.relevant_transcripts(parsed)
        except Exception as exc:                       # noqa: BLE001
            record.notes.append(f"transcript lookup failed: {type(exc).__name__}")
            return
        if not transcripts:
            record.notes.append("no overlapping transcript; genomic level only")
            return

        chosen = self._pick_transcript(transcripts)
        try:
            record.hgvs = str(self._mapper.g_to_c(parsed, chosen))
            record.validation_status = "validated_c"
            record.notes.append(f"transcript={chosen}")
            if len(transcripts) > 1:
                record.notes.append(f"{len(transcripts)} transcripts overlap; "
                                    "selection is not MANE-certified")
        except Exception as exc:                       # noqa: BLE001
            record.notes.append(
                f"g_to_c failed on {chosen}: {type(exc).__name__}: {str(exc)[:80]}")


# ===================================================================
# Orchestration
# ===================================================================
class NormalizationPipeline:
    def __init__(self, parser=None, builder=None, validator=None, reference=None):
        self.parser = parser or VariantParser()
        self.reference = reference
        self.builder = builder or HgvsDraftBuilder(
            fetch=reference.fetch if reference else None)
        self.validator = validator or NullValidator()
        self.failed_records: List[Dict[str, str]] = []

    def run(self, input_path: str) -> List[VariantRecord]:
        records: List[VariantRecord] = []
        self.failed_records = []

        # utf-8-sig strips the Excel BOM; errors=replace keeps a bad byte from
        # killing the whole run - the affected line is reported instead.
        with open(input_path, encoding="utf-8-sig", errors="replace") as handle:
            for line_no, line in enumerate(handle, 1):
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                if "\ufffd" in text:
                    self.failed_records.append(
                        {"line_no": str(line_no), "raw_line": text,
                         "error": "undecodable byte; check the file encoding"})
                    continue
                try:
                    record = self.parser.parse(line, line_no=line_no)
                    if record is None:
                        continue
                    if record.status == "ok":
                        record.hgvs_draft = self.builder.build(record)
                        if record.hgvs_draft is None and record.status == "ok":
                            record.notes.append("no draft could be built")
                            record.status = "needs_review"
                    records.append(record)
                except Exception as exc:               # noqa: BLE001
                    self.failed_records.append(
                        {"line_no": str(line_no), "raw_line": text,
                         "error": f"{type(exc).__name__}: {exc}"})

        logger.info("parsed %d records (%d failed)", len(records), len(self.failed_records))

        for index, record in enumerate(records, 1):
            self.validator.validate(record)
            logger.debug("[%d/%d] %s -> %s", index, len(records),
                         record.raw.strip()[:40], record.validation_status)
        return records

    @staticmethod
    def _header_for(group: Sequence[VariantRecord]) -> List[str]:
        """Union of columns across a bucket - REVIEW holds mixed record types."""
        type_cols: List[str] = []
        for record in group:
            for column in record.TYPE_COLUMNS:
                if column not in type_cols:
                    type_cols.append(column)
        return (list(VariantRecord.COMMON_COLUMNS) + type_cols + ["raw"])

    # Every file this tool may emit. Stale copies from an earlier run must be
    # cleared, otherwise a record reclassified into REVIEW still appears in the
    # previous run's snv.tsv and the output directory contradicts itself.
    OWNED_OUTPUTS: ClassVar[Sequence[str]] = (
        "snv.tsv", "indel.tsv", "sv.tsv", "cnv.tsv", "repeat.tsv",
        "hgvs_input.tsv", "review.tsv", "failed_records.tsv", "run_manifest.txt",
    )

    def _clear_stale_outputs(self, output_dir: str) -> None:
        for name in self.OWNED_OUTPUTS:
            path = os.path.join(output_dir, name)
            if os.path.exists(path):
                os.remove(path)
                logger.debug("removed stale output %s", path)

    def write(self, records, output_dir: str, input_path: str
              ) -> Dict[str, List[VariantRecord]]:
        os.makedirs(output_dir, exist_ok=True)
        self._clear_stale_outputs(output_dir)
        buckets: Dict[str, List[VariantRecord]] = {}
        for record in records:
            key = record.vtype if (record.status == "ok" and record.vtype) else "REVIEW"
            buckets.setdefault(key, []).append(record)

        for key, group in buckets.items():
            header = self._header_for(group)
            path = os.path.join(output_dir, f"{key.lower()}.tsv")
            with open(path, "w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(header)
                for record in group:
                    data = record.as_dict()
                    writer.writerow([data.get(column, "") for column in header])

        if self.failed_records:
            path = os.path.join(output_dir, "failed_records.tsv")
            with open(path, "w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["line_no", "raw_line", "error_message"])
                for failure in self.failed_records:
                    writer.writerow([failure["line_no"], failure["raw_line"],
                                     failure["error"]])

        self._write_manifest(buckets, output_dir, input_path)
        return buckets

    def _write_manifest(self, buckets, output_dir: str, input_path: str) -> None:
        path = os.path.join(output_dir, "run_manifest.txt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(f"tool_version: {__version__}\n")
            handle.write("timestamp: %s\n"
                         % datetime.datetime.now().isoformat(timespec="seconds"))
            handle.write(f"git_commit: {_git_commit()}\n")
            handle.write(f"python: {sys.version.split()[0]}\n")
            handle.write(f"hgvs_library: {_library_version('hgvs')}\n")
            handle.write(f"assembly: {ASSEMBLY}\n")
            handle.write(f"input_file: {os.path.abspath(input_path)}\n")
            handle.write(f"input_sha256: {_sha256(input_path)}\n")
            handle.write(f"sv_min_bp: {SV_MIN_BP}\n")
            handle.write(f"cnv_min_bp: {CNV_MIN_BP}\n")
            handle.write(f"validator: {self.validator.name}\n")
            handle.write(f"validator_data_version: {self.validator.data_version}\n")
            handle.write("reference_fasta: %s\n"
                         % (self.reference.name if self.reference else "none"))
            for key, group in sorted(buckets.items()):
                handle.write(f"{key}: {len(group)}\n")
            handle.write(f"FAILED: {len(self.failed_records)}\n")


def _sha256(path: str) -> str:
    try:
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(65536), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return "unavailable"


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=os.path.dirname(os.path.abspath(__file__)),
            stderr=subprocess.DEVNULL, text=True).strip() or "unavailable"
    except Exception:                                  # noqa: BLE001
        return "unavailable"


def _library_version(name: str) -> str:
    try:
        from importlib.metadata import version
        return version(name)
    except Exception:                                  # noqa: BLE001
        return "not installed"


def load_mane(path: Optional[str]) -> Dict[str, str]:
    """Two-column TSV: gene<TAB>MANE Select transcript accession."""
    if not path:
        return {}
    mapping: Dict[str, str] = {}
    with open(path, encoding="utf-8-sig") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and not line.startswith("#"):
                mapping[parts[0].strip()] = parts[1].strip()
    return mapping


# ===================================================================
# Self-test
# ===================================================================
# (text, expected_vtype, expected_status, expected_draft_or_None)
TEST_CASES = [
    ("chr11-45917701 T>C", "SNV", "ok", "NC_000011.10:g.45917701T>C"),
    ("chrM-8612 T>C", "SNV", "ok", "NC_012920.1:g.8612T>C"),
    # 1-2 digit coordinates used to be dropped by \d{3,}
    ("chrM-73 A>G", "SNV", "ok", "NC_012920.1:g.73A>G"),
    ("chrM-152 T>C", "SNV", "ok", "NC_012920.1:g.152T>C"),
    # left-anchored VCF alleles must become a real deletion / insertion
    ("chr6-43624673 TGAA>T", "INDEL", "ok", "NC_000006.12:g.43624674_43624676del"),
    ("chr22-50221116 T>TCCG", "INDEL", "ok",
     "NC_000022.11:g.50221116_50221117insCCG"),
    ("chr1-100 AT>GC", "INDEL", "ok", "NC_000001.11:g.100_101delinsGC"),
    ("CNTN1:chr12:41.021.576-41.040.185 del 18.6kb", "CNV", "ok",
     "NC_000012.12:g.41021576_41040185del"),
    ("CASK DEL:chrX:41772635-41895004", "CNV", "ok",
     "NC_000023.11:g.41772635_41895004del"),
    ("chr5-126,489,209-126,759,255", None, "needs_review", None),
    ("HTT chr4:3074877 CAG(42)", "REPEAT", "ok",
     "NC_000004.12:g.3074877_3074879CAG[42]"),
    ("chr3:1000-5000 inversion", "SV", "ok", "NC_000003.12:g.1000_5000inv"),
    # prose size disagrees with the coordinates -> flagged, not silently used
    ("chr9:100.000-200.000 del 1,5kb", "CNV", "needs_review", None),
    # word-boundary regressions
    ("Mendelian inheritance pattern in this family", None, "needs_review", None),
    ("model-based caller output", None, "needs_review", None),
    ("chr7:100000-200000 loss of function variant", None, "needs_review", None),
    # bare chromosome, no chr prefix
    ("12:41021576-41040185 del", "CNV", "ok", "NC_000012.12:g.41021576_41040185del"),
    # 1 Mb duplication
    ("chr17:1000000-2000000 duplication 1Mb", "CNV", "ok",
     "NC_000017.11:g.1000000_2000000dup"),
    # already HGVS -> passthrough, never re-parsed
    ("NM_001173464.2:c.1674-1G>A", "HGVS_INPUT", "ok", "NM_001173464.2:c.1674-1G>A"),
    # impossible coordinate for the mitochondrial genome
    ("chrM:45917701 T>C", "SNV", "needs_review", None),
    ("chr45917701 T>C", "SNV", "needs_review", None),
    # accession digits must not become coordinates
    ("variant in exon 5", None, "needs_review", None),
]


def run_self_test() -> bool:
    parser = VariantParser()
    builder = HgvsDraftBuilder()
    passed = failed = 0
    for text, expected_type, expected_status, expected_draft in TEST_CASES:
        record = parser.parse(text, line_no=0)
        actual_type = record.vtype if record else None
        actual_status = record.status if record else None
        actual_draft = None
        if record and record.status == "ok":
            actual_draft = builder.build(record)
            actual_status = record.status

        ok = (actual_type == expected_type and actual_status == expected_status)
        if ok and expected_draft is not None:
            ok = actual_draft == expected_draft

        if ok:
            passed += 1
            print(f"  PASS  {text[:44]:<44} {actual_type}")
        else:
            failed += 1
            print(f"  FAIL  {text[:44]:<44}")
            print(f"        expected {expected_type}/{expected_status}/{expected_draft}")
            print(f"        got      {actual_type}/{actual_status}/{actual_draft}")
            if record and record.notes:
                print(f"        notes    {'; '.join(record.notes)[:110]}")
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
    argp.add_argument("--fasta", help="GRCh38 FASTA for reference checking and "
                                      "3' shifting (offline, recommended)")
    argp.add_argument("--validate", action="store_true",
                      help="validate against UTA (needs network, slow)")
    argp.add_argument("--mane", help="TSV of gene<TAB>MANE Select transcript")
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

    reference = None
    if args.fasta:
        try:
            reference = FastaReference(args.fasta)
            logger.info("reference loaded: %s", reference.name)
        except Exception as exc:                       # noqa: BLE001
            logger.warning("FASTA unavailable (%s); no reference checking",
                           type(exc).__name__)

    validator = NullValidator()
    if args.validate:
        logger.info("connecting to UTA (first connection can be slow)...")
        try:
            validator = UtaValidator(mane=load_mane(args.mane))
            logger.info("UTA ready: %s", validator.data_version)
        except Exception as exc:                       # noqa: BLE001
            logger.warning("UTA unavailable (%s); falling back", type(exc).__name__)
            if reference:
                validator = FastaValidator(reference)
    elif reference:
        validator = FastaValidator(reference)

    pipeline = NormalizationPipeline(validator=validator, reference=reference)
    records = pipeline.run(args.input)
    buckets = pipeline.write(records, args.output_dir, args.input)

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
        print(f"\n  FAILED  ({len(pipeline.failed_records)})  ->  "
              f"{args.output_dir}/failed_records.tsv")

    print("\n" + "=" * 66)
    print(f"  Provenance: {args.output_dir}/run_manifest.txt")


if __name__ == "__main__":
    main()

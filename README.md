# hgvs-normalizer

Turns messy, human-written variant descriptions into valid HGVS.

Clinical and literature sources rarely write variants the way software expects.
`chr12:41.021.576-41.040.185 del 18.6kb`, `chr6-43624673 TGAA>T` and
`HTT chr4:3074877 CAG(42)` all describe real variants, and none of them is
valid HGVS. This tool parses that kind of input, builds a candidate HGVS
description, and optionally verifies it against an external authority.

**Design principle:** an unvalidated string is never written to the `hgvs`
column. Candidates live in `hgvs_draft`, and every record carries an explicit
`validation_status`. Anything uncertain is routed to a review file rather than
guessed at.

## Pipeline

```mermaid
flowchart TD
    A[Free-text input line] --> B{Already valid HGVS?}
    B -- yes --> P[Stage 0: Passthrough]
    B -- no --> C[Stage 1: Parse<br/>locus, alleles, type]
    C --> D[Stage 2: Schema check<br/>required fields, bounds]
    D -- incomplete --> R[review.tsv]
    D -- complete --> E[Stage 3: Draft<br/>allele trimming, optional 3' shift]
    P --> F
    E --> F[Stage 4: Validate - OPTIONAL]
    F --> G[snv.tsv / indel.tsv / sv.tsv<br/>cnv.tsv / repeat.tsv]
    F --> H[run_manifest.txt]
```

Without a validator the pipeline still runs end to end; records stay
`not_validated`. This is deliberate: the core must not depend on a network
service.

## Record types

| Type | Recognised from | Draft form |
|---|---|---|
| `SNV` | single base on both sides | `g.123A>G` |
| `INDEL` | multi-base alleles, `delins`, `ins` | `g.123_125del` |
| `SV` | `del`/`dup`/`inv` below the CNV threshold | `g.1000_5000inv` |
| `CNV` | `del`/`dup` at or above `CNV_MIN_BP` | `g.1000_5000del` |
| `REPEAT` | `CAG(42)`, `CAG[42]`, `(CAG)n` | `g.100_102CAG[42]` |
| `HGVS_INPUT` | already-valid HGVS | passed through |

Mitochondrial accessions use the `m.` prefix rather than `g.`, per HGVS.

## Validation statuses

| Status | Meaning |
|---|---|
| `not_validated` | No validator enabled |
| `no_draft` | Incomplete input, no candidate built |
| `validated_g` | Confirmed at genomic level |
| `validated_c` | Also projected onto a transcript |
| `failed_validation` | Rejected, usually a reference base mismatch |
| `validator_unavailable` | Service unreachable; the draft is preserved |

## Usage

```bash
pip install -r requirements.txt

# core only - no network
python hgvs_normalizer.py --input examples/messy_variants.txt --output-dir output

# canonical form via Mutalyzer 3
python hgvs_normalizer.py --input examples/messy_variants.txt --mutalyzer

# transcript-level projection via UTA
export UTA_DB_URL="postgresql://anonymous:anonymous@uta.biocommons.org/uta/uta_20241220"
python hgvs_normalizer.py --input examples/messy_variants.txt --validate

# known-answer tests
python hgvs_normalizer.py --self-test
```

| Flag | Purpose |
|---|---|
| `--fasta` | GRCh38 FASTA for reference checking and 3' shifting |
| `--mutalyzer [URL]` | Validate with Mutalyzer 3 |
| `--validate` | Validate with UTA and project to `c.` |
| `--mane` | `gene<TAB>transcript` TSV for MANE Select choice |

Validators are chainable; the first conclusive verdict wins.

## Scope

This tool sits between two established layers and replaces neither. Mutalyzer,
VariantValidator and the `hgvs` package validate descriptions that are already
well formed; this tool calls them. tmVar and SETH extract variant mentions from
published literature and normalize to dbSNP identifiers; different input,
different output. The gap addressed here is narrower: in-house coordinate
tables written by people, converted to HGVS with an explicit audit trail.

## Requirements

Python 3.9+. Optional: `pyfaidx` for reference checking, `hgvs` + UTA access
for transcript projection, Docker for the containerised run.

## License

MIT

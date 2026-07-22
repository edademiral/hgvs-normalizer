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
guessed at — in the bundled example, 11 of 23 records go to review, each with
the reason recorded.

## Quick start

The core pipeline depends on **nothing but the Python standard library**.
There is no install step.

```bash
git clone https://github.com/edademiral/hgvs-normalizer.git
cd hgvs-normalizer

python hgvs_normalizer.py --self-test
python hgvs_normalizer.py --input examples/messy_variants.txt --output-dir output
```

`output/` will contain one TSV per variant type, a `review.tsv` holding
everything the tool refused to guess at, and `run_manifest.txt` recording the
input checksum, tool version, git commit and which validator produced the
verdicts.

## Pipeline

```mermaid
flowchart TD
    A[Free-text input line] --> B{Already valid HGVS?}
    B -- yes --> P[Stage 0: Passthrough]
    B -- no --> C[Stage 1: Parse - locus, alleles, type]
    C --> D[Stage 2: Schema check - required fields, bounds]
    D -- incomplete --> R[review.tsv]
    D -- complete --> E[Stage 3: Draft - allele trimming, optional 3-prime shift]
    P --> F
    E --> F[Stage 4: Validate - OPTIONAL]
    F --> G[snv.tsv / indel.tsv / sv.tsv / cnv.tsv / repeat.tsv]
    F --> H[run_manifest.txt]
```

Without a validator the pipeline still runs end to end; records stay
`not_validated`. This is deliberate: the core must not depend on a network
service.

## Optional validators

Each is independent and any subset may be enabled. Install only what you need.

```bash
pip install pyfaidx        # reference base checking and 3-prime shifting
pip install hgvs           # transcript projection via UTA
```

```bash
# offline reference checking against a local GRCh38 FASTA
python hgvs_normalizer.py --input examples/messy_variants.txt --fasta GRCh38.fa

# canonical form via Mutalyzer 3 (standard library HTTP only, nothing to install)
python hgvs_normalizer.py --input examples/messy_variants.txt --mutalyzer

# transcript-level projection via UTA
export UTA_DB_URL="postgresql://anonymous:anonymous@uta.biocommons.org/uta/uta_20241220"
python hgvs_normalizer.py --input examples/messy_variants.txt --validate
```

| Flag | Purpose | Needs |
|---|---|---|
| `--fasta` | Reference base check, 3-prime shifting | `pyfaidx`, a GRCh38 FASTA |
| `--mutalyzer [URL]` | Syntax, reference base, canonical form | network, or a local instance |
| `--validate` | Transcript projection to `c.` | `hgvs`, UTA access |
| `--mane` | MANE Select transcript choice | a `gene<TAB>transcript` TSV |

Validators are chainable and run in order; the first conclusive verdict wins.
An unreachable service is reported per record as `validator_unavailable` and
the draft is preserved, never discarded.

## Docker

```bash
docker build -t hgvs-normalizer .

mkdir -p data && cp examples/messy_variants.txt data/
docker run --rm -v "$(pwd)/data:/data" hgvs-normalizer \
  --input /data/messy_variants.txt --output-dir /data/output
```

Multi-stage build: `psycopg2` needs a compiler to build but not to run, so the
toolchain is discarded after the wheels are made — 692 MB down to 461 MB. The
known-answer suite runs during the build, so a failing test prevents the image
from being produced. The container runs as an unprivileged user, which keeps
bind-mounted output writable from the host.

## Record types

| Type | Recognised from | Draft form |
|---|---|---|
| `SNV` | single base on both sides | `g.123A>G` |
| `INDEL` | multi-base alleles, `delins`, `ins` | `g.123_125del` |
| `SV` | `del`/`dup`/`inv` below the CNV threshold | `g.1000_5000inv` |
| `CNV` | `del`/`dup` at or above `CNV_MIN_BP` | `g.1000_5000del` |
| `REPEAT` | `CAG(42)`, `CAG[42]`, `(CAG)n` | `g.100_102CAG[42]` |
| `HGVS_INPUT` | already-valid HGVS | passed through unchanged |

Mitochondrial accessions take the `m.` prefix rather than `g.`, per HGVS.
VCF-style left-anchored alleles are trimmed before drafting, so
`chr6-43624673 TGAA>T` becomes `g.43624674_43624676del` rather than a delins
that describes no change at all.

## Validation statuses

| Status | Meaning |
|---|---|
| `not_validated` | No validator enabled |
| `no_draft` | Incomplete input, no candidate built |
| `validated_g` | Confirmed at genomic level |
| `validated_c` | Also projected onto a transcript |
| `failed_validation` | Rejected, usually a reference base mismatch |
| `failed_reference_check` | Local FASTA disagrees with the stated reference base |
| `skipped_too_large` | Region exceeds `VALIDATION_MAX_BP` |
| `validator_unavailable` | Service unreachable; the draft is preserved |

## Known limitations

- 3-prime shifting requires reference sequence. Without `--fasta` or
  `--mutalyzer`, drafts are valid but not canonical.
- Transcript selection is not MANE Select unless `--mane` is supplied; when
  several transcripts overlap, the choice is recorded but arbitrary.
- The parser is regex-based and has a ceiling; tools such as tmVar use
  sequence labelling for a reason.
- Throughput is untested beyond a few dozen records; NCBI request limits will
  matter at scale.

## Scope

This tool sits between two established layers and replaces neither. Mutalyzer,
VariantValidator and the `hgvs` package validate descriptions that are already
well formed; this tool calls them. tmVar and SETH extract variant mentions from
published literature and normalize to dbSNP identifiers; different input,
different output. The gap addressed here is narrower: in-house coordinate
tables written by people, converted to HGVS with an explicit audit trail and
with external validation available but never required.

## Requirements

Python 3.9 or newer for the core. Optional: `pyfaidx`, `hgvs` with UTA access,
Docker.

## License

MIT

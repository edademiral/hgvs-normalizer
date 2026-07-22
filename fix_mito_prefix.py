#!/usr/bin/env python3
"""
fix_mito_prefix.py - use m. instead of g. for mitochondrial variants.

HGVS nomenclature uses the m. prefix for variants on the mitochondrial
genome. The draft builder emitted g. for every accession, so chrM drafts
were subtly wrong - Mutalyzer silently corrected them, which is how this
was noticed.

Usage:
    python fix_mito_prefix.py [--dry-run]
"""

import argparse
import ast
import os
import shutil
import sys

HELPER = '''
MITO_ACCESSION = "NC_012920.1"


def hgvs_prefix(accession):
    """Mitochondrial variants take m.; everything else on a chromosome takes g."""
    return "m." if accession == MITO_ACCESSION else "g."

'''

ANCHOR_HELPER = 'def chrom_length(chrom: Optional[str]) -> Optional[int]:\n    return CHROM_LENGTHS_GRCH38.get(norm_chrom(chrom) or "")\n'

ANCHOR_BUILD = '''        vtype = record.vtype
        if vtype in ("SNV", "INDEL"):'''
NEW_BUILD = '''        prefix = hgvs_prefix(acc)
        vtype = record.vtype
        if vtype in ("SNV", "INDEL"):'''

ANCHOR_ALLELE = '''        ref, alt = (record.ref or ""), (record.alt or "")
        pos, ref, alt = trim_alleles(pos, ref, alt)'''
NEW_ALLELE = '''        prefix = hgvs_prefix(acc)
        ref, alt = (record.ref or ""), (record.alt or "")
        pos, ref, alt = trim_alleles(pos, ref, alt)'''

# self-test expectations that change
TEST_FIXES = [
    ('("chrM-8612 T>C", "SNV", "ok", "NC_012920.1:g.8612T>C")',
     '("chrM-8612 T>C", "SNV", "ok", "NC_012920.1:m.8612T>C")'),
    ('("chrM-73 A>G", "SNV", "ok", "NC_012920.1:g.73A>G")',
     '("chrM-73 A>G", "SNV", "ok", "NC_012920.1:m.73A>G")'),
    ('("chrM-152 T>C", "SNV", "ok", "NC_012920.1:g.152T>C")',
     '("chrM-152 T>C", "SNV", "ok", "NC_012920.1:m.152T>C")'),
]


def main():
    argp = argparse.ArgumentParser()
    argp.add_argument("--target", default="hgvs_normalizer.py")
    argp.add_argument("--dry-run", action="store_true")
    args = argp.parse_args()

    if not os.path.exists(args.target):
        print(f"  ERROR  {args.target} not found")
        sys.exit(1)

    with open(args.target, encoding="utf-8") as handle:
        source = handle.read()

    if "def hgvs_prefix" in source:
        print("  SKIP   already fixed; nothing to do")
        return

    for label, anchor in (("helper", ANCHOR_HELPER),
                          ("build()", ANCHOR_BUILD),
                          ("_build_allele()", ANCHOR_ALLELE)):
        if source.count(anchor) != 1:
            print(f"  ERROR  '{label}' anchor found {source.count(anchor)} times "
                  f"(expected 1) - fix by hand")
            sys.exit(1)

    patched = source
    patched = patched.replace(ANCHOR_HELPER, ANCHOR_HELPER + HELPER, 1)
    patched = patched.replace(ANCHOR_BUILD, NEW_BUILD, 1)
    patched = patched.replace(ANCHOR_ALLELE, NEW_ALLELE, 1)

    # every draft f-string hardcodes g. - make it depend on the accession
    swapped = patched.count('f"{acc}:g.')
    patched = patched.replace('f"{acc}:g.', 'f"{acc}:{prefix}')

    for old, new in TEST_FIXES:
        patched = patched.replace(old, new)

    try:
        ast.parse(patched)
    except SyntaxError as exc:
        print(f"  ERROR  result does not parse ({exc}); nothing written")
        sys.exit(1)

    if args.dry_run:
        print(f"  OK     dry run: helper added, {swapped} draft strings would change")
        return

    shutil.copy2(args.target, args.target + ".bak3")
    with open(args.target, "w", encoding="utf-8") as handle:
        handle.write(patched)
    print(f"  OK     fixed {args.target}  ({swapped} draft strings, backup .bak3)")
    print("\n  Try:   python hgvs_normalizer.py --self-test")


if __name__ == "__main__":
    main()

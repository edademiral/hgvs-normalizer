#!/usr/bin/env python3
"""
fix_mutalyzer_endpoint.py - correct the Mutalyzer 3 call shape.

The first version sent the description as a query parameter
(/normalize/?description=...), which returns 404. Mutalyzer 3 takes it as a
PATH segment (/normalize/<description>) and nests errors under "custom".

Usage:
    python fix_mutalyzer_endpoint.py
    python fix_mutalyzer_endpoint.py --dry-run
"""

import argparse
import ast
import os
import shutil
import sys

OLD_GET = '''    def _get(self, path, params):
        self._throttle()
        url = f"{self.base_url}{path}?{urllib.parse.urlencode(params)}"'''

NEW_GET = '''    def _get(self, description):
        self._throttle()
        # Mutalyzer 3 takes the description as a PATH segment, not a query
        # parameter. It contains characters ( : > _ ) that must be encoded.
        url = (f"{self.base_url}{self.normalize_path}"
               f"{urllib.parse.quote(description, safe='')}")'''

OLD_NORMALIZE = '''        if description not in self._cache:
            self._cache[description] = self._get(
                self.normalize_path, {"description": description})
        return self._cache[description]'''

NEW_NORMALIZE = '''        if description not in self._cache:
            self._cache[description] = self._get(description)
        return self._cache[description]'''

OLD_ERRORS = '''        errors = self._first(payload, self.ERROR_KEYS)
        if errors:'''

NEW_ERRORS = '''        # Mutalyzer 3 nests errors under "custom"; older shapes kept as fallback
        custom = payload.get("custom") if isinstance(payload, dict) else None
        errors = ((custom.get("errors") if isinstance(custom, dict) else None)
                  or self._first(payload, self.ERROR_KEYS))
        if errors:'''

OLD_KEYS = '''    NORMALIZED_KEYS = ("normalized_description", "normalized", "description")'''

NEW_KEYS = '''    NORMALIZED_KEYS = ("normalized_description", "corrected_description",
                       "normalized")'''

REPLACEMENTS = [
    ("request URL", OLD_GET, NEW_GET),
    ("cache lookup", OLD_NORMALIZE, NEW_NORMALIZE),
    ("error extraction", OLD_ERRORS, NEW_ERRORS),
    ("response keys", OLD_KEYS, NEW_KEYS),
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

    if "class MutalyzerValidator" not in source:
        print("  ERROR  MutalyzerValidator not found - run apply_mutalyzer_patch.py first")
        sys.exit(1)
    if "quote(description" in source:
        print("  SKIP   endpoint already fixed; nothing to do")
        return

    patched = source
    for label, old, new in REPLACEMENTS:
        if patched.count(old) != 1:
            print(f"  ERROR  '{label}' anchor found {patched.count(old)} times "
                  f"(expected 1) - fix by hand")
            sys.exit(1)
        patched = patched.replace(old, new, 1)

    try:
        ast.parse(patched)
    except SyntaxError as exc:
        print(f"  ERROR  result does not parse ({exc}); nothing written")
        sys.exit(1)

    if args.dry_run:
        print(f"  OK     dry run: {len(REPLACEMENTS)} replacements would be made")
        return

    shutil.copy2(args.target, args.target + ".bak2")
    with open(args.target, "w", encoding="utf-8") as handle:
        handle.write(patched)
    print(f"  OK     fixed {args.target}  (backup: {args.target}.bak2)")
    print("\n  Try:   python hgvs_normalizer.py --input kucuk_test.txt \\")
    print("                --output-dir out_mut --mutalyzer")


if __name__ == "__main__":
    main()

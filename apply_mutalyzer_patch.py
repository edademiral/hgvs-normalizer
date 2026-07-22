#!/usr/bin/env python3
"""
apply_mutalyzer_patch.py - add the optional Mutalyzer validator to
hgvs_normalizer.py without hand-editing.

Usage:
    python apply_mutalyzer_patch.py                     # patches ./hgvs_normalizer.py
    python apply_mutalyzer_patch.py --target other.py
    python apply_mutalyzer_patch.py --dry-run           # show what would change

Safe to run: it writes a .bak first, checks every anchor before touching
anything, and refuses to save if the result does not parse.
"""

import argparse
import ast
import os
import shutil
import sys

BLOCK_IMPORTS = '''import json
import time
import urllib.error
import urllib.parse
import urllib.request
'''

BLOCK_CONSTANTS = '''
MUTALYZER_BASE_URL = "https://mutalyzer.nl/api"
MUTALYZER_NORMALIZE_PATH = "/normalize/"   # verify against your instance
MUTALYZER_MIN_INTERVAL = 0.5               # seconds between calls
'''

BLOCK_CLASSES = '''

class MutalyzerValidator:
    """
    Optional validator backed by Mutalyzer 3.

    Works against the public service or a local instance
    (`pip install mutalyzer-api`, then http://localhost:5000/api).
    A local instance is preferred for batch work: no network dependency and
    no load on a shared public server.

    Mutalyzer checks syntax AND the reference base and returns a normalized
    (3'-shifted) description, so it covers genomic-level validation without
    needing a transcript database.
    """

    name = "mutalyzer"

    # Read defensively: field names are tried in order so a schema change
    # degrades to a clear note instead of a crash.
    NORMALIZED_KEYS = ("normalized_description", "normalized", "description")
    ERROR_KEYS = ("errors", "error")
    INFO_KEYS = ("infos", "warnings")

    def __init__(self, base_url=MUTALYZER_BASE_URL, timeout=20.0,
                 min_interval=MUTALYZER_MIN_INTERVAL,
                 normalize_path=MUTALYZER_NORMALIZE_PATH):
        self.base_url = base_url.rstrip("/")
        self.normalize_path = normalize_path
        self.timeout = timeout
        self.min_interval = min_interval
        self._last_call = 0.0
        self._cache = {}
        self._probe()
        self.data_version = self.base_url

    def _throttle(self):
        wait = self.min_interval - (time.monotonic() - self._last_call)
        if wait > 0:
            time.sleep(wait)
        self._last_call = time.monotonic()

    def _get(self, path, params):
        self._throttle()
        url = f"{self.base_url}{path}?{urllib.parse.urlencode(params)}"
        request = urllib.request.Request(
            url, headers={"Accept": "application/json",
                          "User-Agent": f"hgvs_normalizer/{__version__}"})
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            # An invalid description comes back as 4xx with a JSON body;
            # that body is the answer, not a transport failure.
            body = exc.read().decode("utf-8", errors="replace")
            try:
                return json.loads(body)
            except ValueError:
                raise RuntimeError(f"HTTP {exc.code}: {body[:120]}") from exc

    def _probe(self):
        """Fail once at startup rather than on every record."""
        try:
            urllib.request.urlopen(
                urllib.request.Request(self.base_url + "/",
                                       headers={"Accept": "application/json"}),
                timeout=self.timeout)
        except urllib.error.HTTPError:
            pass            # a 4xx still proves the service answers
        except Exception as exc:
            raise RuntimeError(
                f"Mutalyzer unreachable at {self.base_url}: "
                f"{type(exc).__name__}") from exc

    @staticmethod
    def _first(payload, keys):
        for key in keys:
            if isinstance(payload, dict) and payload.get(key):
                return payload[key]
        return None

    @staticmethod
    def _describe(item):
        if isinstance(item, dict):
            return str(item.get("details") or item.get("message")
                       or item.get("code") or item)
        return str(item)

    def _normalize(self, description):
        if description not in self._cache:
            self._cache[description] = self._get(
                self.normalize_path, {"description": description})
        return self._cache[description]

    def validate(self, record):
        draft = record.hgvs_draft
        if not draft:
            record.validation_status = "no_draft"
            return
        size = getattr(record, "size", None)
        if size and size > VALIDATION_MAX_BP:
            record.validation_status = "skipped_too_large"
            return

        try:
            payload = self._normalize(draft)
        except Exception as exc:                       # noqa: BLE001
            record.validation_status = "validator_unavailable"
            record.notes.append(
                f"mutalyzer error: {type(exc).__name__}: {str(exc)[:100]}")
            return

        errors = self._first(payload, self.ERROR_KEYS)
        if errors:
            items = errors if isinstance(errors, list) else [errors]
            record.validation_status = "failed_validation"
            record.notes.append("mutalyzer: " + "; ".join(
                self._describe(item)[:90] for item in items[:2]))
            return

        normalized = self._first(payload, self.NORMALIZED_KEYS)
        if not normalized:
            record.validation_status = "unparseable"
            keys = sorted(payload) if isinstance(payload, dict) else type(payload).__name__
            record.notes.append(
                f"mutalyzer returned no normalized description (keys: {keys})")
            return

        record.hgvs = str(normalized)
        record.validation_status = "validated_g"
        if str(normalized) != draft:
            record.notes.append(f"mutalyzer rewrote the draft (was {draft})")
        infos = self._first(payload, self.INFO_KEYS) or []
        items = infos if isinstance(infos, list) else [infos]
        for item in items[:2]:
            record.notes.append("mutalyzer: " + self._describe(item)[:90])


class ChainValidator:
    """
    Run several validators in order; the first conclusive verdict wins.

    This is what makes the external services genuinely optional: any subset
    can be supplied, and an inconclusive result (a reference check that
    passed but proves nothing about transcripts) falls through to the next.
    """

    CONCLUSIVE = frozenset({
        "validated_c", "validated_g", "failed_validation",
        "failed_reference_check", "unparseable", "skipped_too_large",
    })

    def __init__(self, validators):
        self.validators = [v for v in validators if v is not None]
        self.name = "+".join(v.name for v in self.validators) or "disabled"
        self.data_version = ("; ".join(f"{v.name}={v.data_version}"
                                       for v in self.validators) or "n/a")

    def validate(self, record):
        record.validation_status = "not_validated"
        for validator in self.validators:
            validator.validate(record)
            if record.validation_status in self.CONCLUSIVE:
                record.notes.append(f"verdict from {validator.name}")
                return

'''

BLOCK_BUILDER = '''
def build_validator(args, reference, logger):
    """Assemble whichever validators are available; none of them is required."""
    stages = []

    if reference:
        stages.append(FastaValidator(reference))

    if getattr(args, "mutalyzer", None):
        try:
            stages.append(MutalyzerValidator(args.mutalyzer))
            logger.info("Mutalyzer ready: %s", args.mutalyzer)
        except Exception as exc:                       # noqa: BLE001
            logger.warning("Mutalyzer unavailable (%s); skipping that stage", exc)

    if args.validate:
        logger.info("connecting to UTA (first connection can be slow)...")
        try:
            stages.append(UtaValidator(mane=load_mane(args.mane)))
            logger.info("UTA ready")
        except Exception as exc:                       # noqa: BLE001
            logger.warning("UTA unavailable (%s); skipping that stage",
                           type(exc).__name__)

    if not stages:
        return NullValidator()
    if len(stages) == 1:
        return stages[0]
    return ChainValidator(stages)


'''

BLOCK_CLI_ARG = '''    argp.add_argument("--mutalyzer", nargs="?", const=MUTALYZER_BASE_URL,
                      default=None, metavar="URL",
                      help="validate with Mutalyzer 3; optionally a base URL "
                           "(local instance: http://localhost:5000/api)")
'''


def fail(message):
    print(f"  ERROR  {message}")
    sys.exit(1)


def main():
    argp = argparse.ArgumentParser()
    argp.add_argument("--target", default="hgvs_normalizer.py")
    argp.add_argument("--dry-run", action="store_true")
    args = argp.parse_args()

    if not os.path.exists(args.target):
        fail(f"{args.target} not found - run this from the project directory")

    with open(args.target, encoding="utf-8") as handle:
        source = handle.read()

    if "class MutalyzerValidator" in source:
        print("  SKIP   MutalyzerValidator is already present; nothing to do")
        return

    # ---- verify every anchor before touching anything ----
    anchors = {
        "imports": "import subprocess\n",
        "classes": "# ===================================================================\n# Orchestration",
        "builder": "def main() -> None:",
        "cli_arg": '    argp.add_argument("--self-test", action="store_true")',
    }
    for label, anchor in anchors.items():
        if source.count(anchor) != 1:
            fail(f"anchor '{label}' found {source.count(anchor)} times "
                 f"(expected exactly 1) - patch by hand instead")

    # the validator-selection block in main() varies if it was edited
    start_marker = "    validator = NullValidator()"
    end_marker = "    pipeline = NormalizationPipeline("
    if start_marker not in source or end_marker not in source:
        fail("could not locate the validator selection block in main()")
    selection = source[source.index(start_marker):source.index(end_marker)]

    # ---- apply ----
    patched = source
    patched = patched.replace(anchors["imports"],
                              anchors["imports"] + BLOCK_IMPORTS, 1)
    patched = patched.replace('NUCLEOTIDES = frozenset("ACGTN")',
                              'NUCLEOTIDES = frozenset("ACGTN")\n' + BLOCK_CONSTANTS, 1)
    patched = patched.replace(anchors["classes"],
                              BLOCK_CLASSES + anchors["classes"], 1)
    patched = patched.replace(anchors["builder"],
                              BLOCK_BUILDER + anchors["builder"], 1)
    patched = patched.replace(anchors["cli_arg"],
                              BLOCK_CLI_ARG + anchors["cli_arg"], 1)
    patched = patched.replace(
        selection, "    validator = build_validator(args, reference, logger)\n\n", 1)

    try:
        ast.parse(patched)
    except SyntaxError as exc:
        fail(f"patched file does not parse ({exc}); nothing was written")

    if args.dry_run:
        print(f"  OK     dry run: {len(patched) - len(source)} characters would be added")
        return

    backup = args.target + ".bak"
    shutil.copy2(args.target, backup)
    with open(args.target, "w", encoding="utf-8") as handle:
        handle.write(patched)

    print(f"  OK     patched {args.target}")
    print(f"  OK     backup saved as {backup}")
    print("\n  Try:   python hgvs_normalizer.py --help | grep mutalyzer")


if __name__ == "__main__":
    main()

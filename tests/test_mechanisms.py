"""Keep the cross-problem taxonomy honest.

mechanisms.json exists so an agent can ask "has this approach family been
tried, and where?" across problems without reading the prose. That only works
if every tag the indexes use is defined there, and every definition is real
prose rather than a placeholder. Same policy as the schema tests: a typo
should fail the PR, not silently degrade the index.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TAXONOMY = json.loads((ROOT / "mechanisms.json").read_text(encoding="utf-8"))
KEBAB = re.compile(r"^[a-z0-9-]+$")


def used_tags() -> dict[str, list[str]]:
    """tag -> ["problem/id", ...] across every prior-art index."""
    out: dict[str, list[str]] = {}
    for index in sorted(ROOT.glob("problems/*/prior-art.json")):
        data = json.loads(index.read_text(encoding="utf-8"))
        for attempt in data["attempts"]:
            for tag in attempt["mechanism"]:
                out.setdefault(tag, []).append(f"{data['problem']}/{attempt['id']}")
    return out


def test_every_used_tag_is_defined():
    """The point of the file: no index tag may be missing from the taxonomy."""
    missing = {
        tag: uses
        for tag, uses in used_tags().items()
        if tag not in TAXONOMY["mechanisms"]
    }
    assert not missing, (
        "mechanism tags used in prior-art.json but not defined in "
        "mechanisms.json (add them with a field and description): "
        + ", ".join(f"{t} ({', '.join(u)})" for t, u in sorted(missing.items()))
    )


def test_every_mechanism_has_field_and_description():
    problems = []
    for tag, meta in TAXONOMY["mechanisms"].items():
        if meta.get("field") not in TAXONOMY["fields"]:
            problems.append(f"{tag}: field {meta.get('field')!r} not in fields")
        if len(meta.get("description", "")) < 20:
            problems.append(f"{tag}: description missing or too short to inform")
        extra = set(meta) - {"field", "description"}
        if extra:
            problems.append(f"{tag}: unknown keys {sorted(extra)}")
    assert not problems, "\n  ".join(["taxonomy defects:"] + problems)


def test_names_are_kebab_case():
    bad = [
        name
        for name in list(TAXONOMY["fields"]) + list(TAXONOMY["mechanisms"])
        if not KEBAB.fullmatch(name)
    ]
    assert not bad, f"non-kebab-case names: {bad}"


def test_entries_are_sorted():
    """Sorted keys keep diffs reviewable when tags are added."""
    for key in ("fields", "mechanisms"):
        names = list(TAXONOMY[key])
        assert names == sorted(names), f"mechanisms.json {key} keys are not sorted"


def test_field_descriptions_exist():
    bad = [f for f, desc in TAXONOMY["fields"].items() if len(desc) < 20]
    assert not bad, f"fields with missing/stub descriptions: {bad}"


def test_query_script_runs():
    """The CLI is part of the contract: gaps/matrix must work on the live repo."""
    for cmd in (["list"], ["matrix"], ["gaps", "collatz"], ["where", "coupling"]):
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "mechanisms.py"), *cmd],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"mechanisms.py {cmd} failed: {result.stderr}"
        assert result.stdout.strip(), f"mechanisms.py {cmd} produced no output"

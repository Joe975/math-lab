"""Validate every prior-art index against docs/prior-art.schema.json.

These indexes are load-bearing twice over: they are the first thing an arriving
agent reads, and CI derives its tier-isolation checks from their `leak_terms`.
A misspelled key therefore does not merely look untidy — it silently removes a
protection. A contributor's typo should fail their PR, not degrade the repo.

Validation is hand-rolled against the subset of JSON Schema the file uses, so
the repo stays dependency-free.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCHEMA = json.loads(
    (ROOT / "docs" / "prior-art.schema.json").read_text(encoding="utf-8")
)
INDEXES = sorted(ROOT.glob("problems/*/prior-art.json"))


def check(value, spec: dict, where: str, errors: list[str]) -> None:
    types = spec.get("type")
    if types:
        wanted = types if isinstance(types, list) else [types]
        ok = any(
            (t == "string" and isinstance(value, str))
            or (t == "array" and isinstance(value, list))
            or (t == "object" and isinstance(value, dict))
            or (t == "null" and value is None)
            for t in wanted
        )
        if not ok:
            errors.append(f"{where}: expected {types}, got {type(value).__name__}")
            return

    if "enum" in spec and value not in spec["enum"]:
        errors.append(f"{where}: {value!r} not one of {spec['enum']}")
    if "pattern" in spec and isinstance(value, str):
        if not re.fullmatch(spec["pattern"], value):
            errors.append(f"{where}: {value!r} does not match {spec['pattern']}")
    if "minLength" in spec and isinstance(value, str):
        if len(value) < spec["minLength"]:
            errors.append(f"{where}: shorter than {spec['minLength']} chars")
    if "minItems" in spec and isinstance(value, list):
        if len(value) < spec["minItems"]:
            errors.append(f"{where}: needs at least {spec['minItems']} item(s)")
    if "items" in spec and isinstance(value, list):
        for i, item in enumerate(value):
            check(item, spec["items"], f"{where}[{i}]", errors)


def check_object(obj: dict, spec: dict, where: str, errors: list[str]) -> None:
    props = spec.get("properties", {})
    for key in spec.get("required", []):
        if key not in obj:
            errors.append(f"{where}: missing required key {key!r}")
    if spec.get("additionalProperties") is False:
        for key in obj:
            if key not in props:
                errors.append(
                    f"{where}: unknown key {key!r} — a typo here silently does nothing"
                )
    for key, value in obj.items():
        if key in props:
            check(value, props[key], f"{where}.{key}", errors)


@pytest.mark.parametrize("path", INDEXES, ids=lambda p: p.parent.name)
def test_index_matches_schema(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    errors: list[str] = []
    check_object(data, SCHEMA, path.parent.name, errors)

    attempt_spec = SCHEMA["$defs"]["attempt"]
    for attempt in data.get("attempts", []):
        ident = attempt.get("id", "?")
        check_object(attempt, attempt_spec, f"{path.parent.name}/{ident}", errors)

    assert not errors, "schema violations:\n  " + "\n  ".join(errors)


@pytest.mark.parametrize("path", INDEXES, ids=lambda p: p.parent.name)
def test_problem_field_matches_directory(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["problem"] == path.parent.name, (
        f"{path}: problem field {data['problem']!r} does not match its directory"
    )


@pytest.mark.parametrize("path", INDEXES, ids=lambda p: p.parent.name)
def test_referenced_records_exist(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    for attempt in data["attempts"]:
        target = path.parent / attempt["file"]
        assert target.exists(), f"{path}: {attempt['file']} does not exist"
        assert attempt["file"].startswith(f"attempts/{attempt['id']}-"), (
            f"{path}: attempt {attempt['id']} points at {attempt['file']}, "
            f"whose number does not match"
        )


@pytest.mark.parametrize("path", INDEXES, ids=lambda p: p.parent.name)
def test_cross_references_resolve(path: Path):
    """verifies/superseded_by must name an attempt that exists."""
    data = json.loads(path.read_text(encoding="utf-8"))
    ids = {a["id"] for a in data["attempts"]}
    for attempt in data["attempts"]:
        for field in ("verifies", "superseded_by"):
            target = attempt.get(field)
            if target and target not in ids:
                pytest.fail(
                    f"{path}: attempt {attempt['id']}.{field} points at {target!r}, "
                    f"which is not an attempt here"
                )


@pytest.mark.parametrize("path", INDEXES, ids=lambda p: p.parent.name)
def test_attempt_ids_are_unique_and_ordered(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    ids = [a["id"] for a in data["attempts"]]
    assert len(ids) == len(set(ids)), f"{path}: duplicate attempt ids: {ids}"
    assert ids == sorted(ids), f"{path}: attempts out of order: {ids}"


def test_scaffolder_produces_a_record_the_schema_rejects_until_filled(tmp_path):
    """new_attempt.py leaves TODO placeholders on purpose.

    A scaffold that validated straight away would let a half-filled index be
    committed. Assert the guard rail is actually there.
    """
    attempt_spec = SCHEMA["$defs"]["attempt"]
    stub = {
        "id": "099",
        "file": "attempts/099-example.md",
        "date": "2026-01-01",
        "mode": "TODO-blind-or-informed",
        "status": "TODO",
        "mechanism": ["TODO-approach-family-tag"],
        "one_line": "TODO: what this established, in one sentence.",
        "leak_terms": [],
    }
    errors: list[str] = []
    check_object(stub, attempt_spec, "stub", errors)
    assert errors, "the scaffolder's placeholders should fail validation"

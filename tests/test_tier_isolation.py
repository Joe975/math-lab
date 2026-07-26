"""Enforce the blind/informed cut.

The repo's whole contribution model rests on one property: a `blind` agent,
given only tier-0 files, learns nothing about what this lab already tried.
That property erodes silently -- someone writes "as we found in 001" into a
PROBLEM.md and nobody notices. These tests make it a red build instead.

The mechanism: every attempt in a `prior-art.json` declares `leak_terms`,
distinctive strings naming this lab's *findings* rather than the published
literature. No tier-0 prose file may contain one.

Scope, stated so the gap is known rather than assumed away: the term check
runs against tier-0 **prose** (PROBLEM.md and the root docs). It does not run
against `harness/` source, because a shared verifier's own identifiers,
flags and default constants are its API -- `singmaster_search.py` cannot avoid
the phrase "canonical cell" and still be readable. Harness code is instead
held to a weaker, structural rule: no narrative back-references to attempts.
Adding a harness file is therefore a human review point.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
TIERS = json.loads((ROOT / "tiers.json").read_text(encoding="utf-8"))

TEXT_SUFFIXES = {".md", ".json", ".py", ".c", ".sh", ".cff", ".txt", ".yml", ".yaml"}


def expand(patterns: list[str]) -> set[Path]:
    """Resolve tier globs to the set of existing files they match."""
    out: set[Path] = set()
    for pattern in patterns:
        for path in ROOT.glob(pattern):
            if path.is_file():
                out.add(path)
    return out


def tier0_files() -> set[Path]:
    return expand(TIERS["tier0"])


def tier0_prose() -> set[Path]:
    """Tier-0 files whose job is to explain, not to compute."""
    return {p for p in tier0_files() if p.suffix.lower() in {".md", ".cff"}}


def bash() -> str:
    """Locate bash: on PATH, or Git for Windows' copy."""
    import shutil

    found = shutil.which("bash")
    if found:
        return found
    for candidate in (
        r"C:\Program Files\Git\bin\bash.exe",
        r"C:\Program Files (x86)\Git\bin\bash.exe",
    ):
        if Path(candidate).exists():
            return candidate
    pytest.skip("bash not available")


def tier1_files() -> set[Path]:
    return expand(TIERS["tier1"])


def leak_terms() -> dict[str, list[str]]:
    """All declared leak terms, keyed by term -> which attempts declared it."""
    terms: dict[str, list[str]] = {}
    for index in sorted(ROOT.glob("problems/*/prior-art.json")):
        data = json.loads(index.read_text(encoding="utf-8"))
        problem = data["problem"]
        for attempt in data["attempts"]:
            for term in attempt.get("leak_terms", []):
                terms.setdefault(term, []).append(f"{problem}/{attempt['id']}")
    return terms


def find_leaks(files: set[Path], terms: dict[str, list[str]]) -> list[str]:
    """Return a human-readable finding for every (file, leaked term) pair."""
    findings = []
    for path in sorted(files):
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            body = path.read_text(encoding="utf-8").lower()
        except UnicodeDecodeError:
            continue
        for term, sources in terms.items():
            if term.lower() in body:
                rel = path.relative_to(ROOT).as_posix()
                findings.append(f"{rel}: leaks {term!r} (from {', '.join(sources)})")
    return findings


# --- the actual guarantees ------------------------------------------------


def test_no_prior_art_leaks_into_tier0_prose():
    """A blind agent must not be able to infer our findings from tier-0 prose."""
    findings = find_leaks(tier0_prose(), leak_terms())
    assert not findings, "prior art leaked into tier 0:\n  " + "\n  ".join(findings)


def test_harness_has_no_narrative_backreferences():
    """Harness code may use its own vocabulary, but must not narrate our history.

    This is the weaker rule that covers what the term check deliberately does
    not: a comment like "as attempt 002 found" turns a shared verifier into
    prior art.
    """
    patterns = [
        "prior-art",
        "PRIOR-ART",
        "attempts/",
        "attempt 00",
        "as we found",
        "we previously",
        "our earlier",
    ]
    findings = []
    for path in sorted(ROOT.glob("harness/**/*")):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        body = path.read_text(encoding="utf-8", errors="ignore")
        for pattern in patterns:
            if pattern in body:
                findings.append(f"{path.relative_to(ROOT).as_posix()}: {pattern!r}")
    assert not findings, "harness narrates prior art:\n  " + "\n  ".join(findings)


def test_leak_detector_actually_detects(tmp_path):
    """Guard the guard: a planted leak must be caught.

    Without this, a bug that made find_leaks always return [] would leave
    test_no_prior_art_leaks_into_tier0 passing forever.
    """
    terms = leak_terms()
    assert terms, "no leak_terms declared anywhere; the leak test would be vacuous"

    planted = sorted(terms)[0]
    decoy = tmp_path / "PROBLEM.md"
    decoy.write_text(f"Statement. Foo.\n\nWe already tried {planted} here.\n", encoding="utf-8")

    # find_leaks reports relative to ROOT, so point it at a file under ROOT.
    scratch = ROOT / "problems" / "_leak_probe.md"
    try:
        scratch.write_text(decoy.read_text(encoding="utf-8"), encoding="utf-8")
        findings = find_leaks({scratch}, terms)
        assert findings, f"planted leak {planted!r} was not detected"
    finally:
        scratch.unlink(missing_ok=True)


def test_tiers_do_not_overlap():
    """A file in both tiers would make the cut meaningless."""
    both = tier0_files() & tier1_files()
    assert not both, "files claimed by both tiers: " + ", ".join(
        sorted(p.relative_to(ROOT).as_posix() for p in both)
    )


def test_every_file_is_tiered():
    """No file may sit outside the cut: unclassified means unreviewed."""
    tracked = subprocess.run(
        ["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    known = {p.relative_to(ROOT).as_posix() for p in tier0_files() | tier1_files()}
    exempt = {".gitignore", ".github/workflows/ci.yml", "scripts/blind.sh"}
    exempt |= {p for p in tracked if p.startswith("tests/")}
    untiered = sorted(set(tracked) - known - exempt)
    assert not untiered, "files in neither tier (add them to tiers.json): " + ", ".join(untiered)


def test_every_attempt_declares_its_mode():
    """`mode` is the dataset: whether prior knowledge helped or anchored."""
    missing = []
    for index in sorted(ROOT.glob("problems/*/prior-art.json")):
        data = json.loads(index.read_text(encoding="utf-8"))
        for attempt in data["attempts"]:
            if attempt.get("mode") not in {"blind", "informed"}:
                missing.append(f"{data['problem']}/{attempt['id']}")
    assert not missing, "attempts with no valid mode: " + ", ".join(missing)


def test_prior_art_index_matches_attempt_files():
    """Every indexed attempt exists, and every attempt file is indexed."""
    problems = [p for p in sorted(ROOT.glob("problems/*")) if p.is_dir()]
    for problem in problems:
        index_path = problem / "prior-art.json"
        assert index_path.exists(), f"{problem.name} has no prior-art.json"
        data = json.loads(index_path.read_text(encoding="utf-8"))

        indexed = {a["file"] for a in data["attempts"]}
        for rel in indexed:
            assert (problem / rel).exists(), f"{problem.name}: indexed but missing: {rel}"

        on_disk = {
            f"attempts/{p.name}" for p in problem.glob("attempts/*.md")
        }
        assert on_disk == indexed, (
            f"{problem.name}: index and attempts/ disagree; "
            f"only on disk: {sorted(on_disk - indexed)}, "
            f"only in index: {sorted(indexed - on_disk)}"
        )


@pytest.mark.parametrize(
    "problem",
    [p.name for p in sorted(ROOT.glob("problems/*")) if p.is_dir()],
)
def test_blind_checkout_contains_no_tier1(problem, tmp_path):
    """The end-to-end guarantee: blind.sh must not emit a single tier-1 file."""
    dest = tmp_path / "work"
    result = subprocess.run(
        [bash(), str(ROOT / "scripts" / "blind.sh"), problem, str(dest)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"blind.sh failed: {result.stderr}"

    forbidden = {"PRIOR-ART.md", "prior-art.json", "STATUS.md", "GUIDANCE.md"}
    forbidden_dirs = {"attempts", "explore", "data"}

    emitted = [p for p in dest.rglob("*") if p.is_file()]
    assert emitted, "blind.sh produced nothing"

    for path in emitted:
        rel = path.relative_to(dest).as_posix()
        assert path.name not in forbidden, f"blind copy contains tier-1 file: {rel}"
        assert not (forbidden_dirs & set(rel.split("/"))), (
            f"blind copy contains tier-1 directory content: {rel}"
        )

    # And the thing it is supposed to contain.
    assert (dest / "problems" / problem / "PROBLEM.md").exists()
    assert (dest / "MODE").read_text(encoding="utf-8").startswith("mode: blind")


def test_blind_checkout_body_is_leak_free(tmp_path):
    """Not just the filenames -- the actual bytes must be clean."""
    dest = tmp_path / "work"
    subprocess.run(
        [bash(), str(ROOT / "scripts" / "blind.sh"), "union-closed", str(dest)],
        capture_output=True,
        text=True,
        check=True,
    )
    terms = leak_terms()
    findings = []
    for path in dest.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        body = path.read_text(encoding="utf-8", errors="ignore").lower()
        for term in terms:
            if term.lower() in body:
                findings.append(f"{path.relative_to(dest).as_posix()}: {term!r}")
    assert not findings, "blind copy leaks prior art:\n  " + "\n  ".join(findings)

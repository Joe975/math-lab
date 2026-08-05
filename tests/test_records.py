"""Hold records and docs to the format CONTRIBUTING describes.

Two failure modes this catches, both of which had already happened:

1. **Docs pointing at files that moved.** `STATUS.md` once told an agent to
   find the cycle prompt in a README section that had been renamed, and the
   prompt did not exist at all. A dangling pointer in the orientation docs is
   worse than a missing one, because the reader trusts it.

2. **A documented record format nothing follows.** A contributor reading
   `CONTRIBUTING.md` should produce a record that looks like the ones already
   here. Where the two diverge, one of them is wrong.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

# Records take one of two shapes. An adversarial review is not an "attempt"
# and forcing it into the attempt sections would make it worse, so it has its
# own contract.
ATTEMPT_SECTIONS = [
    "## Approach",
    "## What was done",
    "## Outcome",
    "what survived",
    "## Leads generated",
]
REVIEW_SECTIONS = [
    "## Claims attacked",
    "## Refutations found",
    "## Claims that survive",
]

# Records written before the format was codified. They are exempt from the
# References requirement only -- never from the section contract. They are NOT
# retro-fitted, because this repo does not edit existing records; the list is
# here rather than hidden so the exception stays visible and finite.
LEGACY_WITHOUT_REFERENCES = {
    "erdos-gyarfas/001-cycle-spectrum-baseline.md",
    "erdos-gyarfas/002-c4free-search-n22plus.md",
    "erdos-straus/001-identity-coverage-map.md",
    "erdos-straus/002-class-601-anomaly.md",
    "graceful-trees/001-labeling-count-statistics.md",
    "lonely-runner/001-k8-near-tight-scan.md",
    "singmaster/001-collision-search.md",
}


def records() -> list[Path]:
    return sorted(ROOT.glob("problems/*/attempts/*.md"))


def rel(path: Path) -> str:
    return f"{path.parent.parent.name}/{path.name}"


@pytest.mark.parametrize("record", records(), ids=rel)
def test_record_follows_one_of_the_two_shapes(record: Path):
    body = record.read_text(encoding="utf-8", errors="ignore")

    missing_attempt = [s for s in ATTEMPT_SECTIONS if s not in body]
    missing_review = [s for s in REVIEW_SECTIONS if s not in body]

    if not missing_attempt or not missing_review:
        return

    pytest.fail(
        f"{rel(record)} matches neither record shape.\n"
        f"  as an attempt, missing: {missing_attempt}\n"
        f"  as a review,  missing: {missing_review}\n"
        f"See CONTRIBUTING.md."
    )


@pytest.mark.parametrize("record", records(), ids=rel)
def test_new_records_cite_their_sources(record: Path):
    if rel(record) in LEGACY_WITHOUT_REFERENCES:
        pytest.skip("predates the format; records are never edited")
    body = record.read_text(encoding="utf-8", errors="ignore")
    assert "## References" in body, (
        f"{rel(record)} has no References section. Mathematical claims need "
        f"their sources, and machine transcriptions need marking as such."
    )


def test_docs_do_not_point_at_missing_files():
    """Every path referenced in the orientation docs must resolve.

    Attempt records are excluded: they refer to the pre-restructure layout and
    are never edited, which AGENTS.md tells the reader.
    """
    pattern = re.compile(
        r"`([A-Za-z0-9_][A-Za-z0-9_./-]*\.(?:md|py|sh|json|c|cff|yml))`"
    )
    broken = []
    for doc in sorted(ROOT.rglob("*.md")):
        if ".git" in doc.parts or "attempts" in doc.parts or ".lake" in doc.parts:
            # .lake/ holds fetched Lean dependencies (e.g. mathlib's own
            # docs); git-ignored, not ours to lint.
            continue
        text = doc.read_text(encoding="utf-8", errors="ignore")
        for match in pattern.finditer(text):
            ref = match.group(1)
            if "<" in ref or "*" in ref or ref.startswith("NNN"):
                continue
            if (ROOT / ref).exists() or (doc.parent / ref).exists():
                continue
            if any(ROOT.rglob(Path(ref).name)):
                continue
            broken.append(f"{doc.relative_to(ROOT).as_posix()} -> {ref}")
    assert not broken, "docs reference files that do not exist:\n  " + "\n  ".join(
        broken
    )


def test_the_orientation_chain_is_intact():
    """A fresh agent lands on CLAUDE.md. Everything else must be reachable.

    Not a style check: if this chain breaks, an agent starting cold silently
    works without the rules rather than failing loudly.
    """
    claude = (ROOT / "CLAUDE.md").read_text(encoding="utf-8")
    for pointer in ("AGENTS.md", "docs/CYCLE.md", "CONTRIBUTING.md", "STATUS.md"):
        assert pointer in claude, f"CLAUDE.md does not point at {pointer}"

    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    for pointer in ("scripts/blind.sh", "prior-art.json", "CONTRIBUTING.md"):
        assert pointer in agents, f"AGENTS.md does not point at {pointer}"

    cycle = (ROOT / "docs" / "CYCLE.md").read_text(encoding="utf-8")
    for pointer in ("STATUS.md", "GUIDANCE.md", "scripts/new_attempt.py"):
        assert pointer in cycle, f"docs/CYCLE.md does not point at {pointer}"

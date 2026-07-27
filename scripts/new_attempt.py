#!/usr/bin/env python3
"""Scaffold a new attempt record and its index entry.

    python scripts/new_attempt.py <problem> <slug> [--mode blind|informed]

Picks the next number by scanning attempts/ rather than trusting the index, so
a half-finished record on disk cannot cause a collision. Writes the record from
docs/attempt-template.md and adds a matching stub to prior-art.json, because a
record without an index entry fails CI and is easy to forget.

Standard library only.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / "docs" / "attempt-template.md"


def next_id(attempts_dir: Path) -> str:
    """Next free number, from what is on disk."""
    used = [
        int(m.group(1))
        for p in attempts_dir.glob("*.md")
        if (m := re.match(r"^(\d{3})-", p.name))
    ]
    return f"{max(used) + 1 if used else 1:03d}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("problem", help="problem slug, e.g. union-closed")
    parser.add_argument("slug", help="short kebab-case name for the approach")
    parser.add_argument(
        "--mode",
        choices=["blind", "informed"],
        help="whether you read the prior art. Omit to fill in by hand.",
    )
    args = parser.parse_args()

    problem_dir = ROOT / "problems" / args.problem
    if not problem_dir.is_dir():
        available = sorted(p.name for p in (ROOT / "problems").iterdir() if p.is_dir())
        print(f"error: no such problem: {args.problem}", file=sys.stderr)
        print(f"available: {', '.join(available)}", file=sys.stderr)
        return 2

    if not re.fullmatch(r"[a-z0-9-]+", args.slug):
        print("error: slug must be lower-case letters, digits and hyphens",
              file=sys.stderr)
        return 2

    attempts = problem_dir / "attempts"
    attempts.mkdir(exist_ok=True)
    ident = next_id(attempts)
    rel = f"attempts/{ident}-{args.slug}.md"
    target = problem_dir / rel
    if target.exists():
        print(f"error: {target} already exists", file=sys.stderr)
        return 1

    body = TEMPLATE.read_text(encoding="utf-8")
    body = body.replace(
        "# NNN — <short title of the approach>",
        f"# {ident} — <short title of the approach>",
    )
    body = body.replace(
        "- **Problem:** <name>, `problems/<slug>/PROBLEM.md`",
        f"- **Problem:** {args.problem}, `problems/{args.problem}/PROBLEM.md`",
    )
    body = body.replace("- **Date:** YYYY-MM-DD", f"- **Date:** {date.today()}")
    if args.mode:
        body = body.replace("- **Mode:** blind | informed", f"- **Mode:** {args.mode}")
    target.write_text(body, encoding="utf-8")

    index_path = problem_dir / "prior-art.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    index["attempts"].append(
        {
            "id": ident,
            "file": rel,
            "date": str(date.today()),
            "mode": args.mode or "TODO-blind-or-informed",
            "status": "TODO",
            "mechanism": ["TODO-approach-family-tag"],
            "one_line": "TODO: what this established, in one sentence.",
            "leak_terms": [],
        }
    )
    index_path.write_text(
        json.dumps(index, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print(f"created {target.relative_to(ROOT)}")
    print(f"added a stub to {index_path.relative_to(ROOT)}")
    print()
    print("Next: fill in both, then run `python -m pytest tests/ -q`.")
    print("Mechanism tags: reuse from mechanisms.json where they fit; a new")
    print("tag needs an entry there (field + description) or CI fails.")
    print("The TODO placeholders will fail the schema test until replaced —")
    print("that is deliberate, so a half-filled index cannot be committed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Query the cross-problem approach taxonomy.

    python scripts/mechanisms.py list [--field FIELD]
    python scripts/mechanisms.py where <tag>
    python scripts/mechanisms.py gaps <problem>
    python scripts/mechanisms.py matrix

`list`   — every tag, grouped by field, with where it has been used.
`where`  — every attempt that used a tag, with status and one-liner.
`gaps`   — the cross-pollination view for one problem: which fields and
           tags have been tried elsewhere but never on it. This is the
           input to an ideation sweep (docs/IDEATE.md).
`matrix` — problems x fields, tag counts per cell.

Reads mechanisms.json (the taxonomy) and problems/*/prior-art.json (the
usage). Tier 1 throughout: do not run this in a blind working copy.
Standard library only.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load():
    taxonomy = json.loads((ROOT / "mechanisms.json").read_text(encoding="utf-8"))
    usage: dict[str, list[dict]] = {}  # tag -> [{problem, id, status, one_line}]
    problems: list[str] = []
    for index in sorted(ROOT.glob("problems/*/prior-art.json")):
        data = json.loads(index.read_text(encoding="utf-8"))
        problems.append(data["problem"])
        for attempt in data["attempts"]:
            for tag in attempt["mechanism"]:
                usage.setdefault(tag, []).append(
                    {
                        "problem": data["problem"],
                        "id": attempt["id"],
                        "status": attempt["status"],
                        "one_line": attempt["one_line"],
                    }
                )
    return taxonomy, usage, problems


def cmd_list(args) -> int:
    taxonomy, usage, _ = load()
    fields = taxonomy["fields"]
    if args.field and args.field not in fields:
        print(f"error: no such field: {args.field}", file=sys.stderr)
        print(f"fields: {', '.join(sorted(fields))}", file=sys.stderr)
        return 2
    for field in sorted(fields):
        if args.field and field != args.field:
            continue
        tags = sorted(
            t for t, m in taxonomy["mechanisms"].items() if m["field"] == field
        )
        print(f"{field} — {fields[field]}")
        if not tags:
            print("  (no tags yet — wide-open lens)")
        for tag in tags:
            used_on = sorted({u["problem"] for u in usage.get(tag, [])})
            where = ", ".join(used_on) if used_on else "unused"
            print(f"  {tag}  [{where}]")
            print(f"    {taxonomy['mechanisms'][tag]['description']}")
        print()
    return 0


def cmd_where(args) -> int:
    taxonomy, usage, _ = load()
    if args.tag not in taxonomy["mechanisms"]:
        print(f"error: unknown tag: {args.tag}", file=sys.stderr)
        near = [t for t in taxonomy["mechanisms"] if args.tag in t]
        if near:
            print(f"did you mean: {', '.join(sorted(near))}", file=sys.stderr)
        return 2
    meta = taxonomy["mechanisms"][args.tag]
    print(f"{args.tag}  (field: {meta['field']})")
    print(f"  {meta['description']}")
    uses = usage.get(args.tag, [])
    if not uses:
        print("  never used in an attempt yet")
        return 0
    for u in uses:
        print(f"  {u['problem']}/{u['id']}  {u['status']}")
        print(f"    {u['one_line']}")
    return 0


def cmd_gaps(args) -> int:
    taxonomy, usage, problems = load()
    if args.problem not in problems:
        print(f"error: no such problem: {args.problem}", file=sys.stderr)
        print(f"available: {', '.join(problems)}", file=sys.stderr)
        return 2

    tried = {t for t, us in usage.items() if any(u["problem"] == args.problem for u in us)}
    tried_fields = {taxonomy["mechanisms"][t]["field"] for t in tried}
    attack_fields = {f for f in taxonomy["fields"] if f != "logic-methodology"}

    print(f"Cross-pollination gaps for {args.problem}")
    print(f"  tried here: {', '.join(sorted(tried)) or '(nothing yet)'}")
    print()
    print("Fields never tried on this problem (ideation lenses, widest first):")
    for field in sorted(attack_fields - tried_fields):
        elsewhere = sorted(
            {
                u["problem"]
                for t, us in usage.items()
                if taxonomy["mechanisms"][t]["field"] == field
                for u in us
            }
        )
        note = f"used on: {', '.join(elsewhere)}" if elsewhere else "unused anywhere"
        print(f"  {field}  ({note})")
        print(f"    {taxonomy['fields'][field]}")
    print()
    print("Tags proven useful elsewhere, never tried here:")
    for tag in sorted(set(taxonomy["mechanisms"]) - tried):
        us = usage.get(tag, [])
        if not us or taxonomy["mechanisms"][tag]["field"] == "logic-methodology":
            continue
        where = ", ".join(sorted({f"{u['problem']}/{u['id']}" for u in us}))
        print(f"  {tag}  [{where}]")
    print()
    print("Caveat: a gap is an untried lens, not a promising one. Check the")
    print("problem's prior-art.json gaps/dead-ends before spending a cycle.")
    return 0


def cmd_matrix(args) -> int:
    taxonomy, usage, problems = load()
    fields = sorted(taxonomy["fields"])
    counts: dict[tuple[str, str], int] = {}
    for tag, us in usage.items():
        field = taxonomy["mechanisms"][tag]["field"]
        for u in us:
            counts[(u["problem"], field)] = counts.get((u["problem"], field), 0) + 1

    width = max(len(p) for p in problems) + 2
    fw = 6
    header = " " * width + "".join(f[:fw - 1].rjust(fw) for f in fields)
    print(header)
    for problem in problems:
        row = problem.ljust(width)
        for field in fields:
            n = counts.get((problem, field), 0)
            row += (str(n) if n else "·").rjust(fw)
        print(row)
    print()
    print("columns: " + ", ".join(f"{f[:fw - 1]}={f}" for f in fields))
    print("cells: attempt-tag uses; · = untried lens (see `gaps <problem>`).")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("list", help="all tags grouped by field")
    p.add_argument("--field", help="restrict to one field")
    p.set_defaults(fn=cmd_list)

    p = sub.add_parser("where", help="attempts that used a tag")
    p.add_argument("tag")
    p.set_defaults(fn=cmd_where)

    p = sub.add_parser("gaps", help="untried fields/tags for one problem")
    p.add_argument("problem")
    p.set_defaults(fn=cmd_gaps)

    p = sub.add_parser("matrix", help="problems x fields tag counts")
    p.set_defaults(fn=cmd_matrix)

    args = parser.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:  # `mechanisms.py list | head` is a normal usage
        raise SystemExit(0)

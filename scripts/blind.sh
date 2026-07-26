#!/usr/bin/env bash
# Materialize a TIER-0 ONLY working copy of one problem.
#
#   scripts/blind.sh <problem> <destdir>
#
# The result contains the problem statement, the shared harness, and the
# contribution rules -- and physically none of this lab's prior art. An agent
# working in <destdir> cannot read what has already been tried, so its attempt
# is independent evidence rather than a variation on ours.
#
# For an INFORMED attempt, just work in a normal clone.

set -euo pipefail

usage() {
  echo "usage: scripts/blind.sh <problem> <destdir>" >&2
  echo >&2
  echo "problems:" >&2
  for d in "$(dirname "$0")/../problems"/*/; do echo "  $(basename "$d")" >&2; done
  exit 2
}

[ $# -eq 2 ] || usage

PROBLEM=$1
DEST=$2
ROOT=$(cd "$(dirname "$0")/.." && pwd)

[ -d "$ROOT/problems/$PROBLEM" ] || { echo "error: no such problem: $PROBLEM" >&2; usage; }
[ -e "$DEST" ] && { echo "error: $DEST already exists; refusing to overwrite" >&2; exit 1; }

mkdir -p "$DEST/problems/$PROBLEM"

# Tier 0: the statement, the harness, and the rules. Nothing else.
cp "$ROOT/problems/$PROBLEM/PROBLEM.md" "$DEST/problems/$PROBLEM/PROBLEM.md"
cp -r "$ROOT/harness" "$DEST/harness"
# NB: README.md is deliberately NOT copied -- it summarizes what the lab found,
# which is the one thing a blind copy must not contain. MODE (below) orients
# instead.
for f in AGENTS.md CONTRIBUTING.md CITATION.cff tiers.json; do
  [ -f "$ROOT/$f" ] && cp "$ROOT/$f" "$DEST/$f"
done
for f in "$ROOT"/LICENSE*; do [ -f "$f" ] && cp "$f" "$DEST/"; done

# Drop harness for other problems: less noise, and it keeps the copy honest
# about what this problem actually needs.
for d in "$DEST/harness"/*/; do
  name=$(basename "$d")
  [ "$name" = "common" ] && continue
  [ "$name" = "$PROBLEM" ] && continue
  rm -rf "$d"
done

cat > "$DEST/MODE" <<EOF
mode: blind
problem: $PROBLEM
source-commit: $(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)

This is a TIER-0 ONLY working copy. It deliberately excludes math-lab's prior
art for this problem: previous attempts, dead ends, route-specific tooling and
the editorial attack surface.

Do not go looking for the excluded material -- the point of this copy is that
your attempt is independent. If you end up consulting it anyway, record
'mode: informed' in your attempt, which is a perfectly good outcome, just a
different one.

When you submit, put 'mode: blind' in your attempt front-matter along with the
source-commit above.
EOF

echo "Blind working copy for '$PROBLEM' at: $DEST"
echo "Contains: PROBLEM.md, harness/$PROBLEM, harness/common, and the rules."
echo "Excludes: attempts/, PRIOR-ART.md, prior-art.json, explore/, data/, STATUS.md, GUIDANCE.md"

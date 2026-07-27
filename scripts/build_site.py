#!/usr/bin/env python3
"""Build the math-lab site.

Pages are generated from the repository itself: the per-problem explainers in
site/explainers/, and the machine-readable prior-art.json indexes. Nothing
about an attempt's status is retyped here, so the site cannot quietly drift
from the records it describes.

    python scripts/build_site.py [--out _site]

Standard library only, deliberately: the site should build from a bare clone
with no install step.
"""

from __future__ import annotations

import argparse
import html
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXPLAINERS = ROOT / "site" / "explainers"
REPO_URL = "https://github.com/Joe975/math-lab"
BLOB = f"{REPO_URL}/blob/main"
ISSUES = f"{REPO_URL}/issues"


def issue_url(template: str, **prefill: str) -> str:
    """A GitHub 'new issue' link with the form fields already filled in.

    Reporting an error should be one click from the thing that is wrong, not a
    scavenger hunt through the repo. Field names must match the ids in
    .github/ISSUE_TEMPLATE/<template>.
    """
    from urllib.parse import urlencode

    query = urlencode({"template": template, **prefill})
    return f"{ISSUES}/new?{query}"

# Presentation for each status term. The vocabulary itself is defined in
# AGENTS.md; this only says how to colour it, and how to explain it.
#
# Every term here is rendered on the Method page as its own badge beside its
# definition, so a reader who meets a badge on a card has somewhere to look it
# up. Adding a term to the vocabulary without a gloss fails the test suite.
STATUS = {
    "VERIFIED": ("verified", "Verified"),
    "VERIFIED_WITH_CORRECTIONS": ("corrected", "Verified, with corrections"),
    "EVIDENCE": ("evidence", "Evidence"),
    "LIVE": ("live", "Live"),
    "REFUTED": ("refuted", "Refuted"),
    "REFUTED_PRIOR_OBSERVATION": ("refuted", "Refuted a prior observation"),
    "MAP": ("map", "Map"),
    "SPECULATION": ("speculation", "Speculation"),
    "NOT_ATTEMPTED": ("none", "Not attempted"),
}

STATUS_GLOSS = {
    "VERIFIED": "Independently re-derived or re-computed, always with a stated "
    "range or scope. For a finite search this describes the range checked, "
    "never the conjecture.",
    "VERIFIED_WITH_CORRECTIONS": "Re-derived independently, and the review that "
    "did it also found errors in the record it was checking. The original is "
    "kept exactly as written; the corrections are listed against it. Read the "
    "corrections before building on the claim.",
    "EVIDENCE": "Computational support, not proof. A search to 10<sup>6</sup> is "
    "evidence about 10<sup>6</sup>.",
    "LIVE": "A route that has survived every adversary tried, with no bound yet "
    "claimed. The most fragile-sounding status and the most interesting one.",
    "MAP": "A survey of the terrain rather than an attack on it: which barriers "
    "sit where, what the small cases look like, which ideas are worth queueing. "
    "Claims nothing about the conjecture itself.",
    "SPECULATION": "Plausible, unproven, and load-bearing. Labelled inline at the "
    "step that needs it, not in a footnote.",
    "REFUTED": "Killed, with the obstruction recorded, so nobody re-treads it "
    "blindly.",
    "REFUTED_PRIOR_OBSERVATION": "Killed one of our own earlier observations "
    "rather than a route. The record that made the observation still stands, "
    "with the refutation recorded against it.",
    "NOT_ATTEMPTED": "Queued and never worked. Recorded so the absence is "
    "visible, instead of looking like a gap in the reporting.",
}


# --- a deliberately small markdown subset ---------------------------------
# Only what the explainers use: ## headings, paragraphs, **bold**, *italic*,
# `code` and [links](url). Attempt records are not rendered here -- they are
# long and GitHub already renders them well -- so nothing more is needed.


def esc(text: str) -> str:
    """Escape for text content. Attribute values use html.escape() directly."""
    return html.escape(text, quote=False)


def inline(text: str) -> str:
    out = html.escape(text, quote=False)
    out = re.sub(r"`([^`]+)`", r"<code>\1</code>", out)
    out = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", out)
    out = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", out)
    out = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', out)
    return out


def markdown(body: str) -> str:
    blocks, buf = [], []

    def flush():
        if buf:
            blocks.append("<p>" + inline(" ".join(buf)) + "</p>")
            buf.clear()

    for line in body.splitlines():
        stripped = line.strip()
        if not stripped:
            flush()
        elif stripped.startswith("## "):
            flush()
            blocks.append(f"<h2>{inline(stripped[3:])}</h2>")
        elif stripped.startswith("### "):
            flush()
            blocks.append(f"<h3>{inline(stripped[4:])}</h3>")
        else:
            buf.append(stripped)
    flush()
    return "\n".join(blocks)


def parse_explainer(path: Path) -> dict:
    raw = path.read_text(encoding="utf-8")
    if not raw.startswith("---"):
        raise ValueError(f"{path.name}: missing front matter")
    _, front, body = raw.split("---", 2)
    meta = {}
    for line in front.strip().splitlines():
        key, _, value = line.partition(":")
        meta[key.strip()] = value.strip()
    for required in ("title", "short", "tagline"):
        if required not in meta:
            raise ValueError(f"{path.name}: front matter missing {required!r}")
    meta["slug"] = path.stem
    meta["body"] = body.strip()
    return meta


# --- page furniture --------------------------------------------------------

CSS = """
:root {
  --bg: #fbfaf8; --fg: #1c1a17; --muted: #6b6560; --rule: #e2ddd6;
  --card: #ffffff; --accent: #8a5a2b; --shadow: 0 1px 2px rgba(0,0,0,.05);
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #16151a; --fg: #e8e4de; --muted: #9c958c; --rule: #2e2c33;
    --card: #1e1d23; --accent: #d3a06a; --shadow: none;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--bg); color: var(--fg);
  font: 17px/1.65 ui-serif, Georgia, "Iowan Old Style", serif;
  -webkit-text-size-adjust: 100%;
}
.wrap { max-width: 46rem; margin: 0 auto; padding: 0 1.25rem; }
header.site { border-bottom: 1px solid var(--rule); margin-bottom: 3rem; }
header.site .wrap { display: flex; gap: 1.5rem; align-items: center;
  padding-top: .35rem; padding-bottom: .35rem; flex-wrap: wrap; }
header.site a { color: var(--fg); text-decoration: none; }
header.site .brand { font-weight: 700; letter-spacing: -.01em;
  display: flex; align-items: center; min-height: 44px; }
header.site nav { margin-left: auto; display: flex; gap: 1.25rem;
  font: 500 .82rem/1 ui-sans-serif, system-ui, sans-serif; }
/* Tap targets, not just text: at .82rem/1 these were 13px tall. */
header.site nav a { color: var(--muted); display: flex; align-items: center;
  min-height: 44px; }
header.site nav a:hover { color: var(--accent); }
h1 { font-size: 2.1rem; line-height: 1.15; letter-spacing: -.02em; margin: 0 0 .6rem; }
h2 { font-size: 1.28rem; letter-spacing: -.01em; margin: 2.6rem 0 .7rem;
     padding-top: 1.1rem; border-top: 1px solid var(--rule); }
h3 { font-size: 1.05rem; margin: 1.8rem 0 .5rem; }
p { margin: 0 0 1.05rem; }
a { color: var(--accent); }
code { font: .86em ui-monospace, SFMono-Regular, Menlo, monospace;
  background: color-mix(in srgb, var(--muted) 13%, transparent);
  padding: .12em .34em; border-radius: 3px; }
.lede { font-size: 1.16rem; color: var(--muted); margin-bottom: 1.8rem; }
.meta { font: .8rem/1.5 ui-sans-serif, system-ui, sans-serif; color: var(--muted);
  margin-bottom: 2.2rem; }
.callout { border-left: 3px solid var(--accent); padding: .1rem 0 .1rem 1.1rem;
  margin: 1.8rem 0; color: var(--muted); }
.callout strong { color: var(--fg); }

/* problem cards on the index */
.cards { display: grid; gap: .9rem; margin: 1.5rem 0 0; padding: 0; list-style: none; }
.cards li { background: var(--card); border: 1px solid var(--rule); border-radius: 9px;
  padding: 1.1rem 1.2rem; box-shadow: var(--shadow); }
.cards a.title { font-weight: 700; text-decoration: none; font-size: 1.06rem; }
.cards a.title:hover { text-decoration: underline; }
.cards p { margin: .4rem 0 0; font-size: .94rem; color: var(--muted); }

/* approach report cards */
.report { background: var(--card); border: 1px solid var(--rule); border-radius: 9px;
  padding: 1.15rem 1.25rem; margin: 1rem 0; box-shadow: var(--shadow); }
.report h3 { margin: 0 0 .1rem; font-size: 1.02rem; }
.report .one-line { margin: .55rem 0 0; }
.report dl { margin: .9rem 0 0; font-size: .9rem; }
.report dt { font: 600 .72rem/1.4 ui-sans-serif, system-ui, sans-serif;
  text-transform: uppercase; letter-spacing: .06em; color: var(--muted); margin-top: .8rem; }
.report dd { margin: .15rem 0 0; }
.report ul { margin: .15rem 0 0; padding-left: 1.1rem; }
.report .actions { margin: 1rem 0 0; padding-top: .8rem; border-top: 1px solid var(--rule);
  font: .84rem/1.5 ui-sans-serif, system-ui, sans-serif; }

/* Colours here are contrast-checked against their own composited background by
   tests/test_site.py -- the text is .68rem, so it needs the 4.5:1 rule. */
.badge { display: inline-block; font: 600 .68rem/1.5 ui-sans-serif, system-ui, sans-serif;
  text-transform: uppercase; letter-spacing: .06em; padding: .12rem .5rem;
  border-radius: 20px; vertical-align: .12em; }
.badge.verified  { background: #1f7a4d1f; color: #176b42; }
/* Verified, but asterisked: the green of a pass, ringed in the amber of the
   corrections. Must not read as an unqualified VERIFIED at a glance. */
.badge.corrected { background: #1f7a4d14; color: #176b42;
  box-shadow: inset 0 0 0 1px #b4790066; }
.badge.evidence  { background: #2563eb1f; color: #2560c8; }
.badge.live      { background: #b4790016; color: #8f6200; }
.badge.refuted   { background: #c0392b1a; color: #b0392e; }
/* Three different meanings that used to share one grey pill. Kept colourless
   on purpose -- they are distinguished by fill, outline and dash instead. */
.badge.map { background: color-mix(in srgb, var(--muted) 14%, transparent);
  color: var(--fg); }
.badge.none { color: var(--muted); box-shadow: inset 0 0 0 1px var(--rule); }
.badge.speculation { color: var(--accent);
  box-shadow: inset 0 0 0 1px color-mix(in srgb, var(--accent) 45%, transparent); }
@media (prefers-color-scheme: dark) {
  .badge.verified  { background: #4ade8018; color: #6ee7a8; }
  .badge.corrected { background: #4ade8010; color: #6ee7a8;
    box-shadow: inset 0 0 0 1px #fbbf2455; }
  .badge.evidence  { background: #60a5fa18; color: #8ebbfb; }
  .badge.live      { background: #fbbf2418; color: #f0c274; }
  .badge.refuted   { background: #f8717118; color: #f79f9f; }
}
/* Status vocabulary on the Method page: the badge itself is the term. */
.glossary { margin: 1.4rem 0 0; }
.glossary dt { margin-top: 1.1rem; }
.glossary dd { margin: .35rem 0 0; }
@media (min-width: 40rem) {
  .glossary { display: grid; grid-template-columns: 11.5rem 1fr; gap: .9rem 1.2rem; }
  .glossary dt { margin-top: .1rem; text-align: right; }
  .glossary dd { margin: 0; }
}

/* Jumping to an attempt from a cross-reference must land somewhere obvious. */
.report:target { border-color: var(--accent); box-shadow: 0 0 0 2px
  color-mix(in srgb, var(--accent) 22%, transparent); scroll-margin-top: 1.5rem; }

.tags { margin: .7rem 0 0; display: flex; flex-wrap: wrap; gap: .35rem; }
.tag { font: .72rem/1.5 ui-monospace, monospace; color: var(--muted);
  border: 1px solid var(--rule); border-radius: 4px; padding: .05rem .38rem; }

footer.site { margin: 4.5rem 0 0; border-top: 1px solid var(--rule); }
footer.site .wrap { padding: 1.6rem 1.25rem 3rem;
  font: .84rem/1.6 ui-sans-serif, system-ui, sans-serif; color: var(--muted); }
.scroll { overflow-x: auto; }
"""


def page(title: str, body: str, description: str) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<meta name="description" content="{html.escape(description)}">
<style>{CSS}</style>
</head>
<body>
<header class="site"><div class="wrap">
  <a class="brand" href="index.html">Math Lab</a>
  <nav>
    <a href="index.html">Problems</a>
    <a href="method.html">Method</a>
    <a href="{REPO_URL}">Repository</a>
  </nav>
</div></header>
<main class="wrap">
{body}
</main>
<footer class="site"><div class="wrap">
  <p>An agent-driven library of attempted approaches to open conjectures.
  No conjecture here is solved. Every record was written and reviewed by AI
  agents and has not been peer-reviewed by a human mathematician.</p>
  <p>MIT licensed. <a href="{REPO_URL}">Source on GitHub</a>.</p>
</div></footer>
</body>
</html>
"""


def badge(status: str) -> str:
    cls, label = STATUS.get(status, ("none", status.replace("_", " ").title()))
    return f'<span class="badge {cls}">{esc(label)}</span>'


# --- rendering -------------------------------------------------------------


def cross_ref(text: str, on_page: set[str]) -> str:
    """Link any attempt id in `text` to that attempt's card on this page.

    The records cite each other constantly -- 004 verifies 003, 003 is
    superseded by 004 -- and following that graph is most of what a reader
    does here, so it has to be clickable rather than merely stated.
    """
    return re.sub(
        r"\b(\d{3})\b",
        lambda m: (
            f'<a href="#attempt-{m.group(1)}">{m.group(1)}</a>'
            if m.group(1) in on_page
            else m.group(1)
        ),
        inline(text),
    )


def render_attempt(attempt: dict, slug: str, on_page: set[str]) -> str:
    ident = attempt["id"]
    href = f"{BLOB}/problems/{slug}/{attempt['file']}"
    parts = [
        f'<div class="report" id="attempt-{html.escape(ident)}">',
        f'<h3><a href="{href}">Attempt {esc(ident)}</a> {badge(attempt["status"])}</h3>',
        f'<p class="one-line">{inline(attempt["one_line"])}</p>',
    ]

    tags = attempt.get("mechanism", [])
    if tags:
        parts.append(
            '<div class="tags">'
            + "".join(f'<span class="tag">{html.escape(t)}</span>' for t in tags)
            + "</div>"
        )

    rows: list[tuple[str, str]] = []
    if attempt.get("range"):
        rows.append(("Range checked", inline(attempt["range"])))
    if attempt.get("refutes"):
        rows.append(("Refutes", cross_ref(attempt["refutes"], on_page)))
    if attempt.get("verifies"):
        rows.append(("Verifies", f"attempt {cross_ref(attempt['verifies'], on_page)}"))
    if attempt.get("superseded_by"):
        rows.append(
            (
                "Superseded by",
                f"attempt {cross_ref(attempt['superseded_by'], on_page)}",
            )
        )
    for key, label in (
        ("corrections", "Corrections"),
        ("gaps", "Open gaps"),
    ):
        items = attempt.get(key) or []
        if items:
            rows.append(
                (
                    label,
                    "<ul>"
                    + "".join(f"<li>{inline(str(i))}</li>" for i in items)
                    + "</ul>",
                )
            )
    for key, label in (
        ("methodological_lesson", "Methodological lesson"),
        ("design_note", "Design note"),
        ("note", "Note"),
    ):
        if attempt.get(key):
            rows.append((label, inline(attempt[key])))
    rows.append(("Mode", esc(attempt.get("mode", "unknown"))))

    if rows:
        parts.append("<dl>")
        for label, value in rows:
            parts.append(f"<dt>{esc(label)}</dt><dd>{value}</dd>")
        parts.append("</dl>")

    challenge = issue_url(
        "challenge.yml",
        record=f"{slug}/{ident}",
        title=f"[challenge] {slug}/{ident}: ",
    )
    parts.append(
        f'<p class="actions"><a href="{href}">Read the full record</a> '
        f'&middot; <a href="{challenge}">Challenge this result</a></p>'
    )
    parts.append("</div>")
    return "\n".join(parts)


def render_problem(meta: dict, index: dict) -> str:
    slug = meta["slug"]
    attempts = index["attempts"]

    body = [
        f"<h1>{esc(meta['title'])}</h1>",
        f'<p class="lede">{inline(meta["tagline"])}</p>',
    ]
    if meta.get("posed"):
        body.append(f'<p class="meta">Posed by {esc(meta["posed"])}.</p>')

    body.append(markdown(meta["body"]))

    body.append("<h2>What this lab tried</h2>")
    body.append(
        f'<p>Standing: {badge(index["route_status"])} '
        f'{inline(index["route_summary"])}</p>'
    )
    if attempts:
        body.append(
            "<p>Every approach below is recorded in full, including the ones that "
            "failed and why. Attempt titles link to the complete record.</p>"
        )
        on_page = {a["id"] for a in attempts}
        body.extend(render_attempt(a, slug, on_page) for a in attempts)
    else:
        body.append(
            "<p>No approaches have been recorded for this problem yet, so there is "
            "nothing here to build on or avoid.</p>"
        )

    body.append("<h2>Go deeper</h2>")
    body.append(
        "<p>"
        f'<a href="{BLOB}/problems/{slug}/PROBLEM.md">Problem statement and '
        "published background</a> &middot; "
        f'<a href="{BLOB}/problems/{slug}/PRIOR-ART.md">Full prior-art notes</a> '
        f'&middot; <a href="{BLOB}/problems/{slug}/prior-art.json">Machine-readable '
        "index</a></p>"
    )
    report = issue_url("attempt.yml", problem=slug, title=f"[attempt] {slug}: ")
    challenge_problem = issue_url("challenge.yml", record=slug, title=f"[challenge] {slug}: ")
    body.append("<h2>Found something, or think we got it wrong?</h2>")
    body.append(
        "<p>Both are welcome, and the second especially. Nothing here has been "
        "peer-reviewed by a human mathematician, so errors are expected rather "
        "than embarrassing — one record in this library was already corrected by "
        "a later review.</p>"
        f'<p><a href="{challenge_problem}">Challenge something on this page</a> '
        f'&middot; <a href="{report}">Report an attempt you ran</a> '
        f'&middot; <a href="{ISSUES}">Browse open issues</a></p>'
    )
    body.append(
        '<div class="callout"><strong>Attacking this yourself?</strong> '
        f'<code>scripts/blind.sh {slug} ../work</code> gives your agent the problem '
        "statement and verification harness with everything on this page removed, so "
        "its attempt is independent. See <a href=\"method.html\">Method</a>.</div>"
    )

    return page(
        f"{meta['title']} — Math Lab",
        "\n".join(body),
        meta["tagline"],
    )


def render_index(problems: list[tuple[dict, dict]]) -> str:
    cards = []
    for meta, index in problems:
        cards.append(
            "<li>"
            f'<a class="title" href="{meta["slug"]}.html">{esc(meta["title"])}</a> '
            f'{badge(index["route_status"])}'
            f'<p>{inline(meta["tagline"])}</p>'
            "</li>"
        )

    attempt_count = sum(len(i["attempts"]) for _, i in problems)
    challenge_any = issue_url("challenge.yml")
    report_any = issue_url("attempt.yml")

    body = f"""
<h1>Math Lab</h1>
<p class="lede">A library of attempted approaches to seven open mathematical
conjectures — including, especially, the ones that failed.</p>

<p>Most research writing records what worked. The dead ends stay in people's
heads, so the same doomed idea gets tried again and again. This is an attempt at
the opposite: an autonomous agent loop worked on these problems and wrote down
every line it pursued, why it was chosen, and precisely where it broke.</p>

<div class="callout"><strong>No conjecture here is solved, and none is close to
solved.</strong> That is the expected outcome. The deliverables are the approach
library, the tooling, and progress on the smaller open problems.</div>

<h2>The problems</h2>
<p>{len(problems)} problems, {attempt_count} recorded approaches. Each page
explains the conjecture in plain terms, says honestly what a breakthrough would
and would not mean, and reports on every approach tried.</p>
<ul class="cards">
{"".join(cards)}
</ul>

<h2>How the results were produced</h2>
<p>An hourly cycle: read the ledger, pick two to four attack lines, fan out
parallel agents, then try hard to <em>refute</em> whatever they claim. Skeptic
agents default to disbelief, and computational claims are re-run by a separate
implementation. Only survivors are recorded as results.</p>

<p>That process visibly caught its own errors, which is the main reason to trust
the rest of it. One skeptic pass found a false proof claim and a wrong headline
constant in the attempt it was reviewing; both records are kept as written, with
the correction recorded against them rather than edited in.
<a href="method.html">More on the method</a>.</p>

<h2>Report an error, or tell us what you tried</h2>
<p>Every record here was written and reviewed by AI agents, never by a human
mathematician. Challenging one is the most useful thing you can do with this
library, and you do not need to be certain to open an issue.</p>
<p><a href="{challenge_any}">Challenge a recorded result</a> &middot;
<a href="{report_any}">Report an attempt you ran</a> &middot;
<a href="{ISSUES}">Browse open issues</a></p>

<h2>Send your own agent</h2>
<p>The library is built to be used by other people's agents, in either of two
modes: <strong>blind</strong>, working from the problem statement alone with our
prior art physically removed, or <strong>informed</strong>, building on what is
here. Which one you pick changes what your result means.
<a href="method.html">How to do it</a>.</p>
"""
    return page(
        "Math Lab — attempted approaches to open conjectures",
        body,
        "A library of attempted approaches to seven open mathematical conjectures, "
        "including the failed ones, with adversarial verification.",
    )


def render_glossary() -> str:
    """The status vocabulary, every term, each shown as the badge you will meet.

    Generated from STATUS rather than written out, because a hand-written
    legend drifts: this one documented five of the eight terms the site was
    actually rendering, and one term it never rendered at all.
    """
    rows = [
        f'<dt>{badge(term)}</dt><dd>{STATUS_GLOSS[term]}</dd>' for term in STATUS
    ]
    return '<dl class="glossary">' + "".join(rows) + "</dl>"


def render_method(problems: list[tuple[dict, dict]]) -> str:
    body = """
<h1>Method</h1>
<p class="lede">How the records were produced, what the status words mean, and
how to point your own agent at the library.</p>

<h2>The loop</h2>
<p>An hourly autonomous cycle. Each pass read the ledger of current state and the
standing human guidance, pulled two to four attack lines from a queue, and fanned
out parallel agents — computational searches, proof attempts from a named angle,
structure surveys. Cycles were kept deliberately small so the loop could run for
a month; depth came from accumulation rather than from any single large run.</p>
<p>No loop is running now, but it is meant to be picked up rather than
reconstructed: <code>docs/CYCLE.md</code> in the repository is the operational
prompt, and the ledger carries a queue of concrete starting points. A cycle can
be run by hand, one line at a time, or restarted on a timer.</p>

<h2>Adversarial verification</h2>
<p>Nothing became a result by being claimed. Every load-bearing claim faced an
independent skeptic agent whose default stance was to <em>refute</em> it, and
every computational claim was re-run by a separate implementation rather than
re-run by the same code.</p>
<p>This is the part that does real work. In the union-closed line, the skeptic
pass confirmed the main interface was sound but found a proof claim that was
false as stated, a headline constant that was wrong because a search grid had
excluded the extremal case, and a data file that could not be reproduced. The
original record was not edited. Both stand, with the relationship between them
recorded, because how often the loop corrects itself is evidence about how much
to trust it.</p>

<h2 id="status-words">What the status words mean</h2>
<p>These are used exactly, and the machine-readable indexes are filtered on
them. Every badge you will meet anywhere on this site is listed here.</p>
__GLOSSARY__

<h2>Blind and informed</h2>
<p>The repository is cut in two. <strong>Tier 0</strong> is published background
only: the problem statement, the literature, and the shared verification harness.
<strong>Tier 1</strong> is everything this lab found — the attempts, the dead
ends, the route-specific tooling, and the editorial view of where to push.</p>

<p>Running <code>scripts/blind.sh &lt;problem&gt; ../work</code> produces a
working copy containing only tier 0. An agent working there cannot re-tread our
dead ends because it cannot see them, and cannot be anchored by our framing. An
informed agent works in a normal clone and reads the machine-readable index
first to check whether its idea's mechanism family is already spent.</p>

<p>The cut is enforced rather than trusted: it is defined in
<code>tiers.json</code>, and a test suite fails the build if a finding of ours
appears in a tier-0 file or in a blind checkout. Each recorded attempt declares
the distinctive terms naming its findings, and those are what the check runs
on.</p>

<h2>Why both modes exist</h2>
<p>Every attempt records which mode produced it. Over time that turns the library
into a dataset on a question nobody has clean data for: <strong>does prior
knowledge help, or does it anchor?</strong> If blind agents keep rediscovering
recorded dead ends, the prior art is load-bearing. If they keep finding things
informed agents were steered away from, the framing here is a local optimum.
Either answer is worth having, which is why a mislabelled mode is worse than no
attempt at all.</p>

<h2>What to trust</h2>
<p>Every record was written by an AI agent and reviewed by other AI agents. None
has been peer-reviewed by a human mathematician. Read <em>verified</em> as "two
independent machine implementations agree" — a real guarantee, and a limited
one — and read anything labelled <em>speculation</em> as exactly that.</p>
"""
    return page(
        "Method — Math Lab",
        body.replace("__GLOSSARY__", render_glossary()),
        "How the records were produced, what the status vocabulary means, and how "
        "to point your own agent at the library in blind or informed mode.",
    )


# --- entry point -----------------------------------------------------------


def load() -> list[tuple[dict, dict]]:
    problems = []
    for path in sorted(EXPLAINERS.glob("*.md")):
        meta = parse_explainer(path)
        meta["order"] = int(meta.get("order", 99))
        index_path = ROOT / "problems" / meta["slug"] / "prior-art.json"
        if not index_path.exists():
            raise SystemExit(f"{path.name}: no problem directory for {meta['slug']!r}")
        problems.append((meta, json.loads(index_path.read_text(encoding="utf-8"))))
    problems.sort(key=lambda pair: (pair[0]["order"], pair[0]["title"]))
    return problems


def build(out: Path) -> list[Path]:
    problems = load()
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    written = []
    for meta, index in problems:
        target = out / f"{meta['slug']}.html"
        target.write_text(render_problem(meta, index), encoding="utf-8")
        written.append(target)

    (out / "index.html").write_text(render_index(problems), encoding="utf-8")
    (out / "method.html").write_text(render_method(problems), encoding="utf-8")
    (out / ".nojekyll").write_text("", encoding="utf-8")
    written += [out / "index.html", out / "method.html"]
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(ROOT / "_site"), help="output directory")
    args = parser.parse_args()
    written = build(Path(args.out))
    print(f"built {len(written)} pages into {args.out}")


if __name__ == "__main__":
    main()

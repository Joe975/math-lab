"""Keep the site honest about the records it describes.

The site is generated from prior-art.json rather than hand-written, so its
statuses cannot drift from the library. These tests cover what generation alone
does not guarantee: that every problem and every attempt actually reaches a
page, that no page links somewhere that does not exist, and that the site never
claims a conjecture has been solved.
"""

from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import build_site  # noqa: E402


def _split_at_rules(css: str) -> tuple[str, list[tuple[str, str]]]:
    """Separate top-level declarations from @media blocks, by brace counting.

    Crude on purpose -- the stylesheet is one literal in build_site.py -- but
    it has to respect nesting: matching @media with a regex silently folded the
    dark palette into the light one and made a contrast test read the wrong
    colours.
    """
    css = re.sub(r"/\*.*?\*/", "", css, flags=re.S)
    top, blocks, i = [], [], 0
    while i < len(css):
        at = css.find("@media", i)
        if at < 0:
            top.append(css[i:])
            break
        top.append(css[i:at])
        open_brace = css.index("{", at)
        depth, j = 1, open_brace + 1
        while depth:
            depth += {"{": 1, "}": -1}.get(css[j], 0)
            j += 1
        blocks.append((css[at:open_brace].strip(), css[open_brace + 1 : j - 1]))
        i = j
    return "".join(top), blocks


def _declarations(css: str) -> dict[str, str]:
    rules: dict[str, str] = {}
    for selectors, body in re.findall(r"([^{}]+)\{([^{}]*)\}", css):
        for selector in selectors.split(","):
            selector = " ".join(selector.split())
            if selector:
                rules[selector] = " ".join(body.split())
    return rules


def _css_rules() -> dict[str, str]:
    """Top-level declarations by selector -- the light palette, no @media."""
    top, _ = _split_at_rules(build_site.CSS)
    return _declarations(top)


def _dark_rules() -> dict[str, str]:
    """Declarations inside the prefers-color-scheme: dark blocks."""
    _, blocks = _split_at_rules(build_site.CSS)
    rules: dict[str, str] = {}
    for condition, body in blocks:
        if "prefers-color-scheme: dark" in condition:
            rules.update(_declarations(body))
    return rules


def _rgba(value: str) -> tuple[float, float, float, float]:
    """#rgb / #rrggbb / #rrggbbaa -> channels in 0..1 plus alpha."""
    h = value.strip().lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)
    alpha = int(h[6:8], 16) / 255 if len(h) == 8 else 1.0
    return (*(int(h[i : i + 2], 16) / 255 for i in (0, 2, 4)), alpha)


def _contrast(fg: str, bg_layers: list[str]) -> float:
    """WCAG contrast of fg over bg_layers composited front-to-back."""

    def luminance(channels):
        srgb = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in channels]
        return 0.2126 * srgb[0] + 0.7152 * srgb[1] + 0.0722 * srgb[2]

    *base, _ = _rgba(bg_layers[-1])
    for layer in reversed(bg_layers[:-1]):
        *top, alpha = _rgba(layer)
        base = [t * alpha + b * (1 - alpha) for t, b in zip(top, base)]
    light, dark = sorted((luminance(_rgba(fg)[:3]), luminance(base)), reverse=True)
    return (light + 0.05) / (dark + 0.05)


@pytest.fixture(scope="module")
def site(tmp_path_factory) -> dict[str, str]:
    out = tmp_path_factory.mktemp("site")
    build_site.build(out)
    return {p.name: p.read_text(encoding="utf-8") for p in out.glob("*.html")}


def problems() -> list[str]:
    return sorted(p.name for p in (ROOT / "problems").iterdir() if p.is_dir())


def test_every_problem_has_an_explainer():
    """A problem with no explainer would be invisible on the site."""
    explained = {p.stem for p in (ROOT / "site" / "explainers").glob("*.md")}
    missing = sorted(set(problems()) - explained)
    assert not missing, f"problems with no site/explainers/<slug>.md: {missing}"


def test_every_problem_has_a_page(site):
    for slug in problems():
        assert f"{slug}.html" in site, f"no page generated for {slug}"


def test_every_attempt_appears_on_its_page(site):
    """The reason the site exists: a report on each approach."""
    missing = []
    for slug in problems():
        index = json.loads(
            (ROOT / "problems" / slug / "prior-art.json").read_text(encoding="utf-8")
        )
        page = site[f"{slug}.html"]
        for attempt in index["attempts"]:
            if f"Attempt {attempt['id']}" not in page:
                missing.append(f"{slug}/{attempt['id']}")
            # Compare against escaped text: summaries contain <= and >=.
            snippet = html.escape(attempt["one_line"])[:40]
            if snippet not in page:
                missing.append(f"{slug}/{attempt['id']} (summary)")
    assert not missing, f"attempts absent from their page: {missing}"


def test_status_badges_match_the_index(site):
    """A status shown on the site must be the status in the record."""
    wrong = []
    for slug in problems():
        index = json.loads(
            (ROOT / "problems" / slug / "prior-art.json").read_text(encoding="utf-8")
        )
        page = site[f"{slug}.html"]
        for attempt in index["attempts"]:
            _, label = build_site.STATUS[attempt["status"]]
            if label not in page:
                wrong.append(f"{slug}/{attempt['id']}: expected badge {label!r}")
    assert not wrong, "\n".join(wrong)


def test_every_status_term_is_renderable():
    """An unmapped status would render as a raw enum on a public page."""
    unmapped = set()
    for slug in problems():
        index = json.loads(
            (ROOT / "problems" / slug / "prior-art.json").read_text(encoding="utf-8")
        )
        for status in [index["route_status"]] + [
            a["status"] for a in index["attempts"]
        ]:
            if status not in build_site.STATUS:
                unmapped.add(status)
    assert not unmapped, f"statuses with no presentation mapping: {sorted(unmapped)}"


def test_no_dead_internal_links(site):
    """Relative links must resolve to a page that was actually generated."""
    dead = []
    for name, body in site.items():
        for href in re.findall(r'href="([^"]+)"', body):
            if href.startswith(("http://", "https://", "#", "mailto:")):
                continue
            if href not in site and not href.endswith(".nojekyll"):
                dead.append(f"{name} -> {href}")
    assert not dead, f"dead internal links: {dead}"


def test_site_does_not_claim_a_solution(site):
    """The repo's credibility rests on not overclaiming. Guard it mechanically."""
    forbidden = [
        "we solved",
        "we have solved",
        "we prove the conjecture",
        "proves the conjecture",
        "solves the conjecture",
        "breakthrough result",
    ]
    hits = []
    for name, body in site.items():
        lowered = body.lower()
        for phrase in forbidden:
            if phrase in lowered:
                hits.append(f"{name}: {phrase!r}")
    assert not hits, f"site overclaims: {hits}"


def test_index_carries_the_honest_framing(site):
    """The 'nothing is solved' disclaimer is load-bearing, not decoration."""
    assert "No conjecture here is solved" in site["index.html"]
    assert "not been peer-reviewed by a human" in site["index.html"]


def test_explainer_front_matter_is_complete():
    for path in sorted((ROOT / "site" / "explainers").glob("*.md")):
        meta = build_site.parse_explainer(path)
        for field in ("title", "short", "tagline"):
            assert meta[field], f"{path.name}: empty {field}"
        assert "## What a breakthrough would mean" in meta["body"], (
            f"{path.name}: missing the implications section, which is the whole "
            f"point of the explainer"
        )


def test_every_page_offers_a_way_to_report_a_problem(site):
    """Reporting an error must be one click from the thing that is wrong."""
    for name, body in site.items():
        if name == "method.html":
            continue
        assert "issues/new?" in body, f"{name} has no issue link"


def test_every_attempt_card_links_to_a_prefilled_challenge(site):
    """The challenge link must name the specific record it is about."""
    import json as _json

    missing = []
    for slug in problems():
        index = _json.loads(
            (ROOT / "problems" / slug / "prior-art.json").read_text(encoding="utf-8")
        )
        page = site[f"{slug}.html"]
        for attempt in index["attempts"]:
            token = f"record={slug}%2F{attempt['id']}"
            if token not in page:
                missing.append(f"{slug}/{attempt['id']}")
    assert not missing, f"attempts with no prefilled challenge link: {missing}"


def test_every_status_term_is_explained_on_the_method_page(site):
    """A badge with nowhere to look it up is a word the reader cannot use.

    The hand-written legend this replaced documented five of the eight terms
    the site was rendering, plus one it never rendered.
    """
    method = site["method.html"]
    missing = [
        term
        for term, (_, label) in build_site.STATUS.items()
        if label not in method or term not in build_site.STATUS_GLOSS
    ]
    assert not missing, f"status terms absent from the Method glossary: {missing}"


def test_the_glossary_shows_the_badge_not_just_the_word(site):
    """A legend with the swatches removed is not a legend."""
    method = site["method.html"]
    for term, (cls, _) in build_site.STATUS.items():
        assert f'<span class="badge {cls}">' in method, (
            f"{term} is glossed on the Method page but its badge is not shown"
        )


def test_neutral_statuses_are_visually_distinguishable(site):
    """MAP, NOT_ATTEMPTED and SPECULATION mean different things.

    They shared one identical grey pill, so three meanings looked the same.
    """
    classes = {cls for cls, _ in build_site.STATUS.values()}
    assert {"map", "none", "speculation"} <= classes
    rules = _css_rules()
    styles = [rules[f".badge.{c}"] for c in ("map", "none", "speculation")]
    assert len(set(styles)) == 3, f"neutral badges share styling: {styles}"


def test_corrections_do_not_look_like_an_unqualified_pass(site):
    """VERIFIED_WITH_CORRECTIONS rendered in the same green as VERIFIED."""
    cls, _ = build_site.STATUS["VERIFIED_WITH_CORRECTIONS"]
    assert cls != build_site.STATUS["VERIFIED"][0]
    rules = _css_rules()
    assert rules[f".badge.{cls}"] != rules[".badge.verified"]


def test_attempt_cards_are_addressable(site):
    """Without an id there is no way to link a reader to the record in question."""
    missing = []
    for slug in problems():
        index = json.loads(
            (ROOT / "problems" / slug / "prior-art.json").read_text(encoding="utf-8")
        )
        page = site[f"{slug}.html"]
        for attempt in index["attempts"]:
            if f'id="attempt-{attempt["id"]}"' not in page:
                missing.append(f"{slug}/{attempt['id']}")
    assert not missing, f"attempt cards with no anchor: {missing}"


def test_cross_references_between_attempts_are_links(site):
    """Following the citation graph is most of what a reader does here."""
    resolvable, unlinked = [], []
    for slug in problems():
        index = json.loads(
            (ROOT / "problems" / slug / "prior-art.json").read_text(encoding="utf-8")
        )
        page = site[f"{slug}.html"]
        on_page = {a["id"] for a in index["attempts"]}
        for attempt in index["attempts"]:
            for field in ("verifies", "superseded_by", "refutes"):
                value = attempt.get(field)
                if not value:
                    continue
                for ident in re.findall(r"\b\d{3}\b", str(value)):
                    if ident not in on_page:
                        continue
                    ref = f"{slug}/{attempt['id']}.{field} -> {ident}"
                    resolvable.append(ref)
                    if f'href="#attempt-{ident}"' not in page:
                        unlinked.append(ref)
    assert resolvable, "no cross-references in the corpus; this test is vacuous"
    assert not unlinked, f"cross-references left as plain text: {unlinked}"


def test_every_card_carries_its_status_on_its_edge(site):
    """Badges start at whatever x the title ends, so a column cannot be scanned.

    The status class on the card is what tints its left edge; without it the
    card falls back to the neutral rule colour.
    """
    missing = []
    for slug in problems():
        index = json.loads(
            (ROOT / "problems" / slug / "prior-art.json").read_text(encoding="utf-8")
        )
        page = site[f"{slug}.html"]
        for attempt in index["attempts"]:
            cls, _ = build_site.STATUS[attempt["status"]]
            card = re.search(
                rf'<div class="report ([^"]*)" id="attempt-{attempt["id"]}"', page
            )
            if not card or f"s-{cls}" not in card.group(1):
                missing.append(f"{slug}/{attempt['id']} ({cls})")
        route, _ = build_site.STATUS[index["route_status"]]
        if f'<li class="s-{route}">' not in site["index.html"]:
            missing.append(f"index card for {slug} ({route})")
    assert not missing, f"cards with no status edge: {missing}"


def test_mode_is_a_chip_not_a_row(site):
    """One categorical word on every card does not deserve a whole dl row."""
    page = site["union-closed.html"]
    assert '<span class="mode">informed</span>' in page
    assert "<dt>Mode</dt>" not in page, "mode is still rendered as a definition row"


def test_report_labels_sit_beside_their_values_when_there_is_room():
    """Stacked, 'MODE / informed' burned two lines on two words."""
    _, blocks = _split_at_rules(build_site.CSS)
    wide = [b for cond, b in blocks if "min-width" in cond]
    assert wide, "no wide-viewport rules at all"
    rules = {}
    for body in wide:
        rules.update(_declarations(body))
    assert "grid" in rules.get(".report dl", ""), (
        "report definition lists never become a label column on wider screens"
    )


@pytest.mark.parametrize("scheme", ["light", "dark"])
def test_badge_text_meets_wcag_aa(scheme):
    """Badges are .68rem uppercase, so they are small text: 4.5:1, not 3:1.

    EVIDENCE sat at 4.00 and LIVE at 3.67 against their own fills.
    """
    card = "#ffffff" if scheme == "light" else "#1e1d23"
    rules = _css_rules() if scheme == "light" else {**_css_rules(), **_dark_rules()}
    failures = []
    for cls in {cls for cls, _ in build_site.STATUS.values()}:
        body = rules.get(f".badge.{cls}", "")
        colour = re.search(r"(?<!-)\bcolor:\s*(#[0-9a-f]{3,8})", body)
        if not colour:  # neutral badges use palette tokens, checked separately
            continue
        fill = re.search(r"background:\s*(#[0-9a-f]{3,8})", body)
        layers = ([fill.group(1)] if fill else []) + [card]
        got = _contrast(colour.group(1), layers)
        if got < 4.5:
            failures.append(f"{scheme} .badge.{cls}: {got:.2f}:1")
    assert not failures, f"badge text below WCAG AA: {failures}"


def test_nav_links_are_reachable_tap_targets():
    """At .82rem/1 with no padding these were 13px tall on a phone."""
    body = _css_rules().get("header.site nav a", "")
    match = re.search(r"min-height:\s*(\d+)px", body)
    assert match, f"nav links declare no min-height: {body!r}"
    assert int(match.group(1)) >= 44, "nav tap target under 44px"


def test_issue_links_name_real_templates(site):
    """A link to a template that does not exist lands the reporter nowhere."""
    available = {
        p.name for p in (ROOT / ".github" / "ISSUE_TEMPLATE").glob("*.yml")
    } - {"config.yml"}
    referenced = set()
    for body in site.values():
        referenced |= set(re.findall(r"issues/new\?template=([a-z_]+\.yml)", body))
    assert referenced, "no issue templates referenced from the site"
    unknown = referenced - available
    assert not unknown, f"site links to missing issue templates: {sorted(unknown)}"

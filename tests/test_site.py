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

"""Keep the swarm tool honest without touching the network.

swarm.py is the one script that talks to an external, paid API, so its
checkpoint/resume logic is what stands between a crashed sweep and a
double-billed or half-lost one. These tests pin the pure parts — template
expansion, job-state classification, payload shape, response parsing — and
drive the CLI far enough to prove that `plan` and `run --dry-run` never
need a key. The network call itself is exercised by actually running a
swarm, not by a mock that would only encode our guesses.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

_spec = importlib.util.spec_from_file_location("swarm", ROOT / "scripts" / "swarm.py")
swarm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(swarm)


def test_expand_template_substitutes_and_numbers():
    jobs = swarm.expand_template(
        "Lens: {{value}}.\nApply {{value}} to the problem.",
        ["Ergodic Theory", "tropical geometry!"],
    )
    assert [name for name, _ in jobs] == [
        "001-ergodic-theory.md",
        "002-tropical-geometry.md",
    ]
    assert jobs[0][1] == "Lens: Ergodic Theory.\nApply Ergodic Theory to the problem."


def test_expand_template_rejects_missing_placeholder():
    """N identical prompts is always a mistake; fail loudly, not expensively."""
    try:
        swarm.expand_template("no placeholder here", ["a"])
    except swarm.SwarmError:
        pass
    else:
        raise AssertionError("template without {{value}} was accepted")


def test_read_values_skips_comments_and_blanks(tmp_path):
    values = tmp_path / "values.txt"
    values.write_text("# lenses\nfirst\n\n  second  \n", encoding="utf-8")
    assert swarm.read_values(values) == ["first", "second"]


def test_job_state_lifecycle(tmp_path):
    """pending -> done -> stale (prompt edited) -> failed each classify right.

    'done' vs 'stale' is the resume guarantee: a re-run must skip paid work
    but must not skip a job whose prompt changed under it.
    """
    outdir = tmp_path / "out"
    outdir.mkdir()
    prompt = tmp_path / "001-job.md"
    prompt.write_text("v1", encoding="utf-8")
    assert swarm.job_state(prompt, outdir) == "pending"

    meta = outdir / "001-job.meta.json"
    meta.write_text(
        json.dumps({"prompt_sha256": swarm.sha256_text("v1")}), encoding="utf-8"
    )
    assert swarm.job_state(prompt, outdir) == "done"

    prompt.write_text("v2", encoding="utf-8")
    assert swarm.job_state(prompt, outdir) == "stale"

    meta.unlink()
    (outdir / "001-job.error.txt").write_text("HTTP 500", encoding="utf-8")
    assert swarm.job_state(prompt, outdir) == "failed"


def test_build_payload_shape():
    payload = swarm.build_payload("gpt-5.6-luna", "hi", "low", 4096)
    assert payload == {
        "model": "gpt-5.6-luna",
        "input": "hi",
        "reasoning": {"effort": "low"},
        "max_output_tokens": 4096,
    }


def test_extract_text_ignores_non_text_items():
    resp = {
        "output": [
            {"type": "reasoning", "content": None},
            {
                "type": "message",
                "content": [
                    {"type": "output_text", "text": "part one "},
                    {"type": "annotation", "text": "IGNORED"},
                    {"type": "output_text", "text": "part two"},
                ],
            },
        ]
    }
    assert swarm.extract_text(resp) == "part one part two"


def test_estimate_cost_known_and_unknown():
    known = {"gpt-5.6-luna": {"input_tokens": 1_000_000, "output_tokens": 1_000_000}}
    assert swarm.estimate_cost(known) == 0.20 + 1.20
    unknown = {"some-future-model": {"input_tokens": 5, "output_tokens": 5}}
    assert swarm.estimate_cost(unknown) is None


def test_cli_plan_then_dry_run_needs_no_key(tmp_path):
    """The full offline path: plan a sweep, dry-run it, with no API key set.

    A director must be able to stage and inspect a sweep before spending
    anything; if dry-run ever demands a key, that property is gone.
    """
    template = tmp_path / "template.md"
    template.write_text("Apply the {{value}} lens.", encoding="utf-8")
    values = tmp_path / "values.txt"
    values.write_text("alpha\nbeta\n", encoding="utf-8")
    jobdir = tmp_path / "jobs"

    env = {k: v for k, v in os.environ.items() if k != "OPENAI_API_KEY"}

    def run(*argv):
        return subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "swarm.py"), *argv],
            capture_output=True,
            text=True,
            env=env,
        )

    planned = run("plan", str(jobdir), "--template", str(template),
                  "--values", str(values))
    assert planned.returncode == 0, planned.stderr
    assert sorted(p.name for p in jobdir.iterdir()) == [
        "001-alpha.md",
        "002-beta.md",
    ]

    dry = run("run", str(jobdir), "--out", str(tmp_path / "out"), "--dry-run")
    assert dry.returncode == 0, dry.stderr
    assert "2 to run" in dry.stdout
    assert "would run 001-alpha.md" in dry.stdout

    # A real run without a key must refuse before any network attempt.
    wet = run("run", str(jobdir), "--out", str(tmp_path / "out"))
    assert wet.returncode == 2
    assert "OPENAI_API_KEY" in wet.stderr


def test_cli_status_on_partial_outdir(tmp_path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "001-a.meta.json").write_text(
        json.dumps(
            {
                "model": "gpt-5.6-luna",
                "usage": {"input_tokens": 100, "output_tokens": 200},
            }
        ),
        encoding="utf-8",
    )
    (outdir / "002-b.error.txt").write_text("HTTP 500: boom", encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "swarm.py"), "status", str(outdir)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "1 done, 1 failed" in result.stdout
    assert "in=100 out=200" in result.stdout
    assert "HTTP 500" in result.stdout

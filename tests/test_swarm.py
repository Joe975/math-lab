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


def test_provider_for_model_routes_by_name():
    assert swarm.provider_for_model("gpt-5.6-luna") == "openai"
    assert swarm.provider_for_model("gemini-3.1-flash-lite") == "gemini"
    assert swarm.provider_for_model("gemma-4-31b-it") == "gemini"
    # An unknown name must not guess Gemini: a proxy's model name belongs on
    # the OpenAI-compatible path.
    assert swarm.provider_for_model("some-proxy-model") == "openai"


def test_resolve_provider_explicit_beats_inference():
    """`--provider` is the escape hatch for an OpenAI-compatible Gemini proxy."""
    assert swarm.resolve_provider("auto", "gemini-3.7-flash").name == "gemini"
    assert swarm.resolve_provider("openai", "gemini-3.7-flash").name == "openai"


def test_resolve_target_picks_each_providers_own_default():
    """`--provider gemini` with no `--model` must not fall to the OpenAI
    default. That failure is silent, and it turns a cross-family skeptic pass
    into a same-family one -- the exact thing the second provider is for."""
    provider, model = swarm.resolve_target("gemini", None)
    assert (provider.name, model) == ("gemini", swarm.DEFAULT_GEMINI_MODEL)
    provider, model = swarm.resolve_target("openai", None)
    assert (provider.name, model) == ("openai", swarm.DEFAULT_MODEL)
    provider, model = swarm.resolve_target("auto", None)
    assert (provider.name, model) == ("openai", swarm.DEFAULT_MODEL)
    # An explicit model always wins over both defaults.
    provider, model = swarm.resolve_target("auto", "gemini-3.1-flash-lite")
    assert (provider.name, model) == ("gemini", "gemini-3.1-flash-lite")


def test_openai_payload_and_endpoint():
    provider = swarm.PROVIDERS["openai"]
    assert provider.payload("gpt-5.6-luna", "hi", "low", 4096) == {
        "model": "gpt-5.6-luna",
        "input": "hi",
        "reasoning": {"effort": "low"},
        "max_output_tokens": 4096,
    }
    assert provider.endpoint("https://x/v1", "gpt-5.6-luna") == "https://x/v1/responses"
    assert "Authorization" in provider.headers("k")


def test_gemini_payload_and_endpoint():
    provider = swarm.PROVIDERS["gemini"]
    assert provider.payload("gemini-3.7-flash", "hi", "high", 4096) == {
        "contents": [{"role": "user", "parts": [{"text": "hi"}]}],
        "generationConfig": {
            "maxOutputTokens": 4096,
            "thinkingConfig": {"thinkingLevel": "high"},
        },
    }
    # The model name rides in the path, not the body — a provider that got
    # this wrong would silently bill the default model.
    assert provider.endpoint("https://x/v1beta", "gemini-3.7-flash") == (
        "https://x/v1beta/models/gemini-3.7-flash:generateContent"
    )
    assert provider.headers("k")["x-goog-api-key"] == "k"


def test_openai_extract_text_ignores_non_text_items():
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
    assert swarm.PROVIDERS["openai"].text(resp) == "part one part two"


def test_gemini_extract_text_drops_thought_parts():
    """Thought summaries are not the answer; a swarm return must be the answer."""
    resp = {
        "candidates": [
            {
                "finishReason": "STOP",
                "content": {
                    "parts": [
                        {"text": "IGNORED", "thought": True},
                        {"text": "part one "},
                        {"text": "part two", "thoughtSignature": "opaque"},
                    ]
                },
            }
        ]
    }
    assert swarm.PROVIDERS["gemini"].text(resp) == "part one part two"


def test_gemini_usage_bills_thoughts_as_output():
    """Gemini's output price includes thinking tokens; leaving thoughtsTokenCount
    out of the output bucket understates every sweep's cost."""
    resp = {
        "usageMetadata": {
            "promptTokenCount": 30,
            "candidatesTokenCount": 10,
            "thoughtsTokenCount": 57,
            "totalTokenCount": 97,
        }
    }
    assert swarm.PROVIDERS["gemini"].usage(resp) == {
        "input_tokens": 30,
        "output_tokens": 67,
    }


def test_incomplete_detection_per_provider():
    """Truncation must fail the job on both dialects, not land as a short draft."""
    openai, gemini = swarm.PROVIDERS["openai"], swarm.PROVIDERS["gemini"]
    assert openai.incomplete({"status": "completed"}) is None
    assert openai.incomplete({"status": "incomplete"}) is not None

    ok = {"candidates": [{"finishReason": "STOP", "content": {"parts": []}}]}
    assert gemini.incomplete(ok) is None
    truncated = {"candidates": [{"finishReason": "MAX_TOKENS"}]}
    assert "MAX_TOKENS" in gemini.incomplete(truncated)
    # A prompt-level block returns no candidate at all.
    assert gemini.incomplete({"promptFeedback": {"blockReason": "SAFETY"}}) is not None


def test_retry_hint_reads_both_channels():
    """OpenAI throttles via a header, Google via RetryInfo in the error body.

    Missing the body channel is not cosmetic: a free-tier Gemini quota window
    is 30s, longer than the whole 2/4/8/16 ladder, so a throttle that should
    cost one wait instead kills the job.
    """
    assert swarm.retry_hint_seconds({"Retry-After": "7"}, "") == 7
    body = (
        '{"error":{"code":429,"details":[{"@type":'
        '"type.googleapis.com/google.rpc.RetryInfo","retryDelay":"30s"}]}}'
    )
    assert swarm.retry_hint_seconds({}, body) == 30
    # Fractional delays round up; waking early just burns an attempt.
    assert swarm.retry_hint_seconds({}, '"retryDelay": "30.5s"') == 31
    assert swarm.retry_hint_seconds({}, "no hint here") is None


def test_gemini_key_env_precedence(monkeypatch):
    """GEMINI_KEY is what the lab environment injects; GEMINI_API_KEY is
    Google's documented name and must win when both are present."""
    provider = swarm.PROVIDERS["gemini"]
    monkeypatch.delenv("GEMINI_API_KEY", raising=False)
    monkeypatch.setenv("GEMINI_KEY", "from-lab")
    assert provider.api_key() == "from-lab"
    monkeypatch.setenv("GEMINI_API_KEY", "from-google")
    assert provider.api_key() == "from-google"
    monkeypatch.delenv("GEMINI_KEY", raising=False)
    monkeypatch.delenv("GEMINI_API_KEY", raising=False)
    assert provider.api_key() is None


def test_estimate_cost_known_and_unknown():
    known = {"gpt-5.6-luna": {"input_tokens": 1_000_000, "output_tokens": 1_000_000}}
    assert swarm.estimate_cost(known) == 0.20 + 1.20
    mixed = {
        "gpt-5.6-luna": {"input_tokens": 1_000_000, "output_tokens": 0},
        "gemini-3.1-flash-lite": {"input_tokens": 0, "output_tokens": 1_000_000},
    }
    assert swarm.estimate_cost(mixed) == 0.20 + 1.50
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

    keys = {"OPENAI_API_KEY", "GEMINI_API_KEY", "GEMINI_KEY"}
    env = {k: v for k, v in os.environ.items() if k not in keys}

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

    # A real run without a key must refuse before any network attempt, and
    # must name the key the chosen provider actually wants.
    wet = run("run", str(jobdir), "--out", str(tmp_path / "out"))
    assert wet.returncode == 2
    assert "OPENAI_API_KEY" in wet.stderr

    wet_gemini = run("run", str(jobdir), "--out", str(tmp_path / "out"),
                     "--model", "gemini-3.1-flash-lite")
    assert wet_gemini.returncode == 2
    assert "GEMINI_API_KEY" in wet_gemini.stderr
    assert "OPENAI_API_KEY" not in wet_gemini.stderr


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

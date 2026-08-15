#!/usr/bin/env python3
"""Fan prompt sets out to a cheap external model, checkpointed and resumable.

    python scripts/swarm.py plan --template FILE --values FILE JOBDIR
    python scripts/swarm.py run JOBDIR [--out DIR] [--model M] [--effort E]
                                [--workers N] [--max-output-tokens K] [--dry-run]
    python scripts/swarm.py one "PROMPT" [--model M] [--effort E]
    python scripts/swarm.py status OUTDIR

`plan`   — expand a template against a values file (one value per line,
           `#` comments and blanks skipped) into numbered prompt files;
           every occurrence of `{{value}}` in the template is replaced.
`run`    — send every prompt file in JOBDIR (*.md, *.txt) to the model,
           WORKERS at a time. Each job writes `<name>.response.md` plus
           `<name>.meta.json` (model, usage, prompt hash) into the out
           dir. Jobs whose meta already matches the prompt hash are
           skipped, so re-running after a crash resumes; a job whose
           prompt changed re-runs. Failures write `<name>.error.txt`
           and are retried on the next run.
`one`    — a single call, response text to stdout. For smoke tests.
`status` — done/failed counts and token/cost totals for an out dir.

Who should use this and how is in docs/SWARM.md: a director model writes
the briefs and filters the returns; the swarm only drafts. Nothing a swarm
returns is a result — the standard adversarial pipeline applies
(docs/CYCLE.md, step 4).

Two provider dialects are built in, picked from the model name unless
`--provider` says otherwise:

    openai   gpt-*      OPENAI_API_KEY   Responses API
    gemini   gemini-*   GEMINI_API_KEY   generateContent  (GEMINI_KEY also read)
                        or gemma-*

`<PROVIDER>_BASE_URL` overrides an endpoint for proxies or compatible
services. Having two families available is not just redundancy: a skeptic
re-implementation written by a *different* model family is stronger
independence than a second worker of the same family (docs/SWARM.md).
Standard library only.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import hashlib
import json
import math
import os
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

# Default worker on each provider as of 2026-08. On the OpenAI side that is
# the cheapest reasoning tier; on the Gemini side it is the current Flash
# release (gemini-3.7-flash, 2026-08-13), which is on introductory pricing at
# roughly three times the flash-lite rate and is the tier tuned for the
# code-writing briefs the swarm is actually given.
DEFAULT_MODEL = "gpt-5.6-luna"
DEFAULT_GEMINI_MODEL = "gemini-3.7-flash"

# $ per 1M tokens (input, output), standard tier, text prompts under 200k,
# read off the OpenAI and Google pricing pages 2026-08. Gemini output prices
# include thinking tokens, which is why usage folds thoughts into output
# below. Only feeds the printed estimate; drift here never affects results.
PRICES = {
    "gpt-5.6-luna": (0.20, 1.20),
    "gpt-5.6-terra": (2.00, 12.00),
    "gpt-5.6-sol": (5.00, 30.00),
    "gemini-3.1-flash-lite": (0.25, 1.50),
    "gemini-3.5-flash-lite": (0.30, 2.50),
    "gemini-3.5-flash": (1.50, 9.00),
    "gemini-3.6-flash": (0.75, 3.75),
    "gemini-3.7-flash": (0.75, 3.75),
    "gemini-3.1-pro-preview": (2.00, 12.00),
}

PROMPT_SUFFIXES = {".md", ".txt"}
PLACEHOLDER = "{{value}}"
# Same retry ladder the repo uses for git pushes.
BACKOFFS = (2, 4, 8, 16)
# Ceiling on a server-requested wait. Gemini free-tier quota windows ask for
# ~30s; anything much past that is better handled by resuming the sweep.
MAX_BACKOFF = 60


class SwarmError(Exception):
    """A job-level failure worth recording, as opposed to a usage error."""


def retry_hint_seconds(headers, detail: str) -> int | None:
    """How long the server asked us to wait, from whichever channel it used.

    OpenAI sends a `Retry-After` header. Google sends nothing in the headers
    and puts a `google.rpc.RetryInfo` in the error body instead — and on a
    free-tier key its quota window (30s) is longer than this module's whole
    backoff ladder, so ignoring it turns a routine throttle into a dead job.
    """
    header = headers.get("Retry-After") if headers else None
    if header and header.strip().isdigit():
        return int(header.strip())
    match = re.search(r'"retryDelay"\s*:\s*"(\d+(?:\.\d+)?)s"', detail)
    if match:
        # Round up: waking a hair early just burns another attempt.
        return math.ceil(float(match.group(1)))
    return None


# --- pure helpers (unit-tested offline) -----------------------------------


def slugify(text: str, maxlen: int = 40) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")
    return slug[:maxlen].rstrip("-") or "job"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def read_values(path: Path) -> list[str]:
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            values.append(line)
    return values


def expand_template(template: str, values: list[str]) -> list[tuple[str, str]]:
    """(filename, prompt) per value. Numbered so name collisions cannot happen."""
    if PLACEHOLDER not in template:
        raise SwarmError(
            f"template contains no {PLACEHOLDER}; every job would be identical"
        )
    return [
        (f"{i:03d}-{slugify(value)}.md", template.replace(PLACEHOLDER, value))
        for i, value in enumerate(values, start=1)
    ]


def job_state(prompt_path: Path, outdir: Path) -> str:
    """'done' | 'stale' (prompt changed since meta) | 'failed' | 'pending'."""
    meta_path = outdir / f"{prompt_path.stem}.meta.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        current = sha256_text(prompt_path.read_text(encoding="utf-8"))
        return "done" if meta.get("prompt_sha256") == current else "stale"
    if (outdir / f"{prompt_path.stem}.error.txt").exists():
        return "failed"
    return "pending"


# --- providers ------------------------------------------------------------
#
# Each provider translates the same four job knobs (model, prompt, effort,
# token cap) into one API's dialect and reads the answer back out in the
# shape the rest of this file expects: text, {input,output}_tokens usage,
# and a reason-string when the response came back short. Adding a third
# provider means adding a class here and a row in PROVIDERS; nothing else
# in the file knows which API it is talking to.


class Provider:
    name = ""
    default_base_url = ""
    key_envs: tuple[str, ...] = ()

    def endpoint(self, base_url: str, model: str) -> str:
        raise NotImplementedError

    def headers(self, api_key: str) -> dict:
        raise NotImplementedError

    def payload(
        self, model: str, prompt: str, effort: str, max_output_tokens: int
    ) -> dict:
        raise NotImplementedError

    def incomplete(self, resp: dict) -> str | None:
        """A reason string if the response is unusable, else None."""
        raise NotImplementedError

    def text(self, resp: dict) -> str:
        raise NotImplementedError

    def usage(self, resp: dict) -> dict:
        raise NotImplementedError

    def model_of(self, resp: dict, fallback: str) -> str:
        return fallback

    def response_id(self, resp: dict) -> str | None:
        return None

    def base_url(self) -> str:
        return os.environ.get(
            f"{self.name.upper()}_BASE_URL", self.default_base_url
        )

    def api_key(self) -> str | None:
        for env in self.key_envs:
            if os.environ.get(env):
                return os.environ[env]
        return None

    def key_hint(self) -> str:
        return " or ".join(self.key_envs)


class OpenAIProvider(Provider):
    name = "openai"
    default_base_url = "https://api.openai.com/v1"
    key_envs = ("OPENAI_API_KEY",)

    def endpoint(self, base_url, model):
        return f"{base_url}/responses"

    def headers(self, api_key):
        return {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        }

    def payload(self, model, prompt, effort, max_output_tokens):
        return {
            "model": model,
            "input": prompt,
            "reasoning": {"effort": effort},
            "max_output_tokens": max_output_tokens,
        }

    def incomplete(self, resp):
        if resp.get("status") == "completed":
            return None
        return json.dumps(
            {
                "status": resp.get("status"),
                "incomplete_details": resp.get("incomplete_details"),
            },
            indent=1,
        )

    def text(self, resp):
        """Concatenated output_text chunks; reasoning items carry no text here."""
        return "".join(
            chunk.get("text", "")
            for item in resp.get("output", [])
            for chunk in item.get("content") or []
            if chunk.get("type") == "output_text"
        )

    def usage(self, resp):
        usage = resp.get("usage", {})
        return {
            "input_tokens": usage.get("input_tokens", 0),
            "output_tokens": usage.get("output_tokens", 0),
        }

    def model_of(self, resp, fallback):
        return resp.get("model", fallback)

    def response_id(self, resp):
        return resp.get("id")


class GeminiProvider(Provider):
    """Google's generateContent. Effort maps straight onto thinkingLevel.

    The level names line up with the OpenAI effort names, so no translation
    table can drift; a level a given model does not support comes back as a
    plain 400 naming it (e.g. `minimal` on gemini-3.7-flash).
    """

    name = "gemini"
    default_base_url = "https://generativelanguage.googleapis.com/v1beta"
    # GEMINI_KEY is what the lab's environment injects; GEMINI_API_KEY is
    # Google's documented name and wins if both are set.
    key_envs = ("GEMINI_API_KEY", "GEMINI_KEY")

    def endpoint(self, base_url, model):
        return f"{base_url}/models/{model}:generateContent"

    def headers(self, api_key):
        return {"x-goog-api-key": api_key, "Content-Type": "application/json"}

    def payload(self, model, prompt, effort, max_output_tokens):
        return {
            "contents": [{"role": "user", "parts": [{"text": prompt}]}],
            "generationConfig": {
                "maxOutputTokens": max_output_tokens,
                "thinkingConfig": {"thinkingLevel": effort},
            },
        }

    def _candidate(self, resp: dict) -> dict:
        candidates = resp.get("candidates") or []
        return candidates[0] if candidates else {}

    def incomplete(self, resp):
        candidate = self._candidate(resp)
        if not candidate:
            # No candidate at all means a prompt-level block; the feedback
            # block is the only place that says why.
            return json.dumps(
                {"promptFeedback": resp.get("promptFeedback")}, indent=1
            )
        reason = candidate.get("finishReason")
        if reason in (None, "STOP"):
            return None
        # MAX_TOKENS is the common one: a truncated draft is worse than a
        # retry, so it fails like OpenAI's `incomplete`.
        return json.dumps({"finishReason": reason}, indent=1)

    def text(self, resp):
        """Answer parts only — thought-summary parts are dropped."""
        return "".join(
            part.get("text", "")
            for part in self._candidate(resp).get("content", {}).get("parts", [])
            if not part.get("thought")
        )

    def usage(self, resp):
        usage = resp.get("usageMetadata", {})
        return {
            "input_tokens": usage.get("promptTokenCount", 0),
            # Thinking tokens bill at the output rate, so they belong in the
            # output bucket or every cost estimate is silently low.
            "output_tokens": (
                usage.get("candidatesTokenCount", 0)
                + usage.get("thoughtsTokenCount", 0)
            ),
        }

    def model_of(self, resp, fallback):
        return resp.get("modelVersion", fallback)

    def response_id(self, resp):
        return resp.get("responseId")


PROVIDERS = {p.name: p for p in (OpenAIProvider(), GeminiProvider())}


def provider_for_model(model: str) -> str:
    """Which dialect a model name speaks. Unknown prefixes fall to openai,
    which is also where an OpenAI-compatible proxy model name belongs."""
    if model.startswith(("gemini-", "gemma-")):
        return "gemini"
    return "openai"


def resolve_provider(name: str, model: str) -> Provider:
    return PROVIDERS[provider_for_model(model) if name == "auto" else name]


def resolve_target(provider_name: str, model: str | None) -> tuple[Provider, str]:
    """(provider, model) from whichever of the two the caller actually gave.

    `--provider gemini` with no `--model` has to land on the Gemini default,
    not on the OpenAI one. Getting this wrong is silent and expensive in the
    way that matters most here: a cross-family skeptic pass that quietly runs
    on the same family it was meant to be independent of.
    """
    if model is None:
        model = DEFAULT_GEMINI_MODEL if provider_name == "gemini" else DEFAULT_MODEL
    return resolve_provider(provider_name, model), model


def estimate_cost(usage_by_model: dict[str, dict[str, int]]) -> float | None:
    """Dollar estimate across models, or None if any model has no price row."""
    total = 0.0
    for model, usage in usage_by_model.items():
        if model not in PRICES:
            return None
        rate_in, rate_out = PRICES[model]
        total += usage.get("input_tokens", 0) * rate_in / 1e6
        total += usage.get("output_tokens", 0) * rate_out / 1e6
    return total


# --- the network part -----------------------------------------------------


def call_api(
    provider: Provider,
    model: str,
    payload: dict,
    api_key: str,
    base_url: str,
    timeout: float,
) -> dict:
    """POST one job with the repo-standard retry ladder.

    Retries transport errors, 429 and 5xx; any other 4xx fails immediately
    because a malformed request will not heal by waiting. 5xx matters more
    on the Gemini side than it looks: `UNAVAILABLE / high demand` is a
    routine 503 there, and without the ladder a sweep would shed jobs to it.
    """
    url = provider.endpoint(base_url, model)
    headers = provider.headers(api_key)
    body = json.dumps(payload).encode("utf-8")
    last = "unknown error"
    for attempt, backoff in enumerate((*BACKOFFS, None)):
        req = urllib.request.Request(url, data=body, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as e:
            detail = e.read().decode("utf-8", errors="replace")[:2000]
            last = f"HTTP {e.code}: {detail}"
            if e.code != 429 and e.code < 500:
                raise SwarmError(last) from None
            hint = retry_hint_seconds(e.headers, detail)
            if hint is not None and backoff is not None:
                # Honour the server over our ladder, but never sleep past the
                # cap: a job stuck behind a multi-minute quota should fail and
                # be resumed later, not hold a worker slot.
                backoff = min(max(backoff, hint), MAX_BACKOFF)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            last = f"{type(e).__name__}: {e}"
        if backoff is None:
            break
        time.sleep(backoff)
    raise SwarmError(f"gave up after {len(BACKOFFS) + 1} attempts; last: {last}")


def run_job(
    prompt_path: Path, outdir: Path, args, provider: Provider, api_key: str
) -> dict:
    name = prompt_path.stem
    prompt = prompt_path.read_text(encoding="utf-8")
    payload = provider.payload(
        args.model, prompt, args.effort, args.max_output_tokens
    )
    started = time.monotonic()
    try:
        resp = call_api(
            provider, args.model, payload, api_key, args.base_url, args.timeout
        )
        reason = provider.incomplete(resp)
        if reason is not None:
            # Keep whatever text came back as evidence, but an incomplete
            # response is a failure: a truncated draft is worse than a retry.
            raise SwarmError(
                f"{reason}\npartial_text: {provider.text(resp)!r}"
            )
    except SwarmError as e:
        (outdir / f"{name}.error.txt").write_text(str(e), encoding="utf-8")
        return {"name": name, "ok": False, "error": str(e).splitlines()[0][:120]}

    (outdir / f"{name}.response.md").write_text(
        provider.text(resp), encoding="utf-8"
    )
    meta = {
        "name": name,
        "provider": provider.name,
        "model": provider.model_of(resp, args.model),
        "effort": args.effort,
        "max_output_tokens": args.max_output_tokens,
        "prompt_file": str(prompt_path),
        "prompt_sha256": sha256_text(prompt),
        "response_id": provider.response_id(resp),
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds"
        ),
        "elapsed_s": round(time.monotonic() - started, 1),
        "usage": provider.usage(resp),
    }
    # Meta lands last: its presence is the done marker, so a crash between
    # the two writes just means the job re-runs.
    (outdir / f"{name}.meta.json").write_text(
        json.dumps(meta, indent=1), encoding="utf-8"
    )
    (outdir / f"{name}.error.txt").unlink(missing_ok=True)
    return {"name": name, "ok": True, "meta": meta}


# --- commands -------------------------------------------------------------


def default_outdir(jobdir: Path) -> Path:
    return Path(os.environ.get("MATHLAB_OUT", "out")) / "swarm" / jobdir.name


def summarize(metas: list[dict]) -> str:
    usage_by_model: dict[str, dict[str, int]] = {}
    for meta in metas:
        bucket = usage_by_model.setdefault(
            meta.get("model", "?"), {"input_tokens": 0, "output_tokens": 0}
        )
        for key in bucket:
            bucket[key] += meta.get("usage", {}).get(key, 0)
    parts = [
        f"{model}: in={u['input_tokens']} out={u['output_tokens']}"
        for model, u in sorted(usage_by_model.items())
    ]
    cost = estimate_cost(usage_by_model)
    if cost is not None:
        parts.append(f"est ${cost:.4f}")
    return "; ".join(parts) if parts else "no usage recorded"


def cmd_plan(args) -> int:
    jobs = expand_template(
        Path(args.template).read_text(encoding="utf-8"),
        read_values(Path(args.values)),
    )
    jobdir = Path(args.jobdir)
    jobdir.mkdir(parents=True, exist_ok=True)
    for filename, prompt in jobs:
        (jobdir / filename).write_text(prompt, encoding="utf-8")
    print(f"planned {len(jobs)} jobs in {jobdir}")
    return 0


def cmd_run(args) -> int:
    jobdir = Path(args.jobdir)
    prompts = sorted(
        p for p in jobdir.iterdir() if p.suffix.lower() in PROMPT_SUFFIXES
    )
    if not prompts:
        print(f"error: no {'/'.join(sorted(PROMPT_SUFFIXES))} files in {jobdir}",
              file=sys.stderr)
        return 2
    outdir = Path(args.out) if args.out else default_outdir(jobdir)
    outdir.mkdir(parents=True, exist_ok=True)

    states = {p: job_state(p, outdir) for p in prompts}
    todo = [p for p, s in states.items() if s != "done"]
    done = len(prompts) - len(todo)
    print(f"{len(prompts)} jobs: {done} done, {len(todo)} to run -> {outdir}")
    if args.dry_run or not todo:
        for p in todo:
            print(f"  would run {p.name} ({states[p]})")
        return 0

    provider, args.model = resolve_target(args.provider, args.model)
    api_key = provider.api_key()
    if not api_key:
        print(f"error: {provider.key_hint()} is not set", file=sys.stderr)
        return 2
    args.base_url = args.base_url or provider.base_url()
    print(f"provider {provider.name}, model {args.model}, effort {args.effort}")

    metas, failures = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(run_job, p, outdir, args, provider, api_key): p
            for p in todo
        }
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result["ok"]:
                metas.append(result["meta"])
                usage = result["meta"]["usage"]
                print(
                    f"  done {result['name']} "
                    f"(in={usage.get('input_tokens', '?')} "
                    f"out={usage.get('output_tokens', '?')} "
                    f"{result['meta']['elapsed_s']}s)"
                )
            else:
                failures.append(result)
                print(f"  FAIL {result['name']}: {result['error']}")

    print(f"run: {len(metas)} ok, {len(failures)} failed; {summarize(metas)}")
    return 1 if failures else 0


def cmd_one(args) -> int:
    prompt = sys.stdin.read() if args.prompt == "-" else args.prompt
    provider, args.model = resolve_target(args.provider, args.model)
    api_key = provider.api_key()
    if not api_key:
        print(f"error: {provider.key_hint()} is not set", file=sys.stderr)
        return 2
    base_url = args.base_url or provider.base_url()
    payload = provider.payload(
        args.model, prompt, args.effort, args.max_output_tokens
    )
    try:
        resp = call_api(
            provider, args.model, payload, api_key, base_url, args.timeout
        )
    except SwarmError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    reason = provider.incomplete(resp)
    print(provider.text(resp))
    usage = provider.usage(resp)
    print(
        f"[{provider.model_of(resp, args.model)} "
        f"in={usage['input_tokens']} out={usage['output_tokens']}"
        + (f" INCOMPLETE {reason}" if reason else "")
        + "]",
        file=sys.stderr,
    )
    return 1 if reason else 0


def cmd_status(args) -> int:
    outdir = Path(args.outdir)
    metas = [
        json.loads(p.read_text(encoding="utf-8"))
        for p in sorted(outdir.glob("*.meta.json"))
    ]
    errors = sorted(outdir.glob("*.error.txt"))
    print(f"{outdir}: {len(metas)} done, {len(errors)} failed; {summarize(metas)}")
    for p in errors:
        first = p.read_text(encoding="utf-8").splitlines()
        print(f"  FAIL {p.stem.removesuffix('.error')}: {first[0][:100] if first else ''}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)

    def add_model_args(p):
        p.add_argument(
            "--model",
            default=None,
            help=f"default: {DEFAULT_MODEL}, or {DEFAULT_GEMINI_MODEL} when "
                 "--provider gemini is given without a model",
        )
        p.add_argument(
            "--provider",
            default="auto",
            choices=["auto", *sorted(PROVIDERS)],
            help="API dialect (default: auto — inferred from the model name)",
        )
        p.add_argument(
            "--effort",
            default="low",
            choices=["minimal", "low", "medium", "high"],
            help="reasoning effort / thinkingLevel (default: low — this is a "
                 "breadth tool). Not every model supports every level.",
        )
        p.add_argument("--max-output-tokens", type=int, default=8192)
        p.add_argument("--timeout", type=float, default=600)
        p.add_argument(
            "--base-url",
            default=None,
            help="default: $OPENAI_BASE_URL / $GEMINI_BASE_URL, else the "
                 "provider's public endpoint",
        )

    p = sub.add_parser("plan", help="expand template x values into prompt files")
    p.add_argument("jobdir")
    p.add_argument("--template", required=True)
    p.add_argument("--values", required=True)
    p.set_defaults(fn=cmd_plan)

    p = sub.add_parser("run", help="send every prompt in JOBDIR, resumable")
    p.add_argument("jobdir")
    p.add_argument("--out", help="default: $MATHLAB_OUT/swarm/<jobdir name>")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--dry-run", action="store_true")
    add_model_args(p)
    p.set_defaults(fn=cmd_run)

    p = sub.add_parser("one", help="single call, text to stdout ('-' = stdin)")
    p.add_argument("prompt")
    add_model_args(p)
    p.set_defaults(fn=cmd_one)

    p = sub.add_parser("status", help="summarize an out dir")
    p.add_argument("outdir")
    p.set_defaults(fn=cmd_status)

    args = parser.parse_args(argv)
    try:
        return args.fn(args)
    except SwarmError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())

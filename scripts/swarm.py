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

Needs OPENAI_API_KEY. Talks to the OpenAI Responses API; OPENAI_BASE_URL
overrides the endpoint for proxies or compatible providers. Standard
library only.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

# Cheapest reasoning tier on the price ladder as of 2026-08.
DEFAULT_MODEL = "gpt-5.6-luna"
DEFAULT_BASE_URL = "https://api.openai.com/v1"

# $ per 1M tokens (input, output), standard tier, read off the OpenAI
# pricing page 2026-08. Only feeds the printed estimate; drift here never
# affects results.
PRICES = {
    "gpt-5.6-luna": (0.20, 1.20),
    "gpt-5.6-terra": (2.00, 12.00),
    "gpt-5.6-sol": (5.00, 30.00),
}

PROMPT_SUFFIXES = {".md", ".txt"}
PLACEHOLDER = "{{value}}"
# Same retry ladder the repo uses for git pushes.
BACKOFFS = (2, 4, 8, 16)


class SwarmError(Exception):
    """A job-level failure worth recording, as opposed to a usage error."""


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


def build_payload(
    model: str, prompt: str, effort: str, max_output_tokens: int
) -> dict:
    return {
        "model": model,
        "input": prompt,
        "reasoning": {"effort": effort},
        "max_output_tokens": max_output_tokens,
    }


def extract_text(resp: dict) -> str:
    """Concatenated output_text chunks; reasoning items carry no text here."""
    return "".join(
        chunk.get("text", "")
        for item in resp.get("output", [])
        for chunk in item.get("content") or []
        if chunk.get("type") == "output_text"
    )


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


def call_api(payload: dict, api_key: str, base_url: str, timeout: float) -> dict:
    """POST to /responses with the repo-standard retry ladder.

    Retries transport errors, 429 and 5xx; any other 4xx fails immediately
    because a malformed request will not heal by waiting.
    """
    body = json.dumps(payload).encode("utf-8")
    last = "unknown error"
    for attempt, backoff in enumerate((*BACKOFFS, None)):
        req = urllib.request.Request(
            f"{base_url}/responses",
            data=body,
            headers={
                "Authorization": f"Bearer {api_key}",
                "Content-Type": "application/json",
            },
        )
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as e:
            detail = e.read().decode("utf-8", errors="replace")[:2000]
            last = f"HTTP {e.code}: {detail}"
            if e.code != 429 and e.code < 500:
                raise SwarmError(last) from None
            retry_after = e.headers.get("Retry-After")
            if retry_after and retry_after.isdigit() and backoff is not None:
                backoff = max(backoff, int(retry_after))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            last = f"{type(e).__name__}: {e}"
        if backoff is None:
            break
        time.sleep(backoff)
    raise SwarmError(f"gave up after {len(BACKOFFS) + 1} attempts; last: {last}")


def run_job(prompt_path: Path, outdir: Path, args, api_key: str) -> dict:
    name = prompt_path.stem
    prompt = prompt_path.read_text(encoding="utf-8")
    payload = build_payload(args.model, prompt, args.effort, args.max_output_tokens)
    started = time.monotonic()
    try:
        resp = call_api(payload, api_key, args.base_url, args.timeout)
        if resp.get("status") != "completed":
            # Keep whatever text came back as evidence, but an incomplete
            # response is a failure: a truncated draft is worse than a retry.
            raise SwarmError(
                json.dumps(
                    {
                        "status": resp.get("status"),
                        "incomplete_details": resp.get("incomplete_details"),
                        "partial_text": extract_text(resp),
                    },
                    indent=1,
                )
            )
    except SwarmError as e:
        (outdir / f"{name}.error.txt").write_text(str(e), encoding="utf-8")
        return {"name": name, "ok": False, "error": str(e).splitlines()[0][:120]}

    (outdir / f"{name}.response.md").write_text(extract_text(resp), encoding="utf-8")
    meta = {
        "name": name,
        "model": resp.get("model", args.model),
        "effort": args.effort,
        "max_output_tokens": args.max_output_tokens,
        "prompt_file": str(prompt_path),
        "prompt_sha256": sha256_text(prompt),
        "response_id": resp.get("id"),
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds"
        ),
        "elapsed_s": round(time.monotonic() - started, 1),
        "usage": resp.get("usage", {}),
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

    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        print("error: OPENAI_API_KEY is not set", file=sys.stderr)
        return 2

    metas, failures = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_job, p, outdir, args, api_key): p for p in todo}
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
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        print("error: OPENAI_API_KEY is not set", file=sys.stderr)
        return 2
    payload = build_payload(args.model, prompt, args.effort, args.max_output_tokens)
    try:
        resp = call_api(payload, api_key, args.base_url, args.timeout)
    except SwarmError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    print(extract_text(resp))
    usage = resp.get("usage", {})
    print(
        f"[{resp.get('model')} in={usage.get('input_tokens')} "
        f"out={usage.get('output_tokens')}]",
        file=sys.stderr,
    )
    return 0


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
        p.add_argument("--model", default=DEFAULT_MODEL)
        p.add_argument(
            "--effort",
            default="low",
            choices=["minimal", "low", "medium", "high"],
            help="reasoning effort (default: low — this is a breadth tool)",
        )
        p.add_argument("--max-output-tokens", type=int, default=8192)
        p.add_argument("--timeout", type=float, default=600)
        p.add_argument(
            "--base-url",
            default=os.environ.get("OPENAI_BASE_URL", DEFAULT_BASE_URL),
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

"""Helpers for working with structured pytest reports."""

from __future__ import annotations

from collections.abc import Iterable
from datetime import UTC, datetime
from typing import Any

DEFAULT_REPO_SLUG = "SatoryKono/ChEMBL_data_acquisition"
REQUIRED_SUMMARY_KEYS = {
    "total",
    "passed",
    "failed",
    "skipped",
    "xfailed",
    "xpassed",
    "error",
    "success_rate",
}


def normalise_message(raw: Any) -> str:
    """Convert structured report message payloads to clean strings."""

    if raw is None:
        return ""
    if isinstance(raw, str):
        return raw.strip()
    if isinstance(raw, Iterable) and not isinstance(raw, (bytes, bytearray)):
        joined = "\n".join(str(part) for part in raw)
        return joined.strip()
    return str(raw).strip()


def validate_summary_report(report: dict[str, Any]) -> None:
    """Validate the subset of a structured report required for Markdown rendering."""

    if not isinstance(report, dict):
        raise ValueError("Report payload must be a dictionary")

    summary = report.get("summary")
    if not isinstance(summary, dict):
        raise ValueError("Report summary is missing or not an object")

    missing = REQUIRED_SUMMARY_KEYS - set(summary)
    if missing:
        missing_keys = ", ".join(sorted(missing))
        raise ValueError(f"Summary missing required keys: {missing_keys}")

    for key in REQUIRED_SUMMARY_KEYS - {"success_rate"}:
        value = summary.get(key)
        if not isinstance(value, int) or value < 0:
            raise ValueError(
                f"Summary field '{key}' must be a non-negative integer, got {type(value).__name__}"
            )

    success_rate = summary.get("success_rate")
    if not isinstance(success_rate, (int, float)):
        raise ValueError("Summary field 'success_rate' must be numeric")

    tests = report.get("tests")
    if not isinstance(tests, list):
        raise ValueError("Report tests section must be a list")
    for index, entry in enumerate(tests):
        if not isinstance(entry, dict):
            raise ValueError(
                f"Test entry at index {index} must be a dictionary, got {type(entry).__name__}"
            )

    meta = report.get("meta")
    if meta is not None and not isinstance(meta, dict):
        raise ValueError("Report meta section must be an object when present")


def build_summary_markdown(report: dict[str, Any]) -> str:
    """Render a Markdown summary document for a structured pytest JSON report."""

    validate_summary_report(report)

    meta = report.get("meta", {})
    summary = report.get("summary", {})
    tests = report.get("tests", [])

    repo = meta.get("repo", DEFAULT_REPO_SLUG)
    commit = meta.get("commit", "unknown")
    branch = meta.get("branch", "unknown")
    timestamp = meta.get("ts_utc", datetime.now(UTC).isoformat())
    duration = float(meta.get("duration_sec", 0.0) or 0.0)
    success_rate = float(summary.get("success_rate", 0.0) or 0.0)
    success_rate_pct = success_rate * 100.0

    lines = [
        "# Test Summary",
        "",
        f"- Repo: `{repo}`",
        f"- Commit: {commit}",
        f"- Branch: {branch}",
        f"- Timestamp (UTC): {timestamp}",
        f"- Duration: {duration:.2f} s",
        f"- Success rate: {success_rate_pct:.2f}%",
        "",
        "| total | passed | failed | skipped | xfailed | xpassed | error |",
        "|------:|-------:|-------:|--------:|--------:|--------:|------:|",
        "| {total:5d} | {passed:5d} | {failed:5d} | {skipped:6d} | {xfailed:6d} | {xpassed:6d} | {error:5d} |".format(
            total=summary.get("total", 0),
            passed=summary.get("passed", 0),
            failed=summary.get("failed", 0),
            skipped=summary.get("skipped", 0),
            xfailed=summary.get("xfailed", 0),
            xpassed=summary.get("xpassed", 0),
            error=summary.get("error", 0),
        ),
        "",
        "## Failed / Error details",
    ]

    failure_rows = []
    for test in tests:
        status = str(test.get("status", "")).lower()
        if status not in {"failed", "error"}:
            continue
        nodeid = str(test.get("nodeid", "<unknown>"))
        message = normalise_message(test.get("error"))
        failure_rows.append((nodeid, status, message))

    if not failure_rows:
        lines.append("- None")
    else:
        for nodeid, status, message in failure_rows:
            lines.append(f"- `{nodeid}` ({status})")
            display_message = message or "<no message>"
            lines.append("  ```")
            lines.extend(f"  {line}" for line in display_message.splitlines())
            lines.append("  ```")

    lines.append("")
    return "\n".join(lines)

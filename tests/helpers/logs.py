from __future__ import annotations

import ast
import re
import shlex
from collections.abc import Iterable
from pathlib import Path
from typing import Any

_NEW_PATTERN = re.compile(
    r"^\[(?P<timestamp>[^\]]+)\]\s\[(?P<level>[^\]]+)\]\s\[(?P<name>[^\]]+)\]\s(?P<message>.*)$"
)
_LEGACY_PATTERN = re.compile(
    r"^(?P<timestamp>.+?)\s\[(?P<level>[^\]]+)\]\s(?P<name>[^\s]+)\s::\s(?P<message>.*)$"
)


def _extract_parts(line: str) -> dict[str, str] | None:
    match = _NEW_PATTERN.match(line)
    if match:
        return match.groupdict()
    legacy = _LEGACY_PATTERN.match(line)
    if legacy:
        return legacy.groupdict()
    return None


def parse_log_lines(text: str) -> list[dict[str, Any]]:
    """Parse log lines in modern or legacy structured formats."""

    entries: list[dict[str, Any]] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        components = _extract_parts(line)
        if not components:
            continue
        timestamp = components["timestamp"].strip()
        level = components["level"].strip()
        name = components["name"].strip()
        message = components["message"].strip()
        tokens = shlex.split(message)
        if not tokens:
            continue
        event = tokens[0]
        data: dict[str, Any] = {}
        for token in tokens[1:]:
            if "=" not in token:
                continue
            key, value_str = token.split("=", 1)
            try:
                value = ast.literal_eval(value_str)
            except Exception:  # pragma: no cover - defensive parsing
                value = value_str
            data[key] = value
        entries.append(
            {
                "timestamp": timestamp,
                "level": level,
                "name": name,
                "event": event,
                "data": data,
                "raw": line,
            }
        )
    return entries


def parse_log_file(path: Path) -> list[dict[str, Any]]:
    """Read ``path`` and parse all log entries."""

    return parse_log_lines(path.read_text(encoding="utf-8"))


def iter_events(entries: Iterable[dict[str, Any]]) -> Iterable[str]:
    for entry in entries:
        yield entry.get("event", "")


__all__ = ["parse_log_lines", "parse_log_file", "iter_events"]

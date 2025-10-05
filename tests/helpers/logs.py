from __future__ import annotations

import ast
import shlex
from pathlib import Path
from typing import Any, Iterable


def parse_log_lines(text: str) -> list[dict[str, Any]]:
    """Parse log lines formatted as ``[ts] [LEVEL] [name] message``."""

    entries: list[dict[str, Any]] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split("] ", 3)
        if len(parts) < 4:
            continue
        timestamp = parts[0].lstrip("[")
        level = parts[1].strip("[]")
        name = parts[2].strip("[]")
        message = parts[3]
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

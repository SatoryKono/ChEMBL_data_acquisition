"""Utilities for emitting per-record JSON logs."""
from __future__ import annotations

import json
from typing import Any, Dict, IO


def emit_log(record: Dict[str, Any], fp: IO[str]) -> None:
    """Write a JSON record followed by a newline to ``fp``.

    Parameters
    ----------
    record:
        Mapping to serialise as JSON.
    fp:
        File-like object opened for text writing.
    """

    fp.write(json.dumps(record, ensure_ascii=False) + "\n")


__all__ = ["emit_log"]

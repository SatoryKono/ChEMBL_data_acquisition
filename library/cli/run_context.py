"""Helpers for deriving deterministic runtime context values."""

from __future__ import annotations

import hashlib
from collections.abc import Sequence
from datetime import UTC, datetime, timedelta

UTC = UTC
_SEED_EPOCH = datetime(2000, 1, 1, tzinfo=UTC)
_MAX_SPAN_SECONDS = 50 * 365 * 24 * 60 * 60


def _parse_date_token(value: str) -> datetime | None:
    """Return a UTC datetime parsed from ``value`` when possible."""

    stripped = value.strip()
    if not stripped:
        return None

    normalised = stripped.replace("Z", "+00:00")
    try:
        parsed = datetime.fromisoformat(normalised)
    except ValueError:
        parsed = None
    if parsed is not None:
        if parsed.tzinfo is None:
            return parsed.replace(tzinfo=UTC)
        return parsed.astimezone(UTC)

    for fmt in ("%Y%m%d", "%Y-%m-%d"):
        try:
            parsed = datetime.strptime(stripped, fmt)
        except ValueError:
            continue
        return parsed.replace(tzinfo=UTC)
    return None


def compute_generated_at(
    *,
    date_token: str | None,
    run_id: str,
    seed_parts: Sequence[str] | None = None,
) -> str:
    """Return a deterministic ``generated_at`` timestamp."""

    if date_token:
        parsed = _parse_date_token(date_token)
        if parsed is not None:
            return parsed.isoformat()

    components: list[str] = []
    if run_id:
        components.append(run_id)
    if seed_parts:
        components.extend(str(part) for part in seed_parts if str(part))
    if not components:
        components.append("chembl-data-acquisition")

    seed = "|".join(components)
    digest = hashlib.sha256(seed.encode("utf-8")).digest()
    seconds = int.from_bytes(digest[:6], "big") % _MAX_SPAN_SECONDS
    generated = _SEED_EPOCH + timedelta(seconds=seconds)
    return generated.isoformat()


__all__ = ["compute_generated_at"]

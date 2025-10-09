"""Runtime context describing the current CLI invocation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(slots=True, frozen=True)
class RunContext:
    """Container storing identifiers shared across the running process."""

    run_id: str
    generated_at: str


_CURRENT: Optional[RunContext] = None


def set_current(context: Optional[RunContext]) -> None:
    """Record ``context`` as the active :class:`RunContext`.

    Parameters
    ----------
    context:
        Context describing the active invocation. ``None`` clears any previous
        value and should be used sparingly (for example in tests).
    """

    global _CURRENT
    _CURRENT = context


def get_current() -> Optional[RunContext]:
    """Return the active :class:`RunContext`, if any."""

    return _CURRENT

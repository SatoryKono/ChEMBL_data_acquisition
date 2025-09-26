"""Proxy module for :mod:`library.utils.cli_tools.check_determinism`."""

from __future__ import annotations

from collections.abc import Sequence
from inspect import signature

from library.utils.cli_tools import check_determinism as _impl

_IMPL_MAIN = _impl.main
_HAS_ARGS = len(signature(_IMPL_MAIN).parameters) > 0


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate execution to :func:`library.utils.cli_tools.check_determinism.main`."""

    if _HAS_ARGS:
        return _IMPL_MAIN(argv)
    return _IMPL_MAIN()


if __name__ == "__main__":  # pragma: no cover - convenience entry point
    raise SystemExit(main())

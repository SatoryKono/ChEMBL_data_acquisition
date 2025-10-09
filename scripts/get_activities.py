"""CLI shim delegating to :mod:`library.utils.cli_tools.get_activities`."""

from __future__ import annotations

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

from library.utils.cli_tools import get_activities


def main() -> int:
    """Execute the `get-activities` helper."""

    return get_activities.main()


if __name__ == "__main__":
    raise SystemExit(main())

"""Thin wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

import sys

from library.cli.entrypoints import activity as _activity

sys.modules[__name__] = _activity
sys.modules.setdefault("scripts.get_activity_data", _activity)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_activity.main())

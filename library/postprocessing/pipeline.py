"""Helpers for declarative post-processing pipeline execution."""

from __future__ import annotations

from typing import Any

import pandas as pd


def passthrough(frame: pd.DataFrame, **_: Any) -> pd.DataFrame:
    """Return *frame* unchanged.

    The helper acts as a placeholder callable for configurable pipeline steps
    that may be overridden by table-specific configurations.  Accepting
    ``**kwargs`` keeps the callable signature compatible with future
    extensions that may inject additional keyword parameters.
    """

    return frame

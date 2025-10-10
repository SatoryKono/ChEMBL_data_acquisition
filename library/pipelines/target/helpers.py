"""High level helpers for the target post-processing pipeline."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from ...common.csv_utils import write_csv_deterministic
from ...config import Config
from ...schemas.targets import TARGETS_COLUMN_ORDER
from . import postprocessing as _postprocessing

__all__ = [
    "read_snapshot",
    "postprocess_target_table",
    "postprocess_target_table_file",
]


def _resolve_sep(cfg: Config | None, sep: str | None) -> str:
    if sep is not None:
        return sep
    if cfg is not None:
        return cfg.io.csv_sep
    return ","


def _resolve_encoding(cfg: Config | None, encoding: str | None) -> str:
    if encoding is not None:
        return encoding
    if cfg is not None:
        return cfg.io.csv_encoding
    return "utf-8"


def read_snapshot(
    path: Path | str,
    *,
    cfg: Config | None = None,
    sep: str | None = None,
    encoding: str | None = None,
    columns: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Load a targets snapshot using deterministic defaults.

    The helper mirrors the behaviour of the legacy Power Query workbook: CSV
    files are read using UTF-8 (falling back to the configured encoding) and all
    values are treated as strings to preserve leading zeros and formatting.
    """

    resolved_sep = _resolve_sep(cfg, sep)
    resolved_encoding = _resolve_encoding(cfg, encoding)
    frame = pd.read_csv(
        path,
        dtype=str,
        sep=resolved_sep,
        encoding=resolved_encoding,
        keep_default_na=False,
    )
    if columns is not None:
        missing = set(columns) - set(frame.columns)
        if missing:
            raise ValueError("missing required columns: " + ", ".join(sorted(missing)))
        frame = frame.loc[:, list(columns)].copy()
    return frame


def postprocess_target_table(df: pd.DataFrame) -> pd.DataFrame:
    """Return the Power Query compatible target table.

    This function is a thin wrapper around
    :func:`library.pipelines.target.postprocessing.finalise_targets`. It is
    exposed as a convenience layer for tests and integration scenarios that need
    to operate on in-memory frames without touching the filesystem.
    """

    result = _postprocessing.finalise_targets(df)
    # ``finalise_targets`` already enforces the canonical ordering via
    # :func:`align_target_columns`. Calling ``loc`` defensively ensures that
    # callers receive a projection that matches ``TARGETS_COLUMN_ORDER`` even
    # when upstream helpers change in the future.
    return result.loc[:, TARGETS_COLUMN_ORDER]


def postprocess_target_table_file(
    input_path: Path | str,
    output_path: Path | str,
    *,
    cfg: Config,
    sep: str | None = None,
    encoding: str | None = None,
) -> Path:
    """Run :func:`postprocess_target_table` on ``input_path`` and write ``output_path``."""

    frame = read_snapshot(input_path, cfg=cfg, sep=sep, encoding=encoding)
    processed = postprocess_target_table(frame)
    resolved_sep = _resolve_sep(cfg, sep)
    resolved_encoding = _resolve_encoding(cfg, encoding)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_csv_deterministic(
        processed,
        output,
        col_order=TARGETS_COLUMN_ORDER,
        key_cols=["target_chembl_id"],
        sep=resolved_sep,
        encoding=resolved_encoding,
        cfg=None,
    )
    return output

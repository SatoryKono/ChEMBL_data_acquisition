"""Path helper functions for pipeline outputs."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from ..config import IoCfg


def default_output_path(
    input_path: str | Path,
    cfg: IoCfg,
    *,
    date: str | None = None,
    stamp_mode: str | None = None,
) -> Path:
    """Construct a deterministic output path based on the input file name.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Location of the CSV supplied to the pipeline. Only the stem is used
        when generating the output file name.
    cfg : IoCfg
        IO configuration containing the destination directory for derived
        artefacts.
    date : str, optional keyword-only
        Date prefix explicitly supplied by the caller. When omitted, the
        function falls back to ``cfg.default_date_prefix`` or generates a
        ``YYYYMMDD`` string for the current UTC date when the effective
        ``stamp_mode`` resolves to ``"date"``.
    stamp_mode : {"omit", "date", "require"}, optional keyword-only
        Overrides the configured strategy controlling how date components are
        appended to the output file name. When ``None`` (the default) the
        function consults ``cfg.output_stamp_mode`` and uses ``"date"`` as the
        ultimate fallback for backwards compatibility.

    Returns
    -------
    pathlib.Path
        Path pointing to ``cfg.output_dir / f"output_<stem>_<YYYYMMDD>.csv"``.

    Raises
    ------
    TypeError
        Raised when ``cfg.output_dir`` cannot be coerced into a
        :class:`pathlib.Path` instance.
    """

    inp = Path(input_path)

    mode = stamp_mode or getattr(cfg, "output_stamp_mode", None) or "date"
    if isinstance(mode, str):
        normalized_mode = mode.strip().lower()
    else:
        normalized_mode = "date"

    if normalized_mode not in {"omit", "date", "require"}:
        raise ValueError(
            "stamp_mode must be 'omit', 'date' or 'require'",
        )

    candidate_date: str | None = None
    if isinstance(date, str):
        candidate_date = date.strip() or None

    if candidate_date is not None:
        date_str = candidate_date
        return Path(cfg.output_dir) / f"output.{inp.stem}_{date_str}.csv"

    if normalized_mode == "omit":
        return Path(cfg.output_dir) / f"output.{inp.stem}.csv"

    if normalized_mode == "require":
        raise ValueError("date must be provided when stamp_mode is 'require'")

    cfg_prefix = getattr(cfg, "default_date_prefix", None)
    if isinstance(cfg_prefix, str):
        cfg_prefix = cfg_prefix.strip() or None
    if cfg_prefix is not None:
        date_str = cfg_prefix
    else:
        date_str = datetime.now(UTC).strftime("%Y%m%d")

    return Path(cfg.output_dir) / f"output.{inp.stem}_{date_str}.csv"

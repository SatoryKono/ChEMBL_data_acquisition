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
        ``YYYYMMDD`` string for the current UTC date.

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

    cfg_prefix = getattr(cfg, "default_date_prefix", None)
    if isinstance(cfg_prefix, str):
        cfg_prefix = cfg_prefix.strip() or None
    if cfg_prefix is not None:
        date_str = cfg_prefix
    else:
        if isinstance(date, str):
            candidate = date.strip()
            date_str = candidate or datetime.now(UTC).strftime("%Y%m%d")
        else:
            date_str = datetime.now(UTC).strftime("%Y%m%d")

    return Path(cfg.output_dir) / f"output.{inp.stem}_{date_str}.csv"

"""Path helper functions for pipeline outputs."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from ..config import IoCfg


def default_output_path(input_path: str | Path, cfg: IoCfg) -> Path:
    """Construct a deterministic output path based on the input file name.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Location of the CSV supplied to the pipeline. Only the stem is used
        when generating the output file name.
    cfg : IoCfg
        IO configuration containing the destination directory for derived
        artefacts.

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
    date_str = datetime.now(UTC).strftime("%Y%m%d")
    return Path(cfg.output_dir) / f"output.{inp.stem}_{date_str}.csv"

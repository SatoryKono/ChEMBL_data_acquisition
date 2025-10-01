"""Path helper functions for pipeline outputs."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from ..config import IoCfg


def default_output_path(input_path: str | Path, cfg: IoCfg) -> Path:
    """Return the default output path for ``input_path``."""

    inp = Path(input_path)
    date_str = datetime.now().strftime("%Y%m%d")
    return Path(cfg.output_dir) / f"output.{inp.stem}_{date_str}.csv"

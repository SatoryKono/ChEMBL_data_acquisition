from __future__ import annotations

from pathlib import Path

import argparse
import pytest

from scripts import get_data


def _args(base: Path, *, limit: int | None) -> argparse.Namespace:
    config_path = base / "config.yaml"
    config_path.write_text("{}")
    input_dir = base / "input"
    input_dir.mkdir()
    args = argparse.Namespace(
        base_path=base,
        input_dir=Path("input"),
        output_dir=Path("output"),
        config=config_path,
        date_prefix="20240101",
        log_level="INFO",
        limit=limit,
        force=False,
        skip_existing=False,
    )
    return args


def test_negative_limit_rejected(tmp_path: Path) -> None:
    args = _args(tmp_path, limit=-1)
    with pytest.raises(ValueError, match="--limit must be zero or a positive integer"):
        get_data._prepare_config(args)


def test_zero_limit_allowed(tmp_path: Path) -> None:
    args = _args(tmp_path, limit=0)
    cfg = get_data._prepare_config(args)
    assert cfg.limit == 0

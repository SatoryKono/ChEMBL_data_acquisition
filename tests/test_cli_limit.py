from __future__ import annotations

import argparse
from pathlib import Path

from library import cli


def test_cli_limit_override(tmp_path: Path) -> None:
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text("")
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=cfg_path, type=Path)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args(["--config", str(cfg_path), "--limit", "5"])
    cfg = cli.apply_config_overrides(
        args, parser, args.config, mapping={"limit": "activity.limit"}
    )
    assert cfg.activity.limit == 5

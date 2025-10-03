import argparse
from pathlib import Path

import pytest

from library import cli
from library.config import ConfigError
from library.utils.config import DEFAULT_CONFIG_PATH


def test_apply_config_overrides_missing_config(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    args = parser.parse_args(["--config", str(tmp_path / "missing.yaml")])
    with pytest.raises(ConfigError, match="configuration file not found"):
        cli.apply_config_overrides(args, parser, args.config)

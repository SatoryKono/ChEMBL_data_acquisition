import argparse
from pathlib import Path

import pytest

from library import cli


def test_apply_config_overrides_missing_config(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args(["--config", str(tmp_path / "missing.yaml")])
    with pytest.raises(SystemExit) as excinfo:
        cli.apply_config_overrides(args, parser, args.config)
    assert excinfo.value.code == 2
    assert "configuration file not found" in capsys.readouterr().err

import argparse
from pathlib import Path

import pytest

from library.cli import add_config_argument, apply_config_overrides


def test_apply_config_overrides_missing_config(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    parser = argparse.ArgumentParser()
    add_config_argument(parser)
    args = parser.parse_args(["--config", str(tmp_path / "missing.yaml")])
    with pytest.raises(SystemExit) as excinfo:
        apply_config_overrides(args, parser, args.config)
    assert excinfo.value.code == 2
    assert "configuration file not found" in capsys.readouterr().err

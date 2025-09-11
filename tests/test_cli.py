import argparse
from pathlib import Path

import pytest

from library.cli import apply_config_overrides


def test_apply_config_overrides_missing_config(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args(["--config", str(tmp_path / "missing.yaml")])
    with pytest.raises(SystemExit) as excinfo:
        apply_config_overrides(args, parser, args.config)
    assert excinfo.value.code == 2
    assert "configuration file not found" in capsys.readouterr().err


def test_no_config_strict_allows_unknown(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg_path = tmp_path / "cfg.yaml"
    cfg_path.write_text("unknown: 1\n", encoding="utf8")
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=cfg_path, type=Path)
    parser.add_argument(
        "--config-strict",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    monkeypatch.setenv(
        "CHEMBL_DA__API__USER_AGENT", "test-agent/1.0 (mailto:test@example.org)"
    )
    args = parser.parse_args(["--config", str(cfg_path), "--no-config-strict"])
    cfg = apply_config_overrides(args, parser, args.config)
    assert cfg.api.rps == 5

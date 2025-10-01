import argparse
from pathlib import Path

from library import cli
from library.cli import add_common_arguments, path_argument


def _write_config(tmp_path: Path) -> Path:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "local:\n  io:\n    csv_sep: '|'\n    csv_encoding: latin1\n"
        "system:\n  log:\n    level: INFO\n"
    )
    return cfg


def test_config_defaults_applied(tmp_path: Path) -> None:
    cfg_path = _write_config(tmp_path)
    parser = argparse.ArgumentParser()
    add_common_arguments(parser)
    parser.add_argument("--config", default=cfg_path, type=path_argument)
    args = parser.parse_args(["--config", str(cfg_path)])
    cfg = cli.apply_config_overrides(args, parser, args.config)
    assert args.sep == "|"
    assert args.encoding == "latin1"
    assert cfg.io.csv_sep == "|"
    assert cfg.io.csv_encoding == "latin1"


def test_cli_overrides_config(tmp_path: Path) -> None:
    cfg_path = _write_config(tmp_path)
    parser = argparse.ArgumentParser()
    add_common_arguments(parser)
    parser.add_argument("--config", default=cfg_path, type=path_argument)
    args = parser.parse_args(
        [
            "--config",
            str(cfg_path),
            "--sep",
            ";",
            "--encoding",
            "utf16",
            "--log-level",
            "DEBUG",
        ]
    )
    cfg = cli.apply_config_overrides(args, parser, args.config)
    assert cfg.io.csv_sep == ";"
    assert cfg.io.csv_encoding == "utf16"
    assert cfg.log.level == "DEBUG"


def test_path_argument_normalises_windows_style(tmp_path: Path) -> None:
    cfg_path = _write_config(tmp_path)
    parser = argparse.ArgumentParser()
    add_common_arguments(parser)
    parser.add_argument("--config", default=cfg_path, type=path_argument)
    args = parser.parse_args(
        ["--config", str(cfg_path), "--input", "data\\input\\file.csv"]
    )
    assert args.input_csv == Path("data/input/file.csv")

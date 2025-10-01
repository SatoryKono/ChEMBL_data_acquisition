from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from library.utils.config import DEFAULT_CONFIG_PATH
from library.utils.cli_tools import mapper_batch_main
from library.utils.cli_tools.table_quality_main import main


def test_table_quality_cli_with_config(tmp_path: Path) -> None:
    csv_path = tmp_path / "data.csv"
    csv_path.write_text("a\n1\n")
    config = Path("tests/fixtures/config.min.yaml")
    rc = main(
        [
            "--config",
            str(config),
            "--input",
            str(csv_path),
            "--table-name",
            "demo",
            "--output",
            str(tmp_path),
        ]
    )
    assert rc == 0


def test_mapper_batch_default_config_outside_repo(tmp_path: Path, monkeypatch) -> None:
    recorded: dict[str, Path] = {}

    def fake_apply_config_overrides(
        args,
        parser,
        config_path,
        mapping=None,
        **kwargs,
    ):
        recorded["config_path"] = Path(config_path)
        return SimpleNamespace(io=SimpleNamespace(exist_ok=True))

    monkeypatch.setattr(
        mapper_batch_main.cli,
        "apply_config_overrides",
        fake_apply_config_overrides,
    )
    monkeypatch.setattr(mapper_batch_main, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(
        mapper_batch_main,
        "run",
        lambda cfg, args: 0,
    )
    monkeypatch.chdir(tmp_path)

    rc = mapper_batch_main.main([])

    assert rc == 0
    assert recorded["config_path"] == DEFAULT_CONFIG_PATH

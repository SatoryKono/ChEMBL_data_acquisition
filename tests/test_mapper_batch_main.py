from __future__ import annotations

from pathlib import Path

import pytest

from library.config import Config
from library.utils.cli_tools import mapper_batch_main as cli


class DummyLogger:
    """Capture log events emitted via :func:`cli.main`."""

    def __init__(self) -> None:
        self.events: list[tuple[str, dict[str, object]]] = []

    def info(self, event: str, **payload: object) -> None:
        self.events.append((event, payload))


def test_main_invokes_run(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """CLI entry point should configure logging and execute ``run``."""

    logger = DummyLogger()
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: logger)

    cfg = Config()
    cfg.api.user_agent = "test@example.org"
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *_, **__: cfg)
    monkeypatch.setattr(cli, "ensure_dirs", lambda *_: None)

    printed: dict[str, Config] = {}
    monkeypatch.setattr(cli, "print_config", lambda c: printed.setdefault("cfg", c))

    called: dict[str, object] = {}

    def fake_run(config: Config, args) -> int:
        called["cfg"] = config
        called["args"] = args
        return 0

    monkeypatch.setattr(cli, "run", fake_run)

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("chembl_id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "out.csv"

    exit_code = cli.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--column",
            "chembl_id",
            "--chunk-size",
            "2",
        ]
    )

    assert exit_code == 0
    assert called["cfg"] is cfg
    assert called["args"].input_csv == input_csv
    assert printed["cfg"] is cfg
    assert any(event == "pipeline_start" for event, _ in logger.events)
    assert any(
        event == "pipeline_end" and payload.get("exit_code") == 0
        for event, payload in logger.events
    )

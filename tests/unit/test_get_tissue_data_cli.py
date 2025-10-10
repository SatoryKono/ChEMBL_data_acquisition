"""Unit tests for the :mod:`scripts.get_tissue_data` CLI helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts import get_tissue_data


class _StubLogger:
    """Collects structured logging events for assertions."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def bind(self, **_: object) -> _StubLogger:  # pragma: no cover - interface parity
        return self

    def info(self, event: str, **data: object) -> None:
        self.events.append(("info", event, dict(data)))

    def warning(
        self, event: str, **data: object
    ) -> None:  # pragma: no cover - defensive
        self.events.append(("warning", event, dict(data)))

    def error(self, event: str, **data: object) -> None:  # pragma: no cover - defensive
        self.events.append(("error", event, dict(data)))


@pytest.mark.unit
def test_get_tissue_data_main__limit_zero_skips_pipeline(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``main()`` should short-circuit when ``--limit 0`` is provided."""

    stub_logger = _StubLogger()
    monkeypatch.setattr(get_tissue_data, "logger", stub_logger)

    def _unexpected_run_cli_command(**_: object) -> None:
        raise AssertionError("run_cli_command should not be invoked when limit is zero")

    monkeypatch.setattr(get_tissue_data, "run_cli_command", _unexpected_run_cli_command)
    monkeypatch.setattr(
        get_tissue_data, "configure_logger", lambda *_args, **_kwargs: None
    )

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("tissue_chembl_id\nCHEMBLT0\n", encoding="utf-8")
    final_out = tmp_path / "output.csv"

    exit_code = get_tissue_data.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(final_out),
            "--limit",
            "0",
        ]
    )

    assert exit_code == 0
    assert ("info", "pipeline_skip_limit", {"limit": 0}) in stub_logger.events
    assert final_out.exists() is False

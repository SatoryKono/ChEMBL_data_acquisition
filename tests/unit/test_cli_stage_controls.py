from __future__ import annotations

import pytest

from library.pipelines.testitem import cli


class DummyLogger:
    def __init__(self) -> None:
        self.messages: list[tuple[str, str, dict[str, object]]] = []

    def info(self, event: str, **payload: object) -> None:
        self.messages.append(("info", event, dict(payload)))

    def warning(self, event: str, **payload: object) -> None:
        self.messages.append(("warning", event, dict(payload)))

    def error(self, event: str, **payload: object) -> None:
        self.messages.append(("error", event, dict(payload)))


@pytest.mark.unit
def test_stage_execution_budget__start_logs_budget(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monotonic_calls = [0.0]

    def fake_monotonic() -> float:
        return monotonic_calls[0]

    monkeypatch.setattr(cli.time, "monotonic", fake_monotonic)
    logger = DummyLogger()

    budget = cli.StageExecutionBudget("fetch", minutes=1, logger=logger)
    budget.start()

    assert (
        "info",
        "fetch_execution_budget_started",
        {"budget_seconds": 60},
    ) in logger.messages


@pytest.mark.unit
def test_stage_execution_budget__enforce_within_budget(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    now = [0.0]

    def fake_monotonic() -> float:
        return now[0]

    monkeypatch.setattr(cli.time, "monotonic", fake_monotonic)
    logger = DummyLogger()

    budget = cli.StageExecutionBudget("stage", minutes=1, logger=logger)
    budget.start()
    now[0] = 30.0
    budget.enforce("chunk")

    assert all(level != "error" for level, _, _ in logger.messages)


@pytest.mark.unit
def test_stage_execution_budget__enforce_after_budget_raises(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    now = [0.0]

    def fake_monotonic() -> float:
        return now[0]

    monkeypatch.setattr(cli.time, "monotonic", fake_monotonic)
    logger = DummyLogger()

    budget = cli.StageExecutionBudget("stage", minutes=1 / 60, logger=logger)
    budget.start()
    now[0] = 5.0

    with pytest.raises(cli.TestitemPipelineStageError):
        budget.enforce("batch")

    assert any(
        event == "stage_execution_budget_exhausted" for _, event, _ in logger.messages
    )


@pytest.mark.unit
def test_stage_watchdog__start_without_timeout(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(cli.time, "monotonic", lambda: 0.0)
    logger = DummyLogger()

    watchdog = cli.StageWatchdog("stage", idle_minutes=0, logger=logger)
    watchdog.start()

    assert (
        watchdog._thread is None
    )  # noqa: SLF001 - internal state asserted for determinism


@pytest.mark.unit
def test_stage_watchdog__ping_logs_progress(monkeypatch: pytest.MonkeyPatch) -> None:
    monotonic_values = [0.0]

    def fake_monotonic() -> float:
        return monotonic_values[0]

    monkeypatch.setattr(cli.time, "monotonic", fake_monotonic)
    logger = DummyLogger()

    class FakeThread:
        def __init__(self, *_, **__):
            self.started = False

        def start(self) -> None:
            self.started = True

        def join(self, timeout: float | None = None) -> None:  # noqa: ARG002
            self.started = False

    monkeypatch.setattr(cli.threading, "Thread", FakeThread)

    watchdog = cli.StageWatchdog("stage", idle_minutes=1, logger=logger)
    watchdog.start()
    monotonic_values[0] = 1.0
    watchdog.ping("chunk", processed=10)

    assert (
        "info",
        "stage_watchdog_progress",
        {"watchdog_event": "chunk", "processed": 10},
    ) in logger.messages


@pytest.mark.unit
def test_stage_watchdog__raise_if_timed_out(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(cli.time, "monotonic", lambda: 0.0)
    logger = DummyLogger()

    watchdog = cli.StageWatchdog("stage", idle_minutes=1, logger=logger)
    watchdog._timed_out.set()

    with pytest.raises(cli.TestitemPipelineStageError):
        watchdog.raise_if_timed_out()


@pytest.mark.unit
def test_stage_watchdog__monitor_emits_timeout(monkeypatch: pytest.MonkeyPatch) -> None:
    monotonic_values = [0.0, 120.0]

    def fake_monotonic() -> float:
        return monotonic_values.pop(0)

    monkeypatch.setattr(cli.time, "monotonic", fake_monotonic)
    logger = DummyLogger()

    watchdog = cli.StageWatchdog(
        "stage", idle_minutes=2, logger=logger, check_interval=1
    )
    watchdog._idle_timeout_seconds = 60
    watchdog._effective_interval = 1
    watchdog._last_activity = 0.0

    class StubEvent:
        def __init__(self) -> None:
            self.calls = 0
            self.flag = False

        def wait(self, _timeout: float) -> bool:
            self.calls += 1
            return self.calls > 1

        def set(self) -> None:
            self.flag = True

        def clear(self) -> None:
            self.flag = False

        def is_set(self) -> bool:
            return self.flag

    watchdog._stop_event = StubEvent()
    watchdog._timed_out = StubEvent()

    watchdog._monitor()

    assert any(event == "stage_watchdog_timeout" for _, event, _ in logger.messages)
    assert watchdog._timed_out.flag is True


@pytest.mark.unit
@pytest.mark.parametrize(
    "items, size, expected",
    [
        ([], 3, []),
        (["a"], 2, [["a"]]),
        (["a", "b", "c", "d"], 2, [["a", "b"], ["c", "d"]]),
        (["a", "b", "c"], 5, [["a", "b", "c"]]),
    ],
)
def test_batched__yields_expected_groups(
    items: list[str], size: int, expected: list[list[str]]
) -> None:
    batches = list(cli._batched(items, size))
    assert batches == expected


@pytest.mark.unit
def test_log_missing_identifier_summary__emits_warning(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = DummyLogger()
    monkeypatch.setattr(cli, "logger", logger)

    cli._log_missing_identifier_summary(["CHEMBL1", "CHEMBL2"])

    assert (
        "warning",
        "chembl_missing_identifiers",
        {"total": 2, "sample": ["CHEMBL1", "CHEMBL2"]},
    ) in logger.messages


@pytest.mark.unit
def test_log_missing_identifier_summary__skips_empty(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = DummyLogger()
    monkeypatch.setattr(cli, "logger", logger)

    cli._log_missing_identifier_summary([])

    assert not logger.messages

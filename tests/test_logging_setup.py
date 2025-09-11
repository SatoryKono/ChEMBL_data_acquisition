"""Tests for :mod:`chembl_da.library.logging_setup`."""

from __future__ import annotations

from io import StringIO
import json
import pytest

from chembl_da.library.logging_setup import LoggerConfig, configure_logger


def parse_lines(buffer: StringIO) -> list[dict[str, object]]:
    buffer.seek(0)
    return [json.loads(line) for line in buffer.getvalue().splitlines() if line]


def test_basic_logging_and_redaction() -> None:
    stream = StringIO()
    cfg = LoggerConfig(level="INFO", run_id="r1", stream=stream)
    logger = configure_logger(cfg)
    logger.info("test_event", api_token="secret", value=1)
    lines = parse_lines(stream)
    assert lines[0]["event"] == "test_event"
    assert lines[0]["api_token"] == "***"
    assert lines[0]["value"] == 1
    assert lines[0]["run_id"] == "r1"


def test_level_filtering() -> None:
    stream = StringIO()
    cfg = LoggerConfig(level="ERROR", run_id="r2", stream=stream)
    logger = configure_logger(cfg)
    logger.info("will_skip")
    logger.error("will_log")
    lines = parse_lines(stream)
    assert len(lines) == 1
    assert lines[0]["event"] == "will_log"


def test_exception_logging() -> None:
    stream = StringIO()
    cfg = LoggerConfig(level="INFO", run_id="r3", stream=stream)
    logger = configure_logger(cfg)
    try:
        raise ValueError("boom")
    except ValueError as exc:
        logger.exception("error_event", exc)
    rec = parse_lines(stream)[0]
    assert rec["event"] == "error_event"
    assert rec["exc_type"] == "ValueError"
    assert rec["exc_message"] == "boom"
    assert "ValueError" in str(rec["traceback"])


def test_stage_context_manager_success() -> None:
    stream = StringIO()
    cfg = LoggerConfig(level="INFO", run_id="r4", stream=stream)
    logger = configure_logger(cfg)
    with logger.stage("stage1") as staged:
        staged.info("inside")
    lines = parse_lines(stream)
    events = [line["event"] for line in lines]
    assert events == ["stage1_start", "inside", "stage1_done"]
    assert "elapsed" in lines[-1]


def test_stage_context_manager_failure() -> None:
    stream = StringIO()
    cfg = LoggerConfig(level="INFO", run_id="r5", stream=stream)
    logger = configure_logger(cfg)
    with pytest.raises(RuntimeError):
        with logger.stage("stage2"):
            raise RuntimeError("fail")
    events = [line["event"] for line in parse_lines(stream)]
    assert events == ["stage2_start", "stage2_fail"]

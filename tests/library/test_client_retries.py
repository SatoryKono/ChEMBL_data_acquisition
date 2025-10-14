from __future__ import annotations

from types import SimpleNamespace

import pytest

from library.clients import pubmed


def test_do_request_retries_until_success(monkeypatch: pytest.MonkeyPatch) -> None:
    responses = [
        (503, "service unavailable", None, "", {}),
        (503, "try again", None, "", {}),
        (200, "{}", {"result": "ok"}, "", {}),
    ]
    calls: list[int] = []

    def _fake_make_request(*args, **kwargs):  # noqa: ANN002, ANN003 - signature dictated by patch
        index = len(calls)
        calls.append(index)
        return responses[index]

    monkeypatch.setattr(pubmed, "_make_request", _fake_make_request)

    data, error = pubmed._do_request(  # pylint: disable=protected-access
        SimpleNamespace(),
        "https://example.test/resource",
        delay=0.0,
        retries=2,
    )

    assert data == {"result": "ok"}
    assert error == ""
    assert calls == [0, 1, 2]


def test_do_request_returns_error_after_retry_exhaustion(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _make_request_fail(*args, **kwargs):  # noqa: ANN002, ANN003 - signature dictated by patch
        return 500, "temporary failure", None, "", {}

    monkeypatch.setattr(pubmed, "_make_request", _make_request_fail)

    data, error = pubmed._do_request(  # pylint: disable=protected-access
        SimpleNamespace(),
        "https://example.test/resource",
        delay=0.0,
        retries=1,
    )

    assert data is None
    assert error.startswith("HTTP 500")

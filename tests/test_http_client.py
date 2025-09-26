"""Tests for HTTP client session utilities."""

from __future__ import annotations

import pytest

responses = pytest.importorskip("responses")

from library.config import ApiCfg, RetryCfg, session_with_retry  # noqa: E402

USER_AGENT = "test-agent/1.0 (mailto:test@example.org)"


@responses.activate
def test_session_with_retry_makes_single_attempt() -> None:
    """The HTTP adapter leaves retry attempts to higher-level clients."""

    url = "http://example.com/post"
    responses.add(responses.POST, url, status=500)
    responses.add(responses.POST, url, json={"ok": True}, status=200)

    session = session_with_retry(
        ApiCfg(user_agent=USER_AGENT),
        RetryCfg(max_attempts=2, backoff_factor=0),
    )
    response = session.post(url)

    assert response.status_code == 500
    assert len(responses.calls) == 1

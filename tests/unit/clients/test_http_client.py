"""Tests for HTTP client session utilities."""

from __future__ import annotations

import pytest

responses = pytest.importorskip("responses")

from library.config import ApiCfg, RetryCfg, session_with_retry  # noqa: E402

USER_AGENT = "test-agent/1.0 (mailto:test@example.org)"


@responses.activate
def test_session_with_retry_retries_configured_attempts() -> None:
    """The session retries failed responses up to the configured attempts."""

    url = "http://example.com/post"
    max_attempts = 3
    for _ in range(max_attempts):
        responses.add(responses.POST, url, status=500)

    session = session_with_retry(
        ApiCfg(user_agent=USER_AGENT),
        RetryCfg(max_attempts=max_attempts, backoff_factor=0),
    )
    response = session.post(url)

    assert response.status_code == 500
    assert len(responses.calls) == max_attempts

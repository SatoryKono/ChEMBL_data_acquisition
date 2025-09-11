"""Tests for HTTP client session utilities."""

from __future__ import annotations

import responses

from library.config import ApiCfg, RetryCfg, session_with_retry


@responses.activate
def test_session_with_retry_retries_post() -> None:
    """Ensure POST requests are retried when configured."""

    url = "http://example.com/post"
    responses.add(responses.POST, url, status=500)
    responses.add(responses.POST, url, json={"ok": True}, status=200)

    session = session_with_retry(ApiCfg(), RetryCfg(max_attempts=2, backoff_factor=0))
    response = session.post(url)

    assert response.status_code == 200
    assert len(responses.calls) == 2

"""Test configuration fixtures.

This module provides common pytest fixtures used across the test suite.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from library.config import Config
import scripts.get_document_data as gdd


@pytest.fixture(autouse=True)
def disable_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Disallow outbound HTTP requests during tests."""

    import requests

    original_request = requests.sessions.Session.request

    def deny(self, method, url, *args, **kwargs):  # type: ignore[override]
        try:
            import responses  # type: ignore[import-not-found]
        except ModuleNotFoundError:  # pragma: no cover - optional dependency
            responses_active = False
        else:
            mock = getattr(responses, "_default_mock", None)
            responses_active = bool(getattr(mock, "_is_mocked", False))

        if responses_active:
            return original_request(self, method, url, *args, **kwargs)

        raise AssertionError("External network access is disabled during tests")

    monkeypatch.setattr(requests.sessions.Session, "request", deny)


def _disable_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Prevent external HTTP requests during tests for determinism."""

    try:
        import requests
    except ModuleNotFoundError:  # pragma: no cover - requests optional in env
        return

    def _deny_request(self, method, url, *args, **kwargs):  # type: ignore[override]
        msg = (
            "External HTTP requests are disabled during tests; "
            f"attempted {method} {url}"
        )
        raise RuntimeError(msg)

    monkeypatch.setattr("requests.sessions.Session.request", _deny_request)


@pytest.fixture()
def cfg() -> Config:
    """Return a baseline :class:`~library.config.Config` instance for tests.

    The configuration requires a valid ``api.user_agent`` value; the fixture
    supplies a deterministic placeholder address so that individual tests do
    not need to construct :class:`Config` instances manually.
    """

    cfg = Config()
    cfg.api.user_agent = "test@example.org"
    return cfg


@pytest.fixture()
def duplicate_document_ids() -> list[str]:
    """Return sample document IDs including duplicates for testing."""

    return ["CHEMBL1", "CHEMBL1", "CHEMBL2"]


@pytest.fixture()
def document_export_postprocess_stub(
    monkeypatch: pytest.MonkeyPatch,
) -> list[Path]:
    """Stub ``postprocess_export_file`` to avoid filesystem I/O in tests."""

    created: list[Path] = []

    def fake_postprocess(
        path: Path,
        *,
        cfg: object,
        output_path: Path | None = None,
        **_: object,
    ) -> Path:
        destination = Path(output_path) if output_path else Path(path).with_name(
            f"preprocessed_{Path(path).name}"
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text("stub")
        created.append(destination)
        return destination

    monkeypatch.setattr(
        gdd.document_export_postprocessing,
        "postprocess_export_file",
        fake_postprocess,
    )
    return created

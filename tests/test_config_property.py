"""Property-based tests for environment configuration overrides."""

from __future__ import annotations

import string
from typing import Any, cast

import pytest
pytest.importorskip("hypothesis")
from hypothesis import assume, given
from hypothesis import strategies as st

from library import config


def _get_by_path(data: dict[str, Any], path: list[str]) -> Any:
    """Return nested value from *data* following *path* components."""

    cur: Any = data
    for part in path:
        cur = cur[part]
    return cur


_VALID_ENV_PATHS: dict[str, list[str]] = {
    "CHEMBL_DA__API__RPS": ["api", "rps"],
    "CHEMBL_DA__OPENALEX__RPS": ["openalex", "rps"],
}


@given(
    item=st.sampled_from(list(_VALID_ENV_PATHS.items())),
    value=st.text(alphabet=string.ascii_letters + string.digits),
)
def test_apply_env_overrides_valid(item: tuple[str, list[str]], value: str) -> None:
    """Valid environment variables override configuration values."""

    env_key, path = item
    data: dict[str, Any] = {}
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setenv(env_key, value)
        config._apply_env_overrides(data)

    assert _get_by_path(data, path) == value


@given(
    parts=st.lists(
        st.text(alphabet=string.ascii_uppercase, min_size=1),
        min_size=1,
        max_size=3,
    ),
    value=st.text(alphabet=string.ascii_letters + string.digits),
)
def test_apply_env_overrides_invalid(parts: list[str], value: str) -> None:
    """Unknown environment variables are ignored and trigger warnings."""

    path = [p.lower() for p in parts]
    assume(not config._is_valid_path(path))
    env_key = "CHEMBL_DA__" + "__".join(parts)
    data: dict[str, Any] = {}
    warnings: list[str] = []

    def fake_warning(msg: str) -> None:
        warnings.append(msg)

    logger = cast(Any, config).logger

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setenv(env_key, value)
        monkeypatch.setattr(logger, "warning", fake_warning)
        config._apply_env_overrides(data)

    assert data == {}
    assert any(env_key in msg for msg in warnings)

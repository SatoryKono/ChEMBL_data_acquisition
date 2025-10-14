"""Simple configuration loader with environment overrides."""

from __future__ import annotations

import os
from collections.abc import Mapping, MutableMapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from .path_utils import ROOT

_ENV_PREFIX = "CHEMBL_DA__"
_DEFAULT_CONFIG_NAME = "config.yaml"
DEFAULT_CONFIG_PATH = ROOT / "config" / _DEFAULT_CONFIG_NAME


@dataclass(frozen=True)
class Config:
    """Dataclass wrapper exposing configuration content."""

    path: Path
    data: Mapping[str, Any]

    def get(self, key: str, default: Any | None = None) -> Any:
        """Return the value stored under ``key`` or ``default`` when missing."""

        return self.data.get(key, default)

    def section(self, key: str) -> Mapping[str, Any]:
        """Return the mapping stored under ``key`` raising ``KeyError`` when missing."""

        value = self.data[key]
        if not isinstance(value, Mapping):
            raise TypeError(f"Section '{key}' is not a mapping")
        return value

    def to_dict(self) -> dict[str, Any]:
        """Return a shallow copy of the configuration mapping."""

        return dict(self.data)


def _coerce_value(value: str) -> Any:
    try:
        return yaml.safe_load(value)
    except yaml.YAMLError:
        return value


def _ensure_mapping(container: MutableMapping[str, Any], key: str) -> MutableMapping[str, Any]:
    existing = container.get(key)
    if isinstance(existing, MutableMapping):
        return existing
    nested: dict[str, Any] = {}
    container[key] = nested
    return nested


def _apply_env_overrides(data: MutableMapping[str, Any]) -> None:
    prefix_len = len(_ENV_PREFIX)
    for env_key, raw_value in os.environ.items():
        if not env_key.startswith(_ENV_PREFIX):
            continue
        tokens = [token for token in env_key[prefix_len:].split("__") if token]
        if not tokens:
            continue
        target = data
        for token in tokens[:-1]:
            target = _ensure_mapping(target, token.lower())
        target[tokens[-1].lower()] = _coerce_value(raw_value)


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise TypeError("Top-level configuration structure must be a mapping")
    return data


def load_config(path: str | Path | None = None) -> Config:
    """Load ``config.yaml`` applying environment overrides."""

    if path is None:
        cfg_path = DEFAULT_CONFIG_PATH
    else:
        candidate = Path(path)
        if not candidate.is_absolute():
            candidate = ROOT / candidate
        cfg_path = candidate
    data = _load_yaml(cfg_path)
    _apply_env_overrides(data)
    return Config(path=cfg_path, data=data)


__all__ = ["Config", "DEFAULT_CONFIG_PATH", "load_config"]

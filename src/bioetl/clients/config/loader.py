from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Mapping

import yaml

from bioetl.clients.config.models import ClientConfig, ClientsCatalog, default_config_path


def load_yaml_config(path: Path) -> Mapping[str, Any]:
    with path.open("r", encoding="utf-8") as fp:
        return yaml.safe_load(fp) or {}


def load_clients_catalog(path: Path) -> ClientsCatalog:
    raw = load_yaml_config(path)
    return ClientsCatalog.load_from_dict(raw)


def load_client_config(source: str, *, path: Path | None = None) -> ClientConfig:
    config_path = path or default_config_path(source)
    payload = load_yaml_config(config_path)
    config = ClientConfig.from_mapping(source, payload)
    _apply_env_overrides(config)
    return config


def _apply_env_overrides(config: ClientConfig) -> None:
    """Переносит секреты из переменных окружения в заголовки/параметры."""

    token_env = os.getenv(f"{config.name.upper()}_TOKEN")
    api_key_env = os.getenv(f"{config.name.upper()}_API_KEY")
    if token_env:
        config.headers = {**config.headers, "Authorization": f"Bearer {token_env}"}
    if api_key_env:
        config.headers = {**config.headers, "X-API-KEY": api_key_env}

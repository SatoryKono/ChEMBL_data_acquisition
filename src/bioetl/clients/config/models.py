from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from pydantic import BaseModel, Field, HttpUrl

from bioetl.clients.base.types import PaginationConfig


class ClientRoute(BaseModel):
    """Описание маршрута клиента."""

    name: str
    path: str


class ClientConfig(BaseModel):
    """Конфигурация клиента внешнего источника."""

    name: str = Field(..., description="Имя источника")
    base_url: HttpUrl | str
    routes: Mapping[str, str]
    headers: Mapping[str, str] = Field(default_factory=dict)
    params: Mapping[str, Any] = Field(default_factory=dict)
    pagination: PaginationConfig | None = None
    timeout: float = 30.0

    @classmethod
    def from_mapping(cls, name: str, payload: Mapping[str, Any]) -> "ClientConfig":
        payload = {"name": name, **payload}
        return cls(**payload)

    @property
    def normalized_base_url(self) -> str:
        return str(self.base_url).rstrip("/")

    class Config:
        arbitrary_types_allowed = True
        validate_assignment = True


class ClientsCatalog(BaseModel):
    """Множество конфигураций клиентов из YAML."""

    clients: Mapping[str, ClientConfig]

    @classmethod
    def load_from_dict(cls, payload: Mapping[str, Any]) -> "ClientsCatalog":
        configs = {name: ClientConfig.from_mapping(name, cfg) for name, cfg in payload.items()}
        return cls(clients=configs)

    def get(self, name: str) -> ClientConfig:
        try:
            return self.clients[name]
        except KeyError as exc:
            raise KeyError(f"Клиент '{name}' не найден в конфигурации") from exc


def default_config_path(source: str) -> Path:
    return Path(__file__).resolve().parent / "yaml" / f"{source}.yml"

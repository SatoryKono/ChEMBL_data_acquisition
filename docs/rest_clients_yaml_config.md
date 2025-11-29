# План перехода REST-клиентов на YAML-конфигурации и фабрику клиентов

## Сводка — цели перехода REST-клиентов на YAML и ожидаемый результат
- Централизация конфигурации REST-источников (базовый URL, аутентификация, ресурсы, поля ответов, пагинация) в расширяемых YAML-файлах вместо захардкоженных параметров в коде клиентов.
- Единый набор pydantic-моделей `clients.config` для безопасной загрузки/валидации конфигураций и их последующего использования клиентами.
- Общая фабрика клиентов создает тонкие реализации REST-клиентов, опираясь только на `SourceConfig` и общий HTTP-бэкенд; клиенты не содержат магических констант.
- Удаление старых способов конфигурирования: все REST-параметры источников переносятся в YAML, клиенты используют только конфигурационные объекты и общий транспорт.

## Структура YAML-конфигов для REST — пример схемы и объяснение ключевых полей
```yaml
source: chembl                # обязательное, уникальное имя источника
protocol: http                # http (по умолчанию) | https
base_url: "https://www.ebi.ac.uk/chembl/api/data"  # обязательное
default_timeout: 30.0         # опциональное, сек
rate_limit:
  requests_per_minute: 60     # опционально; можно расширять (burst/window)
auth:                         # обязательное поле с типом
  type: none                  # none | api_key | bearer | basic | custom
  header_name: null           # для api_key/bearer
  query_param: null           # если ключ передается в query
  token_env: CHEMBL_API_KEY   # имя переменной окружения (опционально)
  value: null                 # фиксированное значение (для non-secret)
  scheme: null                # для bearer/basic, например "Bearer"
headers:                      # опциональные общие заголовки
  User-Agent: "bioetl/1.0"

resources:
  activities:                 # имя ресурса, уникальное в рамках источника
    path: "/activity"        # обязательное
    method: "GET"            # обязательное; GET/POST/PUT/DELETE/HEAD/OPTIONS
    paging:                   # опционально, если ресурс пагинируется
      type: link              # link | offset | page | cursor | token
      page_param: "page"
      page_size_param: "page_size"
      default_page_size: 1000
      max_page_size: 1000
      next_link_path: "page_info.next"   # для link-пагинации
      next_token_path: null              # для cursor/token
      offset_param: null                 # для offset-пагинации
      offset_step: null
      items_path: "activities"          # путь к списку элементов (для валидации)
    query:
      fixed:                   # всегда отправляемые параметры
        format: "json"
      allowed_params:          # whitelisted параметры, пробрасываемые из вызова клиента
        - "target_chembl_id"
        - "assay_chembl_id"
        - "document_chembl_id"
      required_params: []      # опционально
    body: null                 # опционально, для POST/PUT
    headers: {}                # опциональные заголовки на уровне ресурса
    response:
      format: json             # обязательное; json | xml | text | csv
      record_path: "activities"   # обязательное, путь к массиву записей
      id_path: "activity_id"       # опционально, для согласования fetch_one
      fields:                  # обязательное: перечень извлекаемых полей
        - name: "activity_id"
          path: "activity_id"
          type: int
        - name: "standard_value"
          path: "standard_value"
          type: float
        - name: "standard_units"
          path: "standard_units"
          type: str
      extra_metadata:          # опционально, вытаскиваемые метаданные страницы/ответа
        - name: "page"
          path: "page_info.page"
          type: int
      error_path: "error"      # опционально: путь к объекту ошибки в ответе
```

**Общие правила структуры**
- Один YAML на источник (или четко определенную группу) в `config/clients/`.
- Обязательные поля: `source`, `base_url`, `auth.type`, каждый элемент `resources.*` с `path`, `method`, `response.format`, `response.record_path`, `response.fields[*].{name,path,type}`.
- Значения по умолчанию: `protocol=http`, `default_timeout=30`, `auth.type=none`, `response.format=json`, `paging.type=offset` при наличии `page_param`/`page_size_param`.
- Расширяемость: секции `auth`, `paging`, `response` допускают дополнительные ключи (например, `backoff`, `retryable_statuses`, `compression`).

## Python-модели конфигурации — описание и примеры кода
```python
# src/bioetl/clients/config/models.py
from __future__ import annotations
from pathlib import Path
from typing import Literal, Optional
from pydantic import BaseModel, Field, validator

HttpMethod = Literal["GET", "POST", "PUT", "PATCH", "DELETE", "HEAD", "OPTIONS"]
ResponseFormat = Literal["json", "xml", "text", "csv"]
PagingType = Literal["link", "offset", "page", "cursor", "token"]
AuthType = Literal["none", "api_key", "bearer", "basic", "custom"]

class RateLimitConfig(BaseModel):
    requests_per_minute: Optional[int] = Field(None, ge=1)
    burst: Optional[int] = Field(None, ge=1)
    window_seconds: Optional[int] = Field(None, ge=1)

class AuthConfig(BaseModel):
    type: AuthType = "none"
    header_name: Optional[str] = None
    query_param: Optional[str] = None
    token_env: Optional[str] = None
    value: Optional[str] = None
    scheme: Optional[str] = None

class FieldConfig(BaseModel):
    name: str
    path: str
    type: str  # сохранится как строковый тип, маппится далее на python/pandas типы

class ResponseConfig(BaseModel):
    format: ResponseFormat = "json"
    record_path: str
    id_path: Optional[str] = None
    fields: list[FieldConfig]
    extra_metadata: list[FieldConfig] = Field(default_factory=list)
    error_path: Optional[str] = None

class PagingConfig(BaseModel):
    type: PagingType
    page_param: Optional[str] = None
    page_size_param: Optional[str] = None
    default_page_size: Optional[int] = None
    max_page_size: Optional[int] = None
    next_link_path: Optional[str] = None
    next_token_path: Optional[str] = None
    offset_param: Optional[str] = None
    offset_step: Optional[int] = None
    items_path: Optional[str] = None

class QueryConfig(BaseModel):
    fixed: dict[str, str] = Field(default_factory=dict)
    allowed_params: list[str] = Field(default_factory=list)
    required_params: list[str] = Field(default_factory=list)

class ResourceConfig(BaseModel):
    path: str
    method: HttpMethod
    paging: Optional[PagingConfig] = None
    query: QueryConfig = Field(default_factory=QueryConfig)
    body: Optional[dict] = None
    headers: dict[str, str] = Field(default_factory=dict)
    response: ResponseConfig

class SourceConfig(BaseModel):
    source: str
    protocol: str = "http"
    base_url: str
    default_timeout: float = 30.0
    rate_limit: Optional[RateLimitConfig] = None
    auth: AuthConfig = Field(default_factory=AuthConfig)
    headers: dict[str, str] = Field(default_factory=dict)
    resources: dict[str, ResourceConfig]

    @validator("protocol")
    def _validate_protocol(cls, v: str) -> str:
        if v not in {"http", "https"}:
            raise ValueError("protocol must be http or https")
        return v
```

**Загрузка YAML**
```python
# src/bioetl/clients/config/loader.py
import yaml
from pathlib import Path
from .models import SourceConfig

DEFAULT_CONFIG_ROOT = Path(__file__).resolve().parent / "../../config/clients"

class ConfigNotFoundError(FileNotFoundError):
    pass

def load_source_config(name: str, root: Path | None = None) -> SourceConfig:
    root = root or DEFAULT_CONFIG_ROOT
    path = (root / f"{name}.yml").resolve()
    if not path.exists():
        raise ConfigNotFoundError(path)
    with path.open("r", encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}
    cfg = SourceConfig(**raw)
    return cfg

def load_all_sources(root: Path | None = None) -> dict[str, SourceConfig]:
    root = root or DEFAULT_CONFIG_ROOT
    result: dict[str, SourceConfig] = {}
    for file in root.glob("*.yml"):
        cfg = load_source_config(file.stem, root=root)
        result[cfg.source] = cfg
    return result
```

## Фабрика REST-клиентов — интерфейс и пример логики выбора клиента по источнику
```python
# src/bioetl/clients/factory.py
from collections.abc import Mapping
from typing import Any
from bioetl.clients.base.http import HttpBackend
from bioetl.clients.base.client import RestDataClient
from bioetl.clients.config.loader import load_source_config
from bioetl.clients.config.models import SourceConfig

class ClientFactoryContext(BaseModel):
    http_backend: HttpBackend | None = None
    # можно добавить трейсинг/логирование/метрики

_CLIENT_REGISTRY: dict[str, type[RestDataClient]] = {
    "chembl": ChemblClient,
    "pubchem": PubChemClient,
    "pubmed": PubMedClient,
    "openalex": OpenAlexClient,
    "crossref": CrossrefClient,
    "semanticscholar": SemanticScholarClient,
    "uniprot": UniProtClient,
}


def create_client(
    source: str,
    *,
    config: SourceConfig | None = None,
    http_backend: HttpBackend | None = None,
    context: ClientFactoryContext | None = None,
) -> RestDataClient:
    cfg = config or load_source_config(source)
    backend = http_backend or HttpBackend.from_config(cfg)
    cls = _CLIENT_REGISTRY.get(cfg.source.lower())
    if cls is None:
        raise ValueError(f"Unknown client source: {cfg.source}")
    return cls(config=cfg, http_backend=backend, context=context)
```

**Ключевые моменты фабрики**
- Загружает `SourceConfig` из YAML при отсутствии переданного конфига.
- Создает/передает единый `HttpBackend` (сессии, ретраи, таймауты, логирование, трейсинг, пагинация).
- Выбирает класс клиента через реестр по имени источника; клиенты принимают только `SourceConfig` и `HttpBackend`.
- Внутри клиентов запрещены захардкоженные пути/заголовки — используются поля `config.resources[resource_name]`.

## План миграции REST-параметров в YAML — пошаговый план + таблица соответствий
1. **Инвентаризация**: для каждого источника найти в коде базовые URL, маршруты, query/headers, схемы полей, правила пагинации (обычно в `providers/*`, `*client_factory_impl.py`, нормализаторах и стратегиях пагинации).
2. **Создание YAML**: на основе текущих значений создать `config/clients/<source>.yml` по общей схеме (base_url, auth, ресурсы, response.fields, paging).
3. **Добавление моделей**: внедрить `clients.config.models` и `clients.config.loader` для загрузки YAML, настроить валидацию обязательных полей.
4. **Переписывание клиентов**: обновить реализации `clients/<source>/client.py` так, чтобы они использовали только `SourceConfig` и общий `HttpBackend`; удалить локальные константы/старые конфиги.
5. **Обновление фабрик**: заменить разрозненные `*client_factory_impl.py` на единую `clients.factory` с реестром клиентов.
6. **Удаление устаревших конфигов**: убрать старые настройки/константы, адаптеры и фабрики совместимости.
7. **Тестирование**: добавить тесты загрузки YAML и smoke-тесты HTTP-запросов с использованием новых конфигов.

| Источник | Текущий путь параметров (константы/настройки) | Новый YAML-файл | Новый модуль конфигурации |
|----------|-----------------------------------------------|-----------------|---------------------------|
| ChEMBL | `clients/chembl` (будет добавлен), константы в клиенте | `config/clients/chembl.yml` | `clients/config/models.py` |
| PubChem | `providers/pubchem/pubchem_data_client_impl.py`, `pubchem_normalizer_impl.py`, `pubchem_client_factory_impl.py` | `config/clients/pubchem.yml` | `clients/config/models.py` |
| PubMed | `providers/pubmed/*`, `pubmed_client_factory_impl.py` | `config/clients/pubmed.yml` | `clients/config/models.py` |
| Crossref | `providers/crossref/*`, `crossref_client_factory_impl.py` | `config/clients/crossref.yml` | `clients/config/models.py` |
| OpenAlex | `providers/openalex/*`, `openalex_client_factory_impl.py` | `config/clients/openalex.yml` | `clients/config/models.py` |
| Semantic Scholar | `providers/semantic_scholar/*`, `semantic_scholar_client_factory_impl.py` | `config/clients/semantic_scholar.yml` | `clients/config/models.py` |
| UniProt | `providers/uniprot/*`, `uniprot_client_factory_impl.py` | `config/clients/uniprot.yml` | `clients/config/models.py` |
| IUPHAR Target | `providers/iuphar/impl/target_client_impl.py`, `target/factories.py` | `config/clients/iuphar.yml` | `clients/config/models.py` |

## Валидация и health-checks конфигураций — список проверок, которые нужно реализовать
- **Статическая валидация YAML** (pydantic):
  - обязательные поля (`source`, `base_url`, `resources.*.path`, `resources.*.method`, `resources.*.response.record_path`, `fields`).
  - метод ∈ {GET, POST, PUT, PATCH, DELETE, HEAD, OPTIONS}; формат ∈ {json, xml, text, csv}; тип пагинации ∈ {link, offset, page, cursor, token}.
  - при `auth.type` в {api_key, bearer, basic, custom} требовать `header_name` или `query_param` + `token_env`/`value`.
  - при наличии пагинации `page_size_param` ⇒ `default_page_size` задан; `next_link_path`/`next_token_path` обязательны для типов link/cursor/token.
- **Согласованность с контрактом**:
  - `response.record_path` и `paging.items_path` должны быть совместимы с типом ответа (json/xml).
  - `fields.name` уникальны в рамках ресурса.
- **Health-checks (smoke-тесты)**:
  - тестовый запрос к каждому ресурсу с минимальными параметрами (dry-run) через `HttpBackend` с respect rate-limit.
  - валидация, что `record_path` существует в ответе и содержит массив/итерабельную структуру.
  - проверка, что все `fields.path` и `extra_metadata.path` присутствуют хотя бы в первой записи ответа.
  - если есть пагинация, проверка, что `next_link_path`/`next_token_path` возвращают ожидаемый тип (str/URL/token) или отсутствуют на последней странице.
  - логирование/отчет о результатах health-check с указанием, какие поля/пути не найдены.

# Единый контракт и приведение клиентов к нему

## Сводка
- Вводится единый контракт `ExternalDataClient` с согласованными методами (`fetch_one`, `fetch_many`, `iter_pages`, `metadata`, `close`) и унифицированными вспомогательными типами (`ClientRequest`, `Pagination`, `Record`, `Page`, `RequestContext`).
- Все клиенты должны работать через общую инфраструктуру транспорта/пагинации, не раскрывая внешнему коду HTTP/DB детали и не дублируя ретраи/логирование.
- Конкретные клиенты становятся тонкими обертками: выбирают маршрут и параметры из конфигурации, делегируют вызов транспорту и возвращают сырые записи без нормализации.
- Любая бизнес-логика, нормализация, агрегирование и валидация переносится в слои пайплайна/доменных сервисов, что упрощает тестирование и замену источников.
- Старые интерфейсы и адаптеры не используются: каждая реализация напрямую имплементирует новый ABC/протокол, конфигурации хранятся в YAML per-source.

## Целевой контракт клиента
```python
from __future__ import annotations
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from typing import Any, Protocol, runtime_checkable, TypedDict

class Record(TypedDict, total=False):
    """Минимально типизированная запись; поля задаются конфигурацией источника."""

@dataclass(slots=True)
class RequestContext:
    trace_id: str | None = None
    options: Mapping[str, Any] | None = None  # таймауты, лимиты, debug-флаги

@dataclass(slots=True)
class Pagination:
    page_size: int | None = None
    cursor: str | None = None

@dataclass(slots=True)
class ClientRequest:
    route: str  # логическое имя маршрута из YAML
    params: Mapping[str, Any] | None = None  # query/path/body-параметры
    context: RequestContext | None = None
    pagination: Pagination | None = None

@dataclass(slots=True)
class Page:
    items: list[Record]
    next_cursor: str | None = None
    raw: Any | None = None  # исходный ответ транспорта

@runtime_checkable
class ExternalDataClient(Protocol):
    """Единый контракт для всех клиентов внешних источников."""

    def fetch_one(self, request: ClientRequest) -> Record | None: ...

    def fetch_many(self, request: ClientRequest) -> Iterator[Record]: ...

    def iter_pages(self, request: ClientRequest) -> Iterator[Page]: ...

    def metadata(self) -> Mapping[str, Any]: ...  # например, лимиты API, user-agent

    def close(self) -> None: ...
```
Ключевые моменты:
- Входные данные описывает `ClientRequest`: логическое имя маршрута (`route`), параметры (`params`), пагинация и `RequestContext`.
- Возвраты унифицированы: одиночная запись (`Record | None`), поток записей (`Iterator[Record]`), либо поток страниц (`Iterator[Page]`).
- Пагинация/ретраи/логирование инкапсулируются общими компонентами транспорта; клиенты не работают с «сырым» HTTP/DB объектом снаружи.
- Расширяемость достигается через `RequestContext.options`: добавление таймаутов/trace_id не требует изменения сигнатур.

### Общая инфраструктура
- `clients/base/client.py`: реализация общего `ExternalDataClient` на базе транспорта (`HttpTransport`/`DbTransport`), которая знает, как выполнять маршруты и преобразовывать ответы в `Record`/`Page`.
- `clients/base/request.py` и `clients/base/response.py`: описания маршрутов (HTTP-метод, path, query/body, заголовки, auth) и схем ответа (формат JSON/XML/CSV, ключи данных/пагинации). Эти дескрипторы создаются из YAML.
- `clients/base/pagination.py`: единые стратегии cursor/page/offset, принимающие `Pagination` и `ResponseDescriptor`.
- `clients/config/models.py`: pydantic-модели для чтения YAML per-source (endpoint, маршруты, схемы полей, лимиты, пагинация); `clients/config/loader.py` — загрузка и валидация.
- `clients/factory.py`: фабрика, принимающая конфигурационную модель и создающая конкретный клиент с нужным транспортом и маршрутами.

## Сопоставление текущих клиентов с контрактом
| Источник | Текущий класс(ы) | Методы сейчас | Соответствие целевому контракту | Что нужно изменить/удалить |
|---|---|---|---|---|
| Общие ABC | `DataProviderProtocol` и `BaseApiClient` в `clients/base/interfaces.py` | `fetch_one`, `fetch_many`, `iter_pages`, `metadata`, `close` принимают `identifier/params/context`; нет объекта запроса и поток страниц/записей не типизирован под один контракт | Частично совпадает по методам, нет `ClientRequest`, нет TypedDict `Record`, контекст ограничен; транспорт и клиент смешаны | Ввести `ExternalDataClient`, `ClientRequest`, `Pagination`, удалить или заменить старый протокол | 
| PubChem | `PubChemDataClientImpl` | Специализированные `fetch`, `search`, alias-методы с DeprecationWarnings, внутренняя нормализация и пагинация на базе `ProviderDataClientImpl` | Сигнатуры и названия отличаются, возвращаемые типы итераторы Mapping; логика нормализации внутри клиента | Заменить на методы единого интерфейса, убрать нормализацию/варнинги, параметры/роуты/пагинацию вынести в YAML | 
| OpenAlex | `OpenAlexDataClientImpl` | `_iterate_pages_impl`, `_normalize_page_payload`, `search`, алиасы `fetch/fetch_search` для совместимости | Смешение логики пагинации/нормализации, нет `ClientRequest`, методы специфичны | Перейти на единый контракт, удалить кастомную нормализацию/алиасы, маршруты и ключи пагинации в YAML | 
| Crossref | `crossref_data_client_impl`, отдельная стратегия пагинации | Методы специфичны для Crossref, завязаны на свою стратегию и нормализатор | Пагинация и нормализация зашиты, нет объекта запроса | Перевести на общий клиент + конфигурацию пагинации, убрать bespoke стратегию | 
| PubMed | `pubmed_data_client_impl` | Методы efetch/esearch, внутренняя нормализация XML | Смешение логики выборки и нормализации, сигнатуры отличны | Вынести маршруты/форматы в YAML, реализовать `ExternalDataClient` без нормализации | 
| Semantic Scholar | `semantic_scholar_data_client_impl` | Методы поиска с нормализатором и пагинацией | Нет `ClientRequest`, логика нормализации внутри | Переписать на единый контракт, нормализацию вынести | 
| UniProt | `uniprot_data_client_impl` | Методы поиска/детализации + нормализация и пагинация | Специфичные сигнатуры, логика нормализации | Перейти на тонкий клиент с YAML и общим контрактом | 
| IUPHAR Target | `target/factories.py` + `providers/iuphar` | Специальная фабрика и нормализация | Не соответствует: отдельная фабрика, доменная логика в нормализации | Перевести на общий клиент и фабрику, убрать спец. слой | 
| ChEMBL (отсутствует) | — | — | Нет клиента | Добавить тонкий клиент и YAML по единому контракту |

## Пример реализации тонкого клиента
`clients/chembl/client.py` (эскиз):
```python
from bioetl.clients.base.client import BaseExternalDataClient
from bioetl.clients.config.models import ClientConfig

class ChemblClient(BaseExternalDataClient):
    def __init__(self, config: ClientConfig, *, transport_factory: TransportFactory) -> None:
        super().__init__(config=config, transport_factory=transport_factory)

# Использование
config = load_client_config("configs/clients/chembl.yaml")
client = ChemblClient(config, transport_factory=http_transport_factory)
req = ClientRequest(route="fetch_molecule", params={"chembl_id": "CHEMBL25"})
record = client.fetch_one(req)
```
Конкретный клиент не добавляет методов или логики: он просто фиксирует конфигурацию и наследует реализацию из базового класса.

## Запрещённая логика в клиентах
- **Доменные трансформации/нормализация** (маппинг полей в модели, дедупликации, проверка бизнес-правил): переносить в слой трансформаций/доменные сервисы (например, модуль обработки ответов в ETL пайплайне).
- **Агрегирование/обогащение** (джоины, дополнительные запросы к другим источникам, вычисление метрик): выполнять в отдельных шагах пайплайна или сервисах объединения данных.
- **Валидация содержимого** (проверка обязательных полей, схем Pandera, фильтрация по бизнес-правилам): запускать после извлечения в слое валидации/качества данных.
- **Управление ретраями/логированием/трассировкой**: реализуются в общих транспортных компонентах; в клиентах только прокидывается `RequestContext`.
- **Слои совместимости и алиасы методов**: депрецированные методы и адаптеры к старым интерфейсам удаляются, новые клиенты реализуют только единый контракт.

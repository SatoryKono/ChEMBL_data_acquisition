# Архитектурный обзор и целевая модель подсистемы клиентов

## Сводка
- Проведен анализ текущей структуры `src/bioetl/clients/`, включающей базовые абстракции транспорта/провайдеров, фабрики под конкретные источники и слой нормализаторов.
- Текущая реализация смешивает ответственность: часть бизнес-логики и нормализации размещена прямо в клиентах/нормализаторах, конфигурация зашита в коде, нет единого контракта и единой конфигурационной модели.
- Предложена целевая архитектура с едиными ABC/протоколами клиентов данных, унифицированными абстракциями запросов/ответов и вынесением всех параметров источников в YAML-конфигурации.
- В целевой модели клиенты становятся тонкими обертками над инфраструктурными HTTP/DB-компонентами; бизнес-логика, трансформации и валидации выносятся в вышележащие слои (пайплайны/доменные сервисы).
- Введены общие правила: запрет логики в клиентах, обязательное следование единым протоколам, обязательное использование YAML-конфигураций и фабрик без слоев совместимости со старым кодом.

## Текущее состояние `clients/`
```
src/bioetl/clients/
  __init__.py
  enricher_base.py                 # адаптер BaseApiClient -> BaseHttpClientABC
  base/
    http.py                        # ABC для HTTP-транспорта
    interfaces.py                  # протоколы BaseApiClient/DataProviderProtocol/контекст
    pagination.py                  # стратегия пагинации с делегацией транспорту
    normalizers.py                 # INormalizer и IdentityNormalizerImpl
  providers/
    base_provider.py               # BaseDataProviderABC/PagedDataProviderABC с обертками ошибок
    _provider_template.py          # шаблон ProviderDataClientImpl + конфиг маршрутов/нормализация
    pubchem/, pubmed/, crossref/, openalex/, semantic_scholar/, uniprot/, iuphar/
      ..._data_client_impl.py      # реализации провайдеров
      ..._normalizer_impl.py       # нормализаторы (часто бизнес-логика)
      ..._pagination_strategy...py # частные стратегии
    __init__.py
  factories/
    client_factory.py              # протокол фабрик
    helpers.py                     # сборка транспорта/нормализатора из конфигурации
    *client_factory_impl.py        # фабрики под каждый источник (PubChem, PubMed, Crossref, ...)
  target/
    factories.py                   # фабрика Target IUPHAR клиента
```
Ключевые особенности:
- Клиенты реализованы как «провайдеры» с нормализаторами внутри; параметры запросов/роутов и поля ответа зашиты в Python-код.
- Разные фабрики имеют разрозненные контракты и проверки типов; отсутствует единый протокол клиента.
- Нормализаторы и стратегии пагинации реализуют часть бизнес-логики (напр. CrossrefPaginationStrategyImpl, ProviderNormalizer) вместо выделенных доменных слоев.

## Целевая архитектура клиентов
```
src/bioetl/clients/
  base/
    __init__.py
    client.py              # ABC/Protocol: ExternalDataClient (fetch_one/fetch_many/iter_pages/close/metadata)
    transport.py           # HTTP/DB транспортные ABC + retry/rate-limit hooks
    request.py             # абстракции запроса/маршрута (endpoint/path, query, headers, auth)
    response.py            # описатели формата ответа (JSON/XML/CSV), схемы полей
    pagination.py          # общие стратегии пагинации (cursor/page/offset) на базе транспорта
  config/
    __init__.py
    models.py              # pydantic-модели для загрузки YAML конфигов источников (endpoint, auth, поля)
    loader.py              # чтение YAML per-source, валидация, расширение env
  factory.py              # единая фабрика ExternalDataClient из конфигурационных моделей
  chembl/
    client.py             # тонкая обертка: выбор маршрута/endpoint из YAML, вызов базовых абстракций
  pubchem/
    client.py
  pubmed/
    client.py
  openalex/
    client.py
  crossref/
    client.py
  semanticscholar/
    client.py
  uniprot/
    client.py
  ... (прочие источники)
configs/clients/
  chembl.yaml
  pubchem.yaml
  pubmed.yaml
  ...
```
Общие абстракции:
- **ExternalDataClient (Protocol/ABC)**: `fetch_one(identifier, *, params, context)`, `fetch_many(*, params, page_size, context)`, `iter_pages(*, params, page_size, context)`, `metadata()` и `close()`. Возвращаемые записи — сырые структуры (dict/Mapping) без доменных трансформаций.
- **Transport**: `HttpTransport`/`DbTransport` интерфейсы с методами `request`, `get_json`, `paginate`, опционально контекстом для логирования/трассировки.
- **RequestDescriptor**: структура для маршрута (`path/template`), параметров (`query/body`), заголовков, опций аутентификации; связывается с конфигурацией YAML.
- **ResponseDescriptor**: описание формата (`json/xml/text`), схемы извлекаемых полей (ключи/колонки), правил пагинации (`page_key/next_key/page_param`).
- **PaginationStrategy**: унифицированные стратегии (transport-driven cursor/page/offset) подключаемые через конфигурацию.
Конфигурации:
- YAML per-source содержит: базовый endpoint/SQL, параметры аутентификации, схему полей, правила пагинации, ограничения rate-limit/retry, список маршрутов (fetch/search/list) и форматы ответа.
- Pydantic-модели в `config/models.py` валидируют YAML, мапят его на Request/Response/Pagination дескрипторы.
Фабрика:
- `clients/factory.py` принимает pydantic-конфигурации и возвращает конкретный `ExternalDataClient`, создавая транспорт, пагинацию, привязку маршрутов и обертку над общими абстракциями.
Размещение логики:
- Клиенты содержат только выбор маршрута/параметров и вызов транспорта; никакой бизнес-логики, нормализации, валидации доменных правил.
- Вся логика (нормализация, маппинг на доменную модель, валидации Pandera, агрегации) находится в слоях ETL/доменных сервисов, которые используют тонкие клиенты.

## Группы клиентов и соответствия
| Источник | Текущие модули/классы | Целевой модуль/класс | Особенности, которые нужно учесть |
|----------|------------------------|----------------------|------------------------------------|
| PubChem | `providers/pubchem/pubchem_data_client_impl.py`, `pubchem_normalizer_impl.py`, фабрика `pubchem_client_factory_impl.py` | `clients/pubchem/client.py` через общую фабрику | Вынести маршруты поиска/деталей и правила пагинации в YAML; убрать нормализацию из клиента. |
| PubMed | `providers/pubmed/pubmed_data_client_impl.py`, `pubmed_normalizer_impl.py`, фабрика `pubmed_client_factory_impl.py` | `clients/pubmed/client.py` | Описать EFetch/ESearch параметры и формат XML/JSON в YAML; клиент только дергает транспорт. |
| Crossref | `providers/crossref/*`, фабрика `crossref_client_factory_impl.py` | `clients/crossref/client.py` | Унифицировать пагинацию cursor/page через общую стратегию; вынести `items_key`, `next-cursor` в YAML. |
| OpenAlex | `providers/openalex/*`, фабрика `openalex_client_factory_impl.py` | `clients/openalex/client.py` | Маршруты (works/authors/venues) и поля в YAML; нормализация вовне. |
| Semantic Scholar | `providers/semantic_scholar/*`, фабрика `semantic_scholar_client_factory_impl.py` | `clients/semanticscholar/client.py` | Описать поля ответа и лимиты API в YAML; единый клиент-протокол. |
| UniProt | `providers/uniprot/*`, фабрика `uniprot_client_factory_impl.py` | `clients/uniprot/client.py` | Пагинация cursor/`results` ключи в YAML; убрать кастомный configure. |
| IUPHAR Target | `providers/iuphar/impl/target_client_impl.py`, `target/factories.py` | `clients/iuphar/client.py` | Маршруты/поля в YAML; переход на общий ExternalDataClient без спец. фабрики. |
| ChEMBL | (отсутствует в текущем пакете) | `clients/chembl/client.py` | Добавить YAML с SQL/REST маршрутами; клиент только вызывает транспорт. |
| Прочие (Crossref/OpenAlex/etc.) | см. выше | одноименные подкаталоги | Все следуют единому протоколу и YAML-конфигурации. |

## Правила для будущих изменений
1. **Никакой бизнес-логики внутри клиентов**: только вызовы транспорта, параметры запроса/маршрутов и минимальное приведение типов для корректного HTTP/DB-взаимодействия.
2. **Единый протокол клиента**: все клиенты реализуют `ExternalDataClient` (или эквивалентный ABC) с одинаковыми сигнатурами `fetch_one/fetch_many/iter_pages/metadata/close` и едиными форматами возвращаемых данных (сырые `Mapping`).
3. **Конфигурации только в YAML**: все параметры источников (endpoint/SQL, маршруты, схема полей, формат ответа, пагинация, аутентификация, rate-limit/retry) описываются в отдельных YAML-файлах; Python-код читает их через pydantic-модели из `clients/config`.
4. **Тонкие клиенты**: реализации в `clients/<source>/client.py` не содержат нормализаторов, агрегаций, доменной валидации или маппинга на модели; любая такая логика переносится в ETL/доменные сервисы.
5. **Обязательное использование общих абстракций**: транспорт, запрос/ответ, пагинация берутся из `clients/base`; запрещены кастомные стратегии/классы вне общего набора без их формализации в `base/`.
6. **Без слоев совместимости**: удаляем устаревшие фабрики/клиенты, не добавляем адаптеры к старым интерфейсам; миграция выполняется с заменой файлов на новую архитектуру.
7. **Единые фабрики**: создание клиентов происходит через `clients/factory.py`, принимающую конфигурационные модели; точечные фабрики `*client_factory_impl.py` не добавляются.
8. **Минимальные проверки**: внутри клиентов допустимы только проверки корректности конфигурации транспорта/маршрутов; бизнес-валидации и обогащение данных запрещены.

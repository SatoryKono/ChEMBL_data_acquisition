# Руководство по конфигурации

## Обзор

* Все CLI-инструменты читают настройки из [`config/config.yaml`](../config/config.yaml) в корне проекта.
* Значения валидируются в `library.config.load_config`, который вызывает `Config.model_validate` из Pydantic. [`config.schema.json`](../config.schema.json) служит справочной схемой для инструментов и не исполняется при запуске.
* Переопределения применяются в порядке: `config/config.yaml` < переменные окружения < аргументы командной строки.

## Структура `config/config.yaml`

| Раздел | Назначение |
| --- | --- |
| `sources` | Подключения ко внешним сервисам (ChEMBL, UniProt, CrossRef, PubChem, PubMed и др.). |
| `local` | Пути к локальным ресурсам, настройки CSV и рабочие книги инициализации. |
| `system` | Логирование, политика повторов, глобальные лимиты и веса классификации документов. |

Чувствительные данные (токены, персональные e-mail) задавайте через переменные окружения, а не в файле.

## `sources.chembl`

### Клиент API (`sources.chembl.api`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `chembl_base` | `https://www.ebi.ac.uk/chembl/api/data` | Базовый URL ChEMBL REST API. |
| `timeout_connect` | `5` | Таймаут установки соединения (сек.). |
| `timeout_read` | `30` | Таймаут ожидания ответа (сек.). |
| `retries` | `3` | Максимум попыток, выполняемых клиентскими обёртками; общий HTTP-адаптер повторы не делает. |
| `backoff_factor` | `0.5` | Множитель экспоненциального backoff между повторами. |
| `rps` | `20` | Лимит запросов в секунду для rate limiter. |
| `burst` | `20` | Размер «ведра» токенов. |
| `user_agent` | `chembl-da/0.1 (mailto:contact@example.org)` | Заголовок User-Agent; замените e-mail на рабочий. |

### Кэш ответов (`sources.chembl.cache`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `cache_ttl` | `3600` | Время жизни кэша API-ответов (сек.). |
| `cache_maxsize` | `1024` | Максимальное число кэшированных ответов. |

<a id="sources-chembl-molecule-catalog"></a>
### Каталог молекул (`sources.chembl.molecule_catalog`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `cache_path` | `data/cache/molecule_parent_catalog.json` | Путь к локальному JSON с отношениями родитель→потомок, который переиспользуется конвейерами; переопределяется через `CHEMBL_DA_MOLECULE_CATALOG_CACHE` (алиас `CHEMBL_DA__SOURCES__CHEMBL__MOLECULE_CATALOG__CACHE_PATH`). |
| `sqlite_path` | `data/cache/molecule_parent_catalog.sqlite` | Путь к SQLite-кэшу для быстрых запросов по связям; используйте `CHEMBL_DA_SOURCES_CHEMBL_MOLECULE_CATALOG_SQLITE_PATH` или каноничную форму `CHEMBL_DA__SOURCES__CHEMBL__MOLECULE_CATALOG__SQLITE_PATH`. |
| `endpoint` | `molecule` | Ресурс REST API ChEMBL, из которого подкачиваются данные при обновлении кэша. |
| `child_field` | `molecule_chembl_id` | Поле ответа API с идентификатором дочерней молекулы. |
| `parent_field` | `parent_molecule_chembl_id` | Поле ответа API с идентификатором родительской молекулы. |
| `force_refresh_existing` | `false` | При `true` пересобирает связи родитель→потомок даже для записей с уже заполненным родителем, заставляя использовать данные из кэша/ChEMBL. |
| `fields` | `['molecule_chembl_id', 'parent_molecule_chembl_id']` | Список полей, которые запрашиваются у ChEMBL при построении или обновлении каталога; расширяйте его для дополнительных атрибутов. |
| `filters` | `{'parent_molecule_chembl_id__isnull': 'false'}` | Набор фильтров, добавляемых ко всем запросам; по умолчанию выбирает только записи с заполненным родителем в ChEMBL. |
| `hierarchy_lookup_path` | `dictionary/_testitem/molecule_hierarchy.csv` | Необязательный CSV с готовыми связями родитель→потомок, который используется офлайн до обращения к ChEMBL; переопределяйте при распространении собственной витрины или переносе каталога. |
| `hierarchy_lookup_encoding` | `utf-8-sig` | Кодировка, применяемая при чтении CSV иерархии; меняйте, если файл сохранён в другом наборе символов (например, Latin-1 из старых выгрузок). |
| `hierarchy_lookup_delimiter` | `,` | Разделитель столбцов, ожидаемый загрузчиком иерархии; укажите `;` или табуляцию для файлов, подготовленных региональными командами. |
| `page_size` | `500` | Количество записей в одном запросе при перепостроении каталога. |
| `fallback_single_limit` | `null` | Ограничение на число одиночных fallback-запросов после неудачных батчей; `null` оставляет fallback без лимита. |

### Пайплайны (`sources.chembl.pipelines`)

Каждый подпункт задаёт значения по умолчанию для одноимённого CLI-скрипта. Изменения аргументов CLI синхронизируются с конфигурацией перед запуском.

#### Пайплайн активностей (`activity`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `column` | `activity_chembl_id` | Колонка с идентификаторами активностей. |
| `batch_size` | `50` | Размер батча запросов. |
| `timeout` | `30.0` | Таймаут запроса (сек.). |
| `limit` | `null` | Ограничение на число идентификаторов. |
| `dry_run` | `false` | Пропуск сетевых вызовов и записи файлов. |

#### Границы активности (`activity_bounds`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `enable_from_relation` | `true` | Использовать `standard_value` и оператор отношения, если явных границ нет. |
| `enable_from_uncertainty` | `false` | Разбирать выражения `значение ± дельта` из `standard_text_value`. |
| `rounding_digits` | `3` | Количество знаков после запятой при округлении итоговых границ. |
| `clamp_nonnegative` | `true` | Обрезать отрицательные значения для метрик с неотрицательной областью. |
| `log_unknown_relations` | `true` | Логировать предупреждение для нераспознанных отношений. |

#### Обогащение активностей (`activity_enrichment`)

Блок отвечает за дополнительные атрибуты, которые записываются в каждую строку активностей. Он работает совместно с
отдельной конфигурацией `activity_bounds` (описана ниже), отвечающей только за расчёт числовых границ.
`activity_enrichment` состоит из двух независимых подсекций, каждая включается собственным флагом.

##### Маркировка типа действия (`activity_enrichment.action_type`)

* `enabled` — главный переключатель (по умолчанию `true`).
* `column` — название выходной колонки с итоговой меткой (`action_type`).
* Флаги логирования помогают отслеживать диагностические случаи:
  * `log_missing` выводит предупреждение, если метка не определена (`true`).
  * `log_distribution` печатает итоговую статистику распределения (`true`).
* `metrics` сопоставляет типы измерений (IC50, EC50 и т.д.) с базовыми метками (`inhibition`, `activation`, `binding`). Список
  расширяется при необходимости.
* `triages`, `functionality`, `mechanism` — словари переопределений для текстовых совпадений. По умолчанию заполнены только
  типичные функциональные роли (agonist, antagonist и пр.).
* Списки `triage_fields`, `functionality_fields`, `mechanism_fields` задают, какие исходные колонки просматриваются перед применением
  значений по метрикам.
* `allowlist` ограничивает допустимые метки; значения вне списка заменяются на `fallback` после логирования отклонения.
* `positive_label` и `negative_label` задают читаемые ярлыки для положительных/отрицательных модуляторов (`PAM`/`NAM`), а
  `fallback` равен `unknown`.

##### Сводка свойств активности (`activity_enrichment.activity_properties`)

* `enabled` — флаг включения (`true`).
* `column` — исходная колонка со структурированными свойствами (`activity_properties`).
* `summary_column` — зарезервирован под будущий текстовый рендер. Сейчас конвейер сохраняет JSON в `column`,
  вычисляет только детерминированный отпечаток в `hash_column` и не создаёт отдельное поле с резюме.
* `name_field`, `value_field`, `units_field` задают имена ключей внутри записей (`type`, `value`, `units`).
* `separator` и `pair_separator` — устаревшие параметры форматирования, сохранённые для совместимости; текущая
  сериализация JSON их не использует.
* `drop_source_column` удаляет исходную структурированную колонку после агрегации (`true`).
* Флаги логирования по умолчанию выключены (`log_missing=false`, `log_distribution=false`).
* `allowlist` ограничивает перечень сохраняемых групп (measurement, assay, comments, effect_features, triage, mechanism,
  functionality).
* `hash_column` хранит детерминированный отпечаток итоговых свойств (`properties_hash`) для отслеживания изменений далее по цепочке.

##### Границы активности (`activity_bounds`)

Пайплайн активностей дополняет выгрузку нормализованными границами с помощью `compute_activity_bounds` из
`library.processing.activity`. Настройки собраны в отдельном блоке `activity_bounds` (вне `activity_enrichment`) и управляют
последовательностью детерминированных шагов, которые выполняются для каждой строки в следующем порядке：

1. Использовать готовые `standard_lower_value`/`standard_upper_value`, если обе границы заданы.
2. Скомбинировать `standard_value` с противоположной границей и заполнить пропущенное значение.
3. Проанализировать `standard_relation`, если `enable_from_relation = true`, сопоставив операторы (`=`, `≈`, `>=`, `<=`, `between`, `range`) с подходящими границами.
4. Разобрать выражения `значение ± дельта` из `standard_text_value`, когда включён `enable_from_uncertainty`.

Каждый шаг фиксирует найденные значения; отключение этапа просто пропускает его, не затрагивая предыдущие результаты.

Параметры удобно настраивать через YAML или переменные окружения:

```yaml
activity_bounds:
  enable_from_relation: false
  enable_from_uncertainty: true
  rounding_digits: 2
  clamp_nonnegative: true
```

```bash
export CHEMBL_DA__ACTIVITY_BOUNDS__ENABLE_FROM_RELATION=false
export CHEMBL_DA__ACTIVITY_BOUNDS__ROUNDING_DIGITS=2
```

CLI-параметры имеют приоритет над YAML и окружением только для ключей, которые явно прокинуты в парсер (колонка, размер батча, лимит, `--dry-run`). Остальные изменения выполняются через файл настроек или переменные окружения `CHEMBL_DA__ACTIVITY_BOUNDS__*`.

#### Пайплайн ассайев (`assay`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `column` | `assay_chembl_id` | Колонка с идентификаторами ассайев. |
| `batch_size` | `50` | Размер батча запросов. |
| `timeout` | `30.0` | Таймаут запроса (сек.). |
| `limit` | `null` | Ограничение на число идентификаторов. |

#### Пайплайн тест-объектов (`testitem`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `column` | `molecule_chembl_id` | Колонка с идентификаторами соединений. |
| `batch_size` | `1000` | Размер батча запросов. |
| `timeout` | `30.0` | Таймаут запроса (сек.). |
| `limit` | `null` | Ограничение на число идентификаторов. |
| `offset` | `0` | Смещение перед стартом выгрузки, полезно для повторных запусков. |
| `request_limit` | `1000` | Максимальное число объектов в одном запросе API. |
| `retries` | `5` | Количество повторов при временных ошибках ChEMBL. |
| `backoff_factor` | `0.5` | Множитель экспоненциальной задержки между повторами. |
| `fields` | `["molecule_chembl_id", "parent_molecule_chembl_id", "pref_name", "max_phase", "molecule_type", "first_approval", "oral", "parenteral", "topical", "black_box_warning", "structure_type", "molecule_structures.canonical_smiles", "molecule_structures.standard_inchi", "molecule_structures.standard_inchi_key"]` | Список полей, запрашиваемых у ChEMBL. |

Параметры `offset` и `request_limit` помогают управлять пагинацией API: первый задаёт стартовую позицию, второй ограничивает размер одной страницы, который по умолчанию равен максимально допустимым 1000 объектам. Цикл повторов контролируется парами `retries`/`backoff_factor`, позволяя сгладить временные ошибки сервиса. Список `fields` определяет точный набор колонок в ответе и по умолчанию отражает рекомендуемый минимальный профиль тест-объекта.

#### Обогащение тест-объектов (`testitem_molecule_enrichment`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `enable` | `true` | Включает стадию расчёта `salt_chembl_id` и флагов из каталога молекул. |
| `sources.molecule_catalog_path` | `dictionary/_testitem/molecule_catalog.csv` | CSV со столбцами `molecule_chembl_id`, `natural_product`, `prodrug`, `polymer_flag`. |
| `sources.molecule_hierarchy_path` | `dictionary/_testitem/molecule_hierarchy.csv` | CSV с соответствиями дочерней и родительской молекулы. |
| `output.salt_as_null_when_absent` | `true` | При `true` несолевые соединения дают `null`, при `false` — символ `-`. |
| `flags.coerce_to_bool` | `true` | Нормализует значения вида `Y/N`, `1/0`, `yes/no` в булев тип pandas. |
| `flags.parent_fallback` | `true` | Подтягивает флаги из родителя, если у дочерней записи они пусты. |
| `logging.warn_missing_molecule` | `true` | Логирует предупреждение, если молекулы нет в иерархии или каталоге. |
| `logging.warn_inconsistent_flags` | `true` | Сообщает о расхождении флагов между дочерней и родительской записью. |

#### Пайплайн документов (`document`)

| Подсекция | Ключ | Значение | Описание |
| --- | --- | --- | --- |
| `pubmed` | `column` | `PMID` | Колонка с PubMed ID. |
|  | `sleep` | `5.0` | Пауза между циклами опроса (сек.). |
|  | `workers` | `1` | Количество потоков. |
|  | `batch_size` | `5` | Число ID в батче. |
|  | `limit` | `null` | Ограничение на число идентификаторов. |
| `chembl` | `column` | `document_chembl_id` | Колонка с идентификаторами документов ChEMBL. |
|  | `chunk_size` | `50` | Размер батча запросов. |
|  | `timeout` | `30.0` | Таймаут запроса (сек.). |
|  | `limit` | `null` | Ограничение на число идентификаторов. |
| `all` | `column` | `document_chembl_id` | Колонка для объединённого пайплайна. |
|  | `chunk_size` | `5` | Размер выборки ChEMBL при вызове `run_all` → `cl.get_documents`. |
|  | `sleep` | `5.0` | Пауза между циклами опроса PubMed при обогащении. |
|  | `workers` | `1` | Потоки, координирующие загрузку из ChEMBL и обогащение PubMed. |
|  | `batch_size` | `5` | Размер запроса PubMed, передаваемый в `fetch_pubmed_records`. |
|  | `timeout` | `30.0` | Таймаут, применяемый к вызовам ChEMBL и PubMed. |
|  | `limit` | `null` | Ограничение на число идентификаторов в объединённом запуске. |

#### Пайплайн таргетов (`target`)

| Подсекция | Ключ | Значение | Описание |
| --- | --- | --- | --- |
| `uniprot` | `column` | `uniprot_id` | Колонка с UniProt ID. |
|  | `data_dir` | `dictionary/_target/_uniprot` | Каталог с кэшированными JSON UniProt. |
|  | `limit` | `null` | Ограничение на число идентификаторов. |
| `chembl` | `column` | `target_chembl_id` | Колонка с таргетами ChEMBL. |
|  | `chunk_size` | `5` | Размер батча запросов. |
|  | `timeout` | `30.0` | Таймаут запроса (сек.). |
|  | `limit` | `null` | Ограничение на число идентификаторов. |
| `iuphar` | `target_csv` | `dictionary/_target/_IUPHAR/_IUPHAR_target.csv` | Справочник таргетов IUPHAR. |
|  | `family_csv` | `dictionary/_target/_IUPHAR/_IUPHAR_family.csv` | Справочник семейств IUPHAR. |
|  | `limit` | `null` | Ограничение на число идентификаторов. |
| `all` | `data_dir` | `dictionary/_target/_uniprot` | Каталог с данными UniProt. |
|  | `target_csv` | `dictionary/_target/_IUPHAR/_IUPHAR_target.csv` | Таблица таргетов IUPHAR. |
|  | `family_csv` | `dictionary/_target/_IUPHAR/_IUPHAR_family.csv` | Таблица семейств IUPHAR. |
|  | `chunk_size` | `5` | Размер батча при объединении источников. |
|  | `timeout` | `30.0` | Таймаут запроса (сек.). |

|  | `uniprot_column` | `uniprot_id` | Колонка для соединения с UniProt. |
|  | `chembl_out` | `null` | Индивидуальный путь для объединённых данных ChEMBL. |
|  | `uniprot_out` | `null` | Индивидуальный путь для объединённых данных UniProt. |
|  | `iuphar_out` | `null` | Индивидуальный путь для объединённых данных IUPHAR. |
|  | `limit` | `null` | Ограничение на число идентификаторов. |

> При финализации `finalise_targets` вызывает встроенный классификатор таксономии: он использует поля UniProt с родословной (`genus`, `lineage_*`, `taxon_id`, `species_group_flag`) для расчёта итоговой колонки `type` и вспомогательных флагов, поэтому отдельный CSV для `organism` больше не нужен.

## Другие источники (`sources.*`)

| Раздел | Базовый URL | Ключевые параметры |
| --- | --- | --- |
| `openalex` | `https://api.openalex.org` | `timeout_connect=5`, `timeout_read=30`, `retries=3`, `rps=5`, `burst=5`, `mailto="contact@example.org"`. |
| `crossref` | `https://api.crossref.org` | Та же структура, важно указать персональный `mailto`. |
| `uniprot.api` | `https://rest.uniprot.org` | `timeout_connect=5`, `timeout_read=30`, `rps=25`, `burst=25`, `delay=0.25`. |
| `uniprot.mapping` | `https://rest.uniprot.org/idmapping` | `poll_interval=0.5`, `timeout=300.0`, `cache_ttl=null`. |
| `iuphar` | `https://www.guidetopharmacology.org/services` | `timeout_connect=5`, `timeout_read=30`, `rps=5`, `burst=5`. |
| `pubchem` | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | См. [детальные настройки PubChem](#pubchem-sourcespubchem). |
| `pubmed` | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils` | `timeout_connect=5`, `timeout_read=10`, `retries=2`, опциональные `rps`/`burst` для документарных запросов. |
| `semantic_scholar` | `https://api.semanticscholar.org/graph/v1` | `timeout_connect=5`, `timeout_read=10`, `retries=2`, опциональные `rps`/`burst` для документарных запросов. |

Соблюдайте требования сервисов по rate limit и указанию контактной информации.

### PubChem (`sources.pubchem`)

Параметры PubChem применяются прежде всего в пайплайне `testitem`. Каждый ключ отражает содержимое [`config/config.yaml`](../config/config.yaml)
и может быть переопределён переменными окружения (см. [Переменные окружения](#переменные-окружения)). В таблице указаны и
автогенерируемые алиасы `CHEMBL_DA_SOURCES_PUBCHEM_*`, и универсальный формат `CHEMBL_DA__SOURCES__PUBCHEM__*`. Для базового URL
доступен короткий алиас из таблицы в разделе «Переменные окружения».

| Ключ | Значение по умолчанию | Описание | Переменные окружения |
| --- | --- | --- | --- |
| `enable` | `true` | Главный переключатель обогащения данными PubChem. | `CHEMBL_DA_SOURCES_PUBCHEM_ENABLE`, `CHEMBL_DA__SOURCES__PUBCHEM__ENABLE` |
| `base` | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | REST‑эндпоинт PUG, к которому выполняются запросы. | `CHEMBL_DA_PUBCHEM_BASE`, `CHEMBL_DA_SOURCES_PUBCHEM_BASE`, `CHEMBL_DA__SOURCES__PUBCHEM__BASE` |
| `timeout_connect` | `5` | Таймаут установления соединения (сек.). | `CHEMBL_DA_SOURCES_PUBCHEM_TIMEOUT_CONNECT`, `CHEMBL_DA__SOURCES__PUBCHEM__TIMEOUT_CONNECT` |
| `timeout_read` | `60` | Таймаут ожидания ответа (сек.). | `CHEMBL_DA_SOURCES_PUBCHEM_TIMEOUT_READ`, `CHEMBL_DA__SOURCES__PUBCHEM__TIMEOUT_READ` |
| `timeout_seconds` | `30.0` | Максимальная длительность одной попытки разрешения CID с учётом повторов. | `CHEMBL_DA_SOURCES_PUBCHEM_TIMEOUT_SECONDS`, `CHEMBL_DA__SOURCES__PUBCHEM__TIMEOUT_SECONDS` |
| `retries` | `3` | Число попыток, которое выполнит цикл повторов PubChem, прежде чем сдаться. | `CHEMBL_DA_SOURCES_PUBCHEM_RETRIES`, `CHEMBL_DA__SOURCES__PUBCHEM__RETRIES` |
| `rps` | `5` | Локальный лимит запросов в секунду. | `CHEMBL_DA_SOURCES_PUBCHEM_RPS`, `CHEMBL_DA__SOURCES__PUBCHEM__RPS` |
| `burst` | `5` | Ёмкость токен-бакета, связанного с `rps`. | `CHEMBL_DA_SOURCES_PUBCHEM_BURST`, `CHEMBL_DA__SOURCES__PUBCHEM__BURST` |
| `delay` | `0.2` | Пауза между повторами (сек.). | `CHEMBL_DA_SOURCES_PUBCHEM_DELAY`, `CHEMBL_DA__SOURCES__PUBCHEM__DELAY` |
| `backoff_initial_seconds` | `0.5` | Начальная задержка экспоненциального backoff после 429/5xx. | `CHEMBL_DA_SOURCES_PUBCHEM_BACKOFF_INITIAL_SECONDS`, `CHEMBL_DA__SOURCES__PUBCHEM__BACKOFF_INITIAL_SECONDS` |
| `resolve_order` | `cache → smiles → inchikey → inchi → pref_name` | Очерёдность стратегий при поиске PubChem CID. | `CHEMBL_DA_SOURCES_PUBCHEM_RESOLVE_ORDER`, `CHEMBL_DA__SOURCES__PUBCHEM__RESOLVE_ORDER` |
| `cache_ttl` | `3600` | Время жизни in-memory кэша HTTP (сек.). | `CHEMBL_DA_SOURCES_PUBCHEM_CACHE_TTL`, `CHEMBL_DA__SOURCES__PUBCHEM__CACHE_TTL` |
| `cache_ttl_hours` | `null` | TTL (часы) для постоянного CID-кэша; `null` отключает истечение. | `CHEMBL_DA_SOURCES_PUBCHEM_CACHE_TTL_HOURS`, `CHEMBL_DA__SOURCES__PUBCHEM__CACHE_TTL_HOURS` |
| `cid_cache_path` | `"data/cache/pubchem_cid_cache.json"` | Путь к JSON с сохранёнными CID для повторного использования. | `CHEMBL_DA_SOURCES_PUBCHEM_CID_CACHE_PATH`, `CHEMBL_DA__SOURCES__PUBCHEM__CID_CACHE_PATH` |
| `batch_size` | `50` | Размер батча для обработчика PubChem. | `CHEMBL_DA_SOURCES_PUBCHEM_BATCH_SIZE`, `CHEMBL_DA__SOURCES__PUBCHEM__BATCH_SIZE` |
| `prefer_local_smiles` | `false` | Пропускать запросы, если локальные SMILES/InChIKey уже заполнены. | `CHEMBL_DA_SOURCES_PUBCHEM_PREFER_LOCAL_SMILES`, `CHEMBL_DA__SOURCES__PUBCHEM__PREFER_LOCAL_SMILES` |
| `prefer_local_values` | `true` | Сохранять существующие колонки `pubchem_*`, если ответ пуст. | `CHEMBL_DA_SOURCES_PUBCHEM_PREFER_LOCAL_VALUES`, `CHEMBL_DA__SOURCES__PUBCHEM__PREFER_LOCAL_VALUES` |
| `use_parent_for_salts` | `true` | Переходить к родительским структурам, если соль не найдена. | `CHEMBL_DA_SOURCES_PUBCHEM_USE_PARENT_FOR_SALTS`, `CHEMBL_DA__SOURCES__PUBCHEM__USE_PARENT_FOR_SALTS` |
| `allow_polymer` | `false` | Разрешать обработку полимеров и смесей. | `CHEMBL_DA_SOURCES_PUBCHEM_ALLOW_POLYMER`, `CHEMBL_DA__SOURCES__PUBCHEM__ALLOW_POLYMER` |
| `write_not_found_literal` | `false` | Записывать литерал `Not Found`, если CID не найден. | `CHEMBL_DA_SOURCES_PUBCHEM_WRITE_NOT_FOUND_LITERAL`, `CHEMBL_DA__SOURCES__PUBCHEM__WRITE_NOT_FOUND_LITERAL` |

> Совет: поле `resolve_order` принимает любую последовательность поддерживаемых стратегий; чаще всего полезно оставлять `cache`
> первым, чтобы использовать прогретый кэш CID.

## Локальные ресурсы (`local`)

### Справочники (`local.resources`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `dictionary_dir` | `dictionary` | Корневая папка словарей. |
| `iuphar_target_csv` | `dictionary/_target/_IUPHAR/_IUPHAR_target.csv` | Соответствия таргетов IUPHAR. |
| `iuphar_family_csv` | `dictionary/_target/_IUPHAR/_IUPHAR_family.csv` | Справочник семейств IUPHAR. |
| `uniprot_data_dir` | `dictionary/_target/_uniprot` | Кэшированные ответы UniProt. |
| `targets_type_csv` | `dictionary/_Target/targets_type.csv` | Классификация типов таргетов. |


Каталог `dictionary/_target` отражает текущую структуру репозитория; в нём по умолчанию лежат справочники IUPHAR и выгрузки UniProt.


> Таксономическая классификация таргетов вычисляется по данным UniProt в коде пайплайна; отдельный CSV-справочник организмов больше не требуется.

### Настройки ввода-вывода (`local.io`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `output_dir` | `data/output` | Каталог для результирующих наборов данных. |
| `cache_dir` | `.cache` | Каталог HTTP-кэша. |
| `csv_sep` | `,` | Разделитель CSV по умолчанию. |
| `csv_encoding` | `utf-8-sig` | Кодировка экспорта CSV. |
| `csv_chunksize` | `10000` | Размер чанка (строк) при детерминированной записи CSV; значение задано в [`config/config.yaml`](../config/config.yaml). |
| `na_markers` | `["#N/A"]` | Дополнительные маркеры пропусков при чтении CSV. |
| `keep_na_markers` | `false` | Сохранять идентификаторы, совпадающие с `na_markers`, вместо их фильтрации. |
| `exist_ok` | `true` | Создавать каталоги автоматически. |

### Инициализационные книги (`local.init`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `same_doc` | `data/input/ChEMBL/ChEMBL_same_document_20_05.xlsx` | Источник пар «тот же документ». |
| `all_doc` | `data/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx` | Источник пар «разные документы». |
| `output_dir` | `data/output/ChEMBL/processed` | Каталог для подготовленных файлов. |

## Системные настройки (`system`)

| Подсекция | Ключ | Значение | Описание |
| --- | --- | --- | --- |
| `log` | `level` | `INFO` | Уровень логирования. |
|  | `format` | `[%(asctime)s] %(levelname)s %(name)s: %(message)s` | Формат сообщений. |
|  | `datefmt` | `%Y-%m-%d %H:%M:%S` | Формат времени. |
| `rate` | `global_rps` | `100` | Глобальный лимит запросов в секунду. |
|  | `global_burst` | `100` | Ёмкость глобального токен-бакета. |
|  | `limiter_cache_maxsize` | `128` | Максимум кэшированных лимитеров. |
|  | `limiter_cache_ttl` | `600` | TTL записей кэша (сек.). |
| `retry` | `max_attempts` | `3` | Количество повторов при ошибках. |
|  | `backoff_factor` | `0.5` | Базовый множитель экспоненциальной задержки. |
|  | `status_forcelist` | `[429, 500, 502, 503, 504]` | Коды HTTP, инициирующие повтор. |
| `doc_type` | `weights` | `{pubmed: 4, openalex: 3, scholar: 2}` | Веса источников документации. |
|  | `thresholds` | `{review: 1, experimental: 1, unknown: 2}` | Минимальные пороги классификации. |
|  | `limit` | `null` | Ограничение на число классифицированных записей. |

## Переменные окружения

Шаблон переменной: `CHEMBL_DA__РАЗДЕЛ__ПОДРАЗДЕЛ__КЛЮЧ`. Пример:

```bash
export CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT="project/1.0 (mailto:me@example.com)"
export CHEMBL_DA__LOCAL__IO__OUTPUT_DIR=/mnt/datasets
```

Часто используемые алиасы:

| Переменная | Соответствие |
| --- | --- |
| `CHEMBL_DA_BASE` | `sources.chembl.api.chembl_base` |
| `CHEMBL_DA_RPS` | `sources.chembl.api.rps` |
| `CHEMBL_DA_BURST` | `sources.chembl.api.burst` |
| `CHEMBL_DA_TIMEOUT_CONNECT` | `sources.chembl.api.timeout_connect` |
| `CHEMBL_DA_TIMEOUT_READ` | `sources.chembl.api.timeout_read` |
| `CHEMBL_DA_CACHE_TTL` | `sources.chembl.cache.cache_ttl` |
| `CHEMBL_DA_CACHE_MAXSIZE` | `sources.chembl.cache.cache_maxsize` |
| `CHEMBL_DA_OUTDIR` | `local.io.output_dir` |
| `CHEMBL_DA_CACHE_DIR` | `local.io.cache_dir` |
| `CHEMBL_DA_GLOBAL_RPS` | `system.rate.global_rps` |
| `CHEMBL_DA_GLOBAL_BURST` | `system.rate.global_burst` |
| `CHEMBL_DA_LIMITER_CACHE_MAXSIZE` | `system.rate.limiter_cache_maxsize` |
| `CHEMBL_DA_LIMITER_CACHE_TTL` | `system.rate.limiter_cache_ttl` |
| `CHEMBL_DA_LOG_LEVEL` | `system.log.level` |
| `CHEMBL_DA_LOG_FORMAT` | `system.log.format` |
| `CHEMBL_DA_LOG_DATEFMT` | `system.log.datefmt` |
| `CHEMBL_DA_RETRY_MAX_ATTEMPTS` | `system.retry.max_attempts` |
| `CHEMBL_DA_RETRY_BACKOFF_FACTOR` | `system.retry.backoff_factor` |
| `CHEMBL_DA_DICT_DIR` | `local.resources.dictionary_dir` |
| `CHEMBL_DA_UNIPROT_DATA_DIR` | `local.resources.uniprot_data_dir` |
| `CHEMBL_DA_IUPHAR_TARGET_CSV` | `local.resources.iuphar_target_csv` |
| `CHEMBL_DA_IUPHAR_FAMILY_CSV` | `local.resources.iuphar_family_csv` |
| `CHEMBL_DA_TARGETS_TYPE_CSV` | `local.resources.targets_type_csv` |
| `CHEMBL_DA_OPENALEX_BASE` | `sources.openalex.base` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_CONNECT` | `sources.openalex.timeout_connect` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_READ` | `sources.openalex.timeout_read` |
| `CHEMBL_DA_OPENALEX_MAILTO` | `sources.openalex.mailto` |
| `CHEMBL_DA_CROSSREF_BASE` | `sources.crossref.base` |
| `CHEMBL_DA_CROSSREF_TIMEOUT_CONNECT` | `sources.crossref.timeout_connect` |
| `CHEMBL_DA_CROSSREF_TIMEOUT_READ` | `sources.crossref.timeout_read` |
| `CHEMBL_DA_CROSSREF_MAILTO` | `sources.crossref.mailto` |
| `CHEMBL_DA_UNIPROT_BASE` | `sources.uniprot.api.base` |
| `CHEMBL_DA_IUPHAR_BASE` | `sources.iuphar.base` |
| `CHEMBL_DA_PUBCHEM_BASE` | `sources.pubchem.base` |

Любой другой ключ можно задать в длинной форме `CHEMBL_DA__...`.

## CLI-переопределения

* Флаг `--config` позволяет указать альтернативный YAML (по умолчанию `config/config.yaml`).
* `--print-config` выводит итоговую конфигурацию (с учётом окружения и CLI) и завершает работу.
* Аргументы, прокинутые через `apply_config_overrides`, обновляют конфигурацию. Например, `--batch-size 25` задаёт `sources.chembl.pipelines.activity.batch_size` на время запуска.
* Вложенные параметры правятся через `config/config.yaml` или переменные окружения (например, `CHEMBL_DA__SOURCES__CHEMBL__API__RPS=10`). Флаги вида `--sources.…` парсер не поддерживает.

## Процесс валидации

1. Загружается `config/config.yaml`.
2. Применяются переменные окружения и CLI-переопределения.
3. Валидируются значения через `Config.model_validate`: неизвестные ключи и ошибки типов приводят к отказу загрузки.
4. Убеждаемся в наличии каталогов `output_dir` и `cache_dir` (создаются автоматически, если `local.io.exist_ok=true`).

При добавлении новых настроек пересоздайте справочный `config.schema.json` и обновите документацию.

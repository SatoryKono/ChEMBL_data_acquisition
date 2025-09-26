# Руководство по конфигурации

## Обзор

* Все CLI-инструменты читают настройки из [`config.yaml`](../config.yaml) в корне проекта.
* Значения валидируются по [`config.schema.json`](../config.schema.json); неизвестные ключи приводят к ошибке запуска.
* Переопределения применяются в порядке: `config.yaml` < переменные окружения < аргументы командной строки.

## Структура `config.yaml`

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
| `retries` | `3` | Количество HTTP-повторов. |
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
| `cache_path` | `data/cache/molecule_parent_catalog.json` | Путь к локальному JSON с отношениями родитель→потомок, который переиспользуется конвейерами. |
| `endpoint` | `molecule` | Ресурс REST API ChEMBL, из которого подкачиваются данные при обновлении кэша. |
| `child_field` | `molecule_chembl_id` | Поле ответа API с идентификатором дочерней молекулы. |
| `parent_field` | `parent_molecule_chembl_id` | Поле ответа API с идентификатором родительской молекулы. |
| `page_size` | `500` | Количество записей в одном запросе при перепостроении каталога. |

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

Пайплайн активностей дополняет выгрузку нормализованными границами с помощью `compute_activity_bounds` в
`scripts/get_activity_data.py`. Настройки собраны в блоке `activity_bounds` и управляют последовательностью детерминированных
шагов, которые выполняются для каждой строки в следующем порядке：【F:scripts/get_activity_data.py†L212-L353】【F:library/config.py†L371-L388】

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

CLI-параметры имеют приоритет над YAML и окружением только для ключей, которые явно прокинуты в парсер (колонка, размер батча, лимит, `--dry-run`). Остальные изменения выполняются через файл настроек или переменные окружения `CHEMBL_DA__ACTIVITY_BOUNDS__*`.【F:scripts/get_activity_data.py†L536-L603】

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
| `batch_size` | `50` | Размер батча запросов. |
| `timeout` | `30.0` | Таймаут запроса (сек.). |
| `limit` | `null` | Ограничение на число идентификаторов. |

#### Обогащение тест-объектов (`testitem_molecule_enrichment`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `enable` | `true` | Включает стадию расчёта `salt_chembl_id` и флагов из каталога молекул. |
| `sources.molecule_catalog_path` | `dictionary/molecule_catalog.csv` | CSV со столбцами `molecule_chembl_id`, `natural_product`, `prodrug`, `polymer_flag`. |
| `sources.molecule_hierarchy_path` | `dictionary/molecule_hierarchy.csv` | CSV с соответствиями дочерней и родительской молекулы. |
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
|  | `data_dir` | `dictionary/uniprot` | Каталог с кэшированными JSON UniProt. |
|  | `limit` | `null` | Ограничение на число идентификаторов. |
| `chembl` | `column` | `target_chembl_id` | Колонка с таргетами ChEMBL. |
|  | `chunk_size` | `5` | Размер батча запросов. |
|  | `timeout` | `30.0` | Таймаут запроса (сек.). |
|  | `limit` | `null` | Ограничение на число идентификаторов. |
| `iuphar` | `target_csv` | `dictionary/_IUPHAR/_IUPHAR_target.csv` | Справочник таргетов IUPHAR. |
|  | `family_csv` | `dictionary/_IUPHAR/_IUPHAR_family.csv` | Справочник семейств IUPHAR. |
|  | `limit` | `null` | Ограничение на число идентификаторов. |
| `all` | `data_dir` | `dictionary/uniprot` | Каталог с данными UniProt. |
|  | `target_csv` | `dictionary/_IUPHAR/_IUPHAR_target.csv` | Таблица таргетов IUPHAR. |
|  | `family_csv` | `dictionary/_IUPHAR/_IUPHAR_family.csv` | Таблица семейств IUPHAR. |
|  | `chunk_size` | `5` | Размер батча при объединении источников. |
|  | `timeout` | `30.0` | Таймаут запроса (сек.). |
|  | `organism_csv` | `dictionary/_Target/organism.csv` | Классификация организмов и типов таргетов. |
|  | `uniprot_column` | `uniprot_id` | Колонка для соединения с UniProt. |
|  | `chembl_out` | `null` | Индивидуальный путь для объединённых данных ChEMBL. |
|  | `uniprot_out` | `null` | Индивидуальный путь для объединённых данных UniProt. |
|  | `iuphar_out` | `null` | Индивидуальный путь для объединённых данных IUPHAR. |
|  | `limit` | `null` | Ограничение на число идентификаторов. |

## Другие источники (`sources.*`)

| Раздел | Базовый URL | Ключевые параметры |
| --- | --- | --- |
| `openalex` | `https://api.openalex.org` | `timeout_connect=5`, `timeout_read=30`, `retries=3`, `rps=5`, `burst=5`, `mailto="contact@example.org"`. |
| `crossref` | `https://api.crossref.org` | Та же структура, важно указать персональный `mailto`. |
| `uniprot.api` | `https://rest.uniprot.org` | `timeout_connect=5`, `timeout_read=30`, `rps=25`, `burst=25`, `delay=0.25`. |
| `uniprot.mapping` | `https://rest.uniprot.org/idmapping` | `poll_interval=0.5`, `timeout=300.0`, `cache_ttl=null`. |
| `iuphar` | `https://www.guidetopharmacology.org/services` | `timeout_connect=5`, `timeout_read=30`, `rps=5`, `burst=5`. |
| `pubchem` | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | `timeout_connect=5`, `timeout_read=60`, `retries=3`, `rps=5`, `burst=5`, `delay=0.2`, `cache_ttl=3600`, `prefer_local_values=true`, `cache_ttl_hours=null`. |
| `pubmed` | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils` | `timeout_connect=5`, `timeout_read=10`, `retries=2`. |
| `semantic_scholar` | `https://api.semanticscholar.org/graph/v1` | `timeout_connect=5`, `timeout_read=10`, `retries=2`. |

Соблюдайте требования сервисов по rate limit и указанию контактной информации.

## Локальные ресурсы (`local`)

### Справочники (`local.resources`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `dictionary_dir` | `dictionary` | Корневая папка словарей. |
| `iuphar_target_csv` | `dictionary/_IUPHAR/_IUPHAR_target.csv` | Соответствия таргетов IUPHAR. |
| `iuphar_family_csv` | `dictionary/_IUPHAR/_IUPHAR_family.csv` | Справочник семейств IUPHAR. |
| `uniprot_data_dir` | `dictionary/uniprot` | Кэшированные ответы UniProt. |
| `targets_type_csv` | `dictionary/_Target/targets_type.csv` | Классификация типов таргетов. |

### Настройки ввода-вывода (`local.io`)

| Ключ | Значение по умолчанию | Описание |
| --- | --- | --- |
| `output_dir` | `data/output` | Каталог для результирующих наборов данных. |
| `cache_dir` | `.cache` | Каталог HTTP-кэша. |
| `csv_sep` | `,` | Разделитель CSV по умолчанию. |
| `csv_encoding` | `utf-8-sig` | Кодировка экспорта CSV. |
| `na_markers` | `["#N/A"]` | Дополнительные маркеры пропусков при чтении CSV. |
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

* Флаг `--config` позволяет указать альтернативный YAML (по умолчанию `config.yaml`).
* `--print-config` выводит итоговую конфигурацию (с учётом окружения и CLI) и завершает работу.
* Аргументы, прокинутые через `apply_config_overrides`, обновляют конфигурацию. Например, `--batch-size 25` задаёт `sources.chembl.pipelines.activity.batch_size` на время запуска.
* Вложенные параметры правятся через `config.yaml` или переменные окружения (например, `CHEMBL_DA__SOURCES__CHEMBL__API__RPS=10`). Флаги вида `--sources.…` парсер не поддерживает.

## Процесс валидации

1. Загружается `config.yaml`.
2. Применяются переменные окружения и CLI-переопределения.
3. Проверяются типы и неизвестные ключи согласно JSON-схеме.
4. Убеждаемся в наличии каталогов `output_dir` и `cache_dir` (создаются автоматически, если `local.io.exist_ok=true`).

При добавлении новых настроек синхронизируйте схему и документацию.

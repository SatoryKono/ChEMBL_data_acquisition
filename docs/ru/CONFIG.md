# Справочник по конфигурации

Конфигурация загружается из трёх слоёв (по возрастанию приоритета):

1. Значения по умолчанию, определённые в [`library/config/models.py`](../../library/config/models.py).
2. YAML-файлы (`config/config.yaml` и необязательный `config/config.local.yaml`).
3. Переменные окружения и аргументы CLI.

Переменные окружения используют схему `CHEMBL_DA__РАЗДЕЛ__КЛЮЧ`, вложенные ключи
соединяются двойным подчёркиванием. Для распространённых параметров доступны
короткие алиасы (см. `_ALIAS_MAP` в `library/config/models.py`).

```bash
# Пример: увеличить лимит UniProt во время smoke-теста
export CHEMBL_DA__SOURCES__UNIPROT__API__RPS=50
export CHEMBL_DA_UNIPROT_RPS=50  # эквивалентный алиас
```

## Каркас YAML

```yaml
sources:
  chembl:
    api:
      chembl_base: https://www.ebi.ac.uk/chembl/api/data
      user_agent: "chembl-da/1.0 (mailto:team@example.org)"
    pipelines:
      document:
        pubmed:
          column: PMID
local:
  io:
    output_dir: "$CHEMBL_DA_BASE_PATH/output"
  resources:
    dictionary_dir: config/dictionary
system:
  log:
    level: INFO
```

Неуказанные поля наследуют значения из Pydantic-моделей. Храните локальные
переопределения в `config/config.local.yaml`, чтобы не коммитить секреты.

## Верхний уровень

| Раздел | Назначение |
|--------|------------|
| `sources` | Подключения к внешним сервисам, лимиты и дефолты конвейеров. |
| `local` | Пути на файловой системе и вспомогательные скрипты. |
| `activity_enrichment` | Правила классификации активностей. |
| `testitem_molecule_enrichment` | Сопоставления родительских молекул. |
| `activity_bounds` | Вывод нижних/верхних границ по значениям и отношениям. |
| `system` | Логирование, ретраи, лимиты и контроль качества. |

Ниже перечислены все ключи с дефолтами и отметкой обязательности (`Да` — требует
явного значения).

## `sources`

### `sources.chembl.api`

| Ключ | Обязателен | Значение по умолчанию | Описание |
|------|------------|-----------------------|----------|
| `chembl_base` | Нет | `https://www.ebi.ac.uk/chembl/api/data` | Базовый URL API. |
| `timeout_connect` | Нет | `5` | Таймаут установления соединения (сек). |
| `timeout_read` | Нет | `90` | Таймаут чтения (сек). |
| `retries` | Нет | `3` | Попытки повторов. |
| `backoff_factor` | Нет | `0.5` | Множитель экспоненциальной задержки. |
| `rps` | Нет | `20` | Лимит запросов в секунду. |
| `burst` | Нет | `20` | Размер «всплеска» токен-бакета. |
| `user_agent` | **Да** | `chembl-da/1.0 (mailto:chembl-data@ebi.ac.uk)` | Перед продакшеном замените на контакт команды. |

### `sources.chembl.cache`

| Ключ | Значение по умолчанию | Описание |
|------|-----------------------|----------|
| `cache_ttl` | `3600` | TTL кэша (сек). |
| `cache_maxsize` | `1024` | Максимальное число записей. |

### `sources.chembl.molecule_catalog`

| Ключ | Значение по умолчанию | Описание |
|------|-----------------------|----------|
| `cache_path` | `$CHEMBL_DA_BASE_PATH/cache/molecule_parent_catalog.json` | JSON-кэш родительских связей. |
| `sqlite_path` | `$CHEMBL_DA_BASE_PATH/cache/molecule_parent_catalog.sqlite` | SQLite-кэш. |
| `endpoint` | `molecule` | Ресурс API для получения иерархий. |
| `child_field` | `molecule_chembl_id` | Колонка с дочерними идентификаторами. |
| `parent_field` | `parent_molecule_chembl_id` | Колонка с родительскими идентификаторами. |
| `hierarchy_lookup_path` | `config/dictionary/_testitem/molecule_hierarchy.csv` | Локальный CSV с иерархией. |
| `hierarchy_lookup_encoding` | `utf-8-sig` | Кодировка CSV. |
| `hierarchy_lookup_delimiter` | `,` | Разделитель CSV. |
| `page_size` | `500` | Размер страницы API. |

### `sources.chembl.pipelines`

Все значения необязательны и совпадают с дефолтами CLI; их можно переопределять
через YAML, переменные окружения или флаги.

#### `activity`

| Ключ | Значение | Описание |
|------|----------|----------|
| `column` | `activity_chembl_id` | Колонка входного CSV. |
| `batch_size` | `20` | Размер пакета. |
| `workers` | `1` | Количество потоков. |
| `timeout` | `90.0` | Таймаут чтения. |
| `limit` | `null` | Ограничение записей. |
| `dry_run` | `false` | Включает режим сухого прогона. |

#### `assay`

| Ключ | Значение | Описание |
|------|----------|----------|
| `column` | `assay_chembl_id` | Колонка входного CSV. |
| `batch_size` | `50` | Размер пакета. |
| `timeout` | `30.0` | Таймаут. |
| `limit` | `null` | Ограничение записей. |

#### `testitem`

| Ключ | Значение | Описание |
|------|----------|----------|
| `column` | `molecule_chembl_id` | Колонка с идентификаторами молекул. |
| `batch_size` | `250` | Размер батча. |
| `timeout` | `90.0` | Таймаут. |
| `limit` | `null` | Ограничение. |
| `offset` | `0` | Смещение. |
| `request_limit` | `1000` | Ограничение на число запросов. |
| `retries` | `5` | Повторы на уровне батча. |
| `backoff_factor` | `0.5` | Множитель задержки. |
| `batch_retry.enable` | `true` | Включить уменьшение батча при ошибках. |
| `batch_retry.shrink_factor` | `0.5` | Коэффициент уменьшения. |
| `batch_retry.min_size` | `1` | Минимальный размер. |
| `fields` | список | Поля API (см. YAML: идентификаторы родителей, PubChem и т.д.). |

#### `document`

Опции `pubmed`, `chembl`, `all` совпадают с CLI и описаны в
[`USAGE.md`](./USAGE.md).

#### `target`

Подразделы `uniprot`, `chembl`, `iuphar`, `all` задают колонки входа, лимиты,
пути к словарям. Значения совпадают с CLI.

### Внешние сервисы

| Блок | Обязательные поля | Примечание |
|------|--------------------|------------|
| `sources.openalex` | `mailto` (**Да**) | Требуется валидный контакт. |
| `sources.crossref` | `mailto` (**Да**) | CrossRef требует e-mail. |
| `sources.uniprot.api` | нет | Параметры REST API UniProt. |
| `sources.uniprot.mapping` | нет | Настройки ID Mapping (poll interval, timeout). |
| `sources.iuphar` | нет | Базовый URL и лимиты Guide to Pharmacology. |
| `sources.pubchem` | `user_agent` (**Да**) | Также настраиваются RPS, порядок поиска CID, кэши. |
| `sources.pubmed` | нет | Параметры E-utilities (URL, таймауты, повторы). |
| `sources.semantic_scholar` | нет | Базовый URL, таймауты и лимиты Semantic Scholar. |

## `local`

### `local.resources`

| Ключ | Значение | Описание |
|------|----------|----------|
| `dictionary_dir` | `config/dictionary` | Корень словарей. |
| `iuphar_target_csv` | `config/dictionary/_target/_IUPHAR/_IUPHAR_target.csv` | Таблица соответствий. |
| `iuphar_family_csv` | `config/dictionary/_target/_IUPHAR/_IUPHAR_family.csv` | Иерархия семейств. |
| `uniprot_data_dir` | `config/dictionary/_target/_uniprot` | Кэшированные JSON UniProt. |
| `targets_type_csv` | `config/dictionary/_target/targets_type.csv` | Таблица типов таргетов. |

### `local.io`

| Ключ | Значение | Описание |
|------|----------|----------|
| `output_dir` | `$CHEMBL_DA_BASE_PATH/output` | Каталог для экспортов. |
| `cache_dir` | `~/.cache/chembl-da` | Кэш HTTP. |
| `csv_sep` | `,` | Разделитель. |
| `csv_encoding` | `utf-8-sig` | Кодировка. |
| `csv_fallback_separators` | `["\t", ";"]` | Дополнительные разделители при автодетекте. |
| `csv_chunksize` | `10000` | Размер чанка при потоковой записи. |
| `na_markers` | `[#N/A]` | Значения, трактуемые как NA. |
| `keep_na_markers` | `false` | Сохранять NA вместо удаления. |
| `exist_ok` | `true` | Создавать каталоги при отсутствии. |
| `output_stamp_mode` | `omit` | Управляет именами: `omit` оставляет `output.<stem>.csv`, `require` требует явный `--date`. |

### `local.init`

| Ключ | Значение | Описание |
|------|----------|----------|
| `same_doc` | `$CHEMBL_DA_BASE_PATH/input/ChEMBL/ChEMBL_same_document_20_05.xlsx` | Источник для утилиты объединения Excel. |
| `all_doc` | `$CHEMBL_DA_BASE_PATH/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx` | Дополнительная книга. |
| `output_dir` | `$CHEMBL_DA_BASE_PATH/output/ChEMBL/processed` | Куда складываются обработанные файлы. |

## `activity_enrichment`

### `activity_enrichment.action_type`

| Ключ | Значение | Описание |
|------|----------|----------|
| `enabled` | `true` | Включить расчёт action type. |
| `column` | `action_type` | Целевая колонка. |
| `log_missing` | `true` | Логировать пропуски. |
| `log_distribution` | `true` | Выводить распределение. |
| `metrics` | словарь | Соответствие метрик (`ic50`, `ki`, …) action type. |
| `triages` | `{}` | Доп. правила триажа. |
| `functionality` | словарь | Карта функциональных описаний. |
| `mechanism` | `{}` | Зарезервировано под механизмы. |
| `triage_fields` | список | Колонки для поиска подсказок триажа. |
| `functionality_fields` | список | Колонки для функционального анализа. |
| `mechanism_fields` | список | Колонки для механизмов. |
| `allowlist` | список | Допустимые значения. |
| `positive_label` | `PAM` | Метка «положительного» класса. |
| `negative_label` | `NAM` | Метка «отрицательного» класса. |
| `fallback` | `unknown` | Значение по умолчанию. |

### `activity_enrichment.activity_properties`

| Ключ | Значение | Описание |
|------|----------|----------|
| `enabled` | `true` | Включить свёртку свойств. |
| `column` | `activity_properties` | Исходная колонка. |
| `summary_column` | `activity_property_summary` | Целевая колонка. |
| `name_field` | `type` | Поле с названием свойства. |
| `value_field` | `value` | Поле со значением. |
| `units_field` | `units` | Поле с единицами. |
| `separator` | `; ` | Разделитель между свойствами. |
| `pair_separator` | `=` | Разделитель пары имя/значение. |
| `drop_source_column` | `true` | Удалять исходную колонку после свёртки. |
| `log_missing` | `false` | Логировать пропуски. |
| `log_distribution` | `false` | Выводить распределение. |
| `allowlist` | список | Разрешённые категории. |
| `hash_column` | `properties_hash` | Колонка с хэшем структуры. |

## `testitem_molecule_enrichment`

| Ключ | Значение | Описание |
|------|----------|----------|
| `enable` | `true` | Включить обогащение. |
| `sources.molecule_catalog_path` | `config/dictionary/_testitem/molecule_catalog.csv` | Каталог молекул. |
| `sources.molecule_hierarchy_path` | `config/dictionary/_testitem/molecule_hierarchy.csv` | Иерархия молекул. |
| `output.salt_as_null_when_absent` | `true` | Пустые соли записывать как `NULL`. |
| `flags.coerce_to_bool` | `true` | Нормализация булевых значений. |
| `flags.parent_fallback` | `true` | Использовать fallback для родителя. |
| `logging.warn_missing_molecule` | `true` | Логировать отсутствующие молекулы. |
| `logging.warn_inconsistent_flags` | `true` | Предупреждать о конфликтующих флагах. |

## `activity_bounds`

| Ключ | Значение | Описание |
|------|----------|----------|
| `enable_from_relation` | `true` | Выводить границы из пары relation/value. |
| `enable_from_uncertainty` | `false` | Разбирать записи вида `значение ± погрешность`. |
| `rounding_digits` | `3` | Точность округления. |
| `clamp_nonnegative` | `true` | Ограничивать границы снизу нулём. |
| `log_unknown_relations` | `true` | Логировать незнакомые символы отношений. |

## `system`

### `system.log`

| Ключ | Значение | Описание |
|------|----------|----------|
| `level` | `INFO` | Уровень логирования (можно переопределить `--log-level`). |

### `system.rate`

| Ключ | Значение | Описание |
|------|----------|----------|
| `global_rps` | `100` | Глобальный лимит RPS. |
| `global_burst` | `100` | Размер всплеска. |
| `limiter_cache_maxsize` | `128` | Размер кэша лимитеров. |
| `limiter_cache_ttl` | `600` | TTL кэша (сек). |

### `system.retry`

| Ключ | Значение | Описание |
|------|----------|----------|
| `max_attempts` | `3` | Максимум попыток с учётом оригинального запроса. |
| `backoff_factor` | `0.5` | Множитель задержки. |
| `backoff_cap` | `null` | Максимальная задержка (опционально). |
| `status_forcelist` | `[429,500,502,503,504]` | Коды, вызывающие повтор. |

### `system.doc_quality`

| Ключ | Значение | Описание |
|------|----------|----------|
| `enable` | `true` | Включить профили качества. |
| `sample_rows` | `null` | Ограничение по строкам. |
| `include_columns` | `null` | Белый список колонок. |
| `exclude_columns` | `null` | Чёрный список колонок. |
| `fatal_on_error` | `false` | Считать ошибки профилирования фатальными. |

### `system.doc_type`

| Ключ | Значение | Описание |
|------|----------|----------|
| `weights` | `{pubmed:4, openalex:3, scholar:2}` | Веса источников. |
| `thresholds` | `{review:1, experimental:1, unknown:2}` | Пороговые значения. |
| `limit` | `null` | Ограничение по числу строк. |

Подробности и дополнительные поля см. в [`library/config/models.py`](../../library/config/models.py);
данный документ синхронизирован с актуальными моделями.

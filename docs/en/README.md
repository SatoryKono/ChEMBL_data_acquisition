# ChEMBL Data Acquisition Toolkit

Данный документ предоставляет обзор утилит для загрузки, нормализации и экспорта данных, полученных из ChEMBL. Используйте его как основную точку входа для инженеров и аналитиков перед изучением подробных руководств.

## Основные возможности

- **Конвейеры сущностей** для документов, таргетов, ассайев, активностей и тестовых объектов. Каждый конвейер обрабатывает идентификаторы из CSV, обращается к внешним API, выполняет детерминированную валидацию и формирует аудируемые экспорты.
- **Оркестратор** (`get-data`), объединяющий все конвейеры с общей конфигурацией и логированием.
- **Унифицированный CLI-слой** с общими опциями (`--input`, `--final-out`, `--log-level`, `--config` и др.), поддержкой сырых снимков (`--raw-out`, `--raw-format`, `--id-cols`, `--no-reindex-raw`) для таргет-конвейера и точками входа, объявленными в `pyproject.toml`.
- **Триада конфигурации** — значения берутся из `config/config.yaml`, переменных окружения (с префиксом `CHEMBL_DA__`) и флагов командной строки. Локально можно использовать `.env` с помощью `python-dotenv`.
- **Детерминированные выходные данные** — CSV-файлы с фиксированным порядком строк и столбцов, YAML-сайдкары с метаданными, SHA-256-хэши и отчёты о качестве таблиц для каждого запуска.
- **Инструменты для разработчиков** — строгая типизация (`mypy`), линтинг (`ruff`), форматирование (`black`), тестирование (`pytest`), покрытие и проверки детерминизма в CI.

## Структура репозитория

| Путь | Назначение |
|------|------------|
| `scripts/` | CLI-обёртки для каждого конвейера и оркестратора `get-data`. |
| `library/` | Основные модули: HTTP-клиенты, оркестрация, нормализация, валидация, постобработка, IO и вспомогательные функции. |
| `library/cli/commands/` | Точки входа для консольных скриптов, используемых в собранных пакетах. |
| `library/utils/cli_tools/` | Лёгкие утилиты (профилирование таблиц, кэширование таргетов, CSV-хелперы, инструменты сопоставления). |
| `config/` | YAML-конфигурация по умолчанию, схемы и встроенные словари. |
| `config/dictionary/` | Справочные датасеты для конвейеров (UniProt-кэши, таксономии таргетов, QA-фикстуры). |
| `data/` | Входные данные для smoke-тестов и примеры экспортов. |
| `docs/` | Документация проекта (английские `_EN.md` и русские `_RU.md` варианты). |
| `tests/` | Юнит- и интеграционные тесты для конвейеров и CLI-хелперов. |
| `Makefile` | Удобные цели для форматирования, тестов, сборки и линтинга. |

## Поддерживаемые точки входа

Установите проект (`pip install .` или `pip install dist/*.whl`), чтобы зарегистрировать следующие консольные скрипты. Они соответствуют модулям в `scripts/`, `library/cli/commands/` или `library/utils/cli_tools/` и принимают те же аргументы, что и их аналоги `python -m ...`.

| Консольный скрипт | Модуль | Описание |
|-------------------|--------|----------|
| `get-data` | `scripts.get_data:main` | Запуск всех конвейеров последовательно с общей конфигурацией. |
| `get-document-data` | `library.cli.commands.get_document_data:main` | Получение и обогащение документов ChEMBL. |
| `get-target-data` | `library.cli.commands.get_target_data:main` | Агрегация таргетов ChEMBL, UniProt и IUPHAR. |
| `get-assay-data` | `library.cli.commands.get_assay_data:main` | Экспорт метаданных ассайев. |
| `get-activity-data` | `library.cli.commands.get_activity_data:main` | Экспорт нормализованных записей активностей. |
| `get-testitem-data` | `library.cli.commands.get_testitem_data:main` | Обогащение молекул деталями из PubChem. |
| `get-document-type` | `library.utils.cli_tools.get_document_type:main` | Классификация типов публикаций для QA-задач. |
| `csv-utils` | `library.utils.cli_tools.csv_utils_main:main` | Детерминированные утилиты для работы с CSV. |
| `mapper` | `library.utils.cli_tools.mapper_main:main` | Интерактивный сопоставитель UniProt/ChEMBL. |
| `table-quality` | `library.utils.cli_tools.table_quality_main:main` | Генерация профилей качества по столбцам. |
| `chunk-io` | `library.utils.cli_tools.chunk_io_main:main` | Потоковая обработка CSV с сохранением порядка. |
| `get-input-initialisation` | `library.utils.cli_tools.get_input_initialisation:main` | Объединение Excel-файлов инициализации. |
| `get-activities` | `library.utils.cli_tools.get_activities:main` | Генерация синтетических строк активностей для smoke-тестов. |
| `check-determinism` | `library.utils.cli_tools.check_determinism:main` | Сравнение хэшей CSV между запусками. |

Подробное описание аргументов, подкоманд и расширенных сценариев использования смотрите в [`docs/USAGE_EN.md`](./USAGE_EN.md) и [`docs/USAGE_RU.md`](./USAGE_RU.md).

## Document pipeline

Use the document workflow through the unified development entry point `python scripts/get_document_data.py --mode <chembl|pubmed|all>`. The selected mode controls which acquisition stages run while sharing the common CLI flags described in the usage guide.

| `--mode` value | Purpose | Namespaced flags |
|----------------|---------|------------------|
| `chembl` | Retrieve document metadata from the ChEMBL API. | `--chembl-chunk-size`, `--chembl-timeout` (or shared `--chunk-size`, `--timeout`). |
| `pubmed` | Enrich with PubMed, Semantic Scholar, OpenAlex, and CrossRef. | `--pubmed-sleep`, `--pubmed-workers`, `--pubmed-batch-size`, `--openalex-rps`, `--crossref-rps`. |
| `all` | Run the ChEMBL and PubMed stages sequentially before merging the outputs. | Accepts both `chembl` and `pubmed` namespaces plus fallback DOI switches. |

Example invocations:

```bash
# ChEMBL-only export
python scripts/get_document_data.py --mode chembl \
    --input data/input/document.csv \
    --final-out output/documents_chembl.csv \
    --config config/config.yaml

# PubMed + partner services
python scripts/get_document_data.py --mode pubmed \
    --input data/input/document.csv \
    --final-out output/documents_pubmed.csv \
    --config config/config.yaml \
    --openalex-rps 3 --crossref-rps 3

# Full merge with DOI overrides
python scripts/get_document_data.py --mode all \
    --input data/input/document.csv \
    --final-out output/documents_full.csv \
    --config config/config.yaml \
    --fallback-doi-enabled \
    --fallback-doi-path data/input/document_fallback.csv \
    --fallback-doi-overwrite
```

The pipeline writes a deterministic CSV, `<name>.meta.yaml`, `<name>_quality_report_table.csv`, `<name>_data_correlation_report_table.csv`, and `<name>.quality.json` with DOI coverage statistics.

## Target pipeline

The `get-target-data` CLI exposes four acquisition modes (`chembl`, `uniprot`, `iuphar`, `all`). All of them share a consistent set of selector and pagination flags:

- `--column`: input column that stores the identifiers processed by the mode.
- `--chunk-size`: batch size for API calls or cached lookups.
- `--timeout`: network timeout in seconds applied to HTTP requests.
- `--limit`: maximum number of rows processed; omit to handle the entire CSV.
- `--offset`: number of rows skipped before processing begins.

The combined `all` mode forwards these values to its sub-pipelines and supports prefixed overrides such as `--chembl-chunk-size` or `--uniprot-timeout` when a particular source requires a different value. During start-up the script logs the resolved parameter values together with their origin (defaults, YAML configuration, or CLI), simplifying troubleshooting.

Refer to [`docs/en/USAGE.md`](./USAGE.md) for CLI examples and [`docs/en/reference/CONFIG.md`](./reference/CONFIG.md) for matching configuration keys and defaults.

## Требования и установка

| Компонент | Поддерживаемая версия | Последняя протестированная |
|-----------|----------------------|---------------------------|
| Python | 3.11.x | 3.11.12 |
| numpy | >=2.3.3,<3.0 | 2.3.3 |
| pandas | >=2.3.3,<3.0 | 2.3.3 |
| requests | >=2.32.5,<3.0 | 2.32.5 |
| PyYAML | >=6.0.3,<7.0 | 6.0.3 |
| openpyxl | >=3.1.5,<4.0 | 3.1.5 |
| pyarrow | >=17.0.0,<18.0 | 17.0.0 |
| jsonschema | >=4.25.1,<5.0 | 4.25.1 |
| pandera | >=0.26.1,<0.27 | 0.26.1 |
| pydantic | >=2.11.9,<3.0 | 2.11.9 |

Фиксированные зависимости указаны в `requirements-lock.txt`. Перегенерируйте этот файл только после изменения диапазонов в `pyproject.toml`, используя чистое виртуальное окружение и `pip freeze`.

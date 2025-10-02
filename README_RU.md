# Утилиты ChEMBL для получения данных

> **Версия проекта:** 0.2.0 (2025-10-02)
>
> Полный комплект документации на русском и английском языках расположен в
> каталоге [`docs/`](docs/). Каждому файлу соответствует синхронизированная
> версия на втором языке. Консолидированная история изменений ведётся в
> [`CHANGELOG.md`](CHANGELOG.md).

## Краткое описание

Пакет автоматизирует детерминированную выгрузку, нормализацию и экспорта
связанных с ChEMBL датасетов (активности, ассайи, документы, таргеты, тестовые
объекты) с жёсткой валидацией. Ключевые особенности:

- Единый интерфейс командной строки с потоковой обработкой чанков, проверкой
  схем и записью CSV/Parquet, устойчивой к повторным прогонкам.
- Конфигурация через `config/config.yaml`, переменные окружения и ключи CLI,
  валидируемые моделями pydantic и JSON Schema.
- Библиотеки для работы с ChEMBL, PubChem, UniProt, CrossRef, PubMed,
  Semantic Scholar и IUPHAR.
- Полный контур качества: pytest, статическая типизация, проверка
  детерминизма, тесты упаковки.

## Карта документации

| Раздел | Английский | Русский |
|--------|------------|---------|
| Конфигурация | [docs/CONFIG_EN.md](docs/CONFIG_EN.md) | [docs/CONFIG_RU.md](docs/CONFIG_RU.md) |
| CLI-утилиты | [docs/CLI_TOOLS.md](docs/CLI_TOOLS.md) | [docs/CLI_TOOLS_RU.md](docs/CLI_TOOLS_RU.md) |
| Схемы данных | [docs/DATA_SCHEMA_EN.md](docs/DATA_SCHEMA_EN.md) | [docs/DATA_SCHEMA_RU.md](docs/DATA_SCHEMA_RU.md) |
| Обзор ETL | [docs/ETL_PROCESS_EN.md](docs/ETL_PROCESS_EN.md) | [docs/ETL_PROCESS_RU.md](docs/ETL_PROCESS_RU.md) |
| Потоки данных | [docs/ETL_DATA_FLOW_EN.md](docs/ETL_DATA_FLOW_EN.md) | [docs/ETL_DATA_FLOW_RU.md](docs/ETL_DATA_FLOW_RU.md) |
| Результаты | [docs/OUTPUT_EN.md](docs/OUTPUT_EN.md) | [docs/OUTPUT_RU.md](docs/OUTPUT_RU.md) |
| QA и релизы | [docs/QA_PROCESS_EN.md](docs/QA_PROCESS_EN.md) | [docs/QA_PROCESS_RU.md](docs/QA_PROCESS_RU.md) |
| Практические сценарии | [docs/USAGE_EN.md](docs/USAGE_EN.md) | [docs/USAGE_RU.md](docs/USAGE_RU.md) |
| Сводка | [docs/SUMMARY.EN.md](docs/SUMMARY.EN.md) | [docs/SUMMARY.RU.md](docs/SUMMARY.RU.md) |

## Поддерживаемые CLI-пайплайны

При установке командой `pip install .` регистрируются консольные точки входа,
эквивалентные вызовам `python -m …`.

| Консольная команда | Модуль | Назначение |
|--------------------|--------|------------|
| `get-data` | `scripts.get_data` | Сквозной оркестратор всех пайплайнов |
| `get-activity-data` | `scripts.get_activity_data` | Выгрузка активностей с обогащением и валидацией |
| `get-assay-data` | `scripts.get_assay_data` | Получение и постобработка метаданных ассайев |
| `get-document-data` | `scripts.get_document_data` | Сбор документов PubMed, CrossRef, Semantic Scholar |
| `get-target-data` | `scripts.get_target_data` | Агрегация таргетов со стадийными снимками (`--raw-out`) |
| `get-testitem-data` | `scripts.get_testitem_data` | Обогащение тестовых объектов |
| `get-document-type` | `library.utils.cli_tools.get_document_type` | Классификация документов |
| `get-input-initialisation` | `library.utils.cli_tools.get_input_initialisation` | Объединение Excel-инициализаций |
| `csv-utils` | `library.utils.cli_tools.csv_utils_main` | Утилиты обслуживания CSV |
| `mapper` | `library.utils.cli_tools.mapper_main` | Маппинг идентификаторов |
| `table-quality` | `library.utils.cli_tools.table_quality_main` | Профилирование качества |
| `chunk-io` | `library.utils.cli_tools.chunk_io_main` | Проверка потокового IO |
| `get-activities` | `library.utils.cli_tools.get_activities` | Пример минимальной выгрузки активностей |
| `check-determinism` | `library.utils.cli_tools.check_determinism` | Регрессия детерминизма |

Ключи `--raw-out`, `--raw-format`, `--id-cols`, `--no-reindex-raw`,
`--normalize-at-export/--no-normalize-at-export` полностью поддерживаются
конвейером `get-target-data`. Остальные команды принимают параметры, но игнорируют
до внедрения собственных этапов выгрузки «сырых» снапшотов.

## Установка

1. **Подготовка окружения**
   ```bash
   python -m pip install --upgrade pip setuptools wheel
   python -m venv .venv
   source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
   ```

2. **Установка зависимостей**
   ```bash
   pip install -r requirements-lock.txt
   pip install -e .[dev]
   pre-commit install
   ```

3. **Проверка установки**
   ```bash
   pre-commit run --all-files
   pytest
   python -m library.utils.cli_tools.check_determinism --limit 10 --log-level INFO
   ```

Файл `requirements-lock.txt` синхронизирует рабочее окружение с CI. При изменении
диапазонов в `pyproject.toml` создайте новое окружение, выполните
`pip install .[dev]`, затем зафиксируйте версии через `pip freeze > requirements-lock.txt`.

## Конфигурация

- Значения по умолчанию находятся в [`config/config.yaml`](config/config.yaml) и
  описаны схемой [`config.schema.json`](config.schema.json).
- Локальные переопределения: `config/config.local.yaml` или собственный файл,
  переданный через `--config`; переменные окружения используют префикс
  `CHEMBL_DA__...`.
- Основные флаги CLI: `--input`, `--output` (основной результат), `--final-out`
  (стадийный экспорт таргетов), `--config`, `--print-config`, `--sep`,
  `--encoding`, `--chunk-size` / `--batch-size`, `--log-level`.
- Макрос `$CHEMBL_DA_BASE_PATH` подставляется из `--base-path` или переменной
  окружения `CHEMBL_DA_BASE_PATH`.

Полный перечень разделов, параметров ограничения скорости и настройки обогащения
приведён в руководстве по конфигурации.

## Результирующие данные

Документ [docs/OUTPUT_RU.md](docs/OUTPUT_RU.md) описывает структуру итоговых
таблиц по активностям, ассайям, документам, таргетам и тестовым объектам.
Снимки, создаваемые через `--raw-out`, сохраняют порядок колонок API; финальные
таблицы нормализуются и проверяются схемами из каталога [`schemas/`](schemas/).

## Тестирование и QA

- Юнит- и интеграционные тесты: `pytest`, акцент на детерминированность IO,
  проверку схем, пайплайн маппера и упаковку.
- Статический анализ: `pre-commit`, `ruff`, `black`, `mypy`.
- Регрессия: `python -m library.utils.cli_tools.check_determinism`,
  `python -m library.utils.cli_tools.mapper_batch_main` для массовых прогонов.
- Чек-листы QA: [docs/QA_PROCESS_RU.md](docs/QA_PROCESS_RU.md).

## Управление релизами

1. Обновить версии в `pyproject.toml`, README и сопроводительных документах.
2. Зафиксировать изменения в [`CHANGELOG.md`](CHANGELOG.md) и кратко отразить их
   в `docs/RELEASE_NOTES.md`.
3. Запустить полный набор проверок (линтеры, тесты, детерминизм, сборка колеса).
4. Создать тег (`git tag v0.2.0`) и опубликовать артефакты при необходимости.

Историю релизов и синхронизацию локализаций см. в changelog.

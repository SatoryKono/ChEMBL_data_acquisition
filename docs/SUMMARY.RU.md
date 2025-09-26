# Сводка по проекту

Документ описывает назначение утилит ChEMBL Data Acquisition, структуру
репозитория, основные сервисы и вспомогательные процессы.

## Структура репозитория

* `scripts/` — CLI-точки для пайплайнов активностей, ассайев, документов,
  таргетов и тест-объектов, а также обвязка кешированного таргет-пайплайна для
  офлайн-проверок.【F:scripts/get_activity_data.py†L1-L1160】【F:scripts/pipeline_targets_main.py†L1-L141】
* `library/` — переиспользуемые модули: клиенты API, ограничители запросов,
  нормализация, обогащение, валидация, детерминированный ввод-вывод, логирование
  и метаданные.【F:library/__init__.py†L1-L50】【F:library/io.py†L1-L236】
* `schemas/` — `pandera`-схемы и нормализаторы, обеспечивающие стабильные типы,
  порядок колонок и канонические значения.【F:schemas/__init__.py†L1-L16】
* `dictionary/` и `data/` — локальные словари, кэши API и входные Excel/CSV,
  используемые пайплайнами.【F:config.yaml†L96-L154】
* `docs/` — справочная документация.
* `tests/` — модульные и интеграционные тесты конфигурации, обогащений,
  детерминизма и CLI.【F:tests/test_activity_pipeline.py†L1-L220】

## Обзор потока данных

```
[CSV идентификаторов] --read_ids--> [Итератор] --chunk/limit--> [Пакеты]
        │                                           │
        ▼                                           ▼
  ChemblClient / внешние сервисы            Нормализация и обогащение
        │                                           │
        ▼                                           ▼
   Pandera validation  ──►  SidecarErrors (опц. failure CSV)
        │
        ▼
add_pipeline_metadata → write_csv_deterministic →
  <table>.csv + <table>.csv.meta.yaml + отчёты качества
```

* `io.read_ids` построчно считывает идентификаторы, отбрасывает пустые значения
  и проверяет наличие нужной колонки.【F:library/io.py†L87-L160】
* Доступ к API централизован в `ChemblClient` и смежных клиентах, которые
  учитывают лимиты, ретраи и таймауты из `config.yaml`.
* Нормализация и обогащение выполняются в скриптах с использованием модулей
  `document_pipeline`, `target_postprocessing`, `testitem_enrichment`,
  `activity_bounds`. Ошибочные строки собирает `SidecarErrors` и записывает
  рядом с экспортом.【F:library/sidecar.py†L1-L154】
* Детерминированная запись CSV, sidecar-метаданные и отчёты качества обеспечены
  функциями `write_csv_deterministic`, `write_meta_yaml`,
  `add_pipeline_metadata` и `analyze_table_quality`.
  【F:library/csv_utils.py†L451-L603】【F:library/metadata.py†L29-L133】

## Конфигурация

* Базовые значения задаются в [`config.yaml`](../config.yaml) и проверяются по
  [`config.schema.json`](../config.schema.json).
* Ключевые разделы:
  * `sources.*` — базовые URL, политика ретраев, лимиты и настройки пайплайнов
    для ChEMBL, UniProt, IUPHAR, PubMed, Semantic Scholar, OpenAlex, CrossRef и
    PubChem.【F:config.yaml†L11-L258】
  * `local.*` — структура каталогов, параметры CSV и входные рабочие книги.
    【F:config.yaml†L108-L154】
  * `activity_enrichment` / `activity_bounds` — параметры обогащения и расчёта
    границ значений в пайплайне активностей.【F:config.yaml†L155-L238】
  * `testitem_molecule_enrichment` — опциональная логика для солей и флагов
    каталога тест-объектов.【F:config.yaml†L239-L269】
  * `system.*` — логирование, глобальные лимиты, ретраи и веса классификатора
    публикаций.【F:config.yaml†L270-L315】
* Приоритет переопределений: `config.yaml` < переменные окружения < аргументы
  CLI. Доступны короткие алиасы (например, `CHEMBL_DA_RPS` и `CHEMBL_DA_OUTDIR`).
  Полный список приведён в `docs/CONFIG_RU.md`.

## Внешние сервисы

* **ChEMBL REST API** — основное хранилище активностей, ассайев, таргетов,
  документов и молекул; запросы дробятся на чанки и ретраятся согласно
  конфигурации.【F:library/chembl_client.py†L1-L286】
* **PubMed, Semantic Scholar, OpenAlex, CrossRef** — источники библиографических
  данных и DOI для документов.【F:scripts/get_document_data.py†L242-L533】
* **UniProt** — аннотации белков и ID mapping для таргет-пайплайна.
  【F:library/uniprot_library.py†L1-L357】
* **IUPHAR** — классификации рецепторов из локальных CSV-словарей.
  【F:library/target_postprocessing.py†L1-L599】
* **PubChem** — идентификаторы и свойства для обогащения тест-объектов.
  【F:library/testitem_enrichment.py†L17-L216】

## Установка и инструменты

1. Создайте и активируйте виртуальное окружение Python ≥ 3.12.
2. Установите проект с dev-зависимостями: `pip install .[dev]`.
3. Включите контроль качества: `pre-commit install`.
4. Рекомендуемые проверки:
   * `pre-commit run --all-files`
   * `pytest` / `pytest --cov=library --cov=scripts`
   * `ruff check`, `black --check .`, `mypy`

Версии зависимостей описаны в `pyproject.toml` и `requirements-dev.txt`.

## Использование

* Все `scripts/get_*_data.py` принимают стандартные флаги `--config`,
  `--print-config`, `--input`, `--output`, `--log-level`, `--sep`, `--encoding`,
  `--column` и `--batch-size` или `--chunk-size`.
* Дополнительные аргументы (`--timeout`, `--limit`, `--dry-run`, подкоманды для
  документов и таргетов) также прокидываются в конфигурацию через
  `apply_config_overrides` до запуска пайплайна.【F:library/cli.py†L1-L322】
* Логи в формате JSON содержат `run_id`, `event`, `stage` и счётчики по стадиям,
  что упрощает фильтрацию через `jq` и системы сбора логов.

Детальное руководство по запуску приведено в `docs/USAGE_EN.md` и
`docs/USAGE_RU.md`.

## Результаты выгрузки

* Основные CSV сохраняются в `local.io.output_dir` (по умолчанию `data/output`).
* Каждый запуск формирует `<name>.csv.meta.yaml` с конфигурацией, командой,
  статистикой строк и SHA-256 хэшем; ошибки валидации попадают в
  `<name>_failure_cases.csv`, а отчёты качества строятся автоматически.
  【F:library/metadata.py†L29-L133】【F:library/table_quality.py†L1-L192】
* Для документов дополнительно создаётся JSON с качеством, а таргет-пайплайн в
  режиме `all` пишет промежуточные файлы по каждому источнику.

Подробнее см. `docs/OUTPUT_EN.md` и `docs/OUTPUT_RU.md`.

## Тестирование и детерминизм

* Детерминированная запись CSV и контроль хэшей проверяются юнит-тестами и
  утилитами вроде `library.utils.cli_tools.check_determinism`.
  【F:library/utils/cli_tools/check_determinism.py†L1-L145】
* В `tests/data/` хранится набор компактных CSV для smoke-проверок активностей,
  таргетов, документов и тест-объектов.
* Утилита `python -m library.utils.cli_tools.table_quality_main` помогает
  анализировать сторонние данные перед загрузкой.【F:library/utils/cli_tools/table_quality_main.py†L1-L171】

Сводка служит точкой входа в остальную документацию каталога `docs/` и ускоряет
онбординг новых участников команды.

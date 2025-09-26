# Руководство по использованию

## Общие параметры CLI

Все команды `scripts/get_*_data.py` поддерживают единый набор аргументов:

| Опция | Назначение |
| --- | --- |
| `--config` | Путь к YAML-файлу конфигурации (по умолчанию `config.yaml`). |
| `--print-config` | Вывести итоговую конфигурацию после переопределений и завершить работу. |
| `--log-level` | Уровень логирования (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). |
| `--input` | Входной CSV с идентификаторами (по умолчанию `input.csv`). |
| `--output` | Выходной CSV. Если не указан, создаётся `output_<stem>_<YYYYMMDD>.csv` в `local.io.output_dir`. |
| `--sep` | Разделитель CSV; записывается в `cfg.io.csv_sep`. |
| `--encoding` | Кодировка файла; записывается в `cfg.io.csv_encoding`. |
| `--column` | Название колонки с идентификаторами. Значение подтягивается из конфигурации на этапе запуска. |
| `--batch-size` / `--chunk-size` | Максимальное число идентификаторов в одном запросе (название опции зависит от пайплайна). |

Конкретные скрипты добавляют дополнительные ключи (`--timeout`, `--limit`, `--dry-run` и т.п.). После разбора аргументов
`apply_config_overrides` загружает `config.yaml`, применяет переменные окружения, переносит CLI-переопределения в конфигурацию и
возвращает актуальные значения аргументов.

Перед сетевыми вызовами выполняется `library.config.ensure_dirs`, чтобы `local.io.output_dir` и `local.io.cache_dir` существовали
(если `local.io.exist_ok=true`, каталоги создаются автоматически).

## Данные активностей (`get_activity_data.py`)

```bash
python scripts/get_activity_data.py \
  --input data/input-smoke/activity.csv \
  --column activity_chembl_id \
  --batch-size 25 \
  --timeout 45
```

* Использует колонку `sources.chembl.pipelines.activity.column` (по умолчанию `activity_chembl_id`).
* Создаёт основной CSV, sidecar `*.meta.yaml`, при необходимости `*_failure_cases.csv` и отчёты качества.
* Поддерживает `--limit` (ограничение по количеству ID) и `--dry-run` (проверка входных данных без запросов к API).
* Заполняет `lower_value` и `upper_value` на основе канонических полей `standard_*`. Поведение (отношения, разбор `±`, округление, обрезка и логирование) настраивается блоком `activity_bounds.*` в конфигурации.
* Следите за предупреждениями `activity_bounds_unknown_relation` и `activity_bounds_missing_standard_value` в логах — они указывают строки, где границы не удалось восстановить или оператор отношения неизвестен.

## Описания ассайев (`get_assay_data.py`)

```bash
python scripts/get_assay_data.py \
  --input data/input-smoke/assay.csv \
  --column assay_chembl_id \
  --batch-size 25
```

Загружает метаданные ассайев ChEMBL для указанных идентификаторов.

## Метаданные документов (`get_document_data.py`)

```bash
python scripts/get_document_data.py all \
  --input data/input-smoke/documents.csv \
  --column document_chembl_id \
  --batch-size 20
```

Выберите подкоманду `pubmed`, `chembl` или `all` в зависимости от требуемых источников.
Сводку и список ключей смотрите в справке: `python scripts/get_document_data.py --help`
и `python scripts/get_document_data.py <подкоманда> --help`
(например, `--batch-size` управляет размером пакета для PubMed).

## Агрегация таргетов (`get_target_data.py`)

```bash
python scripts/get_target_data.py \
  --input data/input-smoke/targets.csv \
  --column target_chembl_id
```

Комбинирует данные ChEMBL, UniProt и IUPHAR согласно разделу `sources.chembl.pipelines.target.*`.

## Обогащение тест-объектов (`get_testitem_data.py`)

```bash
python scripts/get_testitem_data.py \
  --input data/input-smoke/testitem.csv \
  --column molecule_chembl_id
```

Выгружает дополнительную информацию о соединениях.

### Требования к каталогу родительских молекул

Чтобы получить `parent_molecule_chembl_id` для агрегаций, выгрузку необходимо объединить с каталогом
родителей ChEMBL. Путь к локальному JSON задаётся через
[`sources.chembl.molecule_catalog`](./CONFIG_RU.md#sources-chembl-molecule-catalog) (`cache_path`); убедитесь,
что файл доступен исполнителю, либо переопределите расположение параметрами CLI (`--sources.chembl.molecule_catalog.cache-path`)
или переменными окружения (`CHEMBL_DA_MOLECULE_CATALOG_CACHE`).【F:config.yaml†L25-L33】【F:library/config.py†L487-L551】

Для первичного создания либо обновления файла выполните небольшой Python-скрипт с вызовом
`library.molecule_catalog.load_parent_catalog` — функция считывает готовый кэш и, при его отсутствии,
подкачивает свежие связи ребёнок→родитель из API ChEMBL.【F:library/molecule_catalog.py†L43-L136】

## Инициализация входных данных (`get_input_initialisation.py`)

```bash
python scripts/get_input_initialisation.py \
  --same-doc data/input/ChEMBL/ChEMBL_same_document_20_05.xlsx \
  --all-doc data/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx \
  --out-dir data/output/ChEMBL/processed
```

* Формирует таблицы пар (`pairs_same_document.csv`, `pairs_independent.csv`, `pairs_non_independent.csv`).
* Создаёт срезы по сущностям (`activity_*`, `assay_*`, `document_*`, `target_*`, `testitem_*`, `system_*`).
* Добавляет папку `data_validity_report/` с отчётами качества для каждого файла.

## Профайлер качества таблиц (`table_quality_main.py`)

```bash
python table_quality_main.py \
  --input data/input-smoke/activity.csv \
  --table-name activity
```

Генерирует `<table-name>_quality_report_table.csv` и `<table-name>_data_correlation_report_table.csv`, используя настройки `local.io`.

## Переопределения через CLI

Любой параметр можно задать точечной записью:

```bash
# Временно поднять лимит ChEMBL по RPS
python scripts/get_activity_data.py --sources.chembl.api.rps 10

# Изменить разделитель CSV без правки config.yaml
python scripts/get_assay_data.py --sep ';'
```

## Переменные окружения

Формат `CHEMBL_DA__РАЗДЕЛ__...__КЛЮЧ`. Для популярных путей доступны короткие алиасы (например, `CHEMBL_DA_OUTDIR` →
`local.io.output_dir`). Полный список приведён в `docs/CONFIG_RU.md`.

## Запуск тестов

Создайте виртуальное окружение и установите dev-зависимости:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install .[dev]
```

После этого запустите unit- и smoke-наборы:

```bash
pytest
pytest tests/smoke
```

Для съёма покрытия основных модулей (`library/`, `scripts/`,
`activity_extraction_main.py`) используйте:

```bash
pytest --cov=library --cov=scripts --cov=activity_extraction_main \
       --cov-report=term-missing --cov-report=xml
```

Команда показывает непройденные строки в терминале и сохраняет `coverage.xml`
для CI либо IDE.

Для проверки качества кода можно выполнить:

```bash
ruff check
black --check .
mypy
```

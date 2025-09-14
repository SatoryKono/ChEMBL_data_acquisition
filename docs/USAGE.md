# Usage Examples

Ниже приведены примеры запуска CLI‑утилит на «smoke»‑наборах данных.
Во всех командах поддерживаются флаги:

* `--input` — путь к входному CSV (по умолчанию `input.csv`)
* `--output` — путь к результату  
  (по умолчанию `output_<stem>_YYYYMMDD.csv` в `io.output_dir`)
* `--sep` — разделитель CSV (по умолчанию `,`)
* `--encoding` — кодировка ввода/вывода (по умолчанию `utf-8`)
* `--log-level` — `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`
* `--column` — название столбца с идентификаторами
* `--dictionary` — путь к словарю/справочнику
* `--chunk-size` — размер порции строк при потоковой обработке

Перед началом каждый скрипт вызывает `library.config.ensure_dirs`,
создавая каталоги `io.output_dir` и `io.cache_dir`.

## Activity data

```bash
python scripts/get_activity_data.py \
    --input data/input-smoke/activity.csv \
    --column activity_id
```

Ожидается столбец `activity_id`. Выход — `output_activity_YYYYMMDD.csv`
и метаданные `*.meta.yaml`.

## Assay descriptions

```bash
python scripts/get_assay_data.py \
    --input data/input-smoke/assay.csv \
    --column assay_chembl_id
```

## Document metadata

```bash
python scripts/get_document_data.py \
    --input data/input-smoke/documents.csv \
    --column document_chembl_id
```

## Target data aggregation

```bash
python scripts/get_target_data.py \
    --input data/input-smoke/targets.csv \
    --column target_chembl_id
```

## Test item data enrichment

```bash
python scripts/get_testitem_data.py \
    --input data/input-smoke/testitem.csv \
    --column compound_chembl_id
```

## Input initialisation merging

```bash
python scripts/get_input_initialisation.py \
    --same-doc dictionary/classifications/assay_classification.csv \
    --all-doc dictionary/classifications/target_classification.csv \
    --out-dir data/output-smoke
```

Команда создаёт таблицы соответствий `pairs_*` и
фильтрованные `activity_*`, `assay_*`, `document_*`, `target_*`,
`testitem_*` для каждой категории.

## Table quality profiler

```bash
python table_quality_main.py \
    --input data/input-smoke/activity.csv \
    --table-name activity
```

Генерирует `<table-name>_quality_report_table.csv` и
`<table-name>_data_correlation_report_table.csv`.

## Конфигурационные переопределения

Параметры CLI заменяют значения `config.yaml`.  
Например, `--sep` → `io.csv_sep`, `--encoding` → `io.csv_encoding`,
`--chunk-size` → `jobs.chunk_size`, `--timeout` → `api.timeout_read`,
`--log-level` → `log.level`.

## Переменные окружения

Любое значение из `config.yaml` может быть переопределено переменной
окружения `CHEMBL_DA__SECTION__KEY`. Примеры:

```bash
export CHEMBL_DA__API__TIMEOUT_READ=60
export CHEMBL_DA__LOG__LEVEL=DEBUG
```

## Детерминированный CSV

`library.csv_utils.write_csv_deterministic` сортирует колонки и строки,
что гарантирует одинаковые файлы при повторных запусках.  
`scripts/check_determinism.py` проверяет это поведение.

## Запуск тестов

```bash
pytest
```

## Стиль кода

```bash
black scripts library mapper_main.py table_quality_main.py
ruff check scripts library mapper_main.py table_quality_main.py
mypy scripts library mapper_main.py table_quality_main.py
```


# ChEMBL Data Acquisition Utilities

Набор утилит и библиотек Python 3.12 для скачивания, валидации,
агрегации и экспорта биологических данных из открытых API  
(ChEMBL, PubChem, UniProt, PubMed и др.). Проект демонстрирует
типичный пайплайн обработки табличных данных: от получения
идентификаторов до сериализации нормализованных CSV/Parquet
с сопровождающими метаданными.

## Особенности

* Командные скрипты с унифицированными флагами `--input`, `--output`,
  `--log-level`, `--sep`, `--encoding`, `--column`, `--dictionary`.
* Потоковая обработка больших CSV через чанки, детерминированный вывод.
* Валидаторы схем (`schemas/`) и словари (`dictionary/`) для проверки
  типов, диапазонов и справочников.
* Конфигурация через `config.yaml`, переменные окружения и ключи CLI.
* Логирование через стандартный модуль `logging` с настраиваемым уровнем.
* Полная статическая типизация (PEP 484), линтинг `ruff`, форматирование
  `black`, проверка типов `mypy`, юнит‑тесты `pytest`.

## Требования

| Компонент     | Минимальная версия |
|---------------|-------------------|
| Python        | 3.12              |
| pandas        | 2.1               |
| requests      | 2.31              |
| PyYAML        | 6.0               |

Полный список приведён в `requirements-dev.txt` или `pyproject.toml`.

## Установка

```bash

git clone https://github.com/<org>/ChEMBL_data_acquisition.git

pip install .[dev]
# или установить в editable-режиме
pip install -e .[dev]
```

Install the project together with development tools such as
``black``, ``ruff``, ``mypy`` and ``pytest`` as well as testing utilities like
``hypothesis``, ``responses`` and ``psutil``. Sensitive configuration like API
tokens should be stored in a local ``.env`` file – see
[`Конфигурация через .env`](#конфигурация-через-env) for details.

After installing the dependencies, enable the pre-commit hooks so that
formatting, linting and type checking run automatically before each commit:

```bash
pre-commit install
```

To trigger all checks manually across the repository, execute:

```bash
pre-commit run --all-files
```

## Quick Start
 

1. **Install dependencies** – see [Установка](#установка).

2. **Initialise pre-commit hooks**

   ```bash
   pre-commit install
   pre-commit run --all-files
   ```

   The first command sets up Git hooks to run formatting, linting, static
   type checking and tests before each commit. The second command executes
   the hooks across the entire codebase.

3. **Run a sample script**

   ```bash
   python -m scripts.get_activities --input tests/data/activity_ids_small.csv \
       --output out/activities.csv --limit 10 --log-level INFO
   ```

   The command reads the bundled test identifiers and writes a normalised CSV
   to ``out/activities.csv``. Common CLI flags include ``--input`` and
   ``--output`` for file paths, ``--limit`` to cap processed records,
   ``--log-level`` for verbosity, ``--sep`` for CSV delimiter and
   ``--encoding`` for file encoding. Direct execution via ``python scripts/get_activities.py`` requires installing the project or adding the repository root to ``PYTHONPATH``; using ``-m scripts.get_activities`` avoids this extra setup. Additional examples:

   ```bash
   python mapper_main.py --input tests/data/assays.csv \
       --output out/mapped.csv --log-level DEBUG
   python table_quality_main.py --input tests/data/assays.csv \
       --output out/report.csv --log-level INFO
   ```

4. **Run the tests** – see [Тесты](#тесты).


## Тесты

Запустите линтеры, форматирование и проверки типов:

```bash
pre-commit run --all-files
```

Запустите тесты:

```bash
pytest
```

Для smoke-теста CLI выполните:

```bash
python mapper_batch_main.py --input tests/data/assays.csv \
    --output out/mapped.csv --log-level INFO
```


## Генерация данных

Скрипты из каталога `scripts/` создают CSV-файлы и сохраняют их в `data/output/`. Пример:

```bash
python -m scripts.get_activities --input tests/data/activity_ids_small.csv \
    --output data/output/activities.csv --limit 10 --log-level INFO
```

Результирующие файлы располагаются в `data/output/`. Каталог игнорируется Git и автоматически публикуется как артефакт CI.

## Usage

The examples below illustrate how to run the main CLI tools with common
options like ``--input``, ``--output`` and ``--limit``.

### scripts/get_document_data.py

Retrieve document metadata for a list of PubMed IDs using the bundled
sample file:

```bash
python scripts/get_document_data.py pubmed \
    --input tests/data/pmids.csv \
    --output out/documents.csv \
    --limit 5 \
    --log-level INFO
```

The ``tests/data/pmids.csv`` file contains a small set of PMIDs for
experimentation.

You can also run the PubMed pipeline directly using the library module:

```bash
python -m library.pubmed_library \
    --input-csv tests/data/pmids.csv \
    --output out/documents.csv \
    --log-level INFO
```

### scripts/get_target_data.py

Fetch basic target information from ChEMBL:

```bash
python scripts/get_target_data.py chembl \
    --input path/to/targets.csv \
    --output out/targets.csv \
    --limit 5 \
    --log-level INFO
```

Replace ``path/to/targets.csv`` with a CSV containing a ``target_chembl_id``
column.

The input and output both use ``target_chembl_id`` to align with
validation schemas.
 
### scripts/get_activities.py

Generate dummy activity entries without contacting external services:

```bash
python -m scripts.get_activities --limit 500 --dry-run
```

The command logs that it would generate 500 activity rows and exits without
creating any files.
 

## Updating Dependencies

To keep the environment current, periodically refresh the pinned
libraries and verify that the project remains compatible:

```bash
pip install -U .[dev]
pre-commit run --all-files
```

The first command upgrades runtime and development dependencies to the
newest minor releases permitted by the version ranges. The second command
formats code, lints, runs static type checks and executes the test suite
to confirm nothing broke during the upgrade.

## Конфигурация через `.env`

Часть параметров утилит можно задавать через переменные окружения.
Чтобы не экспортировать их вручную при каждом запуске, поместите пары
``NAME=value`` в файл `.env` и загрузите их с помощью пакета
[`python-dotenv`](https://pypi.org/project/python-dotenv/).

Пример файла:

```dotenv
CHEMBL_DA_LOG_LEVEL=INFO
CHEMBL_API_BASE=https://www.ebi.ac.uk/chembl/api/data
```

Запустить скрипт с автоматической подгрузкой настроек можно так:

```bash
python -m dotenv run -- python scripts/get_assay_data.py --input tests/data/assays.csv \\
    --output out/assays.csv
```

Утилиты читают переменные окружения автоматически, поэтому значения из
`.env` доступны всем CLI без дополнительных аргументов.

## Поле `api.user_agent`

Параметр `api.user_agent` используется для идентификации приложения в
запросах к API и должен содержать контактную информацию. Значение по
умолчанию:

```yaml
api:
  user_agent: "chembl-da/0.1 (mailto:info@example.org)"
```

Параметр можно переопределить в `config.yaml`, через переменную окружения
`CHEMBL_DA__API__USER_AGENT` или флаг CLI `--api.user_agent`.

## Валидация конфигурации

`library.config.load_config` проверяет корректность значений в `config.yaml`.
Некорректный URL приводит к `ValueError` при загрузке:

```yaml
api:
  chembl_base: https://
```

```
ValueError: api.chembl_base must be a valid URL
```

Исправленный вариант задаёт полный адрес службы:

```yaml
api:
  chembl_base: https://www.ebi.ac.uk/chembl/api/data
```

## Ошибки конфигурации

Некорректные значения в `config.yaml` вызывают `ValidationError`. Пример:

```yaml
api:
  rps: -1
```

При загрузке конфигурации:

```python
from library.config import load_config
load_config("config.yaml")
```

Вывод:

```
pydantic_core._pydantic_core.ValidationError: 1 validation error for Config
api.rps
  Input should be greater than or equal to 1 [type=greater_than_equal, input_value=-1, input_type=int]
    For further information visit https://errors.pydantic.dev/2.11/v/greater_than_equal
```

Исправьте значение на положительное число:

```yaml
api:
  rps: 5  # или любое >= 1
```

Диапазоны допустимых значений описаны в [`config.schema.json`](./config.schema.json), где для `api.rps` указан минимум `1`.

## Логирование

Пример включения JSON‑формата через переменную окружения:

```bash
LOG_FORMAT=json python scripts/get_assay_data.py --input tests/data/assays.csv \
    --output out/assays.csv --log-level INFO
```

Уровень логов можно задать флагом `--log-level` или переменной
`CHEMBL_DA_LOG_LEVEL`:

```bash
CHEMBL_DA_LOG_LEVEL=DEBUG python scripts/get_assay_data.py --input tests/data/assays.csv \
    --output out/assays.csv
```

Пример строки лога:

```json
{"ts":"2024-05-01T12:00:00Z","level":"INFO","event":"pipeline_start","run_id":"abc123","stage":"pipeline"}
```

Ключевые поля:

* `ts` – UTC timestamp in ISO 8601 format.
* `level` – severity such as `DEBUG`, `INFO`, `WARN` or `ERROR`.
* `event` – short machine-readable event name.
* `run_id` – unique identifier for the current run.
* `stage` – optional pipeline stage.
* `msg` – optional human-readable message.
* Additional keys – event specific context like `elapsed`, `url` or `rows`.

Dry-run executions emit logs with `event` set to `dry_run`, enabling easy
filtering, for example:

```bash
jq 'select(.event=="dry_run")' log.jsonl
```

The run identifier is generated by the CLI helpers using `uuid.uuid4().hex`
and passed to the logger, which includes it with every record. The value can
be overridden before calling `configure_logger` if a custom identifier is
desired.

Secrets are automatically redacted: values for keys ending in `token`, `key`,
`secret` or `password` are replaced with `"***"`. Log level filtering drops
records below the configured `--log-level` or `CHEMBL_DA_LOG_LEVEL` setting.

Typical log entries look like:

```json
{"ts":"2024-05-01T12:00:00Z","level":"INFO","event":"pipeline_start","run_id":"abc123","stage":"pipeline"}
{"ts":"2024-05-01T12:00:01Z","level":"INFO","event":"request_ok","run_id":"abc123","stage":"fetch","url":"https://api.example.org","status":200}
{"ts":"2024-05-01T12:00:02Z","level":"INFO","event":"validate_done","run_id":"abc123","stage":"validate","rows":42}
{"ts":"2024-05-01T12:00:03Z","level":"INFO","event":"pipeline_done","run_id":"abc123","stage":"pipeline","elapsed":3.2}
```

## Installation

Clone the repository and install the package together with development tools:

```bash
git clone https://example.com/ChEMBL_data_acquisition.git

cd ChEMBL_data_acquisition
pip install .[dev]       # с инструментами разработки
pre-commit install       # git‑хуки: black/ruff/mypy/pytest
pre-commit run --all-files
```

## Быстрый старт

```bash
 
# загрузка активности по идентификаторам из тестового CSV
python -m scripts.get_activity_data \
    --input tests/data/activity_ids_small.csv \
    --output out/activities.csv \
    --limit 10 --log-level INFO
 
python -m scripts.get_activities --limit 10 --dry-run
```

## Reproducibility

The function ``library.csv_utils.write_csv_deterministic`` normalises column
order, row sorting and value serialisation so repeated runs produce identical
files. Every CSV must be stored alongside a ``<name>.meta.yaml`` file capturing
the Git commit, command-line arguments and relevant configuration to allow
others to reproduce the output. Commit both the CSV and its metadata sidecar to
version control.

Verify deterministic behaviour with the helper script ``scripts/check_determinism.py``:

```bash
python scripts/check_determinism.py --log-level INFO
```

The script writes a sample CSV twice using ``write_csv_deterministic`` and
compares SHA-256 hashes. It requires the ``pandas`` package; install it with
``pip install pandas`` if it is not already available in your environment.
This check also runs in the project's CI pipeline and will fail the build
if the hashes differ.

For very large tables, ``write_csv_deterministic`` accepts a ``chunksize``
argument which streams the CSV in smaller pieces to reduce memory usage:

```python
from library.csv_utils import write_csv_deterministic
import pandas as pd

df = pd.read_csv("large.csv")
write_csv_deterministic(df, "out.csv", key_cols=df.columns, chunksize=1000)
```

Rows are still sorted deterministically before writing; ``chunksize`` only
affects how data is flushed to disk.

The higher-level wrapper ``library.io.write_csv`` exposes the same
``chunksize`` argument and additionally writes a metadata sidecar alongside
the CSV:

```python
from library import io, Config
import pandas as pd

cfg = Config()
df = pd.read_csv("large.csv")
io.write_csv(df, "out.csv", cfg=cfg, chunksize=1000)
```

The YAML sidecar records the Git commit and command-line parameters to aid
reproducibility.

Each command-line tool emits a ``<output>.meta.yaml`` sidecar file alongside
every CSV. The YAML document records the SHA-256 checksum, command-line
arguments and timestamps to make results reproducible. The determinism check is
executed in both pre-commit and continuous integration to guarantee that
regressions are detected early.

All commands emit the structured JSON logs described above. Adjust verbosity
with ``--log-level`` or ``CHEMBL_DA_LOG_LEVEL``.

Detailed command line examples using the bundled smoke datasets can be found in
``docs/USAGE.md``.
An overview of the output directory layout and metadata sidecars is available in
``docs/OUTPUT.md``.

### Table quality analysis

``table_quality_main.py`` profiles arbitrary CSV files and reports column
statistics along with correlations between numeric fields. Example usage:

```python
import pandas as pd
from library.table_quality import analyze_table_quality

df = pd.read_csv("data.csv", encoding="utf-8-sig")
quality, corr = analyze_table_quality(df, table_name="data")
```

Running the CLI saves ``data_quality_report_table.csv`` and
``data_data_correlation_report_table.csv`` in the current working directory::

    python table_quality_main.py --input data.csv --table-name data

All scripts share a common set of flags:

## Configuration


Default settings live in ``config.yaml`` and are split into sections for each
API (``api``, ``openalex``, ``crossref``, ``uniprot``, ``iuphar``, ``pubchem``),
I/O and processing (``io``, ``jobs``, ``batch``, ``quality``, ``mapper``) and
general infrastructure (``init``, ``rate``, ``retry``, ``log``). The companion
``config.schema.json`` file documents these fields and is useful for editor
validation, but it must **not** be passed to ``--config`` because it lacks
runtime values such as ``api.user_agent``. A minimal configuration looks like::


    api:
      rps: 5
    io:
      output_dir: data/output
    jobs:
      concurrency: 8

### Переменные окружения

Environment variables override values from the YAML file. Names follow the
``CHEMBL_DA__SECTION__KEY`` pattern with double underscores separating
sections and keys. For example, to enable debug logging:

```bash
export CHEMBL_DA__LOG__LEVEL=DEBUG
```

Most options also provide short aliases. The table lists the supported mappings:

| Alias | Equivalent key |
|-------|----------------|
| `CHEMBL_DA_BASE` | `CHEMBL_DA__API__CHEMBL_BASE` |
| `CHEMBL_DA_TIMEOUT_CONNECT` | `CHEMBL_DA__API__TIMEOUT_CONNECT` |
| `CHEMBL_DA_TIMEOUT_READ` | `CHEMBL_DA__API__TIMEOUT_READ` |
| `CHEMBL_DA_RPS` | `CHEMBL_DA__API__RPS` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_CONNECT` | `CHEMBL_DA__OPENALEX__TIMEOUT_CONNECT` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_READ` | `CHEMBL_DA__OPENALEX__TIMEOUT_READ` |
| `CHEMBL_DA_OPENALEX_RPS` | `CHEMBL_DA__OPENALEX__RPS` |
| `CHEMBL_DA_CROSSREF_TIMEOUT_CONNECT` | `CHEMBL_DA__CROSSREF__TIMEOUT_CONNECT` |
| `CHEMBL_DA_CROSSREF_TIMEOUT_READ` | `CHEMBL_DA__CROSSREF__TIMEOUT_READ` |
| `CHEMBL_DA_CROSSREF_RPS` | `CHEMBL_DA__CROSSREF__RPS` |
| `CHEMBL_DA_UNIPROT_TIMEOUT_CONNECT` | `CHEMBL_DA__UNIPROT__TIMEOUT_CONNECT` |
| `CHEMBL_DA_UNIPROT_TIMEOUT_READ` | `CHEMBL_DA__UNIPROT__TIMEOUT_READ` |
| `CHEMBL_DA_UNIPROT_RPS` | `CHEMBL_DA__UNIPROT__RPS` |
| `CHEMBL_DA_IUPHAR_TIMEOUT_CONNECT` | `CHEMBL_DA__IUPHAR__TIMEOUT_CONNECT` |
| `CHEMBL_DA_IUPHAR_TIMEOUT_READ` | `CHEMBL_DA__IUPHAR__TIMEOUT_READ` |
| `CHEMBL_DA_IUPHAR_RPS` | `CHEMBL_DA__IUPHAR__RPS` |
| `CHEMBL_DA_PUBCHEM_TIMEOUT_CONNECT` | `CHEMBL_DA__PUBCHEM__TIMEOUT_CONNECT` |
| `CHEMBL_DA_PUBCHEM_TIMEOUT_READ` | `CHEMBL_DA__PUBCHEM__TIMEOUT_READ` |
| `CHEMBL_DA_PUBCHEM_RPS` | `CHEMBL_DA__PUBCHEM__RPS` |
| `CHEMBL_DA_OUTDIR` | `CHEMBL_DA__IO__OUTPUT_DIR` |
| `CHEMBL_DA_CONCURRENCY` | `CHEMBL_DA__JOBS__CONCURRENCY` |
| `CHEMBL_DA_CHUNK_SIZE` | `CHEMBL_DA__JOBS__CHUNK_SIZE` |
| `CHEMBL_DA_RETRY_MAX_ATTEMPTS` | `CHEMBL_DA__RETRY__MAX_ATTEMPTS` |
| `CHEMBL_DA_RETRY_BACKOFF_FACTOR` | `CHEMBL_DA__RETRY__BACKOFF_FACTOR` |
| `CHEMBL_DA_LOG_LEVEL` | `CHEMBL_DA__LOG__LEVEL` |
| `CHEMBL_DA_LOG_FORMAT` | `CHEMBL_DA__LOG__FORMAT` |

See ``docs/CONFIG.md`` for a complete overview of all configuration options.

### Schema validation

Configuration values are validated against a JSON Schema via the
``jsonschema`` package. The schema mirrors the dataclass structure and checks
types and value ranges, producing helpful error messages for nested fields.

Command line flags have the highest priority. All utilities accept ``--config``
to point at a configuration file and ``--print-config`` to show the effective
values after all overrides have been applied. The final precedence is::

    YAML < environment variables < CLI options

Only the top-level command line scripts read the configuration file. Modules
under ``library/`` expect a :class:`Config` (or one of its subsections) to be
passed explicitly, making dependencies clear and avoiding hidden global state.
The directories referenced by ``io.output_dir`` and ``io.cache_dir`` are checked
but not created when loading the configuration. Scripts that need these paths
can call :func:`library.config.ensure_dirs` after :func:`load_config` to create
them if they are missing and ``io.exist_ok`` permits it.

Path values such as ``io.output_dir``, ``io.cache_dir`` and the ``init``
workbook paths are exposed as :class:`pathlib.Path` objects. String values in
``config.yaml`` or overrides from the environment and command line are
automatically converted.
 

# профилирование качества таблицы
python table_quality_main.py \
    --input tests/data/activity.csv \
    --table-name activity
```

`--output` по умолчанию формируется как `output_<имя_входа>_YYYYMMDD.csv`
в каталоге, заданном `io.output_dir`.  
Для дополнительных примеров см. [`docs/USAGE.md`](docs/USAGE.md).

## Структура проекта

```
ChEMBL_data_acquisition/
├── config.yaml
├── dictionary/
├── library/
│   ├── __init__.py
│   ├── chembl_client.py
│   ├── csv_utils.py
│   ├── config.py
│   └── ...
├── scripts/
│   ├── get_activity_data.py
│   ├── get_assay_data.py
│   ├── ...
├── tests/
│   └── data/
└── docs/
    ├── CONFIG.md
    ├── OUTPUT.md
    └── USAGE.md
```

## Конфигурация

Параметры читаются из `config.yaml`, переменных окружения
(`CHEMBL_DA__SECTION__KEY`) и ключей CLI.  
Подробности в [`docs/CONFIG.md`](docs/CONFIG.md).

## Вывод и метаданные

Все сгенерированные CSV/Parquet и отчёты сохраняются в `data/output`
(см. [`docs/OUTPUT.md`](docs/OUTPUT.md)).  
Рядом создаются файлы `*.meta.yaml` с коммитом Git, параметрами запуска,
контрольной суммой SHA‑256 и статистикой строк/колонок.

## Разработка и тестирование

```bash
ruff check scripts library mapper_main.py table_quality_main.py
black scripts library mapper_main.py table_quality_main.py
mypy scripts library mapper_main.py table_quality_main.py
pytest
```

Тестовые наборы расположены в `tests/data`.  
Скрипт `scripts/check_determinism.py` проверяет повторяемость CSV‑вывода.

## Лицензия

MIT License. См. файл `LICENSE` (если присутствует).


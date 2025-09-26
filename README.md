# ChEMBL Data Acquisition Utilities

Основная документация и описательные материалы перенесены в каталог [docs/](docs/).

## Особенности

* Командные скрипты с унифицированными флагами `--input`, `--output`,
  `--log-level`, `--sep`, `--encoding`, `--column`, а также `--config` и
  `--print-config` для управления загрузкой настроек. Размер пакетной
  выборки задаётся параметрами `--chunk-size` или `--batch-size` в
  зависимости от конкретного пайплайна.
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
| pandas        | 2.3.2             |
| requests      | 2.32.5            |
| PyYAML        | 6.0.2             |

Полный список приведён в `requirements-dev.txt` или `pyproject.toml`.

### Среда выполнения

* ОС Linux или macOS с доступом к bash/PowerShell. На Windows рекомендуется устанавливать через WSL2.
* Актуальные версии `pip`, `setuptools` и `wheel`, см. шаги в разделе [Installation / Установка](#installation--установка).
* Разрешение на сетевые запросы к API ChEMBL/PubChem/UniProt (порт 443).

## Installation / Установка

### Runtime environment / Среда выполнения

* Linux or macOS shell with Bash or PowerShell support (Windows users should rely on WSL2). / ОС Linux или macOS с доступом к bash/PowerShell (на Windows используйте WSL2).
* Upgrade packaging tools before installing the project. / Перед установкой обновите инструменты распространения пакетов.

  ```bash
  python -m pip install --upgrade pip setuptools wheel
  ```

* Create an isolated virtual environment to keep dependencies local. / Создайте изолированное виртуальное окружение, чтобы зависимости не конфликтовали.

  ```bash
  python -m venv .venv
  source .venv/bin/activate  # Windows: .venv\Scripts\activate
  ```

### Steps / Шаги

**EN.** Clone the repository, switch into it and install the package with development extras. Afterwards enable the pre-commit hooks so formatting, linting, type checking and unit tests run automatically.

**RU.** Клонируйте репозиторий, перейдите в каталог и установите пакет с dev-зависимостями. Затем активируйте pre-commit-хуки для автоматического форматирования, линтинга, проверки типов и тестов.

```bash
git clone https://github.com/<org>/ChEMBL_data_acquisition.git
cd ChEMBL_data_acquisition
pip install .[dev]
pre-commit install
```

Sensitive configuration such as API tokens belongs in a local ``.env`` file – see [`Конфигурация через .env`](#конфигурация-через-env) for usage guidelines.

## Quick Start / Быстрый старт

1. **Install dependencies / Установите зависимости** – follow the steps in [Installation / Установка](#installation--установка).

2. **Initialise pre-commit hooks / Активируйте pre-commit-хуки**

   ```bash
   pre-commit install
   ```

   EN: Git hooks ensure formatting, linting, static type checks and tests run before each commit.

   RU: Git-хуки автоматически запускают форматирование, линтеры, проверки типов и тесты перед каждым коммитом.

3. **Run a sample script / Запустите демонстрационный скрипт**

   ```bash
   python -m library.utils.cli_tools.get_activities --limit 10 --log-level INFO
   ```

   This lightweight helper only emits structured log messages describing the
   dummy activity rows it would generate; it neither reads input files nor
   writes outputs. Use it to verify logging configuration and CLI wiring before
   launching full pipelines. Common CLI flags include ``--limit`` to cap
   processed records, ``--log-level`` for verbosity, ``--sep`` for CSV delimiter
   and ``--encoding`` for file encoding. For end-to-end exports that create
   files, run one of the data pipelines, for example:

  ```bash
  python -m library.utils.cli_tools.mapper_main --input tests/data/chembl_targets_min.csv \
      --column target_chembl_id --output out/targets_mapped.csv --log-level DEBUG
  python -m library.utils.cli_tools.table_quality_main --input tests/data/chembl_targets_min.csv \
      --output out/quality --table-name chembl_targets --log-level INFO
  ```

  Во втором примере аргумент `--output` должен указывать на каталог, в котором
  будут созданы файлы отчёта.

4. **Run the tests / Запустите тесты** – refer to [Tests / Тесты](#tests--тесты).


## Tests / Тесты

**EN.** The `pre-commit` suite runs formatting, linting and static type checks. Execute `pytest` for unit tests and add coverage flags when required. Determinism and smoke checks are available through dedicated CLI helpers.

**RU.** Команда `pre-commit` запускает форматирование, линтеры и проверку типов. Для юнит-тестов используйте `pytest`, при необходимости добавляйте параметры покрытия. Детеминизм и smoke-проверки доступны в отдельных CLI.

```bash
pre-commit run --all-files
pytest
pytest --cov=library --cov=scripts --cov-report=term-missing --cov-report=xml
python -m library.utils.cli_tools.check_determinism --log-level DEBUG
python -m library.utils.cli_tools.mapper_batch_main --input chembl_ids.csv \
    --output out/mapped.csv --log-level INFO
```

Before running the smoke command, create a `chembl_ids.csv` file with a header `chembl_id` and the required identifiers. / Перед запуском smoke-команды создайте `chembl_ids.csv` со столбцом `chembl_id` и нужными идентификаторами.


## Генерация данных

Шесть основных пайплайнов в каталоге `scripts/` формируют CSV-файлы и
размещают их в `data/output/`:

* `get_activity_data.py` — извлекает данные активностей из ChEMBL и
  дополняет их расчётными границами значений.
* `get_assay_data.py` — выгружает описания ассайев.
* `get_document_data.py` — объединяет метаданные публикаций из ChEMBL и
  агрегаторов (PubMed, Semantic Scholar, OpenAlex, Crossref).
* `get_target_data.py` — собирает информацию о таргетах из ChEMBL, UniProt
  и IUPHAR.
* `get_testitem_data.py` — обогащает соединения структурными атрибутами и
  данными PubChem.
* `pipeline_targets_main.py` — лёгкая обвязка для
  `library.pipeline_targets.run_pipeline`, использующая те же параметры CLI,
  что и боевой таргет-пайплайн, но работающая только с локальными файлами и
  подготовленными чанками идентификаторов без сетевых вызовов.

Пример полноценного пайплайна:

```bash
python -m scripts.get_activity_data --input tests/data/activity_ids_small.csv \
    --output data/output/activities.csv --limit 10 --log-level INFO
```

Команда извлекает данные из API ChEMBL, сохраняет таблицу и сопутствующий
`*.meta.yaml`. Утилиты разработки и отладки перенесены в
`library/utils/cli_tools/`, например модуль `get_activities` предназначен
только для демонстрационного логирования и не выполняет файловых операций.
См. [docs/CLI_TOOLS.md](docs/CLI_TOOLS.md) для кратких описаний и типовых
команд. Каталог с результатами игнорируется Git и автоматически публикуется
как артефакт CI.

> **Примечание.** Ранее эта функциональность была доступна через
> `activity_extraction_main.py`. Теперь используйте модульный запуск
> `python -m scripts.get_activity_data`, что упрощает сопровождение и
> работу в виртуальных окружениях.

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

### scripts/pipeline_targets_main.py

Exercise the chunking and batch configuration used by the production
target pipeline without contacting remote services:

```bash
python scripts/pipeline_targets_main.py \
    --input tests/data/chembl_targets_min.csv \
    --output out/targets_cached.csv \
    --chunk-size 25 \
    --batch-size 25 \
    --limit 100
```

The command reads target identifiers from the CSV, chunks them according
to ``--chunk-size`` and ``--limit``, forwards the batch size to
``library.pipeline_targets.run_pipeline`` and writes the cached ChemBL
table with pipeline metadata via ``write_csv``. Use it to verify CLI
overrides, logging and deterministic output before launching the network
backed ``get_target_data`` pipeline.【F:scripts/pipeline_targets_main.py†L1-L141】

### library/utils/cli_tools/get_activities.py

Generate dummy activity entries without contacting external services:

```bash
python -m library.utils.cli_tools.get_activities --limit 500 --dry-run
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

См. также файл `.env.example` с типовыми переменными для контактных e-mail.

Запустить скрипт с автоматической подгрузкой настроек можно так:

```bash
python -m dotenv run -- python scripts/get_assay_data.py --input assay_ids.csv \\
    --output out/assays.csv
```

Файл `assay_ids.csv` должен содержать столбец `assay_chembl_id` с нужными
идентификаторами, например:

```csv
assay_chembl_id
CHEMBL1234567
CHEMBL2345678
```

Утилиты читают переменные окружения автоматически, поэтому значения из
`.env` доступны всем CLI без дополнительных аргументов.

## Поле `api.user_agent`

Параметр `api.user_agent` используется для идентификации приложения в
запросах к API и должен содержать контактную информацию. Значение по
умолчанию:

```yaml
api:
  user_agent: "chembl-da/0.1 (mailto:contact@example.org)"
```

Параметр можно переопределить в `config.yaml`, через переменную окружения
`CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT` или флаг CLI
`--sources.chembl.api.user_agent`.

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

## Logging / Логирование

**EN.** CLI helpers configure structured JSON logging via ``library.logging_setup.configure_logger``. Use environment variables or CLI flags to adjust verbosity and format.

**RU.** CLI-хелперы настраивают структурированное JSON-логирование через ``library.logging_setup.configure_logger``. Управляйте форматом и уровнем логов переменными окружения или ключами CLI.

Пример включения JSON‑формата через переменную окружения:

```bash
LOG_FORMAT=json python scripts/get_assay_data.py --input assay_ids.csv \
    --output out/assays.csv --log-level INFO
```

Уровень логов можно задать флагом `--log-level` или переменной
`CHEMBL_DA_LOG_LEVEL`:

```bash
CHEMBL_DA_LOG_LEVEL=DEBUG python scripts/get_assay_data.py --input assay_ids.csv \
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

## Reproducibility / Воспроизводимость

**EN.** Deterministic CSV writers in ``library.io`` keep outputs and metadata stable across runs.

**RU.** Детерминированные CSV-выгрузки из ``library.io`` обеспечивают повторяемость данных и метаданных между запусками.

The function ``library.csv_utils.write_csv_deterministic`` normalises column
order, row sorting and value serialisation so repeated runs produce identical
files. Every CSV must be stored alongside a ``<name>.meta.yaml`` file capturing
the Git commit, command-line arguments and relevant configuration to allow
others to reproduce the output. Commit both the CSV and its metadata sidecar to
version control.

Verify deterministic behaviour with the helper script ``library.utils.cli_tools.check_determinism``:

```bash
python -m library.utils.cli_tools.check_determinism --log-level INFO
```

The script writes a sample CSV twice using ``write_csv_deterministic`` and
compares SHA-256 hashes. It requires the ``pandas`` package; install it with
``pip install pandas==2.3.2`` if it is not already available in your environment
so the versions stay aligned with ``pyproject.toml``.
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
``docs/USAGE_EN.md`` (русская версия — ``docs/USAGE_RU.md``).
An overview of the output directory layout and metadata sidecars is available in
``docs/OUTPUT_EN.md`` (русская версия — ``docs/OUTPUT_RU.md``).

### Table quality analysis

``library.utils.cli_tools.table_quality_main`` profiles arbitrary CSV files and reports column
statistics along with correlations between numeric fields. Example usage:

```python
import pandas as pd
from library.table_quality import analyze_table_quality

df = pd.read_csv("data.csv", encoding="utf-8-sig")
quality, corr = analyze_table_quality(df, table_name="data")
```

Running the CLI saves ``data_quality_report_table.csv`` and
``data_data_correlation_report_table.csv`` in the current working directory::

    python -m library.utils.cli_tools.table_quality_main --input data.csv --table-name data

All scripts share a common set of flags:

## Configuration


Default settings live in ``config.yaml`` and are grouped into three top-level
sections:

* ``sources`` – external services such as ChEMBL, OpenAlex, Crossref, UniProt,
  IUPHAR and PubChem (including pipeline-specific settings).
* ``local`` – file system inputs and outputs, cached resources and workbook
  paths.
* ``system`` – shared concerns such as logging, retry strategy, rate limiting
  and document-type normalisation.

The companion ``config.schema.json`` file documents these fields and is useful
for editor validation, but it must **not** be passed to ``--config`` because it
lacks runtime values such as ``sources.chembl.api.user_agent``. A minimal
configuration looks like::


    sources:
      chembl:
        api:
          rps: 5
    local:
      io:
        output_dir: data/output

### Переменные окружения

Environment variables override values from the YAML file. Names follow the
``CHEMBL_DA__...`` pattern with double underscores separating each nested
section. For example, to enable debug logging:

```bash
export CHEMBL_DA__LOG__LEVEL=DEBUG
```

Most options also provide short aliases for backwards compatibility. The table
lists every supported alias and the canonical key it maps to:

| Alias | Equivalent key |
|-------|----------------|
| `CHEMBL_DA_BASE` | `CHEMBL_DA__SOURCES__CHEMBL__API__CHEMBL_BASE` |
| `CHEMBL_DA_BURST` | `CHEMBL_DA__SOURCES__CHEMBL__API__BURST` |
| `CHEMBL_DA_CACHE_DIR` | `CHEMBL_DA__LOCAL__IO__CACHE_DIR` |
| `CHEMBL_DA_CACHE_MAXSIZE` | `CHEMBL_DA__SOURCES__CHEMBL__CACHE__CACHE_MAXSIZE` |
| `CHEMBL_DA_CACHE_TTL` | `CHEMBL_DA__SOURCES__CHEMBL__CACHE__CACHE_TTL` |
| `CHEMBL_DA_CROSSREF_BASE` | `CHEMBL_DA__SOURCES__CROSSREF__BASE` |
| `CHEMBL_DA_CROSSREF_MAILTO` | `CHEMBL_DA__SOURCES__CROSSREF__MAILTO` |
| `CHEMBL_DA_CROSSREF_TIMEOUT_CONNECT` | `CHEMBL_DA__SOURCES__CROSSREF__TIMEOUT_CONNECT` |
| `CHEMBL_DA_CROSSREF_TIMEOUT_READ` | `CHEMBL_DA__SOURCES__CROSSREF__TIMEOUT_READ` |
| `CHEMBL_DA_DICT_DIR` | `CHEMBL_DA__LOCAL__RESOURCES__DICTIONARY_DIR` |
| `CHEMBL_DA_GLOBAL_BURST` | `CHEMBL_DA__SYSTEM__RATE__GLOBAL_BURST` |
| `CHEMBL_DA_GLOBAL_RPS` | `CHEMBL_DA__SYSTEM__RATE__GLOBAL_RPS` |
| `CHEMBL_DA_IUPHAR_BASE` | `CHEMBL_DA__SOURCES__IUPHAR__BASE` |
| `CHEMBL_DA_IUPHAR_FAMILY_CSV` | `CHEMBL_DA__LOCAL__RESOURCES__IUPHAR_FAMILY_CSV` |
| `CHEMBL_DA_IUPHAR_TARGET_CSV` | `CHEMBL_DA__LOCAL__RESOURCES__IUPHAR_TARGET_CSV` |
| `CHEMBL_DA_LIMITER_CACHE_MAXSIZE` | `CHEMBL_DA__SYSTEM__RATE__LIMITER_CACHE_MAXSIZE` |
| `CHEMBL_DA_LIMITER_CACHE_TTL` | `CHEMBL_DA__SYSTEM__RATE__LIMITER_CACHE_TTL` |
| `CHEMBL_DA_LOG_DATEFMT` | `CHEMBL_DA__SYSTEM__LOG__DATEFMT` |
| `CHEMBL_DA_LOG_FORMAT` | `CHEMBL_DA__SYSTEM__LOG__FORMAT` |
| `CHEMBL_DA_LOG_LEVEL` | `CHEMBL_DA__SYSTEM__LOG__LEVEL` |
| `CHEMBL_DA_MOLECULE_CATALOG_CACHE` | `CHEMBL_DA__SOURCES__CHEMBL__MOLECULE_CATALOG__CACHE_PATH` |
| `CHEMBL_DA_OPENALEX_BASE` | `CHEMBL_DA__SOURCES__OPENALEX__BASE` |
| `CHEMBL_DA_OPENALEX_MAILTO` | `CHEMBL_DA__SOURCES__OPENALEX__MAILTO` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_CONNECT` | `CHEMBL_DA__SOURCES__OPENALEX__TIMEOUT_CONNECT` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_READ` | `CHEMBL_DA__SOURCES__OPENALEX__TIMEOUT_READ` |
| `CHEMBL_DA_OUTDIR` | `CHEMBL_DA__LOCAL__IO__OUTPUT_DIR` |
| `CHEMBL_DA_PUBCHEM_BASE` | `CHEMBL_DA__SOURCES__PUBCHEM__BASE` |
| `CHEMBL_DA_RETRY_BACKOFF_FACTOR` | `CHEMBL_DA__SYSTEM__RETRY__BACKOFF_FACTOR` |
| `CHEMBL_DA_RETRY_MAX_ATTEMPTS` | `CHEMBL_DA__SYSTEM__RETRY__MAX_ATTEMPTS` |
| `CHEMBL_DA_RPS` | `CHEMBL_DA__SOURCES__CHEMBL__API__RPS` |
| `CHEMBL_DA_TARGETS_TYPE_CSV` | `CHEMBL_DA__LOCAL__RESOURCES__TARGETS_TYPE_CSV` |
| `CHEMBL_DA_TIMEOUT_CONNECT` | `CHEMBL_DA__SOURCES__CHEMBL__API__TIMEOUT_CONNECT` |
| `CHEMBL_DA_TIMEOUT_READ` | `CHEMBL_DA__SOURCES__CHEMBL__API__TIMEOUT_READ` |
| `CHEMBL_DA_UNIPROT_BASE` | `CHEMBL_DA__SOURCES__UNIPROT__API__BASE` |
| `CHEMBL_DA_UNIPROT_DATA_DIR` | `CHEMBL_DA__LOCAL__RESOURCES__UNIPROT_DATA_DIR` |
| `CHEMBL_DA__IO__CACHE_DIR` | `CHEMBL_DA__LOCAL__IO__CACHE_DIR` |
| `CHEMBL_DA__IO__EXIST_OK` | `CHEMBL_DA__LOCAL__IO__EXIST_OK` |

See ``docs/CONFIG_EN.md`` for a complete overview of all configuration options
(русская версия — ``docs/CONFIG_RU.md``).

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
The directories referenced by ``local.io.output_dir`` and ``local.io.cache_dir``
are checked but not created when loading the configuration. Scripts that need
these paths can call :func:`library.config.ensure_dirs` after
:func:`load_config` to create them if they are missing and ``local.io.exist_ok``
permits it.

Path values such as ``local.io.output_dir``, ``local.io.cache_dir`` and the ``local.init``
workbook paths are exposed as :class:`pathlib.Path` objects. String values in
``config.yaml`` or overrides from the environment and command line are
automatically converted.
 
```bash
# профилирование качества таблицы
python -m library.utils.cli_tools.table_quality_main \
    --input tests/data/activity.csv \
    --table-name activity
```

`--output` по умолчанию формируется как `output_<имя_входа>_YYYYMMDD.csv`
в каталоге, заданном `local.io.output_dir`.
Для дополнительных примеров см. [`docs/USAGE_RU.md`](docs/USAGE_RU.md) и английскую версию [`docs/USAGE_EN.md`](docs/USAGE_EN.md).

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
    ├── CONFIG_EN.md / CONFIG_RU.md
    ├── OUTPUT_EN.md / OUTPUT_RU.md
    └── USAGE_EN.md / USAGE_RU.md
```

## Конфигурация

Параметры читаются из `config.yaml`, переменных окружения
(`CHEMBL_DA__...`) и ключей CLI.
Подробности в [`docs/CONFIG_RU.md`](docs/CONFIG_RU.md) и английской версии [`docs/CONFIG_EN.md`](docs/CONFIG_EN.md).

## Output and metadata / Вывод и метаданные

**EN.** Pipelines persist deterministic CSV tables via ``library.io.write_csv`` and place accompanying ``*.meta.yaml`` sidecars in ``data/output``.

**RU.** Пайплайны сохраняют детерминированные CSV с помощью ``library.io.write_csv`` и добавляют рядом файлы ``*.meta.yaml`` в каталоге ``data/output``.

Each sidecar stores the Git commit, launch parameters, SHA‑256 checksum and row/column statistics. See [`docs/OUTPUT_EN.md`](docs/OUTPUT_EN.md) / [`docs/OUTPUT_RU.md`](docs/OUTPUT_RU.md) for layout details.

## Dtype Inspector

The ``dtype_inspector`` utility executes each ``get_*_data`` helper on a small
set of identifiers and logs the resulting pandas dtypes.  Run this periodically
to spot schema drift when upstream services change their responses.

```bash
python -m library.utils.cli_tools.dtype_inspector_main --log-level INFO
```

Consider wiring the script into CI to detect dtype changes early.  The tool is
lightweight and makes only a handful of requests per dataset.

## Development and testing / Разработка и тестирование

**EN.** Individual tools such as ``ruff``, ``black`` and ``mypy`` are wired into ``pre-commit`` but can be executed manually when iterating on specific modules.

**RU.** Отдельные утилиты ``ruff``, ``black`` и ``mypy`` подключены к ``pre-commit``, но их можно запускать вручную при доработке отдельных модулей.

```bash
ruff check scripts library library/utils/cli_tools
black scripts library library/utils/cli_tools
mypy scripts library
pytest
```

Test datasets live in ``tests/data``; ``library.utils.cli_tools.check_determinism`` validates repeatable CSV output. / Тестовые наборы лежат в ``tests/data``; ``library.utils.cli_tools.check_determinism`` проверяет повторяемость CSV-вывода.

## Лицензия

MIT License. См. файл `LICENSE` (если присутствует).


При обновлении справочных материалов добавляйте их непосредственно в папку `docs`.

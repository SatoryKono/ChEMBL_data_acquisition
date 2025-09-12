# ChEMBL Data Acquisition Utilities

Utilities for downloading and processing biological data from public APIs.
The project demonstrates a typical Python 3.12 data pipeline including
parsing, validation, aggregation and export of tabular data.

## Requirements

- Python 3.12 or later. Older interpreters abort with an informative message.

## Установка

```bash
pip install .[dev]
```

Install the project together with development tools such as
``black``, ``ruff``, ``mypy`` and ``pytest`` as well as testing utilities like
``hypothesis``, ``responses`` and ``psutil``. Sensitive configuration like API
tokens should be stored in a local ``.env`` file – see
[`Конфигурация через .env`](#конфигурация-через-env) for details.

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

4. **Run the tests**

   ```bash

   pytest
   ```


   The suite exercises the library modules using fixtures from
   ``tests/data``.


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
pip install .[dev]
```

### Pre-commit hooks

This project uses [pre-commit](https://pre-commit.com/) to run formatting, linting, type checking and tests. Install the hooks and run them across the codebase:

```bash
pre-commit install
pre-commit run --all-files
```

Continuous integration executes the same checks.



## Project structure

```
data/             Example input and output files
dictionary/       Lookup tables used during processing
library/          Reusable data-processing modules
tests/            Pytest suite and sample datasets
scripts/          Command-line utilities and development helpers
mapper_main.py    Mapping CLI
table_quality_main.py  CSV profiling CLI
config.yaml       Global configuration defaults
```

## Command line interface

Individual scripts provide specialised data retrieval utilities:

* ``scripts/get_activities.py`` – fetch ChEMBL activity information.
* ``scripts/get_assay_data.py`` – retrieve assay descriptions from ChEMBL.
* ``scripts/get_document_data.py`` – gather publication metadata.
* ``scripts/get_target_data.py`` – combine ChEMBL, UniProt and IUPHAR target data.
* ``scripts/get_testitem_data.py`` – download compound data and enrich with PubChem.
* ``scripts/get_input_initialisation.py`` – merge ChEMBL initialisation workbooks.

For a quick connectivity check without writing any files, limit the number of
records and enable dry-run mode:

```bash
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

Environment variables override values from the YAML file. Variables use the
``CHEMBL_DA__SECTION__KEY`` pattern and also support short aliases:

* ``CHEMBL_DA__API__CHEMBL_BASE`` / ``CHEMBL_DA_BASE``
* ``CHEMBL_DA__API__TIMEOUT_CONNECT`` / ``CHEMBL_DA_TIMEOUT_CONNECT``
* ``CHEMBL_DA__API__TIMEOUT_READ`` / ``CHEMBL_DA_TIMEOUT_READ``
* ``CHEMBL_DA__API__RPS`` / ``CHEMBL_DA_RPS``

* ``CHEMBL_DA__OPENALEX__TIMEOUT_CONNECT`` / ``CHEMBL_DA_OPENALEX_TIMEOUT_CONNECT``
* ``CHEMBL_DA__OPENALEX__TIMEOUT_READ`` / ``CHEMBL_DA_OPENALEX_TIMEOUT_READ``
* ``CHEMBL_DA__OPENALEX__RPS`` / ``CHEMBL_DA_OPENALEX_RPS``
* ``CHEMBL_DA__CROSSREF__TIMEOUT_CONNECT`` / ``CHEMBL_DA_CROSSREF_TIMEOUT_CONNECT``
* ``CHEMBL_DA__CROSSREF__TIMEOUT_READ`` / ``CHEMBL_DA_CROSSREF_TIMEOUT_READ``
* ``CHEMBL_DA__CROSSREF__RPS`` / ``CHEMBL_DA_CROSSREF_RPS``
* ``CHEMBL_DA__UNIPROT__TIMEOUT_CONNECT`` / ``CHEMBL_DA_UNIPROT_TIMEOUT_CONNECT``
* ``CHEMBL_DA__UNIPROT__TIMEOUT_READ`` / ``CHEMBL_DA_UNIPROT_TIMEOUT_READ``
* ``CHEMBL_DA__UNIPROT__RPS`` / ``CHEMBL_DA_UNIPROT_RPS``
* ``CHEMBL_DA__IUPHAR__TIMEOUT_CONNECT`` / ``CHEMBL_DA_IUPHAR_TIMEOUT_CONNECT``
* ``CHEMBL_DA__IUPHAR__TIMEOUT_READ`` / ``CHEMBL_DA_IUPHAR_TIMEOUT_READ``
* ``CHEMBL_DA__IUPHAR__RPS`` / ``CHEMBL_DA_IUPHAR_RPS``
* ``CHEMBL_DA__PUBCHEM__TIMEOUT_CONNECT`` / ``CHEMBL_DA_PUBCHEM_TIMEOUT_CONNECT``
* ``CHEMBL_DA__PUBCHEM__TIMEOUT_READ`` / ``CHEMBL_DA_PUBCHEM_TIMEOUT_READ``
* ``CHEMBL_DA__PUBCHEM__RPS`` / ``CHEMBL_DA_PUBCHEM_RPS``

* ``CHEMBL_DA__IO__OUTPUT_DIR`` / ``CHEMBL_DA_OUTDIR``
* ``CHEMBL_DA__JOBS__CONCURRENCY`` / ``CHEMBL_DA_CONCURRENCY``
* ``CHEMBL_DA__JOBS__CHUNK_SIZE`` / ``CHEMBL_DA_CHUNK_SIZE``
* ``CHEMBL_DA__RETRY__MAX_ATTEMPTS`` / ``CHEMBL_DA_RETRY_MAX_ATTEMPTS``
* ``CHEMBL_DA__RETRY__BACKOFF_FACTOR`` / ``CHEMBL_DA_RETRY_BACKOFF_FACTOR``
* ``CHEMBL_DA__LOG__LEVEL`` / ``CHEMBL_DA_LOG_LEVEL``
* ``CHEMBL_DA__LOG__FORMAT`` / ``CHEMBL_DA_LOG_FORMAT``

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


Common flags shared by scripts include:

* ``--input`` – input CSV file (default ``input.csv``)
* ``--output`` – destination CSV file (default: auto-generated next to the input)
* ``--log-level`` – logging verbosity (default ``INFO``)
* ``--sep`` – CSV delimiter (default taken from configuration)
* ``--encoding`` – file encoding (default taken from configuration)
* ``--column`` – column containing identifiers (script specific)

Example fetching assay data::

    python scripts/get_assay_data.py --input assays.csv --output assays_out.csv \
        --column assay_chembl_id

Each command validates required columns before querying external APIs and
writes the resulting table to the specified output file.

## Data contracts

Each output table is validated with ``pandera`` to guarantee a consistent
layout. Columns must satisfy the following contracts.

### Activities

Required columns

* ``activity_id`` (*int*, ``>= 0``)
* ``molecule_chembl_id`` (*str*)
* ``standard_value`` (*float*, ``>= 0``)

Optional columns

* ``target_id`` (*str*)
* ``standard_type`` (*str*, one of ``IC50``, ``EC50``, ``Ki``, ``Kd``)
* ``pA_value`` (*float*, ``0–14``)

Valid row

```csv
activity_id,molecule_chembl_id,target_id,standard_type,standard_value,pA_value
1,TST1,TGT1,IC50,50,9
```

Invalid row (``standard_type`` outside enum, ``pA_value`` > 14)

```csv
activity_id,molecule_chembl_id,target_id,standard_type,standard_value,pA_value
2,TST2,TGT2,IC90,100,20
```

### Assays

Required columns

* ``assay_chembl_id`` (*str*)
* ``document_chembl_id`` (*str*)
* ``year`` (*int*, ``1900–2100``)
* ``month`` (*int*, ``1–12``)

Optional columns

* ``target_chembl_id`` (*str*)

Valid row

```csv
assay_chembl_id,document_chembl_id,target_chembl_id,year,month
A1,D1,T1,2023,5
```

Invalid row (``month`` > 12)

```csv
assay_chembl_id,document_chembl_id,target_chembl_id,year,month
A2,D2,T2,2023,13
```

### Documents

Required columns

* ``document_chembl_id`` (*str*)
* ``title`` (*str*)
* ``year`` (*int*, ``1900–2100``)
* ``month`` (*int*, ``1–12``)

Optional columns

* ``doi`` (*str*)
* ``day`` (*int*, ``1–31``)
* ``citation`` (*int*, ``>= 0``)

Valid row

```csv
document_chembl_id,doi,title,year,month,day,citation
D1,10.1000/test,A study,2022,7,15,3
```

Invalid row (``day`` > 31, ``citation`` < 0)

```csv
document_chembl_id,doi,title,year,month,day,citation
D2,10.1000/test2,Another study,2022,7,45,-1
```

### Targets

Required columns

* ``target_chembl_id`` (*str*)
* ``organism`` (*str*)

Optional columns

* ``target_uniprot_id`` (*str*)
* ``pH_dependence`` (*float*, ``0–14``)
* ``isoforms`` (*float*, ``>= 0``)

Valid row

```csv
target_chembl_id,organism,target_uniprot_id,pH_dependence,isoforms
T1,Homo sapiens,P12345,7.4,2
```

Invalid row (``pH_dependence`` > 14)

```csv
target_chembl_id,organism,target_uniprot_id,pH_dependence,isoforms
T2,Mus musculus,P67890,15,1
```

### Testitems

Required columns

* ``salt_chembl_id`` (*str*)
* ``molecule_chembl_id`` (*str*)
* ``molecule_type`` (*str*, ``Small molecule``, ``Biopolymer``,
  ``Oligosaccharide``, ``Unknown``)
* ``mw_freebase`` (*float*, ``0–2000``)

Optional columns

* ``chirality`` (*int*, ``-1``, ``0``, ``1``, ``2``)
* ``num_ro5_violations`` (*float*, ``0–5``)
* ``is_radical`` (*bool*)

Valid row

```csv
salt_chembl_id,molecule_chembl_id,molecule_type,chirality,mw_freebase,num_ro5_violations,is_radical
S1,M1,Small molecule,1,350.5,0,false
```

Invalid row (``molecule_type`` outside enum, ``mw_freebase`` > 2000)

```csv
salt_chembl_id,molecule_chembl_id,molecule_type,chirality,mw_freebase,num_ro5_violations,is_radical
S2,M2,Peptide,0,2500,1,true
```

## Configuration

Default settings such as API endpoints, network timeouts, rate limits and
output directories live in `config.yaml` at the repository root. Each field is
documented inside the file and has a sensible fallback that the utilities use
if the entry is missing. See `docs/CONFIG.md` for a detailed description of all
available options.

Example merging initialisation tables::

    python scripts/get_input_initialisation.py --config config.yaml

The ``same_doc`` and ``all_doc`` workbook paths default to values from
``config.yaml`` but can be overridden on the command line::

    python scripts/get_input_initialisation.py \
      --same-doc path/to/ChEMBL_same_document_20_05.xlsx \
      --all-doc  path/to/ChEMBL_all_10_05_step5.xlsx \
      --out-dir  ./out

The script also profiles each exported table and writes
``<name>_quality_report_table.csv`` and
``<name>_data_correlation_report_table.csv`` alongside the original CSVs
in ``--out-dir``.

## Development

Install the optional developer tools and run the standard quality checks:

* ``pre-commit`` – run formatting, linting, type checking and tests in one go::

      pre-commit run --all-files

* ``black`` – auto-format the code::

      black scripts library mapper_main.py table_quality_main.py

* ``ruff`` – lint the project::

      ruff check scripts library mapper_main.py table_quality_main.py

* ``mypy`` – perform static type checks::

      mypy scripts library mapper_main.py table_quality_main.py

* ``python scripts/check_determinism.py`` – verify deterministic CSV output.

* ``pytest`` – run the unit tests.

Test data live in ``tests/data`` and provide coverage for utility functions in
the library modules.

## FAQ

- Где искать информацию о типичных ошибках конфигурации? См. [Ошибки конфигурации](#ошибки-конфигурации).

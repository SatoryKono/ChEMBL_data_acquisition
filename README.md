# ChEMBL Data Acquisition Utilities

The README is available in multiple languages:

| Language | Link |
|----------|------|
| English  | [README_EN.md](README_EN.md) |
| Русский  | [README_RU.md](README_RU.md) |

 
* Командные скрипты с унифицированными флагами `--input`, `--output`

  (основной путь назначения; `--final-out` пока доступен только в
  `scripts.get_target_data` и `library.utils.cli_tools.pipeline_targets_main`,
  а устаревшие алиасы `--out` сохраняются для совместимости и теперь
  сопровождаются предупреждением),
  `--log-level`, `--sep`, `--encoding`, `--column`, а также `--config` и
  `--print-config` для управления загрузкой настроек. Размер пакетной
  выборки задаётся параметрами `--chunk-size` или `--batch-size` в зависимости
  от конкретного пайплайна. Новые переключатели `--raw-out`, `--raw-format`,
  `--id-cols`, `--no-reindex-raw`, а также пара `--normalize-at-export` /
  `--no-normalize-at-export` уже доступны в пайплайне таргетов; остальные
  команды получат их после расширения общего CLI.

* Потоковая обработка больших CSV через чанки, детерминированный вывод.
* Отдельный «сырой» снимок (`--raw-out`) доступен для пайплайна таргетов;
  поддержка остальных конвейеров запланирована и будет объявлена отдельно.
  По умолчанию «сырой» дамп переупорядочивает столбцы для стабильности
  (`--no-reindex-raw` отключает это), а итоговая таблица нормализуется перед
  записью (`--no-normalize-at-export` сохраняет исходное состояние данных).
* Валидаторы схем (`schemas/`) и словари (`dictionary/`) для проверки
  типов, диапазонов и справочников.
* Конфигурация через `config/config.yaml`, переменные окружения и ключи CLI.
* Логирование через стандартный модуль `logging` с настраиваемым уровнем.
* Полная статическая типизация (PEP 484), линтинг `ruff`, форматирование
  `black`, проверка типов `mypy`, юнит‑тесты `pytest`.

## Console entry points / Консольные точки входа

**EN.** Installing the project via `pip install .` registers dedicated
console scripts for each pipeline. Use them interchangeably with the
`python -m …` invocations shown throughout the documentation.

**RU.** После установки `pip install .` для каждого пайплайна становятся
доступны консольные команды. Их можно использовать вместо вызовов
`python -m …`, приведённых в примерах ниже.

| Console script | Module equivalent | Назначение |
| -------------- | ----------------- | ---------- |
| `get-data` | `python -m scripts.get_data` | Orchestrates the full export / Оркестратор всех этапов |
| `get-activity-data` | `python -m scripts.get_activity_data` | Activity data export / Выгрузка активностей |
| `get-assay-data` | `python -m scripts.get_assay_data` | Assay metadata / Метаданные ассайев |
| `get-document-data` | `python -m scripts.get_document_data` | Document metadata / Метаданные документов |
| `get-document-type` | `python -m library.utils.cli_tools.get_document_type` | Document type classification / Классификация документов по типам публикаций |
| `get-target-data` | `python -m scripts.get_target_data` | Target aggregation / Агрегация таргетов |
| `get-testitem-data` | `python -m scripts.get_testitem_data` | Test item enrichment / Обогащение тест-объектов |
| `csv-utils` | `python -m library.utils.cli_tools.csv_utils_main` | CSV helpers / Утилиты работы с CSV |
| `mapper` | `python -m library.utils.cli_tools.mapper_main` | Identifier mapping / Маппинг идентификаторов |
| `table-quality` | `python -m library.utils.cli_tools.table_quality_main` | Quality reports / Отчёты качества |
| `chunk-io` | `python -m library.utils.cli_tools.chunk_io_main` | Chunked IO harness / Обвязка чтения чанков |
| `get-input-initialisation` | `python -m library.utils.cli_tools.get_input_initialisation` | Combine init workbooks into pair/entity tables / Объединение Excel инициализации в пары и сущности |
| `get-activities` | `python -m library.utils.cli_tools.get_activities` | Smoke logger / Демонстрация логов |
| `check-determinism` | `python -m library.utils.cli_tools.check_determinism` | Determinism checks / Проверка детерминизма |

## Требования

| Компонент     | Минимальная версия | Последняя протестированная |
|---------------|--------------------|-----------------------------|
| Python        | ≥3.11              | 3.12                        |
| numpy         | ≥1.26              | 2.3.3                       |
| pandas        | ≥2.0               | 2.3.3                       |
| requests      | ≥2.31              | 2.32.5                      |
| PyYAML        | ≥6.0               | 6.0.3                       |

**EN.** `requirements-lock.txt` is the single source of truth for pinned
dependencies; `pyproject.toml` documents the declared dependency ranges.
**RU.** Файл `requirements-lock.txt` — единый источник истины с фиксированными
версиями, а `pyproject.toml` описывает поддерживаемые диапазоны зависимостей.

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
pip install -r requirements-lock.txt
pre-commit install
```

**EN.** Installing from the lock file guarantees that local development and
continuous integration use the exact same dependency set. Regenerate the lock
after modifying `pyproject.toml` by creating a fresh virtual environment,
running `pip install .[dev]` and freezing the result with
`pip freeze > requirements-lock.txt`. /
**RU.** Установка по lock-файлу гарантирует идентичные зависимости в локальной
разработке и CI. После изменений в `pyproject.toml` пересоберите lock: создайте
новое виртуальное окружение, выполните `pip install .[dev]`, а затем
`pip freeze > requirements-lock.txt`.

> **EN.** Fresh wheel installs now rely on platform-specific user directories.
> The packaged configuration is copied to the user config home, CSV exports go
> to the user data directory and HTTP caches live in the user cache directory.
> Adjust these paths via the `local.io.*` keys if required or keep the defaults
> listed below. See also the [FAQ entry on wheel vs. source usage](#faq-wheel-vs-source).
>
> **RU.** Новые установки wheel используют каталоги пользователя, подобранные
> `platformdirs`: конфигурация копируется в пользовательский каталог настроек,
> выгрузки — в пользовательский каталог данных, кэши HTTP — в пользовательский
> каталог кэша. При необходимости меняйте их через `local.io.*` или оставляйте
> значения по умолчанию из таблицы. Дополнительные детали собраны в
> [FAQ про wheel и исходники](#faq-wheel-vs-source).

| Platform / Платформа | Config dir / Каталог конфигурации | Output dir / Каталог выгрузки | Cache dir / Каталог кэша |
| --------------------- | --------------------------------- | ----------------------------- | ------------------------ |
| Linux (XDG)           | `~/.config/chembl-data-acquisition/config.yaml` | `~/.local/share/chembl-data-acquisition/output` | `~/.cache/chembl-data-acquisition` |
| macOS                 | `~/Library/Application Support/chembl-data-acquisition/config.yaml` | `~/Library/Application Support/chembl-data-acquisition/output` | `~/Library/Caches/chembl-data-acquisition` |
| Windows               | `%APPDATA%\chembl-data-acquisition\config.yaml` | `%LOCALAPPDATA%\chembl-data-acquisition\output` | `%LOCALAPPDATA%\chembl-data-acquisition\Cache` |

Sensitive configuration such as API tokens belongs in a local ``.env`` file – see [`Конфигурация через .env`](#конфигурация-через-env) for usage guidelines.

## Quick Start / Быстрый старт

1. **Install dependencies / Установите зависимости** – follow the steps in [Installation / Установка](#installation--установка).

2. **Initialise pre-commit hooks / Активируйте pre-commit-хуки**

   ```bash
   pre-commit install
   ```

   EN: Git hooks ensure formatting, linting, static type checks and tests run before each commit.

   RU: Git-хуки автоматически запускают форматирование, линтеры, проверки типов и тесты перед каждым коммитом.

3. **Inspect the unified orchestrator / Изучите унифицированный оркестратор**

   ```bash
   get-data --help
   ```

   EN: The ``get-data`` console script wires every pipeline together and
   exposes common flags for configuration paths, input/output directories and
   logging. Reviewing the help output confirms that the CLI entry point is
   installed correctly. Once comfortable with the options you can launch a full
   export by providing real directories.

   RU: Консольная команда ``get-data`` объединяет все конвейеры и предоставляет
   общие параметры для путей конфигурации, входных/выходных каталогов и
   логирования. Просмотр справки подтверждает корректную установку CLI. После
   ознакомления с опциями можно запускать полный экспорт, указав реальные
   каталоги.

   EN: Without ``--config`` the orchestrator now falls back to the packaged
   ``config/config.yaml`` via ``library.utils.config.DEFAULT_CONFIG_PATH``. Pass
   an explicit path whenever you maintain a local override.

   RU: При запуске без ``--config`` оркестратор использует встроенный
   ``config/config.yaml`` через ``library.utils.config.DEFAULT_CONFIG_PATH``.
   Собственный YAML укажите явным путём.

   For lightweight smoke checks you can still call individual helpers, for
   example:

  ```bash
  python -m library.utils.cli_tools.mapper_main --input tests/data/chembl_targets_min.csv \
      --column target_chembl_id --output out/targets_mapped.csv --log-level DEBUG
  python -m library.utils.cli_tools.table_quality_main --input tests/data/chembl_targets_min.csv \
      --output out/quality --table-name chembl_targets --log-level INFO
  ```

   EN: Maintain a custom configuration by copying the packaged YAML into the
   user config directory and passing it via `--config`. The loader automatically
   merges a sibling `config.local.yaml`, keeping your overrides separate from the
   upstream defaults. Example:

  ```bash
  python - <<'PY'
from importlib import resources
from pathlib import Path

target = Path.home() / ".config" / "chembl-data-acquisition"
target.mkdir(parents=True, exist_ok=True)
source = resources.files("chembl_data_acquisition.config") / "config.yaml"
target_cfg = target / "config.yaml"
target_cfg.write_bytes(source.read_bytes())
PY

  get-data --config ~/.config/chembl-data-acquisition/config.yaml \
      --output-dir ~/.local/share/chembl-data-acquisition/output
  ```

   RU: Чтобы сопровождать собственную конфигурацию, скопируйте штатный YAML в
   пользовательский каталог настроек и передайте путь через `--config`.
   Файл `config.local.yaml`, лежащий рядом, автоматически объединяется с базой,
   поэтому изменения остаются отделены от стандартов. Пример выше можно
   использовать и в Linux/macOS; в PowerShell путь будет
   `%APPDATA%\chembl-data-acquisition\config.yaml`.

  EN: In the reporting example above `--output` sets the destination. `--final-out`
  currently exists only in `scripts.get_target_data` and
  `library.utils.cli_tools.pipeline_targets_main`. The legacy `--out` alias remains
  available but now raises a deprecation warning.

  RU: В примере с отчётностью каталог задаёт `--output`. Флаг `--final-out`
  реализован лишь в `scripts.get_target_data` и
  `library.utils.cli_tools.pipeline_targets_main`. Устаревший алиас `--out`
  сохраняется, но сопровождается предупреждением об удалении.

4. **Run the tests / Запустите тесты** – refer to [Tests / Тесты](#tests--тесты).


## Staged export pipeline / Поэтапный конвейер

```mermaid
flowchart LR
  Fetch --> Raw["Raw CSV / Parquet"] --> Cleanup["Cleanup IDs / Очистка ID"] --> Normalize --> Validate --> Final["Final export / Финальный экспорт"]
```


**EN.** The target pipeline already follows the staged contract with dedicated destinations for raw and cleaned artefacts. Use
`--raw-out` (optionally with `--raw-format parquet`) to capture the raw payload, list composite keys via `--id-cols`, and direct
the cleaned export to `--final-out`. Raw dumps reindex columns alphabetically for deterministic layouts unless
`--no-reindex-raw` keeps the API order. The final CSV is normalized by default; flip the boolean pair
`--normalize-at-export` / `--no-normalize-at-export` when you need the final artefact to mirror the raw payload byte-for-byte.
Placeholder identifiers remain in the raw snapshot and are counted in the metadata (`error_placeholder_counts`), while the
normalized export includes only validated values. Other pipelines will surface `--final-out` once their shared CLI is extended;
until then the deprecated aliases remain available with matching warnings.

**RU.** Пайплайн таргетов уже использует поэтапный контракт с разделением «сырого» и нормализованного вывода. Флаг `--raw-out`
(при необходимости с `--raw-format parquet`) сохраняет исходный ответ, `--id-cols` перечисляет составные ключи, а чистый экспорт
направляется в `--final-out`. По умолчанию «сырой» дамп переупорядочивает столбцы в алфавитном порядке — флаг
`--no-reindex-raw` сохраняет исходную раскладку. Финальный CSV нормализуется автоматически; переключение пары
`--normalize-at-export` / `--no-normalize-at-export` позволяет либо выполнить нормализацию непосредственно перед записью, либо
сохранить артефакт идентичным «сырому» снимку. Временные идентификаторы остаются в «сыром» файле и учитываются в метаданных
(`error_placeholder_counts`), тогда как нормализованный экспорт содержит только прошедшие валидацию значения. Прочие пайплайны
получат `--final-out` после расширения общего CLI; до тех пор устаревшие алиасы остаются доступными и сопровождаются теми же
предупреждениями.

> **EN.** `--raw-out`, `--final-out`, `--raw-format`, `--id-cols`, `--no-reindex-raw`, and the boolean pair
> `--normalize-at-export` / `--no-normalize-at-export` are currently exposed via `get-target-data` and
> `library.utils.cli_tools.pipeline_targets_main`. Other commands will adopt these switches once the shared parser lands.
> **RU.** Флаги `--raw-out`, `--final-out`, `--raw-format`, `--id-cols`, `--no-reindex-raw`, а также пара
> `--normalize-at-export` / `--no-normalize-at-export` уже доступны в `get-target-data` и
> `library.utils.cli_tools.pipeline_targets_main`. Остальные команды получат их после доработки общего парсера.



## Tests / Тесты

**EN.** The `pre-commit` suite runs formatting, linting and static type checks. Execute `pytest` for unit tests and add coverage flags when required. Determinism and smoke checks are available through dedicated CLI helpers. The canonical checklist lives in the QA process documents listed below.

**RU.** Команда `pre-commit` запускает форматирование, линтеры и проверку типов. Для юнит-тестов используйте `pytest`, при необходимости добавляйте параметры покрытия. Детеминизм и smoke-проверки доступны в отдельных CLI. Актуальный список проверок приведён в документах по процессу обеспечения качества ниже.

| Language / Язык | Checklist / Чек-лист |
|-----------------|----------------------|
| English         | [docs/QA_PROCESS_EN.md](docs/QA_PROCESS_EN.md) |
| Русский         | [docs/QA_PROCESS_RU.md](docs/QA_PROCESS_RU.md) |


```bash
pre-commit run --all-files
pytest
pytest --cov=library --cov=scripts --cov-report=term-missing --cov-report=xml
python -m scripts.get_data --help
tmp_dir=$(mktemp -d) && python -m library.utils.cli_tools.pipeline_targets_main \
    --input tests/data/chembl_targets_min.csv \
    --final-out "${tmp_dir}/targets.csv" --log-level INFO --limit 2
python -m library.utils.cli_tools.check_determinism --log-level DEBUG
python -m library.utils.cli_tools.mapper_batch_main --input chembl_ids.csv \
    --output out/mapped.csv --log-level INFO
```

Before running the smoke command, create a `chembl_ids.csv` file with a header `chembl_id` and the required identifiers. / Перед запуском smoke-команды создайте `chembl_ids.csv` со столбцом `chembl_id` и нужными идентификаторами.


## Генерация данных

**EN.** Five production pipelines live in `scripts/` and write CSV outputs to
`~/.local/share/chembl-da/output` by default (honouring `--base-path`). /
**RU.** Пять рабочих пайплайнов расположены в каталоге `scripts/` и по
умолчанию сохраняют CSV-файлы в `~/.local/share/chembl-da/output`
(учитывается `--base-path`):

* **EN.** `get_activity_data.py` extracts activity records from ChEMBL and
  adds calculated bounds. / **RU.** `get_activity_data.py` извлекает данные
  активностей из ChEMBL и дополняет их расчётными границами значений.
* **EN.** `get_assay_data.py` exports assay descriptions. / **RU.**
  `get_assay_data.py` выгружает описания ассайев.
* **EN.** `get_document_data.py` merges publication metadata from ChEMBL and
  partner aggregators (PubMed, Semantic Scholar, OpenAlex, Crossref). /
  **RU.** `get_document_data.py` объединяет метаданные публикаций из ChEMBL и
  агрегаторов (PubMed, Semantic Scholar, OpenAlex, Crossref).
* **EN.** `get_target_data.py` collects target information from ChEMBL,
  UniProt and IUPHAR. / **RU.** `get_target_data.py` собирает информацию о
  таргетах из ChEMBL, UniProt и IUPHAR.
* **EN.** `get_testitem_data.py` enriches compounds with structural attributes
  and PubChem data. / **RU.** `get_testitem_data.py` обогащает соединения
  структурными атрибутами и данными PubChem.

### Управляющий скрипт get_data

**EN.** `get_data.py` orchestrates the five production pipelines, forwarding a
single set of CLI options to each step. / **RU.** `get_data.py` запускает все
пять пайплайнов последовательно, прокидывая единые параметры командной строки
в каждый модуль.

Common options / Общие параметры:

* `--base-path` — **EN.** base directory for inputs and outputs. / **RU.**
  корневой каталог данных.
* `--input-dir`, `--output-dir` — **EN.** sub-directories for source CSV files
  and generated artefacts. / **RU.** подкаталоги с исходными CSV и
  выгруженными файлами.
* `--config` — **EN.** path to the shared YAML configuration. / **RU.** путь к
  общей YAML-конфигурации.
* `--date` — **EN.** date or prefix used in output filenames. / **RU.** дата
  или префикс, добавляемый к именам выгрузок.
* `--log-level`, `--force`, `--skip-existing` — **EN.** logging verbosity and
  overwrite policy. / **RU.** уровень логирования и политика перезаписи.
* `--limit` — **EN.** cap forwarded to every pipeline; ``0`` skips execution
  without touching inputs or outputs. / **RU.** ограничение для всех шагов;
  значение ``0`` пропускает выполнение без обращений к входным данным и
  файловой системе.

Setting ``--limit 0`` logs ``pipeline_skip_limit`` for every stage and leaves
existing artefacts untouched.

Example / Пример запуска:

```bash
python -m scripts.get_data \
    --base-path data \
    --input-dir input \
    --output-dir output \
    --config config/config.yaml \
    --date 20240101 \
    --log-level INFO
```

The command reads files such as `data/input/document.csv`, writes outputs like
`data/output/output.documents_20240101.csv` and stops if any step fails. / Команда
использует входы вида `data/input/document.csv`, сохраняет результаты в
файлы `data/output/output.documents_20240101.csv` и прерывает выполнение при ошибке
любого шага.

**EN.** `library.utils.cli_tools.pipeline_targets_main` is a cached harness
that reuses the CLI contract of the production target pipeline while working
solely with local files and prepared identifier batches. / **RU.**
`library.utils.cli_tools.pipeline_targets_main` — кешируемая обвязка, которая
использует те же CLI-параметры, что и боевой таргет-пайплайн, но работает
только с локальными файлами и подготовленными чанками идентификаторов без
сетевых вызовов.

Пример полноценного пайплайна:

```bash
python -m scripts.get_activity_data --input tests/data/activity_ids_small.csv \
    --output data/output/activities.csv --limit 10 --log-level INFO
```

Команда извлекает данные из API ChEMBL, сохраняет таблицу и сопутствующий
`*.meta.yaml`. Утилиты разработки и отладки перенесены в
`library/utils/cli_tools/`, например модуль `get_activities` предназначен
только для демонстрационного логирования и не выполняет файловых операций.
См. [docs/CLI_TOOLS.md](docs/CLI_TOOLS.md) (English) и
[docs/CLI_TOOLS_RU.md](docs/CLI_TOOLS_RU.md) (Русский) для кратких описаний и
типовых команд. Каталог с результатами игнорируется Git и автоматически публикуется
как артефакт CI.

> **Примечание.** Ранее эта функциональность была доступна через
> `activity_extraction_main.py`. Теперь используйте модульный запуск
> `python -m scripts.get_activity_data`, что упрощает сопровождение и
> работу в виртуальных окружениях.

## Usage

The examples below illustrate how to run the main CLI tools with common
options like ``--input``, ``--output`` (primary destination flag for most
commands; ``--final-out`` is currently restricted to ``scripts.get_target_data``
and ``library.utils.cli_tools.pipeline_targets_main``) and ``--limit``. The
compatibility alias ``--out`` stays available for now but triggers explicit
deprecation warnings. Passing ``--limit 0`` short-circuits processing before any
network or filesystem access, which is handy for configuration smoke tests. The
target pipeline already exposes ``--raw-out``, ``--final-out``, ``--raw-format``
and ``--id-cols``; other commands will gain the staging switches once the shared
CLI is extended.

### scripts/get_document_data.py

Retrieve document metadata for a list of PubMed IDs using the bundled
sample file:

```bash
python -m scripts.get_document_data pubmed \
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
python -m scripts.get_target_data chembl \
    --input path/to/targets.csv \
    --final-out out/targets.csv \
    --limit 5 \
    --log-level INFO
```

Replace ``path/to/targets.csv`` with a CSV containing a ``target_chembl_id``
column.

The input and output both use ``target_chembl_id`` to align with
validation schemas.

### library.utils.cli_tools.pipeline_targets_main

Exercise the chunking and batch configuration used by the production
target pipeline without contacting remote services:

```bash
python -m library.utils.cli_tools.pipeline_targets_main \
    --input tests/data/chembl_targets_min.csv \
    --final-out out/targets_cached.csv \
    --chunk-size 25 \
    --batch-size 25 \
    --limit 100
```

The command reads target identifiers from the CSV, chunks them according
to ``--chunk-size`` and ``--limit``, forwards the batch size to
``library.pipeline_targets.run_pipeline`` and writes the cached ChemBL
table with pipeline metadata via ``write_csv``. Use it to verify CLI
overrides, logging and deterministic output before launching the network
backed ``get_target_data`` pipeline.

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
pip install -r requirements-lock.txt --upgrade
pre-commit run --all-files
```

The upgrade reinstalls every pinned dependency at the versions recorded in the
lock file and surfaces conflicts immediately. When intentionally moving to
newer releases, update `pyproject.toml`, recreate the lock as described in the
installation section and commit the refreshed `requirements-lock.txt` together
with the source changes.

## Конфигурация через `.env`

Часть параметров утилит можно задавать через переменные окружения.
Чтобы не экспортировать их вручную при каждом запуске, поместите пары
``NAME=value`` в файл `.env` и загрузите их с помощью пакета
[`python-dotenv`](https://pypi.org/project/python-dotenv/).

Пример файла:

```dotenv
CHEMBL_DA_LOG_LEVEL=INFO
CHEMBL_DA_BASE=https://www.ebi.ac.uk/chembl/api/data
```

**EN.** Use either the short alias `CHEMBL_DA_BASE` or the fully qualified
`CHEMBL_DA__SOURCES__CHEMBL__API__CHEMBL_BASE`; both expand to the same setting.
Refer to the [alias table](library/config.py#L1531-L1634) in
`library/config.py` for all supported mappings. / **RU.** Допустимо указывать
как короткий алиас `CHEMBL_DA_BASE`, так и полное имя
`CHEMBL_DA__SOURCES__CHEMBL__API__CHEMBL_BASE` — обе переменные настраивают
одно и то же значение. Полный список алиасов приведён в
[`library/config.py`](library/config.py#L1531-L1634).

См. также файл `.env.example` с типовыми переменными для контактных e-mail.

Запустить скрипт с автоматической подгрузкой настроек можно так:

```bash
python -m dotenv run -- python -m scripts.get_assay_data --input assay_ids.csv \\
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
запросах к API и должен содержать **реальные** контактные данные. Если
оставить шаблон `contact@example.org`, загрузка конфигурации завершится
ошибкой валидации. Значение по умолчанию:

```yaml
api:
  user_agent: "chembl-da/1.0 (mailto:chembl-data@ebi.ac.uk)"
```

Параметр можно переопределить в `config/config.yaml` или через переменную
окружения `CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT`. Перед выводом решения в
продакшен замените `chembl-data@ebi.ac.uk` на собственный рабочий адрес — это
лишь документированное значение по умолчанию. Валидатор по-прежнему отвергает
заполнитель `contact@example.org`, поэтому наличие шаблона блокирует запуск.
Отдельного CLI-флага для `api.user_agent` не предусмотрено (см.
`library/cli/parser.py`), поэтому значение задаётся только через
конфигурационный файл или окружение.


## Валидация конфигурации

`library.config.load_config` проверяет корректность значений в `config/config.yaml`.
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

Некорректные значения в `config/config.yaml` вызывают `ValidationError`. Пример:

```yaml
api:
  rps: -1
```

При загрузке конфигурации:

```python
from library.config import load_config
load_config("config/config.yaml")
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

Диапазоны допустимых значений описаны в [`config.schema.json`](./config.schema.json) — этот файл экспортирован из Pydantic-модели и служит справочным артефактом, где для `api.rps` указан минимум `1`.

## Logging / Логирование

**EN.** CLI helpers configure structured JSON logging via ``library.logging_setup.configure_logger``. Use environment variables or CLI flags to adjust verbosity. The JSON layout is fixed and now stamps the staging phase (`fetch`, `raw`, `cleanup`, `normalize`, `validate`, `final_export`).

**RU.** CLI-хелперы настраивают структурированное JSON-логирование через ``library.logging_setup.configure_logger``. Управляйте уровнем логов переменными окружения или ключами CLI. Формат JSON фиксирован и теперь дополнительно фиксирует стадию (`fetch`, `raw`, `cleanup`, `normalize`, `validate`, `final_export`).


Уровень логов можно задать флагом `--log-level` или переменной
`CHEMBL_DA_LOG_LEVEL`:

```bash
CHEMBL_DA_LOG_LEVEL=DEBUG python -m scripts.get_assay_data --input assay_ids.csv \
    --output out/assays.final.csv
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

Smoke fixtures for full orchestration live in ``tests/data/input-smoke/``. The expected JSON structure (including stage names and placeholder counters) is validated by ``tests/test_logging.py``, ``tests/test_logging_setup.py`` and the smoke harness ``tests/smoke/test_get_data_scripts.py``.

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

Use ``--output`` to redirect these artefacts to another directory. The value must be a directory path
(do not include the final file name). The compatibility alias ``--out`` remains wired in but now raises a
deprecation warning whenever it is invoked. ``--final-out`` is currently exclusive to
``scripts.get_target_data`` and ``library.utils.cli_tools.pipeline_targets_main`` until the shared CLI is
updated for the remaining commands. When ``local.io.exist_ok`` is set to ``false`` the directory has to
exist beforehand; otherwise it is created automatically. The target pipeline uses the additional
``--raw-out``/``--final-out`` switches to separate staging outputs.

All scripts share a common set of flags:

## Configuration


Default settings live in ``config/config.yaml`` and are grouped into three top-level
sections:

* ``sources`` – external services such as ChEMBL, OpenAlex, Crossref, UniProt,
  IUPHAR and PubChem (including pipeline-specific settings).
* ``local`` – file system inputs and outputs, cached resources and workbook
  paths.
* ``system`` – shared concerns such as logging, retry strategy, rate limiting
  and document-type normalisation.

Because the wheel bundles ``config/config.yaml``, any custom configuration file
can simply extend it by overriding the keys you need while omitting the rest.
Create a lightweight YAML such as:

```
sources:
  chembl:
    api:
      rps: 10
system:
  log:
    level: DEBUG
```

Pass this file via ``--config`` (or ``CHEMBL_DA__...`` environment variables)
and the loader will merge your overrides with the packaged defaults.

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
        output_dir: "$CHEMBL_DA_BASE_PATH/output"

### Переменные окружения

Environment variables override values from the YAML file. Names follow the
``CHEMBL_DA__...`` pattern with double underscores separating each nested
section. For example, to enable debug logging:

```bash
export CHEMBL_DA__LOG__LEVEL=DEBUG
```

Most options also provide short aliases for backwards compatibility. The table
lists every supported alias and the canonical key it maps to. See
[`_ALIAS_OVERRIDES`](library/config.py#L1569-L1634) and
[`_ALIAS_MAP`](library/config.py#L1637-L1640) for the authoritative source:

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
| `CHEMBL_DA_LOG_LEVEL` | `CHEMBL_DA__SYSTEM__LOG__LEVEL` |
| `CHEMBL_DA_MOLECULE_CATALOG_CACHE` | `CHEMBL_DA__SOURCES__CHEMBL__MOLECULE_CATALOG__CACHE_PATH` |
| `CHEMBL_DA_OPENALEX_BASE` | `CHEMBL_DA__SOURCES__OPENALEX__BASE` |
| `CHEMBL_DA_OPENALEX_MAILTO` | `CHEMBL_DA__SOURCES__OPENALEX__MAILTO` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_CONNECT` | `CHEMBL_DA__SOURCES__OPENALEX__TIMEOUT_CONNECT` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_READ` | `CHEMBL_DA__SOURCES__OPENALEX__TIMEOUT_READ` |
| `CHEMBL_DA_OUTDIR` | `CHEMBL_DA__LOCAL__IO__OUTPUT_DIR` |
| `CHEMBL_DA_PUBCHEM_BASE` | `CHEMBL_DA__SOURCES__PUBCHEM__BASE` |
| `CHEMBL_DA_PUBCHEM_USER_AGENT` | `CHEMBL_DA__SOURCES__PUBCHEM__USER_AGENT` |
| `CHEMBL_DA_RETRY_BACKOFF_FACTOR` | `CHEMBL_DA__SYSTEM__RETRY__BACKOFF_FACTOR` |
| `CHEMBL_DA_RETRY_MAX_ATTEMPTS` | `CHEMBL_DA__SYSTEM__RETRY__MAX_ATTEMPTS` |
| `CHEMBL_DA_RPS` | `CHEMBL_DA__SOURCES__CHEMBL__API__RPS` |
| `CHEMBL_DA_PUBMED_RPS` | `CHEMBL_DA__SOURCES__PUBMED__RPS` |
| `CHEMBL_DA_PUBMED_BURST` | `CHEMBL_DA__SOURCES__PUBMED__BURST` |
| `CHEMBL_DA_SEMANTIC_SCHOLAR_RPS` | `CHEMBL_DA__SOURCES__SEMANTIC_SCHOLAR__RPS` |
| `CHEMBL_DA_SEMANTIC_SCHOLAR_BURST` | `CHEMBL_DA__SOURCES__SEMANTIC_SCHOLAR__BURST` |
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

Configuration values are validated by :func:`library.config.load_config`, which
calls :meth:`Config.model_validate <pydantic.BaseModel.model_validate>` from
Pydantic. Validation follows the model definitions and produces detailed error
messages for nested fields. The accompanying ``config.schema.json`` is generated
from the same Pydantic model for IDE hints and documentation; it is not
executed at runtime and should not be passed to ``jsonschema``.

Command line flags have the highest priority. All utilities accept ``--config``
to point at a configuration file and ``--print-config`` to show the effective
values after all overrides have been applied. When a ``config.local.yaml`` file
is present next to the primary configuration (including custom paths provided to
``--config``) it is deep-merged after the base YAML to allow per-environment
defaults. The final precedence is::

    YAML < config.local.yaml < environment variables < CLI options

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
``config/config.yaml`` or overrides from the environment and command line are
automatically converted.
 
```bash
# профилирование качества таблицы
python -m library.utils.cli_tools.table_quality_main \
    --input tests/data/activity.csv \
    --table-name activity
```

`--final-out` по умолчанию формируется как `output.<имя_входа>_YYYYMMDD.csv`
в каталоге, заданном `local.io.output_dir`. Устаревшие алиасы
`--output`/`--out` по-прежнему работают, но сопровождаются предупреждением и будут
удалены после миграции всех пайплайнов. Пайплайн таргетов может использовать
`--raw-out` и `--final-out` (с `--raw-format`) для разведения «сырого» снимка и
финального экспорта. Для дополнительных примеров см. [`docs/USAGE_RU.md`](docs/USAGE_RU.md)
и английскую версию [`docs/USAGE_EN.md`](docs/USAGE_EN.md).

## Структура проекта

```
ChEMBL_data_acquisition/
├── config/config.yaml
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

Параметры читаются из `config/config.yaml`, переменных окружения
(`CHEMBL_DA__...`) и ключей CLI.
Подробности в [`docs/CONFIG_RU.md`](docs/CONFIG_RU.md) и английской версии [`docs/CONFIG_EN.md`](docs/CONFIG_EN.md).

## Output and metadata / Вывод и метаданные

**EN.** Pipelines persist deterministic CSV tables via ``library.io.write_csv`` and place accompanying ``*.meta.yaml`` sidecars in ``~/.local/share/chembl-da/output``.

**RU.** Пайплайны сохраняют детерминированные CSV с помощью ``library.io.write_csv`` и добавляют рядом файлы ``*.meta.yaml`` в каталоге ``~/.local/share/chembl-da/output``.

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

## FAQ

<a id="faq-wheel-vs-source"></a>
### Wheel vs. source installation / Установка wheel против исходников

**EN.** Install the published wheel (`pip install chembl-data-acquisition`) when you only need the pipelines: the package ships with the default configuration, dictionary resources and schema validators, and now places writable artefacts under platform-specific user directories (see the table above). Keep long-term tweaks in a sibling `config.local.yaml` or point CLI commands at a cloned copy via `--config`.

**RU.** Используйте опубликованный wheel (`pip install chembl-data-acquisition`), если достаточно готовых пайплайнов: пакет включает конфигурацию по умолчанию, словари и схемы, а артефакты записывает в пользовательские каталоги (см. таблицу выше). Долгосрочные правки сохраняйте в `config.local.yaml` рядом с файлом или передавайте путь к собственной копии через `--config`.

**EN.** Clone the repository when developing new features or debugging: editable installs keep tests, documentation sources and tooling such as `pre-commit` at hand, and you can pin dependencies via `requirements-lock.txt`.

**RU.** Клонируйте репозиторий для разработки и отладки: так остаются доступны тесты, исходники документации и инструменты (`pre-commit` и др.), а зависимости фиксируются через `requirements-lock.txt`.

## Лицензия

MIT License. См. файл `LICENSE` (если присутствует).


При обновлении справочных материалов добавляйте их непосредственно в папку `docs`.
 
Additional reference materials live in the [docs/](docs/) directory.
 

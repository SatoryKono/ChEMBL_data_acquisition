# Утилиты ChEMBL Data Acquisition

Основная справка и расширенная документация расположены в каталоге [docs/](docs/).

## Особенности

* Унифицированные CLI-флаги `--input`, `--output`, `--log-level`, `--sep`, `--encoding`, `--column`, а также `--config` и
  `--print-config` для управления конфигурацией. Размер партий задаётся параметрами `--chunk-size` или `--batch-size`
  в зависимости от пайплайна.
* Потоковая обработка крупных CSV-файлов с детерминированным выводом.
* Валидаторы схем в [`schemas/`](schemas/) и словари в [`dictionary/`](dictionary/) для проверки типов, диапазонов и справочных данных.
* Конфигурация через `config/config.yaml`, переменные окружения и CLI-переопределения.
* Логирование на базе стандартного модуля `logging` с настраиваемой детализацией.
* Полная статическая типизация (PEP 484), линтинг `ruff`, форматирование `black`, проверка типов `mypy` и тесты `pytest`.

## Требования

| Компонент | Минимально поддерживаемая | Последняя протестированная |
|-----------|---------------------------|----------------------------|
| Python    | 3.11                      | 3.12                       |
| numpy     | 2.3.3                     | 2.3.3                      |
| pandas    | 2.3.2                     | 2.3.2                      |
| requests  | 2.32.5                    | 2.32.5                     |
| PyYAML    | 6.0.2                     | 6.0.2                      |

Полный список доступен в `requirements-dev.txt` или `pyproject.toml`. Рабочие зависимости зафиксированы с помощью совместимых
диапазонов, поэтому обновления исправлений внутри минорной версии поддерживаются. CI прогоняет проверки на минимальной и
актуальной версиях из таблицы.

### Среда выполнения

* ОС Linux или macOS с доступом к Bash/PowerShell (для Windows — WSL2).
* Актуальные версии `pip`, `setuptools` и `wheel`; следуйте шагам из раздела [Установка](#установка).
* Доступ к сетевым API ChEMBL/PubChem/UniProt (порт 443).

## Установка

### Подготовка среды

* Обновите инструменты распространения пакетов перед установкой проекта.

  ```bash
  python -m pip install --upgrade pip setuptools wheel
  ```

* Создайте изолированное виртуальное окружение, чтобы избежать конфликтов зависимостей.

  ```bash
  python -m venv .venv
  source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
  ```

### Шаги

Клонируйте репозиторий, перейдите в каталог и установите пакет с dev-зависимостями. Затем активируйте pre-commit-хуки, чтобы
форматирование, линтинг, проверки типов и юнит-тесты запускались автоматически.

```bash
git clone https://github.com/<org>/ChEMBL_data_acquisition.git
cd ChEMBL_data_acquisition
pip install .[dev]
pre-commit install
```

Чувствительные настройки (например, API-токены) храните в локальном `.env` — подробности см. в разделе [Конфигурация через `.env`](#конфигурация-через-env).

## Быстрый старт

1. **Установите зависимости** — выполните шаги из раздела [Установка](#установка).
2. **Активируйте pre-commit-хуки**

   ```bash
   pre-commit install
   ```

   Git-хуки гарантируют запуск форматирования, линтеров, статических проверок и тестов перед каждым коммитом.

3. **Запустите демонстрационный скрипт**

   ```bash
   python -m library.utils.cli_tools.get_activities --limit 10 --log-level INFO
   ```

   Лёгкий хелпер выводит структурированные логи о фиктивных строках активностей и не выполняет чтение/запись файлов. Используйте его для проверки настроек логирования и обвязки CLI перед запуском полноценных пайплайнов. Распространённые флаги: `--limit` для ограничения записей, `--log-level` для детализации, `--sep` для разделителей CSV и `--encoding` для кодировок. Для полной выгрузки запустите один из рабочих пайплайнов, например:

   ```bash
   python -m library.utils.cli_tools.mapper_main --input tests/data/chembl_targets_min.csv \
       --column target_chembl_id --output out/targets_mapped.csv --log-level DEBUG
   python -m library.utils.cli_tools.table_quality_main --input tests/data/chembl_targets_min.csv \
       --output out/quality --table-name chembl_targets --log-level INFO
   ```

   Во втором примере аргумент `--output` должен указывать на каталог, куда будут сохранены файлы отчёта.

4. **Запустите тесты** — см. раздел [Тесты](#тесты).

## Тесты

`pre-commit` выполняет форматирование, линтинг и статические проверки. Для юнит-тестов используйте `pytest`, при необходимости добавляйте параметры покрытия. Дополнительные smoke- и детерминизм-проверки доступны через отдельные CLI.

```bash
pre-commit run --all-files
pip check
pytest
pytest --cov=library --cov=scripts --cov-report=term-missing --cov-report=xml
python -m library.utils.cli_tools.check_determinism --log-level DEBUG
python -m library.utils.cli_tools.mapper_batch_main --input chembl_ids.csv \
    --output out/mapped.csv --log-level INFO
```

Перед запуском smoke-команды создайте `chembl_ids.csv` с заголовком `chembl_id` и нужными идентификаторами.

## Генерация данных

Шесть основных пайплайнов из каталога [`scripts/`](scripts/) формируют CSV-файлы и сохраняют их в `data/output/`:

* `get_activity_data.py` — выгружает активности из ChEMBL и обогащает расчётными границами значений.
* `get_assay_data.py` — загружает описания ассайев.
* `get_document_data.py` — объединяет метаданные публикаций из ChEMBL и агрегаторов (PubMed, Semantic Scholar, OpenAlex, Crossref).
* `get_target_data.py` — собирает информацию о таргетах из ChEMBL, UniProt и IUPHAR.
* `get_testitem_data.py` — дополняет соединения структурными атрибутами и данными PubChem.
* `library.utils.cli_tools.pipeline_targets_main` — облегчённая обёртка над `library.pipeline_targets.run_pipeline`, использующая те же CLI-параметры, что и боевой таргет-пайплайн, но работающая только с локальными файлами и подготовленными чанками идентификаторов без сетевых вызовов.

Пример запуска полного пайплайна:

```bash
python -m scripts.get_activity_data --input tests/data/activity_ids_small.csv \
    --output data/output/activities.csv --limit 10 --log-level INFO
```

Команда обращается к API ChEMBL, сохраняет таблицу и сопровождающий `*.meta.yaml`. Утилиты разработки находятся в `library/utils/cli_tools/`; например, модуль `get_activities` предназначен лишь для демонстрационного логирования и не выполняет файловых операций. См. [`docs/CLI_TOOLS.md`](docs/CLI_TOOLS.md) для кратких описаний и типовых команд. Каталог результатов игнорируется Git и публикуется как артефакт CI.

> **Примечание.** Ранее использовался входной скрипт `activity_extraction_main.py`. Теперь применяйте модульный запуск `python -m scripts.get_activity_data`, что упрощает сопровождение и работу в виртуальных окружениях.

## Использование

Ниже приведены примеры запуска основных CLI-инструментов с типовыми флагами (`--input`, `--output`, `--limit`).

### `scripts/get_document_data.py`

Получите метаданные публикаций по списку PubMed ID с использованием тестового файла:

```bash
python -m scripts.get_document_data pubmed \
    --input tests/data/pmids.csv \
    --output out/documents.csv \
    --limit 5 \
    --log-level INFO
```

Файл `tests/data/pmids.csv` содержит небольшой набор PMIDs для экспериментов.

Запустить PubMed-пайплайн можно и напрямую через модуль библиотеки:

```bash
python -m library.pubmed_library \
    --input-csv tests/data/pmids.csv \
    --output out/documents.csv \
    --log-level INFO
```

### `scripts/get_target_data.py`

Получите базовую информацию о таргетах из ChEMBL:

```bash
python -m scripts.get_target_data chembl \
    --input path/to/targets.csv \
    --output out/targets.csv \
    --limit 5 \
    --log-level INFO
```

Замените `path/to/targets.csv` на CSV со столбцом `target_chembl_id`. Вход и выход используют одноимённый столбец для соответствия схемам валидации.

### `library.utils.cli_tools.pipeline_targets_main`

Проверьте настройку чанков и батчей, используемую боевым пайплайном таргетов, без обращений к внешним сервисам:

```bash
python -m library.utils.cli_tools.pipeline_targets_main \
    --input tests/data/chembl_targets_min.csv \
    --output out/targets_cached.csv \
    --chunk-size 25 \
    --batch-size 25 \
    --limit 100
```

Команда считывает идентификаторы таргетов из CSV, разбивает их по `--chunk-size` и `--limit`, передаёт размер батча в `library.pipeline_targets.run_pipeline` и сохраняет кешированную таблицу ChemBL вместе с метаданными через `write_csv`. Используйте её для проверки переопределений CLI, логирования и детерминированного вывода перед запуском сетевого пайплайна `get_target_data`.

### `library/utils/cli_tools/get_activities.py`

Сгенерируйте фиктивные активности без обращений к внешним службам:

```bash
python -m library.utils.cli_tools.get_activities --limit 500 --dry-run
```

Команда лишь логирует факт генерации 500 строк и завершает работу без создания файлов.

## Обновление зависимостей

Периодически обновляйте закреплённые зависимости и проверяйте совместимость:

```bash
pip install -U .[dev]
pre-commit run --all-files
```

Первая команда обновляет рабочие и dev-зависимости до новых минорных версий в допустимых диапазонах. Вторая запускает форматирование, линтинг, статические проверки и тесты, чтобы убедиться в корректности после обновления.

## Конфигурация через `.env`

Часть параметров CLI можно задавать через переменные окружения. Чтобы не экспортировать их вручную, сохраните пары `NAME=value` в файле `.env` и загрузите их с помощью [`python-dotenv`](https://pypi.org/project/python-dotenv/).

Пример файла:

```dotenv
CHEMBL_DA_LOG_LEVEL=INFO
CHEMBL_API_BASE=https://www.ebi.ac.uk/chembl/api/data
```

Типовые переменные с контактными e-mail приведены в `.env.example`.

Запуск скрипта с автоматической подгрузкой настроек:

```bash
python -m dotenv run -- python -m scripts.get_assay_data --input assay_ids.csv \
    --output out/assays.csv
```

Файл `assay_ids.csv` должен содержать столбец `assay_chembl_id` с нужными идентификаторами, например:

```csv
assay_chembl_id
CHEMBL1234567
CHEMBL2345678
```

Утилиты читают переменные окружения автоматически, поэтому значения из `.env` доступны всем CLI без дополнительных аргументов.

## Параметр `api.user_agent`

`api.user_agent` идентифицирует приложение в API-запросах и должен включать контактные данные. Значение по умолчанию:

```yaml
api:
  user_agent: "chembl-da/0.1 (mailto:contact@example.org)"
```

Переопределите параметр в `config/config.yaml` или переменной окружения `CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT`. Отдельного CLI-флага нет (см. `library/cli/parser.py`), поэтому используются только конфигурационные файлы или окружение.

## Валидация конфигурации

`library.config.load_config` проверяет значения в `config/config.yaml`. Некорректный URL вызовет `ValueError` при загрузке:

```yaml
api:
  chembl_base: https://
```

```
ValueError: api.chembl_base must be a valid URL
```

Используйте полный адрес службы:

```yaml
api:
  chembl_base: https://www.ebi.ac.uk/chembl/api/data
```

## Ошибки конфигурации

Некорректные значения в `config/config.yaml` приводят к `ValidationError`. Пример:

```yaml
api:
  rps: -1
```

```python
from library.config import load_config
load_config("config/config.yaml")
```

Результат:

```
pydantic_core._pydantic_core.ValidationError: 1 validation error for Config
api.rps
  Input should be greater than or equal to 1 [type=greater_than_equal, input_value=-1, input_type=int]
    For further information visit https://errors.pydantic.dev/2.11/v/greater_than_equal
```

Замените значение на положительное число:

```yaml
api:
  rps: 5  # любое >= 1
```

Диапазоны допустимых значений описаны в [`config.schema.json`](./config.schema.json), экспортированном из Pydantic-моделей — там для `api.rps` указан минимум `1`.

## Логирование

CLI-хелперы настраивают структурированное JSON-логирование через `library.logging_setup.configure_logger`. Управляйте уровнем и форматом логов переменными окружения или флагами CLI.

Пример включения JSON-формата через переменную окружения:

```bash
LOG_FORMAT=json python -m scripts.get_assay_data --input assay_ids.csv \
    --output out/assays.csv --log-level INFO
```

Уровень логов можно задать флагом `--log-level` или переменной `CHEMBL_DA_LOG_LEVEL`:

```bash
CHEMBL_DA_LOG_LEVEL=DEBUG python -m scripts.get_assay_data --input assay_ids.csv \
    --output out/assays.csv
```

Пример строки лога:

```json
{"ts":"2024-05-01T12:00:00Z","level":"INFO","event":"pipeline_start","run_id":"abc123","stage":"pipeline"}
```

Ключевые поля:

* `ts` — метка времени в формате ISO 8601 (UTC).
* `level` — уровень (`DEBUG`, `INFO`, `WARN`, `ERROR`).
* `event` — краткое машинно-читаемое имя события.
* `run_id` — уникальный идентификатор запуска.
* `stage` — опциональный этап пайплайна.
* `msg` — опциональное человекочитаемое сообщение.
* Дополнительные ключи — контекст события (`elapsed`, `url`, `rows`).

Запуски с `--dry-run` пишут логи с `event = "dry_run"`, что упрощает фильтрацию, например:

```bash
jq 'select(.event=="dry_run")' log.jsonl
```

Идентификатор запуска генерируется CLI-хелперами через `uuid.uuid4().hex` и передаётся логгеру, который включает его в каждую запись. При необходимости подставьте собственное значение перед вызовом `configure_logger`.

Секреты маскируются автоматически: значения ключей, оканчивающихся на `token`, `key`, `secret` или `password`, заменяются на `"***"`. Фильтрация по уровню логов отбрасывает записи ниже заданного `--log-level` или `CHEMBL_DA_LOG_LEVEL`.

Типичные записи:

```json
{"ts":"2024-05-01T12:00:00Z","level":"INFO","event":"pipeline_start","run_id":"abc123","stage":"pipeline"}
{"ts":"2024-05-01T12:00:01Z","level":"INFO","event":"request_ok","run_id":"abc123","stage":"fetch","url":"https://api.example.org","status":200}
{"ts":"2024-05-01T12:00:02Z","level":"INFO","event":"validate_done","run_id":"abc123","stage":"validate","rows":42}
{"ts":"2024-05-01T12:00:03Z","level":"INFO","event":"pipeline_done","run_id":"abc123","stage":"pipeline","elapsed":3.2}
```

## Воспроизводимость

Детерминированные писатели CSV в `library.io` обеспечивают стабильность данных и метаданных между запусками. Функция `library.csv_utils.write_csv_deterministic` нормализует порядок колонок и метаданные, чтобы повторные запуски давали идентичные файлы. Все CLI-скрипты поддерживают общий набор флагов:

```bash
python -m library.utils.cli_tools.table_quality_main --input tests/data/activity.csv \
    --table-name activity
```

По умолчанию `--output` формируется как `output_<имя_входа>_YYYYMMDD.csv` в каталоге, указанном в `local.io.output_dir`. Дополнительные примеры приведены в [`docs/USAGE_RU.md`](docs/USAGE_RU.md) (английская версия — [`docs/USAGE_EN.md`](docs/USAGE_EN.md)).

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

Параметры читаются из `config/config.yaml`, переменных окружения (`CHEMBL_DA__...`) и CLI-флагов. Подробности смотрите в [`docs/CONFIG_RU.md`](docs/CONFIG_RU.md) (английская версия — [`docs/CONFIG_EN.md`](docs/CONFIG_EN.md)).

### Переменные окружения

Переменные окружения переопределяют значения из YAML. Имена следуют шаблону `CHEMBL_DA__...`, где двойное подчёркивание разделяет уровни. Например, для включения debug-логов:

```bash
export CHEMBL_DA__LOG__LEVEL=DEBUG
```

Большинство опций имеют короткие алиасы для обратной совместимости. Таблица перечисляет все поддерживаемые алиасы и их канонические ключи:

| Alias | Эквивалентный ключ |
|-------|--------------------|
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

Полное описание параметров приведено в [`docs/CONFIG_RU.md`](docs/CONFIG_RU.md) (английская версия — [`docs/CONFIG_EN.md`](docs/CONFIG_EN.md)).

### Проверка схемы

Значения конфигурации валидируются функцией :func:`library.config.load_config`, которая вызывает :meth:`Config.model_validate <pydantic.BaseModel.model_validate>`. Проверка следует определениям моделей и формирует подробные сообщения об ошибках для вложенных полей. Файл `config.schema.json` сгенерирован из той же Pydantic-модели для подсказок IDE; он не исполняется во время работы и не должен передаваться в `jsonschema`.

CLI-флаги имеют высший приоритет. Все утилиты принимают `--config` для указания файла и `--print-config` для вывода итоговых значений после всех переопределений. Итоговый порядок приоритетов:

```
YAML < переменные окружения < CLI-флаги
```

Только верхнеуровневые скрипты читают конфигурационный файл. Модули в `library/` ожидают явную передачу `Config` (или его частей), что делает зависимости прозрачными и исключает скрытое глобальное состояние.

Каталоги, указанные в `local.io.output_dir` и `local.io.cache_dir`, проверяются, но не создаются при загрузке конфигурации. Скрипты, которым нужны эти пути, могут вызвать `library.config.ensure_dirs` после `load_config`, чтобы создать каталоги, если они отсутствуют и флаг `local.io.exist_ok` это разрешает.

Пути вроде `local.io.output_dir`, `local.io.cache_dir` и рабочие книги `local.init` возвращаются как объекты `pathlib.Path`. Строки из `config/config.yaml` или переопределения из окружения и CLI автоматически преобразуются.

```bash
python -m library.utils.cli_tools.table_quality_main \
    --input tests/data/activity.csv \
    --table-name activity
```

По умолчанию `--output` формируется как `output_<имя_входа>_YYYYMMDD.csv` в каталоге `local.io.output_dir`. Дополнительные примеры см. в [`docs/USAGE_RU.md`](docs/USAGE_RU.md) (английская версия — [`docs/USAGE_EN.md`](docs/USAGE_EN.md)).

## Вывод и метаданные

Пайплайны записывают детерминированные CSV через `library.io.write_csv` и сохраняют рядом `*.meta.yaml` в `data/output`. Файлы метаданных содержат Git-коммит, параметры запуска, SHA-256 и статистику по строкам/колонкам. Подробности описаны в [`docs/OUTPUT_RU.md`](docs/OUTPUT_RU.md) (английская версия — [`docs/OUTPUT_EN.md`](docs/OUTPUT_EN.md)).

## Dtype Inspector

Утилита `dtype_inspector` запускает каждый `get_*_data` на небольшом наборе идентификаторов и логирует полученные типы pandas. Запускайте её периодически, чтобы выявлять дрейф схем при изменениях внешних сервисов.

```bash
python -m library.utils.cli_tools.dtype_inspector_main --log-level INFO
```

Интегрируйте скрипт в CI для раннего обнаружения изменений — инструмент лёгкий и делает лишь несколько запросов на набор данных.

## Разработка и тестирование

Индивидуальные инструменты `ruff`, `black` и `mypy` подключены к `pre-commit`, но их можно запускать вручную при доработке отдельных модулей.

```bash
ruff check scripts library library/utils/cli_tools
black scripts library library/utils/cli_tools
mypy scripts library
pytest
```

Тестовые наборы лежат в `tests/data`; `library.utils.cli_tools.check_determinism` проверяет повторяемость CSV-вывода.

## Лицензия

MIT License. См. файл `LICENSE` (если присутствует).

При обновлении справочных материалов добавляйте их напрямую в каталог `docs`.

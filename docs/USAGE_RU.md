# Руководство по использованию

Документ описывает стандартные аргументы CLI, типовые сценарии запуска и
вспомогательные утилиты. Для английской версии обратитесь к `docs/USAGE_EN.md`.

## Общие параметры CLI

Все команды `scripts/get_*_data.py` поддерживают единый набор аргументов.
После установки пакета командой `pip install .` каждый пайплайн доступен
и как консольный скрипт:

| Консольная команда | Запуск через `python -m` |
| ------------------ | ------------------------ |
| `get-data` | `python -m scripts.get_data` |
| `get-activity-data` | `python -m scripts.get_activity_data` |
| `get-assay-data` | `python -m scripts.get_assay_data` |
| `get-document-data` | `python -m scripts.get_document_data` |
| `get-target-data` | `python -m scripts.get_target_data` |
| `get-testitem-data` | `python -m scripts.get_testitem_data` |

Выбирайте подходящий формат — оба варианта принимают идентичные аргументы.
Подкоманды для выбора источника (`chembl`, `uniprot`, `iuphar`, `all`)
доступны только у `scripts/get_document_data.py` и
`scripts/get_target_data.py`; остальные скрипты (`get_activity_data.py`,
`get_assay_data.py`, `get_testitem_data.py`) представляют собой одиночные
CLI без вложенных подкоманд.

| Опция | Назначение |
| --- | --- |
| `--config` | Путь к YAML-файлу конфигурации (по умолчанию `config/config.yaml`). |
| `--print-config` | Вывести итоговую конфигурацию после переопределений и завершить работу. |
| `--log-level` | Уровень логирования (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). |
| `--input` | Входной CSV с идентификаторами (по умолчанию `input.csv`). |
| `--output` | Выходной CSV. Если не указан, создаётся `output.<stem>_<YYYYMMDD>.csv` в `local.io.output_dir`. |
| `--sep` | Разделитель CSV; записывается в `cfg.io.csv_sep`. |
| `--encoding` | Кодировка файла; записывается в `cfg.io.csv_encoding`. |
| `--column` | Название колонки с идентификаторами. Значение подтягивается из конфигурации на этапе запуска. |
| `--batch-size` / `--chunk-size` | Максимальное число идентификаторов в одном запросе (название опции зависит от пайплайна). |

Конкретные скрипты добавляют дополнительные ключи (`--timeout`, `--limit`, `--openalex-rps` и т.п.). После разбора аргументов
`apply_config_overrides` загружает `config/config.yaml`, применяет переменные окружения, переносит CLI-переопределения в конфигурацию и
возвращает актуальные значения аргументов.

Перед сетевыми вызовами выполняется `library.config.ensure_dirs`, чтобы `local.io.output_dir` и `local.io.cache_dir` существовали
(если `local.io.exist_ok=true`, каталоги создаются автоматически).

### Примерные входные файлы

В репозитории есть несколько компактных наборов в `tests/data/` для быстрой проверки пайплайнов (например, `tests/data/activity_ids_small.csv` и `tests/data/input-smoke/testitem.csv`). Пайплайны без готовых примеров требуют пользовательских CSV с колонкой идентификатора из `config/config.yaml` или переданной через `--column`.

### Мониторинг структурированных логов

Все утилиты используют `library.logging_setup.Logger` и пишут JSON-строки с уникальным `run_id` и дополнительными полями (`status`, `rps` и др.). Основные события:

| Событие | Когда появляется |
| --- | --- |
| `pipeline_start` | После настройки логирования и перед валидацией конфигурации. |
| `documents_processed` | Счётчик прогресса документного пайплайна после каждого обработанного батча. |
| `process_limit` | Фиксируется при обрезке списка идентификаторов опцией `--limit` или её эквивалентами в конфигурации. |
| `pipeline_skip_limit` | Логируется, когда `--limit 0` завершает работу до сетевых и файловых операций. |
| `write_done` | Успешная запись CSV с указанием пути и числа строк. |
| `pipeline_done` / `pipeline_fail` | Финальный статус перед выходом из программы. |

Предупреждения `pubmed_batch_request_failed` и `pubmed_batch_unexpected_error`
теперь содержат `pmids_count` и укороченный `pmids_sample` вместо полного
списка идентификаторов, поэтому записи остаются компактными и при этом дают
контекст для разбора проблемных батчей.

Для онлайн-контроля направляйте вывод через `jq` или аналогичный инструмент:

```bash
get-document-data all --input documents.csv --column document_chembl_id \
  | tee run.log | jq -r '"\(.level) \(.event) :: \(.msg // "")"'
```

При необходимости повышайте детализацию ключом `--log-level DEBUG`. JSON-структура совместима с системами сбора логов без дополнительных форматтеров.

## Данные активностей (`get_activity_data.py`)

Вариант консольного скрипта:

```bash
get-activity-data --input tests/data/activity_ids_small.csv \
  --column activity_id \
  --batch-size 25 \
  --timeout 45
```

Запуск через модуль:

```bash
python -m scripts.get_activity_data \
  --input tests/data/activity_ids_small.csv \
  --column activity_id \
  --batch-size 25 \
  --timeout 45
```

* В репозитории доступен файл `tests/data/activity_ids_small.csv`; задайте `--column activity_id` (либо переименуйте колонку в `activity_chembl_id`, чтобы использовать настройки по умолчанию).
* Использует колонку `sources.chembl.pipelines.activity.column` (по умолчанию `activity_chembl_id`).
* Создаёт основной CSV, sidecar `*.meta.yaml`, при необходимости `*_failure_cases.csv` и отчёты качества.
* Поддерживает `--limit` (ограничение по количеству ID) и `--dry-run` (проверка входных данных без запросов к API). Значение `--limit 0` мгновенно завершает выполнение без обращений к внешним сервисам и файловой системе.
* Заполняет `lower_value` и `upper_value` на основе канонических полей `standard_*`. Поведение (отношения, разбор `±`, округление, обрезка и логирование) настраивается блоком `activity_bounds.*` в конфигурации.
* Следите за предупреждениями `activity_bounds_unknown_relation` и `activity_bounds_missing_standard_value` в логах — они указывают строки, где границы не удалось восстановить или оператор отношения неизвестен.

## Описания ассайев (`get_assay_data.py`)

Вариант консольного скрипта:

```bash
get-assay-data --input path/to/assay_ids.csv \
  --column assay_chembl_id \
  --batch-size 25
```

Запуск через модуль:

```bash
python -m scripts.get_assay_data \
  --input path/to/assay_ids.csv \
  --column assay_chembl_id \
  --batch-size 25
```

Загружает метаданные ассайев ChEMBL для указанных идентификаторов. Подготовьте CSV с заголовком `assay_chembl_id` (либо передайте своё имя колонки через `--column`) — готового smoke-файла в репозитории нет.

## Метаданные документов (`get_document_data.py`)

Вариант консольного скрипта:

```bash
get-document-data all --input path/to/documents.csv \
  --column document_chembl_id \
  --batch-size 20
```

Запуск через модуль:

```bash
python -m scripts.get_document_data all \
  --input path/to/documents.csv \
  --column document_chembl_id \
  --batch-size 20
```


Команда объединяет данные ChEMBL и PubMed. Подготовьте CSV с колонкой `document_chembl_id` или передайте имя своей колонки через
`--column` — smoke-набор для пайплайна в репозитории отсутствует. Для разовых корректировок используйте доступные флаги (`--chunk-size`, `--sleep`, `--workers`, `--batch-size`, `--timeout`, `--limit`, `--offset`, `--openalex-rps`, `--crossref-rps`) в зависимости от выбранной подкоманды. Подкоманда PubMed дополнительно поддерживает аргументы `--fallback-doi-*` для подстановки DOI. Вложенные параметры (например `sources.chembl.pipelines.document.pubmed.batch_size`) меняются через конфигурацию или переменные окружения, например:

```bash
CHEMBL_DA__SOURCES__CHEMBL__PIPELINES__DOCUMENT__PUBMED__BATCH_SIZE=20 \
  get-document-data all --input path/to/documents.csv \
    --column document_chembl_id
```

Альтернативно обновите значение в `config/config.yaml`.

Выберите подкоманду `pubmed`, `chembl` или `all` в зависимости от требуемых источников.
Сводку и список ключей смотрите в справке: `get-document-data --help`
и `get-document-data <подкоманда> --help`
(например, `--batch-size` управляет размером пакета для PubMed).

Вариант консольного скрипта:

```bash
get-document-data pubmed --input path/to/documents.csv \
  --column PMID \
  --openalex-rps 2.5 \
  --crossref-rps 1.5 \
  --fallback-doi-csv path/to/doi_overrides.csv \
  --fallback-doi-pmid-column pmid_override \
  --fallback-doi-value-column doi_override
```

Запуск через модуль:

```bash
python -m scripts.get_document_data pubmed \
  --input path/to/documents.csv \
  --column PMID \
  --openalex-rps 2.5 \
  --crossref-rps 1.5 \
  --fallback-doi-csv path/to/doi_overrides.csv \
  --fallback-doi-pmid-column pmid_override \
  --fallback-doi-value-column doi_override
```

Флаги `--openalex-rps` и `--crossref-rps` позволяют временно изменить лимиты без правки YAML, а параметры `--fallback-doi-*`
подключают лёгкий CSV с соответствиями PMID→DOI до обращения к внешним сервисам. Подготовьте файл с колонками из аргументов `--fallback-doi-pmid-column` и `--fallback-doi-value-column`; если не задавать их явно, CLI ожидает заголовки `PMID` и `DOI`, а в примере выше показано переименование через явные параметры. Подкоманда PubMed по умолчанию ожидает колонку `PMID`, поэтому флаг `--column` можно опустить, если CSV уже использует это имя.


## Агрегация таргетов (`get_target_data.py`)

Вариант консольного скрипта:

```bash
get-target-data chembl --input path/to/targets.csv \
  --column target_chembl_id
```

Запуск через модуль:

```bash
python -m scripts.get_target_data chembl \
  --input path/to/targets.csv \
  --column target_chembl_id
```

Комбинирует данные ChEMBL, UniProt и IUPHAR согласно разделу `sources.chembl.pipelines.target.*`. Соберите CSV с колонкой `target_chembl_id` (по одной записи в строке); готовый smoke-набор отсутствует.

## Обвязка таргет-пайплайна (`library.utils.cli_tools.pipeline_targets_main`)

```bash
python -m library.utils.cli_tools.pipeline_targets_main \
  --input tests/data/chembl_targets_min.csv \
  --output out/targets_cached.csv \
  --chunk-size 50 \
  --batch-size 50 \
  --limit 200
```

Облегчённая CLI-команда повторяет интерфейс `get_target_data.py`, но запускает
`library.pipeline_targets.run_pipeline` на подготовленных чанках без сетевых
запросов. Идентификаторы читаются через `read_ids` с учётом `--chunk-size`,
`--limit`, разделителя и кодировки, размер батча прокидывается в пайплайн, а
результат записывается после `add_pipeline_metadata` и `write_csv`, что
гарантирует ту же детерминированность, что и у основного пайплайна. Утилита
подходит для проверки переопределений конфигурации, логирования и параметров
батчирования до запуска `get_target_data` с обращениями к внешним API.

## Обогащение тест-объектов (`get_testitem_data.py`)

Вариант консольного скрипта:

```bash
get-testitem-data --input tests/data/input-smoke/testitem.csv \
  --column molecule_chembl_id
```

Запуск через модуль:

```bash
python -m scripts.get_testitem_data \
  --input tests/data/input-smoke/testitem.csv \
  --column molecule_chembl_id
```

Выгружает дополнительную информацию о соединениях. Можно использовать комплект `tests/data/input-smoke/testitem.csv` или собственный CSV с нужной колонкой идентификаторов.

Параметры пагинации теперь берутся из секции `sources.chembl.pipelines.testitem` файла `config.yaml`. По умолчанию пайплайн запрашивает 1000 молекул за вызов (`batch_size`/`request_limit`) и ограничивает ответ полями из `testitem.fields`, чтобы не грузить лишние данные. При необходимости уменьшить батч или добавить дополнительные колонки достаточно изменить конфигурацию или задать переменные окружения:

```yaml
sources:
  chembl:
    pipelines:
      testitem:
        offset: 500        # пропустить первые 500 идентификаторов
        request_limit: 750 # зафиксировать предел страницы ниже максимального значения API
        fields:
          - molecule_chembl_id
          - pref_name
          - structure_type
```

Пример запуска с кастомными настройками:

Вариант консольного скрипта:

```bash
get-testitem-data --config config/config.yaml \
  --input data/input/testitem_ids.csv \
  --batch-size 750 \
  --offset 500
```

Запуск через модуль:

```bash
python -m scripts.get_testitem_data \
  --config config/config.yaml \
  --input data/input/testitem_ids.csv \
  --batch-size 750 \
  --offset 500
```

### Контроль `properties_hash`

PubChem-дополнение добавляет детерминированные свойства (`pubchem_cid`, `pubchem_iupac_name`, `pubchem_molecular_formula`,
`pubchem_isomeric_smiles`, `pubchem_canonical_smiles`, `pubchem_inchi`, `pubchem_inchikey`). Чтобы
отслеживать изменения между выгрузками, выгрузите только эти колонки во временный файл и посчитайте SHA-256 с помощью
`library.metadata.file_sha256` или `library.csv_utils.sha256_file`. Полученное значение `properties_hash` удобно сохранять в журнале релиза или sidecar, чтобы фиксировать сдвиги в данных PubChem даже при неизменном количестве строк.

### Требования к каталогу родительских молекул

Чтобы получить `parent_molecule_chembl_id` для агрегаций, выгрузку необходимо объединить с каталогом
родителей ChEMBL. Путь к локальному JSON задаётся через
[`sources.chembl.molecule_catalog.cache_path`](./CONFIG_RU.md#sources-chembl-molecule-catalog); убедитесь,
что файл доступен исполнителю, либо задайте новое расположение переменной окружения
`CHEMBL_DA_MOLECULE_CATALOG_CACHE` (алиас для `CHEMBL_DA__SOURCES__CHEMBL__MOLECULE_CATALOG__CACHE_PATH`) или правкой

`config/config.yaml`.


Для первичного создания либо обновления файла выполните небольшой Python-скрипт с вызовом
`library.molecule_catalog.load_parent_catalog` — функция считывает готовый кэш и, при его отсутствии,
подкачивает свежие связи ребёнок→родитель из API ChEMBL.

### Обогащение солей и флагов каталога

Опциональный блок `testitem_molecule_enrichment` добавляет к `testitem.csv`
столбцы `salt_chembl_id`, `natural_product`, `prodrug`, `polymer_flag` на
основе двух CSV-словарей:

* `tests/data/input-smoke/molecule_hierarchy.csv` со столбцами `molecule_chembl_id`,
  `parent_molecule_chembl_id` описывает связи соли и родителя.
* `tests/data/input-smoke/molecule_catalog.csv` со столбцами `molecule_chembl_id`,
  `natural_product`, `prodrug`, `polymer_flag` содержит булевы признаки.

Если молекулы нет в словарях, в лог попадают предупреждения
`testitem_enrichment_missing_child_flags` или
`testitem_enrichment_missing_parent_flags` — проверьте данные или временно
отключите блок через `testitem_molecule_enrichment.enable=false`. Сообщение
`testitem_enrichment_inconsistent_flag` сигнализирует о расхождении флагов
между дочерней и родительской записью до применения фолбэка. Поведение можно
тонко настроить через `testitem_molecule_enrichment.flags.*`, отключив
фолбэк или приведение к булевому типу, если downstream-потребители требуют
исходные текстовые значения.

## Инициализация входных данных (`library/utils/cli_tools/get_input_initialisation.py`)

```bash
python -m library.utils.cli_tools.get_input_initialisation \
  --same-doc path/to/ChEMBL_same_document.xlsx \
  --all-doc path/to/ChEMBL_all.xlsx \
  --out-dir path/to/output
```

* Формирует таблицы пар (`pairs_same_document.csv`, `pairs_independent.csv`, `pairs_non_independent.csv`).
* Создаёт срезы по сущностям (`activity_*`, `assay_*`, `document_*`, `target_*`, `testitem_*`, `system_*`).
* Добавляет папку `data_validity_report/` с отчётами качества для каждого файла. Используйте исходные Excel-выгрузки вашего процесса ChEMBL — демонстрационные файлы в репозитории не поставляются.

## Профайлер качества таблиц (`library/utils/cli_tools/table_quality_main.py`)

```bash
python -m library.utils.cli_tools.table_quality_main \
  --input tests/data/activities_valid.csv \
  --table-name activity
```

Генерирует `<table-name>_quality_report_table.csv` и `<table-name>_data_correlation_report_table.csv`, используя настройки `local.io`. При необходимости подставьте собственные CSV.

## Переопределения конфигурации во время запуска

CLI-флаги покрывают задокументированные аргументы каждого скрипта. Например, пайплайн активностей поддерживает
`--batch-size`, `--timeout`, `--limit`, `--dry-run`:

```bash
get-activity-data --batch-size 25 --timeout 45
```

Вложенные параметры меняются через `config/config.yaml` или переменные окружения. Чтобы временно увеличить лимит запросов
к ChEMBL без правки файла, экспортируйте переменную и запустите команду в том же сеансе:

```bash
export CHEMBL_DA__SOURCES__CHEMBL__API__RPS=10
get-activity-data
```

При необходимости проверьте итоговую конфигурацию флагом `--print-config` до запуска пайплайна.

## Мониторинг структурированных логов

Все CLI пишут JSON-логи через `library.logging_setup`. Запись включает отметку времени (`ts`), уровень (`level`), событие
(`event`) и `run_id`, унаследованный от параметров CLI; дополнительные поля добавляются после автоматической маскировки
секретов. Применяйте `jq` или подобные инструменты, чтобы фильтровать события по `event`, `stage` или кодам предупреждений
(`activity_bounds_*`, `parent_lookup_*` и т.д.). Меняйте уровень детализации флагом `--log-level` или переменными окружения без
правок `config/config.yaml`.

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

Для съёма покрытия основных модулей (`library/`, `scripts/`) используйте:

```bash
pytest --cov=library --cov=scripts \
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

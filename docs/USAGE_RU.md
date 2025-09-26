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

### Мониторинг структурированных логов

Все утилиты используют `library.logging_setup.Logger` и пишут JSON-строки с уникальным `run_id` и дополнительными полями (`status`, `rps` и др.). Основные события:

| Событие | Когда появляется |
| --- | --- |
| `pipeline_start` | После настройки логирования и перед валидацией конфигурации.【F:scripts/get_activity_data.py†L558-L565】 |
| `documents_processed` / `activities_processed` | Периодические счётчики прогресса внутри рабочих циклов.【F:scripts/get_document_data.py†L399-L407】【F:scripts/get_activity_data.py†L451-L459】 |
| `write_done` | Успешная запись CSV с указанием пути и числа строк.【F:scripts/get_document_data.py†L640-L642】【F:scripts/get_activity_data.py†L506-L522】 |
| `pipeline_done` / `pipeline_fail` | Финальный статус перед выходом из программы.【F:scripts/get_activity_data.py†L571-L579】【F:scripts/get_document_data.py†L1191-L1208】 |

Для онлайн-контроля направляйте вывод через `jq` или аналогичный инструмент:

```bash
python scripts/get_document_data.py all --input documents.csv --column document_chembl_id \
  | tee run.log | jq -r '"\(.level) \(.event) :: \(.msg // "")"'
```

При необходимости повышайте детализацию ключом `--log-level DEBUG`. JSON-структура совместима с системами сбора логов без дополнительных форматтеров.【F:library/logging_setup.py†L1-L120】

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


Команда объединяет данные ChEMBL и PubMed. Для разовых корректировок доступны флаги `--batch-size`, `--timeout`, `--limit`,
`--dry-run`. Вложенные параметры (например `sources.chembl.pipelines.document.pubmed.batch_size`) меняются через конфигурацию
или переменные окружения, например:

```bash
CHEMBL_DA__SOURCES__CHEMBL__PIPELINES__DOCUMENT__PUBMED__BATCH_SIZE=20 \
  python scripts/get_document_data.py \
    --input data/input-smoke/documents.csv \
    --column document_chembl_id
```

Альтернативно обновите значение в `config.yaml`.

Выберите подкоманду `pubmed`, `chembl` или `all` в зависимости от требуемых источников.
Сводку и список ключей смотрите в справке: `python scripts/get_document_data.py --help`
и `python scripts/get_document_data.py <подкоманда> --help`
(например, `--batch-size` управляет размером пакета для PubMed).

```bash
python scripts/get_document_data.py pubmed \
  --input data/input-smoke/documents.csv \
  --column document_chembl_id \
  --openalex-rps 2.5 \
  --crossref-rps 1.5 \
  --fallback-doi-csv data/input-smoke/doi_overrides.csv \
  --fallback-doi-pmid-column pmid_override \
  --fallback-doi-value-column doi_override
```

Флаги `--openalex-rps` и `--crossref-rps` позволяют временно изменить лимиты без правки YAML, а параметры `--fallback-doi-*`
подключают лёгкий CSV с соответствиями PMID→DOI до обращения к внешним сервисам.【F:scripts/get_document_data.py†L989-L1041】


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

### Контроль `properties_hash`

PubChem-дополнение добавляет детерминированные свойства (`pubchem_cid`, `pubchem_iupac_name`, `pubchem_molecular_formula`,
`pubchem_isomeric_smiles`, `pubchem_canonical_smiles`, `pubchem_inchi`, `pubchem_inchikey`).【F:schemas/testitems.py†L30-L38】 Чтобы
отслеживать изменения между выгрузками, выгрузите только эти колонки во временный файл и посчитайте SHA-256 с помощью
`library.metadata.file_sha256` или `library.csv_utils.sha256_file`.【F:library/metadata.py†L29-L70】【F:library/csv_utils.py†L530-L560】 Полученное значение `properties_hash` удобно сохранять в журнале релиза или sidecar, чтобы фиксировать сдвиги в данных PubChem даже при неизменном количестве строк.

### Требования к каталогу родительских молекул

Чтобы получить `parent_molecule_chembl_id` для агрегаций, выгрузку необходимо объединить с каталогом
родителей ChEMBL. Путь к локальному JSON задаётся через
[`sources.chembl.molecule_catalog.cache_path`](./CONFIG_RU.md#sources-chembl-molecule-catalog); убедитесь,
что файл доступен исполнителю, либо задайте новое расположение переменной окружения
`CHEMBL_DA_MOLECULE_CATALOG_CACHE` (алиас для `CHEMBL_DA__SOURCES__CHEMBL__MOLECULE_CATALOG__CACHE_PATH`) или правкой

`config.yaml`.【F:config.yaml†L25-L33】【F:library/config.py†L487-L551】


Для первичного создания либо обновления файла выполните небольшой Python-скрипт с вызовом
`library.molecule_catalog.load_parent_catalog` — функция считывает готовый кэш и, при его отсутствии,
подкачивает свежие связи ребёнок→родитель из API ChEMBL.【F:library/molecule_catalog.py†L43-L136】

### Обогащение солей и флагов каталога

Опциональный блок `testitem_molecule_enrichment` добавляет к `testitem.csv`
столбцы `salt_chembl_id`, `natural_product`, `prodrug`, `polymer_flag` на
основе двух CSV-словарей:

* `dictionary/_testitem/molecule_hierarchy.csv` со столбцами `molecule_chembl_id`,
  `parent_molecule_chembl_id` описывает связи соли и родителя.
* `dictionary/_testitem/molecule_catalog.csv` со столбцами `molecule_chembl_id`,
  `natural_product`, `prodrug`, `polymer_flag` содержит булевы признаки.

Если молекулы нет в словарях, в лог попадают предупреждения
`testitem_enrichment_missing_child_flags` или
`testitem_enrichment_missing_parent_flags` — проверьте данные или временно
отключите блок через `testitem_molecule_enrichment.enable=false`. Сообщение
`testitem_enrichment_inconsistent_flag` сигнализирует о расхождении флагов
между дочерней и родительской записью до применения фолбэка. Поведение можно
тонко настроить через `testitem_molecule_enrichment.flags.*`, отключив
фолбэк или приведение к булевому типу, если downstream-потребители требуют
исходные текстовые значения.【F:library/testitem_enrichment.py†L17-L216】

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

## Профайлер качества таблиц (`scripts/table_quality_main.py`)

```bash
python scripts/table_quality_main.py \
  --input data/input-smoke/activity.csv \
  --table-name activity
```

Генерирует `<table-name>_quality_report_table.csv` и `<table-name>_data_correlation_report_table.csv`, используя настройки `local.io`.

## Переопределения конфигурации во время запуска

CLI-флаги покрывают задокументированные аргументы каждого скрипта. Например, пайплайн активностей поддерживает
`--batch-size`, `--timeout`, `--limit`, `--dry-run`:

```bash
python scripts/get_activity_data.py --batch-size 25 --timeout 45
```

Вложенные параметры меняются через `config.yaml` или переменные окружения. Чтобы временно увеличить лимит запросов
к ChEMBL без правки файла, экспортируйте переменную и запустите команду в том же сеансе:

```bash
export CHEMBL_DA__SOURCES__CHEMBL__API__RPS=10
python scripts/get_activity_data.py
```

При необходимости проверьте итоговую конфигурацию флагом `--print-config` до запуска пайплайна.

## Мониторинг структурированных логов

Все CLI пишут JSON-логи через `library.logging_setup`. Запись включает отметку времени (`ts`), уровень (`level`), событие
(`event`) и `run_id`, унаследованный от параметров CLI; дополнительные поля добавляются после автоматической маскировки
секретов. Применяйте `jq` или подобные инструменты, чтобы фильтровать события по `event`, `stage` или кодам предупреждений
(`activity_bounds_*`, `parent_lookup_*` и т.д.). Меняйте уровень детализации флагом `--log-level` или переменными окружения без
правок `config.yaml`.【F:library/logging_setup.py†L65-L205】

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

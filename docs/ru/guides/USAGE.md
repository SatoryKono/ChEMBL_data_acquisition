# Руководство по использованию

Русская версия описывает запуск пайплайнов ChEMBL и вспомогательных утилит. За
английским вариантом обращайтесь к [`../../en/guides/USAGE.md`](../../en/guides/USAGE.md).

## Общий шаблон CLI

Каждый пайплайн доступен как консольная команда (после установки пакета) и как
вызов `python -m scripts.<name>` в режиме разработки. Параметры делятся на три
группы:

1. **Общие флаги** из `library.cli.parser.add_common_arguments`:
   - `--input / --final-out` — путь к входному CSV и результирующему файлу.
     Сохраняется только скрытый алиас `--out` для обратной совместимости,
     который выводит предупреждение; документы и новые скрипты должны
     использовать `--final-out`.
   - `--log-level` — уровень логирования (`DEBUG`, `INFO`, `WARNING`, `ERROR`).
   - `--verbose` — включает детальный (`DEBUG`) вывод без изменения YAML-конфига.
   - `--run-id` — явный идентификатор запуска для логов и sidecar-файлов; по
     умолчанию вычисляется детерминированно или берётся из `CHEMBL_DA_RUN_ID`,
     если переменная задана.
   - `--sep`, `--encoding` — разделитель и кодировка CSV (`utf-8-sig` по умолчанию).
   - `--base-path`, `--input-dir`, `--output-dir`, `--date` — удобные шорткаты,
     которые оркестратор применяет для построения путей и имён файлов.
   - `--force`, `--skip-existing` — перезаписать существующий файл или пропустить
     этап, если выгрузка уже создана.
   - `--config` — альтернативный YAML-конфиг (по умолчанию `config/config.yaml`).
   - `--print-config` — вывести эффективную конфигурацию (с учётом CLI и
     переменных окружения) и завершить работу.
2. **Общие опции пагинации**, задаваемые в конкретных пайплайнах:
   `--column`, `--batch-size` / `--chunk-size`, `--timeout`, `--limit`,
   `--offset`, `--workers` (для параллельных запросов), `--dry-run`
   (только для активностей).
3. **Специфические параметры**, управляющие поведением отдельных пайплайнов:
   стадийные флаги таргет-пайплайна, настройки DOI-фолбэков для документов и т.п.

Любая команда возвращает ненулевой код при ошибках валидации, проблемах с IO или
неудачных запросах к внешним сервисам. Для каждого запуска создаётся текстовый лог
`logs/<script>_<YYYYMMDD>.log` в корне репозитория. Переменная
`CHEMBL_DA_BASE_PATH` переносит директорию в `<base>/logs`. Формат записей
остается `[время] [уровень] [логгер] сообщение`, что упрощает последующую
диагностику.

### Шаблоны входных файлов

В репозитории лежат минимальные CSV-шаблоны в `data/input` (по одному на каждый
пайплайн: `document.csv`, `target.csv`, `assay.csv`, `activity.csv`,
`testitem.csv`) и расширенные примеры в `data/input/full`. Скопируйте нужный
файл, заполните идентификаторы и передайте путь через `--input`. Если вы
используете собственные списки, разместите их в любом доступном каталоге и
убедитесь, что заголовки соответствуют [`../DATA_SCHEMA.md`](../DATA_SCHEMA.md).

## Оркестратор (`get-data`)

```
get-data --base-path /data/chembl \
    --input-dir input --output-dir exports \
    --config /data/chembl/config.yaml \
    --date 20250101 --limit 100 --log-level INFO
```

Оркестратор готовит аргументы и запускает пайплайны в порядке: документы (`all`),
таргеты (`all`), ассайи, тест-объекты, активности. В примере выше входные CSV
читаются из каталога `/data/chembl/input` — скопируйте туда шаблоны или укажите
собственный путь. Каждая команда получает параметр `--final-out`, поэтому
выгрузки создаются в целевых каталогах без устаревших алиасов. Значение
`--limit 0` пропускает выполнение, а `--dry-run` выводит план без записи файлов.

Полезные переключатели при отладке:

| Флаг | Действие | Типичный сценарий |
|------|----------|-------------------|
| `--rerun-postprocess` | Пересобирает выравненные выгрузки даже при наличии кеша. | Обновить результаты после правок нормализации. |
| `--debug` | Включает `DEBUG` и сохраняет промежуточные артефакты. | Разобрать падение шага локально. |
| `--keep-intermediate` | Сохраняет промежуточные каталоги без включения `DEBUG`. | Собрать диагностические файлы «на потом». |
| `--dry-run` | Формирует и логирует план без записи файлов. | Проверить конфигурацию в CI или ноутбуке. |
| `--disable-pubchem` | Отключает обогащение PubChem на шаге testitem. | Сравнить с устаревшими выгрузками или изолировать отклонения enrichment. |
| `--force-pubchem` | Принудительно включает PubChem даже при отключении в YAML. | Временный обход экспериментов с отключённым обогащением. |
| `--print-config` | Выводит итоговую конфигурацию и завершает работу. | Зафиксировать параметры запуска для отчёта. |
| `--run-id` | Переопределяет вычисленный идентификатор в логах и манифестах. | Связать прогон с внешним планировщиком или тикетом. |

Приоритет PubChem-флагов следующий: `--disable-pubchem` всегда выигрывает, затем учитывается `--force-pubchem`, а при отсутствии обоих признаков используется значение `sources.pubchem.enable` из конфигурации.【F:library/cli/commands/get_data.py†L1060-L1108】【F:library/cli/commands/get_data.py†L1427-L1471】

Полный перечень базовых и расширенных флагов приведён в [USAGE.md](../USAGE.md); таблица выше выделяет наиболее востребованные при ежедневной отладке опции.【F:library/cli/commands/get_data.py†L949-L1108】

### Примеры частичного запуска документа

Для повторных прогонов часто достаточно перезапустить только документный шаг в
режиме `chembl` или `pubmed`, не трогая остальные пайплайны. Это делается через
`--override-subcommand`, которое меняет подкоманду шага и прокидывается дальше по
всему плану.

```bash
python scripts/get_data.py \
  --base-path /data/chembl \
  --input-dir input --output-dir exports \
  --config /data/chembl/config.yaml \
  --date 20250101 \
  --override-subcommand document=chembl
```

Команда выше выполнит только ChEMBL-часть конвейера документов и завершит
оркестрацию, как только документный шаг закончится. Аналогично можно собрать
данные PubMed, сохранив настройки лимитов и директорий:

```bash
python scripts/get_data.py \
  --base-path /data/chembl \
  --input-dir input --output-dir exports \
  --config /data/chembl/config.yaml \
  --date 20250101 \
  --override-subcommand document=pubmed
```

Если нужно переименовать результат (например, чтобы отличать `chembl` и
`pubmed`-варианты), добавьте `--override-output-stem document=document_pubmed` —
этот флаг воздействует только на выбранный шаг и не требует правок реестра.

## Пайплайн документов (`python scripts/get_document_data.py`)

Конвейер документов запускается единым скриптом `python scripts/get_document_data.py`.
Если режим не указан, по умолчанию выполняется объединённый прогон `all`, но флаг
`--mode <chembl|pubmed|all>` остаётся доступен для фокусировки на конкретной
стадии. Он заменяет прежние позиционные подкоманды, сохраняя общие аргументы CLI
в духе остальных пайплайнов.

### Краткая шпаргалка

| `--mode` | Назначение | Колонка по умолчанию | Пространственные флаги |
|----------|------------|----------------------|------------------------|
| `chembl` | Выгружает документные метаданные из ChEMBL. | `document_chembl_id` | `--chembl-chunk-size`, `--chembl-timeout` (алиасы: `--chunk-size`, `--timeout`). |
| `pubmed` | Обогащает данными PubMed, Semantic Scholar, OpenAlex и CrossRef. | `PMID` | `--pubmed-sleep`, `--pubmed-workers`, `--pubmed-batch-size`, `--openalex-rps`, `--crossref-rps`. |
| `all` | Последовательно выполняет стадии ChEMBL и PubMed и объединяет результат. | `document_chembl_id` | Принимает оба пространства имён и поддерживает параметры DOI-фолбэков. |

### Фрагмент справки

```
$ python scripts/get_document_data.py --mode all --help
...
  --chembl-chunk-size CHEMBL_CHUNK_SIZE, --chunk-size CHEMBL_CHUNK_SIZE
                        Максимум идентификаторов в одном запросе ChEMBL
  --pubmed-sleep PUBMED_SLEEP, --sleep PUBMED_SLEEP
                        Пауза между запросами к PubMed в секундах
  --pubmed-workers PUBMED_WORKERS, --workers PUBMED_WORKERS
                        Количество параллельных запросов PubMed
  --pubmed-batch-size PUBMED_BATCH_SIZE, --batch-size PUBMED_BATCH_SIZE
                        Максимум PMID в одном запросе PubMed
  --chembl-timeout CHEMBL_TIMEOUT, --timeout CHEMBL_TIMEOUT
                        Таймаут HTTP-запросов к ChEMBL в секундах
  --openalex-rps OPENALEX_RPS
                        Ограничение RPS для OpenAlex
  --crossref-rps CROSSREF_RPS
                        Ограничение RPS для CrossRef

Fallback DOI overrides:
  --fallback-doi-enabled
                        Включить чтение DOI из CSV-файла
  --fallback-doi-path FALLBACK_DOI_PATH
                        CSV с DOI, сопоставленными по PMID
  --fallback-doi-col-pmid FALLBACK_DOI_COL_PMID
                        Колонка с PubMed ID во fallback-файле
  --fallback-doi-col-doi FALLBACK_DOI_COL_DOI
                        Колонка с DOI во fallback-файле
  --fallback-doi-delimiter FALLBACK_DOI_DELIMITER
                        Разделитель CSV (по умолчанию: io.csv_sep)
  --fallback-doi-encoding FALLBACK_DOI_ENCODING
                        Кодировка CSV (по умолчанию: io.csv_encoding)
  --fallback-doi-overwrite
                        Разрешить перезапись существующих DOI значениями из CSV
```

### Флаги DOI-фолбэков

| Флаг | Значение по умолчанию | Описание |
|------|-----------------------|----------|
| `--fallback-doi-enabled` | Выключен | Активирует чтение переопределений из CSV. |
| `--fallback-doi-path` | _обязателен при включении_ | Путь к CSV с парами PMID → DOI. |
| `--fallback-doi-col-pmid` | `PMID` | Имя колонки с PubMed ID во fallback-файле. |
| `--fallback-doi-col-doi` | `DOI` | Имя колонки с DOI во fallback-файле. |
| `--fallback-doi-delimiter` | `local.io.csv_sep` (по умолчанию `,`) | Разделитель при чтении CSV. |
| `--fallback-doi-encoding` | `local.io.csv_encoding` (по умолчанию `utf-8-sig`) | Кодировка fallback-файла. |
| `--fallback-doi-overwrite` | Выключен | Разрешает замещать существующие DOI значениями из CSV. |

### Примеры запуска

Ниже режимы указаны явно для наглядности. Без `--mode` команда выполняет
дефолтный объединённый прогон `all`.

```bash
# Только ChEMBL
python scripts/get_document_data.py --mode chembl \
    --input data/input/document.csv \
    --final-out output/documents_chembl.csv \
    --config config/config.yaml

# PubMed и партнёры с ограничением RPS
python scripts/get_document_data.py --mode pubmed \
    --input data/input/document.csv \
    --final-out output/documents_pubmed.csv \
    --config config/config.yaml \
    --openalex-rps 3 --crossref-rps 3

# Полный прогон с ручными DOI
python scripts/get_document_data.py --mode all \
    --input data/input/document.csv \
    --final-out output/documents_full.csv \
    --config config/config.yaml \
    --fallback-doi-enabled \
    --fallback-doi-path data/input/document_fallback.csv \
    --fallback-doi-overwrite
```

Выгрузка включает детерминированный CSV, `<имя>.meta.yaml`, отчёты `<имя>_quality_report_table.csv`, `<имя>_data_correlation_report_table.csv` и `<имя>.quality.json` со статистикой по DOI.

### Повторный запуск постобработки

Если нужно обновить итоговые выгрузки, сохранив подготовленные стадии,
добавьте `--rerun-postprocess`. Скрипт пересоберёт согласованные с этапами
артефакты даже при наличии результатов прошлых запусков. В паре с
`--emit-legacy-artifacts` обновятся и стандартные каталоги вывода, и legacy-CSV.
Для диагностики можно дополнительно указать `--keep-intermediate` или запустить
оркестратор с `--debug`, чтобы сохранить временные каталоги:

```bash
python scripts/get_document_data.py --mode all \
    --rerun-postprocess --emit-legacy-artifacts \
    --final-out output/documents_full.csv
```

## Пайплайн таргетов (`get-target-data`)

Поддерживаются латинские и кириллические алиасы (`chembl`/`сруьид`, `uniprot`/`гтшзкще`, `iuphar`/`шгзрфк`, `all`/`фдд`).

> **Требование по UniProt**: во входном CSV обязательно должна быть колонка с
> UniProt-идентификаторами. По умолчанию ожидается `uniprot_id`. Можно указать
> альтернативное имя через `target.all.uniprot_column`, либо положиться на
> автоматический поиск по ключевым словам вроде `uniprot` или `accession`.

### Стадийные флаги

Только таргет-пайплайн реализует следующие опции:

- `--raw-out` — путь к объединённому набору до финальной нормализации.
- `--raw-format` — формат «сырого» снимка (`csv` или `parquet`).
- `--id-cols` — колонки, используемые для детерминированной сортировки при
  записи «сырого» файла.
- `--no-reindex-raw` — сохранить порядок колонок, полученный от внешних API.
- `--normalize-at-export / --no-normalize-at-export` — применять ли нормализацию
  перед записью финального CSV (`--no-normalize-at-export` оставляет сырые данные).

### Режимы

| Режим | Описание | Основные параметры |
|-------|----------|--------------------|
| `chembl` | Получает таргеты из ChEMBL, нормализует, валидирует и экспортирует. | `--column`, `--chunk-size`, `--timeout`, `--limit`, `--offset`. |
| `uniprot` | Разрешает UniProt-идентификаторы через локальный кеш или REST API. | `--column`, `--limit`. |
| `iuphar` | Сопоставляет UniProt с семействами IUPHAR по встроенным словарям. | `--limit`. |
| `all` | Последовательно запускает `chembl`, `uniprot`, `iuphar`, объединяет и экспортирует результат. | Все параметры выше плюс стадийные флаги. |

Пример со «сбором сырья»:

```
get-target-data all \
    --input data/input/target.csv \
    --final-out output/targets_$(date +%Y%m%d).csv \
    --raw-out output/targets_raw_$(date +%Y%m%d).parquet \
    --raw-format parquet --id-cols target_chembl_id uniprot_id \
    --config config/config.yaml --log-level INFO
```

Для офлайн-проверок используйте `python -m library.utils.cli_tools.pipeline_targets_main` —
он повторяет боевую логику, читая кешированные ответы.

## Пайплайн ассайев (`get-assay-data`)

```
get-assay-data --input data/input/assay.csv \
    --final-out output/assay_$(date +%Y%m%d).csv \
    --batch-size 50 --timeout 60 --limit 200 --postprocess
```

Скрипт выгружает данные из ChEMBL, вычисляет счётчики по таргетам, нормализует и
валидирует таблицу, затем формирует стандартный набор артефактов.
Флаг `--postprocess` запускает этап очистки, который удаляет устаревшие
колонки и пересоздаёт QA-артефакты тем же кодом, что и оркестратор.【F:library/cli/commands/get_assay_data.py†L864-L908】

## Пайплайн активностей (`get-activity-data`)

```
get-activity-data --input data/input/activity.csv \
    --final-out output/activities_$(date +%Y%m%d).csv \
    --batch-size 20 --workers 4 --timeout 90 --limit 500
```

Флаг `--dry-run` уникален для данного пайплайна: он проверяет аргументы и входной
файл и завершает работу без сетевых запросов. Пост-обработка вычисляет
`lower_value`/`upper_value` по алгоритму, описанному в [`../OUTPUT.md`](../OUTPUT.md).

## Пайплайн тест-объектов (`get-testitem-data`)

```
get-testitem-data --input data/input/testitem.csv \
    --final-out output/testitems_$(date +%Y%m%d).csv \
    --batch-size 250 --timeout 90 --limit 400 --pubchem-enable --postprocess
```

Пайплайн объединяет молекулы из ChEMBL с данными PubChem, нормализует результат и
по умолчанию сохраняет три CSV: датасет, `*_quality_report_table.csv` и
`*_data_correlation_report_table.csv`. Используйте `--pubchem-enable`, чтобы
принудительно включить обогащение (а `--no-pubchem-enable` — для тестовых прогонов)
и `--postprocess`, чтобы запустить детерминированную постобработку с иерархиями и
QA. В связке с `--emit-legacy-artifacts` это возвращает исторические диагностические
файлы (failure cases, `.meta.yaml`, отчёты постобработки).【F:library/cli/commands/get_testitem_data.py†L896-L934】【F:library/pipelines/testitem/cli.py†L812-L902】

## Вспомогательные утилиты

Каждый модуль экспортирует функцию `main(argv)` и подходит для использования как
через `python -m`, так и через зарегистрированные консольные команды.

| Модуль | Команда | Назначение |
|--------|---------|------------|
| `library.utils.cli_tools.check_determinism` | `check-determinism --log-level INFO` | Проверка детерминированности CSV через сравнение хэшей тестовых файлов. |
| `library.utils.cli_tools.chunk_io_main` | `chunk-io --input data.csv --final-out copy.csv` | Потоковое чтение/запись CSV с сохранением порядка и Unix-разделителями строк. |
| `library.utils.cli_tools.csv_utils_main` | `csv-utils --input data.csv --final-out clean.csv --sep ,` | Нормализация разделителей, кавычек и порядка колонок. |
| `library.utils.cli_tools.dtype_inspector_main` | `python -m library.utils.cli_tools.dtype_inspector_main` | Диагностика типов pandas-таблиц. |
| `library.utils.cli_tools.get_activities` | `python scripts/get_activities.py --limit 10 --dry-run` | Генерация тестовых активностей для проверки логов и CLI; по умолчанию используется колонка `activity_id`. |
| `library.utils.cli_tools.get_activities` | `get-activities --limit 10 --final-out output/activities_smoke.csv` | Генерация тестовых активностей с записью детерминированного CSV и `.meta.yaml` для smoke-проверок; по умолчанию используется колонка `activity_id`. |
| `library.utils.cli_tools.get_document_type` | `get-document-type --input docs.csv` | Применение встроенных правил классификации публикаций. |
| `library.utils.cli_tools.get_input_initialisation` | `get-input-initialisation --same-doc init.xlsx --all-doc pairs.xlsx` | Объединение Excel-книг инициализации. |
| `library.utils.cli_tools.mapper_main` | `mapper --input ids.csv --final-out mapped.csv --column target_chembl_id` | Интерактивный маппер UniProt/ChEMBL. |
| `library.utils.cli_tools.mapper_batch_main` | `python -m library.utils.cli_tools.mapper_batch_main --input ids.csv` | Пакетный маппер для автоматизации и QA. |
| `library.utils.cli_tools.pipeline_targets_main` | `python -m library.utils.cli_tools.pipeline_targets_main --input targets.csv` | Повтор таргет-пайплайна на кешированных ответах с поддержкой стадийных флагов. |
| `library.utils.cli_tools.table_quality_main` | `table-quality --input data.csv --table-name chembl_targets --final-out reports/` | Профилирование таблиц и отчёты качества. |

## Советы для масштабных запусков

- Явно задавайте пути `--final-out` и `--raw-out`, если в течение дня запускается
  несколько выгрузок; дефолтное имя включает `output.<stem>` и дату из `--date`.
- Для автоматизации задайте `CHEMBL_DA_BASE_PATH`, чтобы контролировать каталоги
  вывода (см. секцию `local.io` в `config/config.yaml`).
- Документный пайплайн поддерживает CSV с руками проставленными DOI через
  `--fallback-doi-csv`, `--fallback-doi-pmid-column`, `--fallback-doi-value-column`.
- Все команды читают `CHEMBL_DA_LOG_LEVEL`. При необходимости объединяйте его с
  `--verbose` (или `--log-level DEBUG`) для подробной диагностики без правки YAML.

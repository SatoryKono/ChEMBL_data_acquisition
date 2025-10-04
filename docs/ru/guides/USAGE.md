# Руководство по использованию

Русская версия описывает запуск пайплайнов ChEMBL и вспомогательных утилит. За
английским вариантом обращайтесь к [`docs/USAGE_EN.md`](./USAGE_EN.md).

## Общий шаблон CLI

Каждый пайплайн доступен как консольная команда (после установки пакета) и как
вызов `python -m scripts.<name>` в режиме разработки. Параметры делятся на три
группы:

1. **Общие флаги** из `library.cli.parser.add_common_arguments`:
   - `--input / --final-out` — путь к входному CSV и результирующему файлу.
     Алиасы `--output` и `--out` оставлены для совместимости, но сопровождаются
     предупреждением; используйте `--final-out`.
   - `--log-level` — уровень логирования (`DEBUG`, `INFO`, `WARNING`, `ERROR`).
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

При запуске документного пайплайна через `python scripts/get_document_data.py`
используйте `--mode <chembl|pubmed|all>`, чтобы выбрать нужную стадию. У
установленной консольной команды сохраняются позиционные подкоманды для
обратной совместимости.

Любая команда возвращает ненулевой код при ошибках валидации, проблемах с IO или
неудачных запросах к внешним сервисам.

### Шаблоны входных файлов

В репозитории лежат минимальные CSV-шаблоны в `data/input` (по одному на каждый
пайплайн: `document.csv`, `target.csv`, `assay.csv`, `activity.csv`,
`testitem.csv`) и расширенные примеры в `data/input/full`. Скопируйте нужный
файл, заполните идентификаторы и передайте путь через `--input`. Если вы
используете собственные списки, разместите их в любом доступном каталоге и
убедитесь, что заголовки соответствуют [`docs/DATA_SCHEMA_RU.md`](./DATA_SCHEMA_RU.md).

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

## Пайплайн документов (`get-document-data`)

Пайплайн доступен через установленную команду `get-document-data` и напрямую
через скрипт. Во втором случае используйте алиас `--mode`:

```bash
python scripts/get_document_data.py --mode <chembl|pubmed|all> [опции]
```

Запуск справки подчёркивает новый флаг:

```text
$ python scripts/get_document_data.py --help
usage: get_document_data.py [-h] ... [--mode {pubmed,chembl,all}] {pubmed,chembl,all} ...

optional arguments:
  --mode {pubmed,chembl,all}
                        Acquisition mode alias for the positional sub-command.
```

Общие селекторы и их значения по умолчанию (эм-дэш означает отсутствие опции в
режиме):

| Опция | `chembl` | `pubmed` | `all` | Описание |
|-------|---------:|---------:|------:|----------|
| `--column` | `document_chembl_id` | `PMID` | `document_chembl_id` | Колонка с идентификаторами во входном CSV. |
| `--chunk-size` | `5` | — | `5` | Размер пачки для запросов ChEMBL. |
| `--batch-size` | — | `100` | `50` | Количество PMID в одном запросе PubMed. |
| `--sleep` | — | `5.0` | `5.0` | Пауза между запросами PubMed (секунды). |
| `--workers` | — | `1` | `1` | Число параллельных потоков для PubMed. |
| `--timeout` | `30.0` | — | `30.0` | Таймаут HTTP-запросов к ChEMBL (секунды). |
| `--limit` | `None` | `None` | `None` | Максимум идентификаторов; `0` пропускает шаг. |
| `--offset` | `0` | `0` | `0` | Количество строк, пропускаемых перед обработкой. |
| `--openalex-rps` | — | `None` | `None` | Override RPS для OpenAlex. |
| `--crossref-rps` | — | `None` | `None` | Override RPS для CrossRef. |

Фолбэки DOI работают в режимах обогащения и сохраняют идемпотентность запусков:

| Опция | По умолчанию | Режимы | Описание |
|-------|--------------|--------|----------|
| `--fallback-doi-csv` | `None` | `pubmed`, `all` | CSV с ручными соответствиями PMID → DOI. |
| `--fallback-doi-pmid-column` | `PMID` | `pubmed`, `all` | Колонка с PMID во фолбэк-файле. |
| `--fallback-doi-value-column` | `DOI` | `pubmed`, `all` | Колонка со значениями DOI во фолбэк-файле. |
| `--fallback-overwrite` / `--fallback-doi-overwrite` | `False` | `pubmed`, `all` | Принудительно заменять DOI значением из CSV. |

Примеры запуска:

```bash
# Минимальная выгрузка из ChEMBL
python scripts/get_document_data.py --mode chembl \
    --input data/input/document.csv \
    --final-out output/documents_chembl.csv \
    --config config/config.yaml

# Обогащение PubMed с ручными DOI
python scripts/get_document_data.py --mode pubmed \
    --input data/input/document.csv \
    --final-out output/documents_pubmed.csv \
    --fallback-doi-csv data/input/manual_doi.csv \
    --fallback-overwrite

# Полный цикл с настройкой namespace-ограничений
python scripts/get_document_data.py --mode all \
    --input data/input/document.csv \
    --final-out output/documents_full.csv \
    --openalex-rps 2 --crossref-rps 3
```

На выходе формируется детерминированный CSV, `<имя>.meta.yaml`, отчёты
`<имя>_quality_report_table.csv`, `<имя>_data_correlation_report_table.csv` и
`<имя>.quality.json` с показателями покрытия DOI.

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
    --final-out output/assays_$(date +%Y%m%d).csv \
    --batch-size 100 --timeout 60 --limit 200
```

Скрипт выгружает данные из ChEMBL, вычисляет счётчики по таргетам, нормализует и
валидирует таблицу, затем формирует стандартный набор артефактов.

## Пайплайн активностей (`get-activity-data`)

```
get-activity-data --input data/input/activity.csv \
    --final-out output/activities_$(date +%Y%m%d).csv \
    --batch-size 50 --workers 4 --timeout 60 --limit 500
```

Флаг `--dry-run` уникален для данного пайплайна: он проверяет аргументы и входной
файл и завершает работу без сетевых запросов. Пост-обработка вычисляет
`lower_value`/`upper_value` по алгоритму, описанному в [`docs/OUTPUT_RU.md`](./OUTPUT_RU.md).

## Пайплайн тест-объектов (`get-testitem-data`)

```
get-testitem-data --input data/input/testitem.csv \
    --final-out output/testitems_$(date +%Y%m%d).csv \
    --batch-size 1000 --timeout 60 --limit 400
```

Пайплайн объединяет молекулы из ChEMBL с данными PubChem, нормализует результат и
записывает стандартный комплект файлов.

## Вспомогательные утилиты

Каждый модуль экспортирует функцию `main(argv)` и подходит для использования как
через `python -m`, так и через зарегистрированные консольные команды.

| Модуль | Команда | Назначение |
|--------|---------|------------|
| `library.utils.cli_tools.check_determinism` | `check-determinism --log-level INFO` | Проверка детерминированности CSV через сравнение хэшей тестовых файлов. |
| `library.utils.cli_tools.chunk_io_main` | `chunk-io --input data.csv --final-out copy.csv` | Потоковое чтение/запись CSV с сохранением порядка. |
| `library.utils.cli_tools.csv_utils_main` | `csv-utils --input data.csv --final-out clean.csv --sep ,` | Нормализация разделителей, кавычек и порядка колонок. |
| `library.utils.cli_tools.dtype_inspector_main` | `python -m library.utils.cli_tools.dtype_inspector_main` | Диагностика типов pandas-таблиц. |
| `library.utils.cli_tools.get_activities` | `get-activities --limit 10` | Генерация тестовых активностей для проверки логов и CLI. |
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
  `--log-level DEBUG` для подробной диагностики без правки YAML.

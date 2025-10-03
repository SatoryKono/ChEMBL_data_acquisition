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

Любая команда возвращает ненулевой код при ошибках валидации, проблемах с IO или
неудачных запросах к внешним сервисам.

## Оркестратор (`get-data`)

```
get-data --base-path /data/chembl \
    --input-dir seeds --output-dir exports \
    --config /data/chembl/config.yaml \
    --date 20250101 --limit 100 --log-level INFO
```

Оркестратор готовит аргументы и запускает пайплайны в порядке: документы (`all`),
таргеты (`all`), ассайи, тест-объекты, активности. Значение `--limit 0` пропускает
выполнение, а `--dry-run` выводит план без записи файлов.

## Пайплайн документов (`get-document-data`)

Подкоманды:

| Режим | Описание | Основные параметры |
|-------|----------|--------------------|
| `chembl` | Выгружает метаданные документов из ChEMBL. | `--column`, `--chunk-size`, `--timeout`, `--limit`, `--offset`. |
| `pubmed` | Обогащает данными PubMed, Semantic Scholar, OpenAlex и CrossRef. | `--column`, `--sleep`, `--workers`, `--batch-size`, `--limit`, `--offset`, `--openalex-rps`, `--crossref-rps`, `--fallback-doi-*`. |
| `all` | Запускает `chembl`, объединяет внешние источники и формирует итоговую таблицу. | Параметры режима `pubmed` плюс `--fallback-doi-*` для CSV с ручными DOI. |

Пример:

```
get-document-data all \
    --input seeds/document_ids.csv \
    --final-out output/documents_$(date +%Y%m%d).csv \
    --config config/config.yaml \
    --limit 500 --log-level INFO
```

На выходе — детерминированный CSV, файл `<имя>.meta.yaml`, отчёты
`<имя>_quality_report_table.csv`, `<имя>_data_correlation_report_table.csv` и
`<имя>.quality.json` с покрытием DOI.

## Пайплайн таргетов (`get-target-data`)

Поддерживаются латинские и кириллические алиасы (`chembl`/`сруьид`, `uniprot`/`гтшзкще`, `iuphar`/`шгзрфк`, `all`/`фдд`).

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
    --input seeds/target_ids.csv \
    --final-out output/targets_$(date +%Y%m%d).csv \
    --raw-out output/targets_raw_$(date +%Y%m%d).parquet \
    --raw-format parquet --id-cols target_chembl_id uniprot_id \
    --config config/config.yaml --log-level INFO
```

Для офлайн-проверок используйте `python -m library.utils.cli_tools.pipeline_targets_main` —
он повторяет боевую логику, читая кешированные ответы.

## Пайплайн ассайев (`get-assay-data`)

```
get-assay-data --input seeds/assay_ids.csv \
    --final-out output/assays_$(date +%Y%m%d).csv \
    --batch-size 100 --timeout 60 --limit 200
```

Скрипт выгружает данные из ChEMBL, вычисляет счётчики по таргетам, нормализует и
валидирует таблицу, затем формирует стандартный набор артефактов.

## Пайплайн активностей (`get-activity-data`)

```
get-activity-data --input seeds/activity_ids.csv \
    --final-out output/activities_$(date +%Y%m%d).csv \
    --batch-size 50 --workers 4 --timeout 60 --limit 500
```

Флаг `--dry-run` уникален для данного пайплайна: он проверяет аргументы и входной
файл и завершает работу без сетевых запросов. Пост-обработка вычисляет
`lower_value`/`upper_value` по алгоритму, описанному в [`docs/OUTPUT_RU.md`](./OUTPUT_RU.md).

## Пайплайн тест-объектов (`get-testitem-data`)

```
get-testitem-data --input seeds/molecule_ids.csv \
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
| `library.utils.cli_tools.check_determinism` | `check-determinism --input a.csv --previous b.csv` | Сравнение SHA-256 и метаданных разных запусков. |
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

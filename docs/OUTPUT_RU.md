# Выходные артефакты

Документ описывает состав файлов, которые формируют пайплайны ChEMBL Data
Acquisition, а также вспомогательные модули, отвечающие за их создание.

## Структура каталогов

Экспорты сохраняются в `local.io.output_dir` (по умолчанию `data/output`).
Запись автоматически создаёт родительские каталоги, если `local.io.exist_ok`
равно `true`.

```
data/output/
└── ChEMBL/
    └── processed/
        ├── activity.csv
        ├── activity.csv.meta.yaml
        ├── activity_failure_cases.csv
        ├── activity_quality_report_table.csv
        ├── activity_data_correlation_report_table.csv
        └── ...
```

Промежуточные файлы таргет-пайплайна в режиме `all`
(`*_chembl.csv`, `*_uniprot.csv`, `*_iuphar.csv`) сохраняются в той же
директории, если CLI-аргументы не задают другие пути.

## Sidecar с метаданными (`*.csv.meta.yaml`)

Каждый CSV сопровождается `<name>.csv.meta.yaml`, который создаёт
`library.metadata.write_meta_yaml`. Файл содержит:

* `generated_at` — отметку времени в формате ISO 8601 (UTC).
* `git_sha` — хэш коммита на момент запуска.
* `python_version`, `platform` — сведения о среде выполнения.
* `command` — точную команду CLI.
* `config` — применённые настройки (секреты автоматически маскируются).
* `inputs` — описание входных файлов и аргументов.
* `stats` — `rows_total`, `rows_kept`, `rows_dropped` и контрольную сумму
  `output_sha256`.
* `schema` — имя схемы валидации.

Если sidecar уже существует, содержимое аккуратно объединяется, чтобы сохранить
ручные пометки.

## Артефакты валидации

* При ошибках Pandera строки попадают в `<stem>_failure_cases.csv` через
  `SidecarErrors`, что позволяет изучить проблемы без остановки основного
  пайплайна.
* `library.table_quality.analyze_table_quality` формирует
  `<stem>_quality_report_table.csv` и
  `<stem>_data_correlation_report_table.csv`. CLI-утилиты сохраняют отчёты рядом с
  выгрузкой, а `library.utils.cli_tools.get_input_initialisation` использует
  подкаталог `<output>/data_validity_report/`.

Все файлы пишутся в UTF-8 и соблюдают детерминированный порядок строк и колонок,
что упрощает сравнение версий.

## Детерминированная запись CSV

`library.io.write_csv` вызывает `library.io.writers.write_csv_deterministic`,
который сортирует строки и колонки по явным ключам и учитывает параметры
`cfg.io.csv_sep`, `cfg.io.csv_encoding`, а также опциональные `key_cols` и
`col_order`, переданные пайплайном.

## Технические колонки пайплайна

Перед валидацией в таблицы добавляются `pipeline_version` и `timestamp_utc`
(функция `library.pipeline_metadata.add_pipeline_metadata`). Колонки описаны в
схемах активностей, документов, тест-объектов и других сущностей, поэтому они
присутствуют даже при пустой выгрузке.

* `pipeline_version` — версия установленного пакета или значение из
  `pyproject.toml`; единое для всех таблиц в рамках запуска.
* `timestamp_utc` — время старта пайплайна (ISO 8601). Используйте как
  служебный маркер, а не временную метку строк.

## Классификация публикаций

`scripts/get_document_data.py` обогащает документные выгрузки баллами и метками,
которые рассчитывает `library.processing.document.merge_metadata`.
В таблице и схеме появляются поля:

| Колонка | Описание |
| --- | --- |
| `publication_types_normalised` | Упорядоченный через точку с запятой список типов публикаций из ChEMBL, PubMed, Semantic Scholar, OpenAlex и CrossRef. |
| `publication_type_score_review` | Целочисленный вес признаков обзора. |
| `publication_type_score_experimental` | Вес экспериментальных признаков. |
| `publication_type_score_unknown` | Вес неоднозначных или неизвестных типов. |
| `publication_class` | Итоговая метка (`review`, `experimental`, `unknown`), вычисленная классификатором с учётом порогов. |

При отсутствии известных токенов баллы равны нулю, а итоговая метка — `unknown`.

## Границы активности (`lower_value`, `upper_value`)

`activity.csv` включает диапазоны значений, рассчитанные из канонических полей
ChEMBL `standard_*`. Приоритет действий:

1. Использовать `standard_lower_value` и `standard_upper_value`, если оба
   присутствуют.
2. Дополнить отсутствующую границу из пары `standard_value` + явная нижняя/верхняя граница.
3. При `activity_bounds.enable_from_relation = true` интерпретировать символы
   отношений (`=`, `≈`, `>=`, `<=`, `between`, `range`).
4. При включённом `activity_bounds.enable_from_uncertainty` разобрать записи вида
   `значение ± дельта` из `standard_text_value`.

Результаты округляются до `activity_bounds.rounding_digits` знаков (по умолчанию
`3`) и обрезаются до нуля для концентрационных метрик, если
`activity_bounds.clamp_nonnegative = true`.

## JSON-отчёт о качестве документов

Документный пайплайн записывает `<stem>.quality.json`. `library.processing.document.build_quality_report`
формирует сводку, а `save_quality_report` сохраняет её в стабильном формате для
сравнения между запусками. Структура отчёта:

| Ключ | Тип | Описание |
| --- | --- | --- |
| `rows_total` | Целое | Количество строк в выгрузке. |
| `doi_coverage` | Число | Доля документов с DOI (0.0–1.0). |
| `publication_class_counts` | Объект | Количество строк по меткам `review`, `experimental`, `unknown`; отсутствующие категории представляются нулём. |
| `error_counts` | Объект | Число ошибок по источникам (`pubmed`, `semantic_scholar`, `openalex`, `crossref`). |

## Экспорт тест-объектов

`scripts/get_testitem_data.py` формирует `testitem.csv` с метаданными и отчётами
качества. Каждая строка объединяет данные ChEMBL, обогащение PubChem и служебные
колонки, обеспечивая опорное измерение соединений. Для работы требуется каталог
родительских молекул (`sources.chembl.molecule_catalog.cache_path`) и, по
возможности, CSV-словари с иерархией и признаками солей.

## Кешированный таргет-пайплайн (`library.utils.cli_tools.pipeline_targets_main`)

Утилита `library.utils.cli_tools.pipeline_targets_main` повторяет аргументы боевого таргет-пайплайна,
но работает только с кешированными чанками ChEMBL. Идентификаторы читаются через
`read_ids`, далее вызывается `library.pipeline_targets.run_pipeline`, добавляются
метаданные и выполняется детерминированная запись CSV/sidecar. Инструмент полезен
для проверки конфигурации и параметров батчирования без обращений к внешним
сервисам.

## Рекомендации по сопровождению

* Храните результаты в папках с датой (`YYYYMMDD/`) — это облегчает сравнение.
* Архивируйте устаревшие выгрузки: sidecar сохраняет достаточно информации для
  воспроизводимости.
* Контролируйте свободное место при длительных или параллельных запусках.

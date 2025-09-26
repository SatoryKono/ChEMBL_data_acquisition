# Структура выходных данных

## Базовый каталог

Все наборы сохраняются в `local.io.output_dir` (по умолчанию `data/output`). Если аргумент `--output` не задан, функции
используют `library.io.default_output_path` и формируют имя `output_<stem входного файла>_<YYYYMMDD>.csv` в указанном каталоге.
При `local.io.exist_ok=true` необходимые папки создаются автоматически.

Пример:

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

## Метаданные sidecar

Каждый CSV сопровождается `<name>.csv.meta.yaml`, который записывает `library.metadata.write_meta_yaml`. Файл включает:

* `generated_at` — отметку времени в формате ISO 8601 (UTC).
* `git_sha` — хэш коммита репозитория на момент запуска.
* `python_version`, `platform` — сведения о среде выполнения.
* `command` — точную команду запуска.
* `config` — релевантные настройки (секреты маскируются).
* `inputs` — описание исходных файлов и аргументов.
* `stats` — показатели `rows_total`, `rows_kept`, `rows_dropped`, а также контрольную сумму `output_sha256`.
* `schema` — имя схемы валидации, применённой к данным.

Если sidecar уже существует, новый набор полей аккуратно объединяется с текущим содержимым.

## Артефакты валидации

* При ошибках проверки Pandera проблемные строки сохраняются в `<stem>_failure_cases.csv` рядом с основным файлом.
* `library.table_quality.analyze_table_quality` создаёт отчёты `<stem>_quality_report_table.csv` и
  `<stem>_data_correlation_report_table.csv`. CLI-утилиты размещают их рядом с экспортом, а
  `scripts/get_input_initialisation.py` — в подкаталоге `<output>/data_validity_report/`.

Все отчёты записываются в кодировке UTF-8 и наследуют детерминированный порядок строк и колонок.

## Детерминированный экспорт

`library.io.write_csv` вызывает `library.csv_utils.write_csv_deterministic`, сортируя колонки и строки по ключевым столбцам.
Поведение следует настройкам `cfg.io.csv_sep`, `cfg.io.csv_encoding` и учитывает аргументы `key_cols`/`col_order`.

## Метаданные пайплайна

Все выгрузки дополняются двумя служебными колонками, которые добавляет `library.pipeline_metadata.add_pipeline_metadata` до валидации и записи CSV: `pipeline_version` фиксирует версию установленного пакета (или значение из `pyproject.toml`), а `timestamp_utc` содержит время запуска в формате ISO 8601.【F:library/pipeline_metadata.py†L24-L84】 Колонки описаны в схемах валидации для активностей, документов, тест-объектов и других сущностей, поэтому они присутствуют даже при пустом наборе данных.【F:schemas/activities.py†L52-L55】【F:schemas/documents.py†L111-L112】【F:schemas/testitems.py†L41-L42】

`pipeline_version` неизменен в пределах одного запуска и подходит для соединений между таблицами, полученными в той же сессии. `timestamp_utc` отражает время оркестратора, поэтому его следует трактовать как служебный маркер, а не временную метку строк.

## Колонки классификации публикаций

`scripts/get_document_data.py` дополняет объединённые метаданные детерминированными баллами и метками, которые возвращает `library.document_pipeline.merge_metadata`.【F:scripts/get_document_data.py†L607-L676】【F:library/document_pipeline.py†L160-L208】 В итоговом `document.csv` и схеме валидации присутствуют следующие поля:

| Колонка | Описание |
| --- | --- |
| `publication_types_normalised` | Список уникальных типов публикаций, собранных из ChEMBL, PubMed, Semantic Scholar, OpenAlex и CrossRef. Значения сортируются и разделяются точкой с запятой для устойчивых диффов. |
| `publication_type_score_review` | Целочисленный вес, накопленный при голосовании за обзоры. Баллы неотрицательны и растут при подтверждении из нескольких источников. |
| `publication_type_score_experimental` | Целочисленный вес для экспериментальных материалов; вычисляется тем же алгоритмом, что и для обзоров. |
| `publication_type_score_unknown` | Целочисленный вес, отражающий неявные или явно неизвестные типы. |
| `publication_class` | Итоговая метка (`review`, `experimental` или `unknown`), выбранная по результатам сравнения весов с порогами классификатора.【F:library/document_type_classifier.py†L7-L74】 |

При отсутствии известных терминов все баллы равны `0`. Классификатор выбирает максимальный балл, превышающий минимальные пороги, и возвращает `unknown`, если сигнал неоднозначен — не стоит трактовать колонку как показатель качества ручной модерации.

## Границы активности (`lower_value`, `upper_value`)

В `activity.csv` появились поля `lower_value` и `upper_value`, вычисляемые после нормализации в `scripts/get_activity_data.py`. Алгоритм использует только канонические поля `standard_*`, уже приведённые ChEMBL к общим единицам. Очерёдность источников для каждой строки:

1. Явные границы `standard_lower_value` и `standard_upper_value`.
2. Пары `standard_value` + `standard_upper_value`: минимальное значение заполняет `lower_value`, максимальное — `upper_value`, если соответствующая граница ещё пуста.
3. Интерпретация отношений при `activity_bounds.enable_from_relation = true`: `=`/`≈`/`~` задают обе границы, `>=` — только `lower_value`, `<=` — только `upper_value`, `between`/`range` требуют второго нормализованного числа. Неизвестные маркеры оставляют поля пустыми и логируются. Если доступно только исходное `value` без нормализованного аналога, выводится предупреждение `activity_bounds_missing_standard_value`.
4. При включённом `activity_bounds.enable_from_uncertainty` разбираются записи вида `значение ± дельта` из `standard_text_value` при наличии канонического значения.

Результат округляется до `activity_bounds.rounding_digits` знаков (по умолчанию `3`) и обрезается до нуля для концентрационных метрик, когда `activity_bounds.clamp_nonnegative = true`. Определение таких метрик опирается на эвристики по `standard_type` и `standard_units`. Остальные столбцы не меняются, порядок колонок остаётся детерминированным согласно схеме.【F:scripts/get_activity_data.py†L1-L234】【F:config.yaml†L108-L147】【F:library/config.py†L358-L420】【F:schemas/activities.py†L32-L64】

## JSON-отчёт о качестве документов

`scripts/get_document_data.py` дополнительно формирует файл `<stem>.quality.json`. Хелпер
`library.document_pipeline.build_quality_report` собирает общую численность строк, долю заполненных DOI, распределение классов
публикаций и счётчики ошибок по источникам; `save_quality_report` записывает структуру в UTF-8 JSON со стабильным форматированием для удобства сравнения.【F:scripts/get_document_data.py†L636-L674】【F:library.document_pipeline.py†L300-L356】 Отчёт позволяет контролировать полноту метаданных без чтения основного CSV, например отслеживать падение покрытия DOI.

JSON-файл содержит:

| Ключ | Тип | Описание |
| --- | --- | --- |
| `rows_total` | Целое | Общее число строк в выгрузке. |
| `doi_coverage` | Число с плавающей точкой | Доля документов с заполненным DOI (в диапазоне 0.0–1.0). |
| `publication_class_counts` | Объект | Словарь с количествами по меткам `review` / `experimental` / `unknown`; при отсутствии значения используется `unknown`. |
| `error_counts` | Объект | Словарь с числом ошибок обогащения по источникам (`pubmed`, `semantic_scholar`, `openalex`, `crossref`). |

Все ключи присутствуют всегда, отсутствующие категории представлены нулевыми счётчиками, поэтому мониторинги могут опираться на стабильную схему между запусками.

## Рекомендации по ведению архива

* Выделяйте под каждый запуск отдельную папку с датой (`YYYYMMDD/`) для удобного сравнения.
* Архивируйте устаревшие результаты — метаданные сохраняют достаточно информации для воспроизводимости.
* Следите за свободным местом при длительных или параллельных выгрузках.

## Экспорт тест-айтемов

`scripts/get_testitem_data.py` формирует `testitem.csv`, дополняя его стандартными `*.meta.yaml`,
опциональными `*_failure_cases.csv` и отчётами качества с детерминированным порядком, описанным выше.【F:scripts/get_testitem_data.py†L151-L299】
Каждая строка объединяет поля ChEMBL (`molecule_chembl_id`, структурные дескрипторы, статусные признаки),
обогащение PubChem и метаданные запуска, превращая выгрузку в опорное измерение соединений.【F:scripts/get_testitem_data.py†L36-L193】【F:schemas/testitems.py†L12-L31】

Перед распространением данных выполните объединение с каталогом родительских молекул, чтобы добавить
`parent_molecule_chembl_id` для агрегаций. Отображение хранится в JSON по пути
`sources.chembl.molecule_catalog.cache_path` и загружается функцией
`library.molecule_catalog.load_parent_catalog`, которая при необходимости обновляет кэш через API ChEMBL.【F:config.yaml†L25-L33】【F:library/molecule_catalog.py†L43-L136】

Дополнительный этап обогащения использует `dictionary/molecule_hierarchy.csv` и
`dictionary/molecule_catalog.csv`, чтобы до валидации заполнить `salt_chembl_id`
и булевы признаки `natural_product`, `prodrug`, `polymer_flag`. Соль
распознаётся, когда `parent_molecule_chembl_id` отличается от
`molecule_chembl_id`, после чего в `salt_chembl_id` попадает идентификатор
дочерней молекулы. Флаги нормализуются в nullable-тип pandas `boolean`, а
пропуски у дочерних записей при необходимости подтягиваются из родителя. Если
молекула отсутствует в словарях, лог пишется с предупреждением для диагностики.【F:scripts/get_testitem_data.py†L205-L233】【F:library/testitem_enrichment.py†L17-L216】

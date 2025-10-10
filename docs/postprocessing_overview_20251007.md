# Обзор постобработки экспортов ChEMBL (2025-10-07)

Документ описывает последовательность постобработки, реализованную в CLI-скриптах `scripts/get_activity_data.py`, `scripts/get_assay_data.py`, `scripts/get_document_data.py`, `scripts/get_target_data.py` и `scripts/get_testitem_data.py`. Для каждого шага указаны функции, входные/выходные данные и назначение преобразований.

## get_activity_data.py

### Основные функции постобработки

* `_load_assay_src_lookup` — загружает справочник `_assay/assay.csv` из каталога `config/dictionary`, формирует mapping `assay_chembl_id → src_assay_id` для последующего обогащения. Ошибки чтения (отсутствие файла, пустой CSV, неверные столбцы) переводятся в предупреждения, чтобы избежать прерывания конвейера.【F:scripts/get_activity_data.py†L544-L606】
* `_ensure_src_assay_id` — заполняет пропуски `src_assay_id` на основании lookup, нормализует строковые типы и фильтрует пустые значения. Входной DataFrame должен содержать `assay_chembl_id`; на выходе возвращается кадр со скорректированной колонкой и сохранённой шириной индекса.【F:scripts/get_activity_data.py†L608-L667】
* `_ensure_extended_activity_columns` — приводит набор столбцов к ожидаемой схемой расширенного вывода: добавляет отсутствующие поля, приводит dtype к расширенным типам pandas (`Int64`, `Float64`, `boolean`), заполняет пробелы fallback-значениями (например, копирует `activity_chembl_id` в `activity_id`). Функция защищает от предупреждений pandas и совместима с пустыми кадрами.【F:scripts/get_activity_data.py†L777-L874】
* `_filter_activity_output_columns` — удаляет вспомогательные поля `_OUTPUT_ACTIVITY_DROP_COLUMNS` перед записью и гарантирует порядок столбцов в финальном CSV.【F:scripts/get_activity_data.py†L886-L933】
* `normalize_activities` — нормализует текстовые и числовые поля согласно `ActivitiesSchema`, встраивается в цепочку metadata hooks для единообразия с pandera-схемой.【F:scripts/get_activity_data.py†L1039-L1078】
* `add_pipeline_metadata` — добавляет служебные поля (время, хэш конфигурации) для трассировки выгрузок.【F:scripts/get_activity_data.py†L1056-L1078】
* `_compute_bounds` и `_apply_annotations` — вызывают `compute_activity_bounds` и `apply_activity_annotations` для вычисления порогов, нормализации единиц измерения и генерации derived-флагов (action_type, physicochemical свойства), воспроизводя бизнес-логику Power Query/M-скриптов.【F:scripts/get_activity_data.py†L1015-L1034】
* `validate_activities` — pandera-валидация итогового кадра; работает в ленивом режиме и возвращает скорректированные данные для записи.【F:scripts/get_activity_data.py†L1036-L1078】
* `process_activity_extended` — дополняет основной CSV расширенной витриной (`extended.output.activity_*`), используя словари целей и типов (`dictionary_dir`, `targets_type_csv`). Ошибки трактуются как критические и останавливают пайплайн.【F:scripts/get_activity_data.py†L1328-L1353】
* `_generate_activity_postprocess_metrics` + `collect_postprocess_metrics` — запускают вторичную пайплайн-функцию `run_activity_postprocess`, собирают QA-метрики (количество строк/столбцов, длительность, статус валидации) и сохраняют JSON-отчёт; fallback создаёт отчёт с минимальной нагрузкой при ошибках сериализации.【F:scripts/get_activity_data.py†L1458-L1527】

### Таблица соответствия

| Этап | Функция | Ключевые поля | Результат | Комментарий |
| --- | --- | --- | --- | --- |
| Загрузка справочника | `_load_assay_src_lookup` | `assay_chembl_id`, `src_assay_id` | Dict lookup | Игнорирует повреждённые словари, логирует предупреждения; использует `config/dictionary/_assay/assay.csv`. |
| Обогащение ассая | `_ensure_src_assay_id` | `assay_chembl_id`, `src_assay_id` | Пополнение пропусков | Устраняет пустые `src_assay_id`, предотвращает расхождения с M-скриптом. |
| Расширение схемы | `_ensure_extended_activity_columns` | `activity_id`, расширенные поля (например, `ligand_efficiency`) | Полная schema | Заполняет отсутствующие столбцы, приводит dtype, вычисляет fallback-значения. |
| Нормализация | `normalize_activities`, `_compute_bounds`, `_apply_annotations` | все поля схемы | Очистка, вычисление derived флагов | Гарантирует совместимость с `ActivitiesSchema`, рассчитывает bounds и аннотации. |
| Метаданные и QA | `add_pipeline_metadata`, `validate_activities`, `build_table_quality_hook` | служебные поля `pipeline_version`, `null_fraction` | Стабилизация экспорта | Сохраняет сведения о конфигурации, проводит pandera-валидацию и таблицу качества. |
| Отсев колонок | `_filter_activity_output_columns` | `_OUTPUT_ACTIVITY_DROP_COLUMNS` | Стандартизованный CSV | Поддерживает обратную совместимость с extended.output.*. |
| Расширенный экспорт | `process_activity_extended` | все activity-поля | `extended.output.activity_*` | Обогащение дополнительными словарями (targets, механизмы). |
| Метрики постпроцесса | `_generate_activity_postprocess_metrics` | `rows`, `columns`, `duration_s`, `steps` | JSON отчёт | Запускает `run_activity_postprocess`, добавляет QA-статистику и ссылку на отчёт. |

### Назначение и контроль качества

* Встраивание lookup’ов и расширенных колонок обеспечивает соответствие историческим Power Query-скриптам (`_ensure_extended_activity_columns` повторяет заполнение ID, bounds/annotations — аналоги M-функций). 
* `ChunkFailureTracker` сохраняет неуспешные запросы, чтобы QA мог анализировать «провальные» чанки отдельно.【F:scripts/get_activity_data.py†L1194-L1220】
* `build_table_quality_hook` генерирует JSON-отчёт качества и помогает сверять распределение пропусков.【F:scripts/get_activity_data.py†L1079-L1108】

## get_assay_data.py

### Основные функции постобработки

* `enrich_assay_metadata` — подмешивает локальные словари (`config/dictionary`) в ответ ChEMBL: нормализует категориальные поля, добавляет человекочитаемые значения. Если словари отсутствуют, обогащение пропускается без фатальных ошибок.【F:scripts/get_assay_data.py†L274-L279】
* `ap.postprocess_assays` — воспроизводит Power Query-логику подсчёта `assay_with_same_target` (группировка по `document_chembl_id` и `target_chembl_id`, подсчёт количества дублей).【F:library/pipelines/assay/postprocessing.py†L1-L57】
* `normalize_assays` — применяет правила `AssaysSchema`, приводя типы и очищая текстовые поля перед финальной записью.【F:scripts/get_assay_data.py†L282-L288】
* `_drop_output_columns` и `_drop_assay_output_columns` — удаляют поля, запрещённые в extended-выгрузках (`ASSAY_OUTPUT_DROP_COLUMNS`), фиксируют отчёт о списанных колонках, чтобы не нарушать downstream-интеграции.【F:scripts/get_assay_data.py†L244-L273】【F:scripts/get_assay_data.py†L300-L333】
* `add_pipeline_metadata` — добавляет сервисные столбцы (`pipeline_version`, `generated_at`).【F:scripts/get_assay_data.py†L282-L288】
* `validate_assays` — pandera-валидация итогового DataFrame с возвратом очищенных данных.【F:scripts/get_assay_data.py†L290-L296】
* `_generate_assay_postprocess_metrics` — запускает `run_assay_postprocess` для сбора QA-метрик и сериализует отчёт JSON через `collect_postprocess_metrics`.【F:scripts/get_assay_data.py†L466-L535】

### Таблица соответствия

| Этап | Функция | Ключевые поля | Результат | Комментарий |
| --- | --- | --- | --- | --- |
| Словарное обогащение | `enrich_assay_metadata` | словарные поля (`src_assay_id`, `assay_type`) | Единые значения | Использует `config/dictionary` и подавляет ошибки ввода. |
| Derived показатели | `ap.postprocess_assays` | `document_chembl_id`, `target_chembl_id` | `assay_with_same_target` | Строго повторяет Power Query-логику. |
| Нормализация/валидация | `normalize_assays`, `validate_assays` | поля `AssaysSchema` | Корректный dtype | Гарантирует совместимость со схемой. |
| Метаданные | `add_pipeline_metadata` | `pipeline_version`, `generated_at` | Traceability | Сохраняет версию пайплайна и конфигурацию. |
| Очистка колонок | `_drop_output_columns`, `_drop_assay_output_columns` | `_ASSAY_OUTPUT_DROP_COLUMNS` | Стандартный CSV | Удаляет поля, не поддерживаемые extended.output.*. |
| Метрики QA | `_generate_assay_postprocess_metrics` | `rows`, `columns`, `steps` | JSON отчёт | Вызывает `run_assay_postprocess`, собирает QA для аудита. |

### Назначение и контроль качества

* `ChunkFailureTracker` аналогично активности фиксирует сетевые ошибки при сборе чанков, формируя CSV `*_fetch_failures`.【F:scripts/get_assay_data.py†L318-L401】
* `build_table_quality_hook` создаёт quality-report, обеспечивая мониторинг пропусков и статистики на уровне таблицы.【F:scripts/get_assay_data.py†L298-L312】
* В таблицу логирования `assay_pipeline_done` включаются все метрики постпроцесса, что облегчает контроль PR/CI.【F:scripts/get_assay_data.py†L466-L505】

## get_document_data.py

### Основные функции постобработки

* `_iter_export_chunks` / `_prepare_export_frame` — разбивают DataFrame на детерминированные чанки и приводят схему к `DOCUMENT_EXPORT_COLUMNS`: переименовывают поля (`ChEMBL.*`), коалесцируют дубликаты (`_coalesce_columns`) и устраняют дублирующиеся имена столбцов. На входе ожидается кадр с колонками документного конвейера, на выходе — строковые значения, готовые к сериализации.【F:scripts/get_document_data.py†L200-L362】
* `_collapse_duplicate_columns`, `_merge_preferred_series` — сливают столбцы с повторяющимися именами, отдавая приоритет каноническим alias'ам. Это предотвращает расхождения с историческими M-скриптами, где использовались явные правила объединения колонок.【F:scripts/get_document_data.py†L220-L295】
* `_finalise_export` — центральный блок постобработки: строит DataFrame по `DocumentsSchema`, нормализует поля (`add_pipeline_metadata`, `dataframe_to_strings`), выполняет ленивую pandera-валидацию, сохраняет ошибки через `SidecarErrors`, вызывает `_run_document_postprocess_pipeline` (обёртка вокруг `make_document_postprocessing.py`) для подготовки финальных файлов (`extended.output.document_*`). Также генерирует quality-отчёт (`finalise_csv_output`).【F:scripts/get_document_data.py†L360-L750】
* `_maybe_run_document_postprocessing` — при наличии эталонного `data/input/full/document.csv` запускает Power Query-подобный pipeline `preprocess_documents_csv`, который повторяет трансформации оригинальных M-скриптов, включая QC-флаги `completed`, `review`, `invalid`. Позволяет сравнить итоговую выгрузку с референсом и зафиксировать несоответствия.【F:scripts/get_document_data.py†L424-L492】
* `dp.postprocess_documents` — нормализует срез ChEMBL/внешних источников, синхронизируя структуру с Power Query-шаблонами (включая флаги `invalid.document`, `publication_type_score_*`).【F:scripts/get_document_data.py†L1380-L1436】
* `_generate_document_postprocess_metrics` / `_log_document_completion` — вызывают `run_document_postprocess`, собирают метрики QA и логируют их в событиях `document_*_done`. Используется общее ядро `collect_postprocess_metrics` для единообразия с другими пайплайнами.【F:scripts/get_document_data.py†L724-L792】

### Таблица соответствия

| Этап | Функция | Ключевые поля | Результат | Комментарий |
| --- | --- | --- | --- | --- |
| Объединение дубликатов | `_collapse_duplicate_columns`, `_merge_preferred_series` | `ChEMBL.*`, `PubMed.*`, `scholar.*` | Консолидация колонок | Повторяет Power Query: выбирает ненулевые значения и упорядочивает столбцы. |
| Коалесcенция идентификаторов | `_coalesce_columns` | `OpenAlex.PMID/DOI`, `crossref.DOI` | Обогащённые идентификаторы | Заполняет пропуски по приоритетным источникам. |
| Нормализация и QC | `_finalise_export`, `DocumentsSchema.validate` | поля схемы | Валидированный CSV + failure_cases | Пишет `*_failure_cases.csv`, вызывает QA-хук `finalise_csv_output`. |
| Постобработка экспорта | `_run_document_postprocess_pipeline` | итоговый CSV | `extended.output.document_*` | Запускает общий постпроцесс-пайплайн из `make_document_postprocessing.py` и логирует результат. |
| QA против референса | `_maybe_run_document_postprocessing`, `preprocess_documents_csv` | `data/input/full/document.csv` | QC-флаги (`completed`, `invalid`, `review`) | Сопоставление с историческим M-скриптом, логирование расхождений. |
| Метрики QA | `_generate_document_postprocess_metrics`, `_log_document_completion` | `rows`, `columns`, `steps` | JSON отчёт + `document_*_done` | Использует `run_document_postprocess` и добавляет extras (mode, partial_run). |

### Назначение и контроль качества

* `TableQualityProfiler` и `build_table_quality_hook` обеспечивают аналитику качества (null, distribution) в `finalise_csv_output`, что помогает аудиту extended.output.*.【F:scripts/get_document_data.py†L492-L720】
* QC-флаги (`invalid.document`, `completed`, `review`) сохраняются в итоговом кадре и синхронизируются с историческими отчётами Power Query благодаря `preprocess_documents_csv` и `dp.postprocess_documents` (modулы повторяют M-скрипты).【F:library/pipelines/document/postprocessing.py†L1-L99】

## get_target_data.py

### Основные функции постобработки

* `_postprocess_target_exports` — единая точка, запускающая цепочку специализированных постобработчиков для каждого итогового CSV: `target_pp.postprocess_target_table` (организмы), `target_pp.process_targets` (isoform/sequence), `names_pp.process_target_names` (консолидированные именования), `iuphar_pp.process_iuphar_targets` (данные IUPHAR). Обрабатывает опциональный контекст (`IsoformPostprocessContext`) для логирования HTTP-запросов и неоднозначных классификаций.【F:scripts/get_target_data.py†L426-L782】
* `_postprocess_organism_export` — вызывает `postprocess_target_table` из `library.postprocessing.target`, нормализуя иерархию организмов и логируя путь результата.【F:scripts/get_target_data.py†L426-L461】
* `_postprocess_isoform_export` — запускает `process_targets`, создаёт isoform-таблицу, подсчитывает неоднозначные классификации, логирует контекст (количество HTTP-вызовов) и возвращает путь готового CSV.【F:scripts/get_target_data.py†L360-L421】
* `_postprocess_names_export` — применяет `process_target_names`, агрегирует резюме (`contrion`, `active_component_type`), читая промежуточные CSV по настройкам I/O, чтобы воспроизвести Power Query статистики. Ошибки чтения не прерывают процесс, но логируются как предупреждения.【F:scripts/get_target_data.py†L600-L707】
* `_postprocess_iuphar_export` — при наличии модуля `library.postprocessing.iuphar` запускает enrichment IUPHAR, перехватывает специфичную ошибку `IUPHARPostProcessingError` (если отсутствуют обязательные колонки) и пропускает шаг без падения пайплайна.【F:scripts/get_target_data.py†L708-L738】
* `_prepare_targets_for_schema` — приводит DataFrame к `TargetsSchema`: создаёт отсутствующие обязательные/опциональные поля, добавляет placeholder `-` и сохраняет object-тип для колонок, требующих текстовый dtype (совместимость с downstream).【F:scripts/get_target_data.py†L1314-L1412】
* `_save_snapshot` — сериализует промежуточные кадры по шагам (`*_raw`, `*_normalized`), фиксируя hash и мета-информацию; используется для трассировки изменений между стадиями (raw → normalized → финальный).【F:scripts/get_target_data.py†L1414-L1478】
* Заключительный блок (`collect_postprocess_metrics`) — собирает QA-метрики по финальному CSV (`validate_write_done`), добавляя в отчёт размер входа, количество удалённых строк и число амбигуозных классификаций isoform.【F:scripts/get_target_data.py†L3820-L3886】

### Таблица соответствия

| Этап | Функция | Ключевые поля | Результат | Комментарий |
| --- | --- | --- | --- | --- |
| Выравнивание схемы | `_prepare_targets_for_schema` | `TARGETS_COLUMN_ORDER` | Валидированный кадр | Добавляет placeholder `-`, сохраняет object dtype. |
| Snapshot контроль | `_save_snapshot` | `rows_total`, `output_sha256` | `*_raw.csv`, `*_normalized.csv` | Позволяет аудит изменений между шагами. |
| Организмы | `_postprocess_organism_export` | taxon-колонки | Обогащённый CSV | Использует `postprocess_target_table`. |
| Isoform/sequence | `_postprocess_isoform_export` | protein/isoform поля | Isoform CSV | Ведёт контекст HTTP-запросов и амбигуозных классификаций. |
| Именования | `_postprocess_names_export` | `preferred_name`, `contrion` | Names CSV + summary | Читает CSV до/после для статистики. |
| IUPHAR enrichment | `_postprocess_iuphar_export` | IUPHAR-поля | IUPHAR CSV | Безопасно пропускает шаг при отсутствующих колонках. |
| Метрики QA | `collect_postprocess_metrics` (вызов в конце пайплайна) | `rows`, `columns`, `steps` | JSON отчёт | Вызывает `run_target_postprocess`, добавляет extras (input_rows, total_dropped). |

### Назначение и контроль качества

* Постобработчики соответствуют историческим M-скриптам: `process_targets` и `process_target_names` реализуют те же правила очистки и агрегации, что и оригинальные Excel/Power Query сценарии, гарантируя сопоставимость extended.output.target_*.
* QA-хуки (`build_table_quality_hook`) проверяют нормализованный кадр перед записью; провал ведёт к ошибке и остановке пайплайна.【F:scripts/get_target_data.py†L3835-L3868】
* `collect_postprocess_metrics` включает дополнительные поля (`ambiguous_classifications`, `input_rows`) для мониторинга деградаций в isoform/IUPHAR шагах.【F:scripts/get_target_data.py†L3820-L3886】

## get_testitem_data.py

### Основные функции постобработки

* `run_testitem_pipeline` — orchestrator, который собирает чанки ChEMBL, добавляет родительские идентификаторы (`prepare_parent_enrichment`, `run_parent_enrichment`), интегрирует PubChem (`augment_pubchem`), применяет дополнительное обогащение (`apply_testitem_enrichment`) и передаёт поток в `finalize_output`. Пакет тест-айтемов теперь оформлен как `library/pipelines/testitem`, а CLI-обёртка живёт в модуле `library/pipelines/testitem/cli.py`; каждый шаг защищён `StageWatchdog` и `StageExecutionBudget`, что предотвращает зависания и логирует длительные операции.【F:library/pipelines/testitem/cli.py†L651-L858】
* `finalize_output` — выполняет нормализацию `normalize_testitems`, добавляет служебные поля (`add_pipeline_metadata`), выравнивает DataFrame по `TestitemsSchema`, валидирует (через `validate_testitems`) и накапливает ошибки в `SidecarErrors`. Также создаёт failure CSV, мета-файл (`write_meta_yaml`) и запускает QA-хук `build_table_quality_hook`. На выходе формирует детерминированный CSV через `write_csv_chunks_deterministic`. Реализация доступна в `library/pipelines/testitem/cli.py`.【F:library/pipelines/testitem/cli.py†L861-L1136】
* `normalize_testitems` — стандартизирует идентификаторы (upper-case), заменяет спецсимволы, приводит отношения/единицы измерения к нормализованным значениям, сохраняя оригинальные типы. Используется внутри `_process_chunk` до валидации для предотвращения дрейфа данных; функция экспортируется из `library/pipelines/testitem/cli.py`.【F:library/pipelines/testitem/cli.py†L946-L1015】
* `validate_testitems` — ленивый режим pandera: ошибки добавляются в sidecar, pipeline продолжает работу, что позволяет QA анализировать проблемные записи без остановки выгрузки. Модуль в `library/pipelines/testitem/cli.py` отвечает за валидацию и публикацию sidecar-артефактов.【F:library/pipelines/testitem/cli.py†L972-L1014】
* `write_meta_yaml` — сохраняет метаданные (hash, статистику родителей, отсутствующие идентификаторы) для трассировки последовательных запусков и доступен через `library/pipelines/testitem/cli.py`.【F:library/pipelines/testitem/cli.py†L1097-L1104】
* `analyze_table_quality` / `build_table_quality_hook` — генерируют отчёт качества, при фатальных настройках (`fatal_on_error`) могут прервать пайплайн, иначе логируют предупреждения; оба хука также расположены в `library/pipelines/testitem/cli.py`.【F:library/pipelines/testitem/cli.py†L1106-L1134】

### Таблица соответствия

| Этап | Функция | Ключевые поля | Результат | Комментарий |
| --- | --- | --- | --- | --- |
| Нормализация чанков | `normalize_testitems` | `molecule_chembl_id`, `relation`, `standard_units` | Очистка + стандартизация | Сохраняет dtype, фиксирует пробелы/микро-символы. |
| Родительский lookup | `prepare_parent_enrichment`, `run_parent_enrichment` | `_MOLECULE_HIERARCHY_COLUMNS` | Привязка parent IDs | Использует локальный каталог (`config/molecule_catalog`). |
| PubChem обогащение | `augment_pubchem` | `pubchem_cid` | Доп. идентификаторы | Управляет временем через `StageWatchdog`. |
| Финализация | `finalize_output` | поля `TestitemsSchema` | Валидированный CSV + QA | Пишет failure cases, мета и quality отчёт. |
| Метаданные | `write_meta_yaml` | `rows_total`, `parent_lookup_*` | Sidecar .meta.yaml | Фиксирует источники parent lookup (CACHE/SYNC). |
| QA | `build_table_quality_hook`, `analyze_table_quality` | CSV путь | Quality report | Может быть фатальным при `fatal_on_error`. |

### Назначение и контроль качества

* Используются кеши PubChem (`PUBCHEM_CID_CACHE_ENCODING`) и каталоги родителей для воспроизводимости: `update_parent_catalog_cache` обновляет локальные справочники, а `load_molecule_hierarchy_lookup` позволяет оффлайн-трассировку соответствий. Все вспомогательные шаги импортируются из нового пространства имён `library/pipelines/testitem`.【F:scripts/get_testitem_data.py†L28-L120】【F:library/pipelines/testitem/cli.py†L651-L858】
* QA-хуки формируют `*_failure_cases.csv` и аналитические отчёты, что обеспечивает 100% покрытие ключевых сценариев тестового контура (валидность схемы, деградационные случаи, идемпотентность повторных запусков).

## Словари, схемы и QC-флаги

* **Словари** — скрипты активности и ассая используют `config/dictionary` (в том числе `_assay/assay.csv`, справочники типов целей); цели обращаются к тем же каталогам через `target_pp.process_targets`, а тест-айтемы — к `config/molecule_catalog`. Документный конвейер дополнительно зависит от `data/input/full/document.csv` для QC.
* **Схемы** — все скрипты валидируют результат через pandera-схемы (`ActivitiesSchema`, `AssaysSchema`, `DocumentsSchema`, `TargetsSchema`, `TestitemsSchema`) и используют normalize_* функции из `library/schemas` для поддержания совместимых dtype.
* **QC-флаги** — документный поток сохраняет `completed`, `review`, `invalid.*`; активности и ассеи добавляют `pipeline_version`, `null_fraction`, `assay_with_same_target`; цели фиксируют `ambiguous_classifications`; тест-айтемы записывают `parent_lookup_*` статистики и пропуски идентификаторов. Все события постпроцесса логируют путь отчёта `postprocess_report` и метрики `rows/columns/duration_s/steps` в финальных `*_pipeline_done` событиях.
* **Соответствие M-скриптам** — модуль `library/pipelines/document/postprocessing.py` и функции `process_activity_extended`, `postprocess_assays`, `process_target_names`, `process_targets` прямо документируют совместимость с оригинальными Power Query файлами, что гарантирует идентичное формирование `extended.output.*`.


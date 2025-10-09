# Выходные данные

Каждый конвейер записывает детерминированный CSV, файл метаданных (`.meta.yaml`),
профиль качества (`_quality_report_table.csv`) и JSON-итог (`.quality.json`).
Ниже описаны схемы CSV, типы данных и проверки Pandera (`library/schemas/*.py`).

## Метаданные и QC

Для файла `output.targets_20250228.csv` дополнительно создаются:

- `output.targets_20250228.meta.yaml` — название конвейера, параметры запуска,
  итоговая конфигурация, версия схемы и SHA-256 хэш CSV.
- `output.targets_20250228_quality_report_table.csv` — табличный профиль (число
  строк, пропуски, числовая статистика, покрытия регулярных выражений).
- `output.targets_20250228.quality.json` — JSON-сводка с теми же метриками и
  уровнями серьёзности (`info`/`warn`/`error`).

## Таблица документов (`documents`)

Схема: [`library/schemas/documents.py`](../../library/schemas/documents.py).

### Базовые поля ChEMBL

| Колонка | Тип | Обязательная | Описание |
|---------|-----|--------------|----------|
| `document_chembl_id` | string | Да | Основной идентификатор. |
| `title` | string | Нет | Название публикации. |
| `abstract` | string | Нет | Аннотация. |
| `doi` | string | Нет | DOI из ChEMBL. |
| `year` | string/int | Нет | Год публикации (приводится к числу при возможности). |
| `journal` | string | Нет | Название журнала. |
| `journal_abbrev` | string | Нет | Аббревиатура журнала. |
| `volume` | string/int | Нет | Том. |
| `issue` | string/int | Нет | Номер выпуска. |
| `first_page` | string/int | Нет | Первая страница. |
| `last_page` | string/int | Нет | Последняя страница. |
| `pubmed_id` | string/int | Нет | PubMed ID из ChEMBL. |
| `authors` | string | Нет | Список авторов. |
| `source` | string | Нет | Источник записи в ChEMBL. |

### Производные поля

| Колонка | Тип | Описание |
|---------|-----|----------|
| `doi_normalised` | string | Нормализованный DOI (нижний регистр, без пробелов). |
| `publication_types_normalised` | string | Канонизированные типы публикаций. |
| `publication_type_score_review` / `publication_type_score_experimental` / `publication_type_score_unknown` | int | Баллы для определения класса публикации. |
| `publication_class` | string | Итоговая классификация (`review`, `experimental`, `unknown`). |
| `fetch_status` | string | Статус загрузки. |
| `error_source` | string | Источник ошибки (PubMed, CrossRef и т.п.). |
| `date_code` | string | Код периода (год-месяц). |
| `Index` | string/int | Исходный индекс строки. |
| `pipeline_version` | string | Версия пакета. |
| `timestamp_utc` | string | Время выгрузки (UTC, ISO 8601). |

### Поля `PubMed.*`

Все значения строковые (кроме дат, которые приводятся к строкам). Содержат
неструктурированные или списочные данные из E-utilities.

| Колонка | Описание |
|---------|----------|
| `PubMed.PMID` | Идентификатор PubMed. |
| `PubMed.DOI` | DOI из PubMed. |
| `PubMed.ArticleTitle` | Заголовок статьи. |
| `PubMed.Abstract` | Аннотация. |
| `PubMed.JournalTitle` | Название журнала. |
| `PubMed.Volume` / `PubMed.Issue` | Том / номер. |
| `PubMed.StartPage` / `PubMed.EndPage` | Диапазон страниц. |
| `PubMed.PublicationType` | Типы публикаций. |
| `PubMed.MeSH_Descriptors` / `PubMed.MeSH_Qualifiers` | MeSH-термины. |
| `PubMed.ChemicalList` | Список химических веществ. |
| `PubMed.DayRevised` / `MonthRevised` / `YearRevised` | Дата пересмотра. |
| `PubMed.YearCompleted` / `MonthCompleted` / `DayCompleted` | Дата завершения индексации. |
| `PubMed.Error` | Сообщение об ошибке загрузки. |
| `PubMed.ISSN` | ISSN. |
| `PubMed.is_review` | Флаг «обзор» из постобработки. |

### Поля `scholar.*`

| Колонка | Описание |
|---------|----------|
| `scholar.PMID` | Связанный PubMed ID. |
| `scholar.Venue` | Место публикации. |
| `scholar.PublicationTypes` | Типы публикаций. |
| `scholar.SemanticScholarId` | Идентификатор Semantic Scholar. |
| `scholar.ExternalIds` | JSON-подобный набор внешних ID. |
| `scholar.DOI` | DOI от Semantic Scholar. |
| `scholar.Error` | Ошибка вызова API. |
| `scholar.is_review` | Флаг обзора. |

### Поля `OpenAlex.*`

| Колонка | Описание |
|---------|----------|
| `OpenAlex.PublicationTypes` | Типы публикаций OpenAlex. |
| `OpenAlex.TypeCrossref` | Тип из CrossRef. |
| `OpenAlex.Genre` | Жанровая классификация. |
| `OpenAlex.Id` | Идентификатор OpenAlex. |
| `OpenAlex.Venue` | Место публикации. |
| `OpenAlex.MeshDescriptors` / `OpenAlex.MeshQualifiers` | MeSH-термины. |
| `OpenAlex.Error` | Сообщение об ошибке. |
| `OpenAlex.is_review` | Флаг обзора. |

### Поля `crossref.*`

| Колонка | Описание |
|---------|----------|
| `crossref.Type` | Основной тип CrossRef. |
| `crossref.Subtype` | Подтип. |
| `crossref.Title` | Заголовок. |
| `crossref.Subtitle` | Подзаголовок. |
| `crossref.Subject` | Тематики. |
| `crossref.Error` | Ошибка при обращении. |

## Таблица таргетов (`targets`)

Схема: [`library/schemas/targets.py`](../../library/schemas/targets.py). Порядок
колонок соответствует `TARGETS_COLUMN_ORDER`.

### Идентификаторы и наименования

| Колонка | Тип | Описание |
|---------|-----|----------|
| `target_chembl_id` | string | Основной ID ChEMBL. |
| `pref_name` | string | Предпочтительное название. |
| `target_type` | string | Тип таргета. |
| `gene_symbol` | string | Основной символ гена. |
| `gene_symbol_list` | string | Дополнительные символы. |
| `protein_name_canonical` | string | Каноническое название белка. |
| `protein_name_alt` | string | Альтернативные названия. |
| `protein_synonym_list` | string | Список синонимов. |
| `species_group_flag` | string | Флаг группировки по видам. |
| `organism` | string | Организм. |
| `tax_id` | string | Taxonomy ID из ChEMBL. |
| `taxon_id` | string/int | Taxonomy ID из UniProt. |

### Дополнительная классификация (пост-обработка)

Модульная пост-обработка таргетов
(`library.postprocess.targets.run_target_pipeline`) формирует агрегированное
представление для QA-метрик и сервисов, которым нужны компактные
классификационные поля. Хелпер запускается автоматически после записи основного
CSV и проверяет схему из `library/postprocess/targets/schema.py`.

| Колонка | Описание |
|---------|----------|
| `target_class` | Основная классификация, извлечённая из иерархии классов ChEMBL. |
| `protein_family` | Описание белкового семейства верхнего уровня. |
| `synonyms` | Детерминированный список синонимов (предпочтительные имена, описания компонентов, альтернативные названия). |
| `pipeline_version` | Версия, зафиксированная пост-обработчиком (может отличаться при отдельном запуске). |

### Данные UniProt

| Колонка | Описание |
|---------|----------|
| `uniprot_id_primary` | Основной акцессий. |
| `uniprot_ids_all` | Все акцессии. |
| `isoform_ids` / `isoform_names` / `isoform_synonyms` | Сведения по изоформам. |
| `uniprot_last_update` | Дата обновления записи UniProt. |
| `uniprot_version` | Версия последовательности. |
| `secondaryAccessions` | Вторичные акцессии. |
| `recommendedName` / `geneName` / `secondaryAccessionNames` | Называния из записи UniProt. |
| `uniProtkbId` | Идентификатор UniProtKB. |

### Функциональные аннотации

| Колонка | Описание |
|---------|----------|
| `molecular_function`, `cellular_component`, `subcellular_location` | Аннотации GO и текстовые описания. |
| `features_signal_peptide`, `features_transmembrane`, `features_topology` | Фичи из UniProt. |
| `ptm_glycosylation`, `ptm_lipidation`, `ptm_disulfide_bond`, `ptm_modified_residue` | Посттрансляционные модификации. |
| `glycosylation`, `lipidation`, `disulfide_bond`, `modified_residue`, `phosphorylation`, `acetylation`, `ubiquitination` | Развёрнутые списки модификаций. |
| `topology`, `transmembrane`, `intramembrane`, `signal_peptide`, `propeptide` | Структурные характеристики. |

### Кросс-ссылки

| Колонка | Описание |
|---------|----------|
| `xref_chembl`, `xref_uniprot`, `xref_ensembl`, `xref_iuphar` | Кросс-ссылки на источники. |
| `xref_pdb`, `xref_alphafold` | Структурные базы. |
| `GuidetoPHARMACOLOGY`, `gtop_target_id`, `gtop_synonyms`, `gtop_natural_ligands_n`, `gtop_interactions_n`, `gtop_function_text_short` | Данные Guide to PHARMACOLOGY. |
| `cross_references` | Консолидированные кросс-ссылки UniProt. |

### Классификации и домены

| Колонка | Описание |
|---------|----------|
| `pfam`, `Pfam` | Домены Pfam (из разных источников). |
| `interpro`, `InterPro` | Аннотации InterPro. |
| `SUPFAM`, `PROSITE`, `PRINTS`, `TCDB` | Дополнительные классификации. |
| `protein_classifications` | Иерархия ChEMBL. |
| `protein_class_pred_L1` … `protein_class_pred_confidence` | Предсказанные классы, правило, доказательства, уверенность. |

### Реакции и компоненты

| Колонка | Описание |
|---------|----------|
| `target_components` | Состав таргета. |
| `reactions` | Перечень реакций. |
| `reaction_ec_numbers` | EC-номера. |

### IUPHAR

| Колонка | Описание |
|---------|----------|
| `iuphar_target_id` | ID из Guide to PHARMACOLOGY. |
| `iuphar_family_id` | Идентификатор семейства. |
| `iuphar_type`, `iuphar_class`, `iuphar_subclass` | Классификация. |
| `iuphar_chain` | Цепь/субъединица. |
| `iuphar_name` | Полное название. |
| `iuphar_full_id_path`, `iuphar_full_name_path` | Полные пути через `>`.

### Аудит

| Колонка | Описание |
|---------|----------|
| `pipeline_version` | Версия пакета. |
| `timestamp_utc` | Время выгрузки. |

## Таблица ассайев (`assays`)

Схема: [`library/schemas/assays.py`](../../library/schemas/assays.py).

| Колонка | Тип | Описание |
|---------|-----|----------|
| `assay_chembl_id` | string | Основной ID ассая. |
| `accession` | string | UniProt-акцессий целевого белка. |
| `assay_cell_type` | string | Тип клеток. |
| `assay_subcellular_fraction` | string | Субклеточная фракция. |
| `assay_group` | string | Группировка, передаваемая вместе с выгрузкой. |
| `assay_tissue` | string | Ткань/орган. |
| `assay_strain` | string | Штамм (если указан в ChEMBL). |
| `bao_format` | string | Код BAO. |
| `description` | string | Описание. |
| `document_chembl_id` | string | Связанный документ. |
| `isoform` | string/nullable | Изоформа; тип не фиксируется, чтобы принимать строки и числа. |
| `mutation` | string | Детали мутации. |
| `target_chembl_id` | string | ID связанного таргета. |
| `year` | numeric/string | Год публикации (строка или число). |
| `pipeline_version` | string | Версия пакета. |
| `timestamp_utc` | string | Время выгрузки (ISO 8601). |

> **Примечание пост-обработки.** Колонки из списка
> [`_ASSAY_OUTPUT_DROP_COLUMNS`](../../scripts/get_assay_data.py) удаляются
> после выгрузки, чтобы устаревшие поля (например, `ASSAY_ID`, `Target TYPE` или
> `acts_per_assay_step5`) не возвращались в опубликованные CSV.

## Таблица активностей (`activities`)

Схема: [`library/schemas/activities.py`](../../library/schemas/activities.py).

| Колонка | Тип | Проверка |
|---------|-----|----------|
| `activity_id` | string/int | Идентификатор записи. |
| `molecule_chembl_id` | string | ID молекулы. |
| `assay_chembl_id` | string | ID ассая. |
| `activity_comment` | string | Комментарий. |
| `assay_description` | string | Описание ассая. |
| `assay_variant_accession` | string | Акцессий варианта. |
| `assay_variant_mutation` | string | Мутация. |
| `bao_format` / `bao_label` | string | Аннотации BAO. |
| `data_validity_comment` / `data_validity_description` | string | Замечания по валидности. |
| `document_chembl_id` | string | Исходный документ. |
| `pchembl_value` | numeric/string | pChEMBL. |
| `potential_duplicate` | bool/string | Флаг потенциального дубля. |
| `qudt_units` | string | Единицы QUDT. |
| `record_id` | string/int | Идентификатор источника. |
| `relation` | string | Сырой оператор отношения. |
| `src_assay_id` / `src_id` | string/int | Идентификаторы источника. |
| `standard_relation` | string | Нормализованное отношение. |
| `standard_units` | string | Нормализованные единицы. |
| `type` / `units` / `value` | string / string / numeric | Исходные значения измерения. |
| `standard_type` | string | **Проверка:** значение должно входить в список из конфигурации (`metrics`). |
| `standard_value` | float | **Проверка:** ≥ 0, жёсткое приведение к float. |
| `lower_value` / `upper_value` | float | Нижняя/верхняя границы. |
| `activity_properties` | string | Исходный JSON-подобный блок свойств. |
| `action_type` | string | **Проверка:** одно из `PAM`, `NAM`, `activation`, `inhibition`, `binding`, `triaged`, `unknown`. |
| `properties_hash` | string | Хэш структуры свойств. |
| `pipeline_version` | string | Версия пакета. |
| `timestamp_utc` | string | Время выгрузки. |

## Таблица тестовых объектов (`testitems`)

Схема: [`library/schemas/testitems.py`](../../library/schemas/testitems.py).

| Колонка | Тип | Описание |
|---------|-----|----------|
| `molecule_chembl_id` | string | Основной идентификатор. |
| `parent_molecule_chembl_id` | string | Родительская молекула. |
| `salt_chembl_id` | string | ID соли. |
| `natural_product` / `prodrug` / `polymer_flag` | boolean | Флаги ChEMBL. |
| `black_box_warning` | string | Предупреждение FDA. |
| `first_approval` | string | Год первого одобрения. |
| `max_phase` | string | Максимальная стадия разработок. |
| `canonical_smiles` | string | Канонический SMILES. |
| `standard_inchi` | string | Стандартный InChI. |
| `standard_inchi_key` | string | Стандартный InChIKey. |
| `molecule_type` | string | Тип молекулы. |
| `oral` / `parenteral` / `topical` | string/bool | Флаги способов введения. |
| `pref_name` | string | Предпочтительное название. |
| `pubchem_canonical_smiles` | string | SMILES из PubChem. |
| `pubchem_cid` | string/int | CID PubChem. |
| `pubchem_inchi` / `pubchem_inchikey` / `pubchem_isomeric_smiles` / `pubchem_iupac_name` / `pubchem_molecular_formula` | string | Поля, полученные из PubChem. |
| `structure_type` | string | Тип структуры. |
| `pipeline_version` | string | Версия пакета. |
| `timestamp_utc` | string | Время выгрузки. |

## Таблица тканей (`tissues`)

Схема: [`library/schemas/tissues.py`](../../library/schemas/tissues.py). Порядок
колонок задаётся `TISSUES_COLUMN_ORDER`.

### Колонки

| Колонка | Тип | Описание |
|---------|-----|----------|
| `tissue_chembl_id` | string | Основной идентификатор, возвращаемый эндпоинтом `/tissue`. Для отсутствующих записей добавляется строка-заглушка с исходным идентификатором. |
| `pref_name` | string | Предпочтительное название ткани, если присутствует. |
| `uberon_id` | string | Кросс-ссылка UBERON, передаваемая ChEMBL. |
| `efo_id` | string | Кросс-ссылка EFO, передаваемая ChEMBL. |
| `bto_id` | string | Кросс-ссылка BRENDA Tissue Ontology из ChEMBL. |
| `caloha_id` | string | Идентификатор анатомической онтологии Caloha. |
| `pipeline_version` | string | Версия пакета, добавляемая :func:`library.pipelines.common.add_pipeline_metadata`. |
| `timestamp_utc` | string | Временная метка запуска пайплайна в формате ISO 8601 (UTC). |

### Сортировка, пропуски и дополнительные файлы

- Строки сортируются по `tissue_chembl_id` (по возрастанию). При совпадении
  идентификаторов детерминированный порядок обеспечивает CSV-писатель.
- Кросс-ссылки сохраняются в виде nullable-строк pandas и при повторной загрузке
  читаются как `<NA>` (`dtype="string"`).
- Метаданные проставляются для каждой строки, включая заглушки для отсутствующих
  записей.
- Если `--final-out` не указан, CLI создаёт
  `output.tissue_<YYYYMMDD>.csv`, `output.tissue_<YYYYMMDD>.meta.yaml`,
  `output.tissue_<YYYYMMDD>_quality_report_table.csv` и
  `output.tissue_<YYYYMMDD>.quality.json`.

## Экспорт клеточных линий (`cellline`)

Схема описана в [`library/schemas/celllines.py`](../../library/schemas/celllines.py);
порядок колонок соответствует `CELL_LINE_COLUMN_ORDER` и включает служебные
поля `pipeline_version`, `timestamp_utc`.

| Колонка | Тип | Описание |
|---------|-----|----------|
| `cell_chembl_id` | string | Основной идентификатор ChEMBL (обязательный, уникальный). |
| `cell_name` | string | Предпочтительное название клеточной линии. |
| `cell_description` | string | Текстовое описание из ChEMBL. |
| `cell_id` | integer | Внутренний числовой идентификатор ChEMBL (может отсутствовать). |
| `cell_source_organism` | string | Организм-источник клеточной линии. |
| `cell_source_tax_id` | integer | Таксон NCBI для организма-источника (nullable). |
| `cell_source_tissue` | string | Исходная ткань/орган. |
| `cellosaurus_id` | string | Акцессия Cellosaurus при наличии. |
| `cl_lincs_id` | string | Идентификатор LINCS (nullable). |
| `clo_id` | string | Идентификатор в Cell Line Ontology. |
| `efo_id` | string | Ссылка на EFO (nullable). |
| `pipeline_version` | string | Версия пакета на момент экспорта. |
| `timestamp_utc` | string | Временная метка экспорта (UTC). |

Нормализация (`normalize_cell_lines`) приводит идентификаторы к типу `string`,
числовые поля — к nullable-типу pandas `Int64`. Пропуски выгружаются как пустые
строки, что обеспечивает совместимость с CSV-пайплайнами. Строки сортируются по
`cell_chembl_id`, обеспечивая стабильный порядок между запусками.

## Профили качества

CSV-отчёт содержит:

| Колонка | Описание |
|---------|----------|
| `column` | Имя анализируемой колонки. |
| `row_count` | Количество строк. |
| `non_null` | Непустые значения (pandas notnull). |
| `non_empty_ratio` | Доля значений, прошедших `_non_empty_mask`. |
| `distinct_count` | Число уникальных значений. |
| `numeric_min` / `numeric_mean` / `numeric_max` | Статистика по числовым данным. |
| `pattern_email_ratio`, `pattern_doi_ratio`, `pattern_url_ratio`, `pattern_issn_ratio` | Покрытие регулярных выражений. |
| `bool_like_ratio` | Доля значений, похожих на булевы. |

JSON-документ содержит те же метрики с уровнями серьёзности, настроенными через
`system.doc_quality`.

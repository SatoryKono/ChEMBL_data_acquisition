# Выходные таблицы и файлы

Основные скрипты (`scripts/get_*_data.py`) формируют нормализованные CSV и YAML
sidecar-файлы в каталоге `data/output` (путь настраивается через `io.output_dir`).
Ниже перечислены стандартные выгрузки, структура колонок и правила формирования
ключей. Все таблицы сопровождаются `pipeline_version` и `timestamp_utc`, что
позволяет отслеживать источник данных.

## Обзор файлов

| Скрипт | Основной файл | Дополнительно |
|--------|---------------|---------------|
| `get_activity_data.py` | `activities.csv` | `activities.meta.yaml`, `activities_failure_cases.csv` (при ошибках) |
| `get_assay_data.py` | `assays.csv` | `assays.meta.yaml`, контроль качества `assays_quality.csv` |
| `get_document_data.py` | `documents.csv` | `documents.meta.yaml`, отчёты классификации/ошибок |
| `get_target_data.py` | `targets.csv` | `targets.meta.yaml`, при включении дополнительных режимов — `targets_chembl.csv`, `targets_uniprot.csv`, `targets_iuphar.csv` |
| `get_testitem_data.py` | `testitems.csv` | `testitems.meta.yaml` |
| `check_determinism.py` | `determinism_report.csv` | Используется для сравнения эталонов |
| `table_quality_main.py` | `*_quality.csv` | Метрики качества входных таблиц |

Каждый CSV сопровождается YAML-файлом (`*.meta.yaml`), описанным в разделе
[Метаданные](#метаданные-sidecar). При валидационных ошибках дополнительно
создаётся `*_failure_cases.csv` с описанием строк, нарушивших схему Pandera.

## Таблица `activities.csv`

* **Назначение:** витрина активностей молекула-мишень с основными показателями
  ChEMBL и вычисленным значением `standard_value`.
* **Источники:** ChEMBL Activities API (`/activity`), справочники BAO, данные
  документов.
* **Ключи:**
  * `activity_id` — исходный идентификатор активности (может быть `null`).
  * `molecule_chembl_id` + `assay_chembl_id` — составной ключ для дедупликации.

| Колонка | Тип | Источник/логика |
|---------|-----|-----------------|
| `activity_id` | any | Идентификатор активности из ChEMBL. |
| `molecule_chembl_id` | `str` | Ключ молекулы. |
| `assay_chembl_id` | `str` | Ключ ассея, связывает с `assays.csv`. |
| `activity_comment` | `str` | Комментарий ChEMBL. |
| `assay_description` | `str` | Описание ассея. |
| `assay_variant_accession` | `str` | UniProt ID варианта. |
| `assay_variant_mutation` | `str` | Мутация в мишени. |
| `bao_format` / `bao_label` | `str` | Теги BAO, используются при агрегации. |
| `data_validity_comment` / `data_validity_description` | `str` | Проверки качества ChEMBL. |
| `document_chembl_id` | `str` | Связь с `documents.csv`. |
| `pchembl_value`, `relation`, `standard_relation` | `object/str` | Метаданные количественных значений. |
| `standard_value` | `float` | Нормализованное значение, coercion к `float` с проверкой >= 0. |
| `standard_type`, `standard_units` | `str` | Тип и единицы измерения. |
| `type`, `units`, `value` | `object/str` | Исходные измерения. |
| `potential_duplicate` | any | Флаг дубля. |
| `qudt_units` | `str` | Нормализованные единицы из QUDT. |
| `record_id`, `src_assay_id`, `src_id` | any | Исходные ID в источнике. |
| `pipeline_version`, `timestamp_utc` | `str` | Метаданные пайплайна. |

## Таблица `assays.csv`

* **Назначение:** справочник биологических тестов, связанных с активностями.
* **Источники:** ChEMBL Assays API и локальные словари (`dictionary/_Target`).
* **Ключи:** `assay_chembl_id` (основной), `ASSAY_ID` (внутренний идентификатор).

| Колонка | Тип | Описание |
|---------|-----|----------|
| `assay_chembl_id` | `str` | Primary key ChEMBL. |
| `ASSAY_ID` | `str` | Внутренний ID из исходного CSV. |
| `Target TYPE`, `target_name`, `target_chembl_id` | `str` | Тип и идентификаторы мишени. |
| `accession`, `isoform`, `mutation` | `str` | Данные о белке/изоформе. |
| `assay_cell_type`, `assay_tissue`, `assay_subcellular_fraction` | `str` | Контекст эксперимента. |
| `bao_format` | `str` | Типология BAO. |
| `description` | `str` | Текстовое описание. |
| `document_chembl_id` | `str` | Связь с публикацией. |
| `acts_per_assay_step5`, `month`, `year`, `version` | any | Статистические поля (тип гибкий). |
| `cited_assay_corr`, `error_assay_corr`, `higly_correlated_cit`, `shuffled_cit`, `shuffled_target_assay` | any | Флаги качества. |
| `substrate_name` | `str` | Используемый субстрат. |
| `pipeline_version`, `timestamp_utc` | `str` | Метаданные. |

## Таблица `documents.csv`

* **Назначение:** агрегация публикаций из ChEMBL, PubMed, Semantic Scholar,
  OpenAlex, CrossRef.
* **Источники:** ChEMBL Document API, PubMed E-utilities, OpenAlex API,
  CrossRef Works, Semantic Scholar Graph API.
* **Ключи:**
  * `document_chembl_id` — основной идентификатор.
  * `pubmed_id` (`PubMed.PMID`) и `doi`/`doi_normalised` — альтернативные ключи
    для объединения данных из внешних сервисов.

Из-за большого числа колонок таблица разбита на блоки. Типы колонок соответствуют
`schemas/documents.py`.

### Блок ChEMBL

`document_chembl_id`, `title`, `abstract`, `doi`, `year`, `journal`,
`journal_abbrev`, `volume`, `issue`, `first_page`, `last_page`, `pubmed_id`,
`authors`, `source` — строковые/числовые поля с метаданными публикации.

### Derived-поля

`doi_normalised`, `publication_types_normalised`,
`publication_type_score_review`, `publication_type_score_experimental`,
`publication_type_score_unknown`, `publication_class`, `date_code` — используются
для классификации и фильтрации документов.

### PubMed

Колонки с префиксом `PubMed.` содержат значения из E-utilities: `PMID`, `DOI`,
`ArticleTitle`, `Abstract`, `JournalTitle`, сведения о ревизиях (`YearRevised` и
т.д.), `PublicationType`, `MeSH_*`, `ChemicalList`, `ISSN`. Все поля допускают
`null` и приводятся к строке/объекту с `coerce=True`.

### Semantic Scholar

`scholar.PMID`, `scholar.Venue`, `scholar.PublicationTypes`,
`scholar.SemanticScholarId`, `scholar.ExternalIds`, `scholar.DOI`,
`scholar.Error`, `scholar.is_review` — метаданные и статус обработки в API.

### OpenAlex

`OpenAlex.PublicationTypes`, `OpenAlex.TypeCrossref`, `OpenAlex.Genre`,
`OpenAlex.Id`, `OpenAlex.Venue`, `OpenAlex.MeshDescriptors`,
`OpenAlex.MeshQualifiers`, `OpenAlex.Error`, `OpenAlex.is_review` — сведения из
OpenAlex и результирующая классификация.

### CrossRef

`crossref.Type`, `crossref.Subtype`, `crossref.Title`, `crossref.Subtitle`,
`crossref.Subject`, `crossref.Error` — атрибуты записи из CrossRef.

### Метаданные пайплайна

`Index`, `pipeline_version`, `timestamp_utc` — служебная информация для аудита.

## Таблица `targets.csv`

* **Назначение:** объединённый справочник таргетов с полями из ChEMBL, UniProt и
  IUPHAR.
* **Источники:** ChEMBL Target API, UniProt REST/ID Mapping, словари IUPHAR,
  локальные классификаторы (`dictionary/_Target`).
* **Ключи:**
  * `target_chembl_id` — основной ключ.
  * `uniprot_id_primary` и `iuphar_target_id` — альтернативные идентификаторы
    для сопоставления с внешними системами.

Колонки следуют порядку `TARGETS_COLUMN_ORDER`.

### Идентификаторы и базовые поля

`target_chembl_id`, `pref_name`, `target_type`, `organism`, `taxon_id`,
`species_group_flag`, `target_components`, `cross_references`, `gene_symbol` и др.
Содержат основную информацию из ChEMBL.

### UniProt

`uniprot_id_primary`, `uniprot_ids_all`, `isoform_ids`, `isoform_names`,
`isoform_synonyms`, `protein_name_canonical`, `protein_name_alt`, `sequence_length`,
`uniprot_last_update`, `uniprot_version`, множество столбцов по аннотациям белка
(`molecular_function`, `cellular_component`, `subcellular_location`,
`topology`, `transmembrane`, `intramembrane`, `glycosylation`, `lipidation`,
`disulfide_bond`, `modified_residue`, `phosphorylation`, `acetylation`,
`ubiquitination`, `signal_peptide`, `propeptide`). Типы — строки или `object`
(для массивов/списков в сериализованном виде).

### Доменные базы

`pfam`, `interpro`, `SUPFAM`, `PROSITE`, `InterPro`, `Pfam`, `PRINTS`, `TCDB`,
`GuidetoPHARMACOLOGY`, `xref_pdb`, `xref_alphafold` — ссылки на внешние базы.

### IUPHAR

`iuphar_target_id`, `iuphar_family_id`, `iuphar_type`, `iuphar_class`,
`iuphar_subclass`, `iuphar_chain`, `iuphar_name`, `iuphar_full_id_path`,
`iuphar_full_name_path`, `gtop_*` — поля из словарей IUPHAR/Guide to
PHARMACOLOGY.

### Классификации

`protein_classifications`, `protein_class_pred_L1/L2/L3`,
`protein_class_pred_rule_id`, `protein_class_pred_evidence`,
`protein_class_pred_confidence`, а также `gene_symbol_list`,
`protein_synonym_list`, `reactions`, `reaction_ec_numbers` — результаты
постобработки.

### Метаданные

`pipeline_version`, `timestamp_utc` — отметки генерации.

## Таблица `testitems.csv`

* **Назначение:** витрина тест-айтемов (препаратов и молекул) для клинических и
  доклинических исследований.
* **Источники:** ChEMBL Molecule API, PubChem PUG REST.
* **Ключи:** `molecule_chembl_id` (основной), `pubchem_cid` (внешний идентификатор).

| Колонка | Тип | Источник |
|---------|-----|----------|
| `molecule_chembl_id` | `str` | ChEMBL идентификатор молекулы. |
| `black_box_warning`, `oral`, `parenteral`, `topical` | any | Флаги допуска/предупреждений. |
| `first_approval`, `max_phase` | any/`str` | Регуляторные сведения. |
| `canonical_smiles`, `standard_inchi`, `standard_inchi_key`, `molecule_type` | `str` | Структурные данные. |
| `pref_name` | `str` | Название молекулы. |
| `pubchem_*` | `str`/`object` | Поля, полученные из PubChem (SMILES, InChI, IUPAC, формула). |
| `structure_type` | `str` | Тип структуры. |
| `pipeline_version`, `timestamp_utc` | `str` | Метаданные. |

## Метаданные sidecar

Каждый CSV сопровождается YAML-файлом с моделью `CsvMetaSchema`:

| Поле | Тип | Описание |
|------|-----|----------|
| `generated_at` | `datetime` | UTC-время генерации. |
| `git_sha` | `str` | Коммит репозитория. |
| `columns` | `list[str]` | Итоговый порядок колонок в CSV. |
| `dtypes` | `dict[str, str]` | Типы столбцов по отчёту pandas. |
| `command` | `str` | CLI-команда, запустившая выгрузку (если доступно). |
| `config` | `dict[str, Any]` | Снимок конфигурации, применённой при запуске. |

## Failure cases и отчёты качества

* `*_failure_cases.csv` — строки, не прошедшие валидацию Pandera. Содержат номер
  индекса, имя поля, ожидаемый и фактический тип/значение.
* `*_quality.csv` (из `table_quality_main.py`) — агрегированные метрики, включая
  количество пропусков, распределение типов и дубликатов.
* `determinism_report.csv` — сравнительный отчёт `check_determinism.py`,
  используемый для мониторинга изменений в витринах.

## Форматы данных

* Основные файлы — CSV с разделителем `io.csv_sep` и кодировкой `io.csv_encoding`.
* Дополнительно возможно сохранение в Parquet (при активации в скриптах
  постобработки).
* Все идентификаторы и внешние ключи представлены в текстовом виде, чтобы
  сохранялась совместимость между системами хранения.

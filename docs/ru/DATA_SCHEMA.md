# Схема данных

Документ описывает входные и выходные структуры, которые контролируются пайплайнами. Выходные таблицы валидируются Pandera-схемами из `library/schemas`, поэтому описание синхронизировано с кодом.

## Входные таблицы

| Файл | Назначение | Обязательные колонки |
|------|------------|----------------------|
| `activity.csv` | Идентификаторы для `get_activity_data`. | `activity_chembl_id` |
| `assay.csv` | Идентификаторы для `get_assay_data`. | `assay_chembl_id` |
| `documents.csv` | Идентификаторы для `get_document_data` (`chembl`/`all`). | `document_chembl_id` |
| `targets.csv` | Идентификаторы для `get_target_data`. | `target_chembl_id` |
| `testitem.csv` | Идентификаторы для `get_testitem_data`. | `molecule_chembl_id` |

## Выходные таблицы

Ниже приведены таблицы, сформированные напрямую из Pandera-схем, поэтому они отражают текущую реализацию.

### activity.csv (финальный экспорт)
- **Назначение:** нормализованные измерения активности с расчётом границ, аннотациями и служебными полями.
- **Схема:** `library/schemas/activities.py` (`ActivitiesSchema`).

| Колонка | Pandera dtype | Обязательная | Допустимы NULL |
| --- | --- | --- | --- |
| `activity_id` | `None` | Да | Да |
| `molecule_chembl_id` | `str` | Да | Да |
| `assay_chembl_id` | `str` | Да | Да |
| `activity_comment` | `str` | Нет | Да |
| `assay_description` | `str` | Нет | Да |
| `assay_variant_accession` | `str` | Нет | Да |
| `assay_variant_mutation` | `str` | Нет | Да |
| `bao_format` | `str` | Нет | Да |
| `bao_label` | `str` | Нет | Да |
| `data_validity_comment` | `str` | Нет | Да |
| `data_validity_description` | `str` | Нет | Да |
| `document_chembl_id` | `str` | Нет | Да |
| `pchembl_value` | `object` | Нет | Да |
| `potential_duplicate` | `None` | Нет | Да |
| `qudt_units` | `str` | Нет | Да |
| `record_id` | `None` | Нет | Да |
| `relation` | `object` | Нет | Да |
| `src_assay_id` | `None` | Нет | Да |
| `src_id` | `None` | Нет | Да |
| `standard_relation` | `str` | Нет | Да |
| `standard_units` | `str` | Нет | Да |
| `type` | `str` | Нет | Да |
| `units` | `str` | Нет | Да |
| `value` | `object` | Нет | Да |
| `standard_type` | `str` | Нет | Да |
| `standard_value` | `float64` | Да | Да |
| `lower_value` | `float64` | Нет | Да |
| `upper_value` | `float64` | Нет | Да |
| `activity_properties` | `str` | Нет | Да |
| `action_type` | `str` | Нет | Да |
| `properties_hash` | `str` | Нет | Да |
| `pipeline_version` | `str` | Нет | Да |
| `timestamp_utc` | `str` | Нет | Да |

> CLI удаляет промежуточные `standard_lower_value`/`standard_upper_value` после расчёта границ, чтобы зафиксировать заголовок CSV.【F:library/cli/entrypoints/activity.py†L1239-L1360】

### assay.csv (финальный экспорт)
- **Назначение:** консолидация метаданных ассев и служебных полей.
- **Схема:** `library/schemas/assays.py` (`AssaysSchema`).

| Колонка | Pandera dtype | Обязательная | Допустимы NULL |
| --- | --- | --- | --- |
| `assay_chembl_id` | `str` | Да | Да |
| `accession` | `str` | Нет | Да |
| `assay_cell_type` | `str` | Нет | Да |
| `assay_subcellular_fraction` | `str` | Нет | Да |
| `assay_group` | `str` | Нет | Да |
| `assay_tissue` | `str` | Нет | Да |
| `assay_strain` | `str` | Нет | Да |
| `bao_format` | `str` | Нет | Да |
| `description` | `str` | Нет | Да |
| `document_chembl_id` | `str` | Нет | Да |
| `isoform` | `None` | Нет | Да |
| `mutation` | `str` | Нет | Да |
| `target_chembl_id` | `str` | Нет | Да |
| `year` | `None` | Нет | Да |
| `pipeline_version` | `str` | Нет | Да |
| `timestamp_utc` | `str` | Нет | Да |

### target.csv (финальный экспорт)
- **Назначение:** объединённые сведения из ChEMBL, UniProt и IUPHAR с фиксированным порядком колонок.
- **Схема:** `library/schemas/targets.py` (`TargetsSchema`).

| Колонка | Pandera dtype | Обязательная | Допустимы NULL |
| --- | --- | --- | --- |
| `target_chembl_id` | `str` | Да | Да |
| `uniprot_id_primary` | `str` | Нет | Да |
| `uniprot_ids_all` | `str` | Нет | Да |
| `isoform_ids` | `str` | Нет | Да |
| `isoform_names` | `str` | Нет | Да |
| `isoform_synonyms` | `str` | Нет | Да |
| `hgnc_id` | `str` | Нет | Да |
| `gene_symbol` | `str` | Нет | Да |
| `protein_name_canonical` | `str` | Нет | Да |
| `protein_name_alt` | `str` | Нет | Да |
| `organism` | `str` | Нет | Да |
| `taxon_id` | `object` | Нет | Да |
| `lineage_superkingdom` | `str` | Нет | Да |
| `lineage_phylum` | `str` | Нет | Да |
| `lineage_class` | `str` | Нет | Да |
| `sequence_length` | `str` | Нет | Да |
| `features_signal_peptide` | `object` | Нет | Да |
| `features_transmembrane` | `object` | Нет | Да |
| `features_topology` | `str` | Нет | Да |
| `ptm_glycosylation` | `object` | Нет | Да |
| `ptm_lipidation` | `object` | Нет | Да |
| `ptm_disulfide_bond` | `object` | Нет | Да |
| `ptm_modified_residue` | `object` | Нет | Да |
| `xref_chembl` | `str` | Нет | Да |
| `xref_uniprot` | `str` | Нет | Да |
| `xref_ensembl` | `str` | Нет | Да |
| `xref_iuphar` | `str` | Нет | Да |
| `gtop_target_id` | `str` | Нет | Да |
| `gtop_synonyms` | `str` | Нет | Да |
| `gtop_natural_ligands_n` | `str` | Нет | Да |
| `gtop_interactions_n` | `str` | Нет | Да |
| `gtop_function_text_short` | `str` | Нет | Да |
| `uniprot_last_update` | `str` | Нет | Да |
| `uniprot_version` | `str` | Нет | Да |
| `pipeline_version` | `str` | Нет | Да |
| `timestamp_utc` | `str` | Нет | Да |
| `pfam` | `str` | Нет | Да |
| `interpro` | `str` | Нет | Да |
| `xref_pdb` | `str` | Нет | Да |
| `xref_alphafold` | `str` | Нет | Да |
| `hgnc_name` | `str` | Нет | Да |
| `uniProtkbId` | `str` | Нет | Да |
| `secondaryAccessions` | `str` | Нет | Да |
| `recommendedName` | `str` | Нет | Да |
| `geneName` | `str` | Нет | Да |
| `secondaryAccessionNames` | `str` | Нет | Да |
| `molecular_function` | `str` | Нет | Да |
| `cellular_component` | `str` | Нет | Да |
| `subcellular_location` | `str` | Нет | Да |
| `topology` | `str` | Нет | Да |
| `transmembrane` | `object` | Нет | Да |
| `intramembrane` | `object` | Нет | Да |
| `glycosylation` | `object` | Нет | Да |
| `lipidation` | `object` | Нет | Да |
| `disulfide_bond` | `object` | Нет | Да |
| `modified_residue` | `object` | Нет | Да |
| `phosphorylation` | `object` | Нет | Да |
| `acetylation` | `object` | Нет | Да |
| `ubiquitination` | `object` | Нет | Да |
| `signal_peptide` | `object` | Нет | Да |
| `propeptide` | `object` | Нет | Да |
| `GuidetoPHARMACOLOGY` | `str` | Нет | Да |
| `family` | `str` | Нет | Да |
| `SUPFAM` | `str` | Нет | Да |
| `PROSITE` | `str` | Нет | Да |
| `InterPro` | `str` | Нет | Да |
| `Pfam` | `str` | Нет | Да |
| `PRINTS` | `str` | Нет | Да |
| `TCDB` | `str` | Нет | Да |
| `pref_name` | `str` | Нет | Да |
| `AddCellularitySmart ` | `str` | Нет | Да |
| `tax_id` | `str` | Нет | Да |
| `species_group_flag` | `str` | Нет | Да |
| `target_components` | `str` | Нет | Да |
| `protein_classifications` | `str` | Нет | Да |
| `cross_references` | `str` | Нет | Да |
| `gene_symbol_list` | `str` | Нет | Да |
| `protein_synonym_list` | `str` | Нет | Да |
| `reactions` | `str` | Нет | Да |
| `reaction_ec_numbers` | `str` | Нет | Да |
| `protein_class_pred_L1` | `str` | Нет | Да |
| `protein_class_pred_L2` | `str` | Нет | Да |
| `protein_class_pred_L3` | `str` | Нет | Да |
| `protein_class_pred_rule_id` | `str` | Нет | Да |
| `protein_class_pred_evidence` | `str` | Нет | Да |
| `protein_class_pred_confidence` | `str` | Нет | Да |
| `iuphar_target_id` | `str` | Нет | Да |
| `iuphar_family_id` | `str` | Нет | Да |
| `iuphar_type` | `str` | Нет | Да |
| `iuphar_class` | `str` | Нет | Да |
| `iuphar_subclass` | `str` | Нет | Да |
| `iuphar_chain` | `str` | Нет | Да |
| `iuphar_name` | `str` | Нет | Да |
| `iuphar_full_id_path` | `str` | Нет | Да |
| `iuphar_full_name_path` | `str` | Нет | Да |

### testitems.csv (финальный экспорт)
- **Назначение:** расширенный каталог молекул с привязкой к родителям и обогащением PubChem.
- **Схема:** `library/schemas/testitems.py` (`TestitemsSchema`).

| Колонка | Pandera dtype | Обязательная | Допустимы NULL |
| --- | --- | --- | --- |
| `molecule_chembl_id` | `str` | Да | Да |
| `parent_molecule_chembl_id` | `str` | Нет | Да |
| `salt_chembl_id` | `str` | Нет | Да |
| `natural_product` | `boolean` | Нет | Да |
| `prodrug` | `boolean` | Нет | Да |
| `polymer_flag` | `boolean` | Нет | Да |
| `black_box_warning` | `None` | Нет | Да |
| `first_approval` | `None` | Нет | Да |
| `max_phase` | `str` | Нет | Да |
| `canonical_smiles` | `str` | Нет | Да |
| `standard_inchi` | `str` | Нет | Да |
| `standard_inchi_key` | `str` | Нет | Да |
| `molecule_type` | `str` | Нет | Да |
| `oral` | `None` | Нет | Да |
| `parenteral` | `None` | Нет | Да |
| `pref_name` | `str` | Нет | Да |
| `pubchem_canonical_smiles` | `str` | Нет | Да |
| `pubchem_cid` | `object` | Нет | Да |
| `pubchem_inchi` | `str` | Нет | Да |
| `pubchem_inchikey` | `str` | Нет | Да |
| `pubchem_isomeric_smiles` | `str` | Нет | Да |
| `pubchem_iupac_name` | `str` | Нет | Да |
| `pubchem_molecular_formula` | `str` | Нет | Да |
| `structure_type` | `str` | Нет | Да |
| `topical` | `None` | Нет | Да |
| `pipeline_version` | `str` | Нет | Да |
| `timestamp_utc` | `str` | Нет | Да |

### documents.csv (финальный экспорт)
- **Назначение:** единая библиографическая запись на основе ChEMBL, PubMed, Semantic Scholar, OpenAlex и CrossRef.
- **Схема:** декларативный YAML `config/schema/document.yaml`, материализованный через `DOCUMENT_DECLARATION` (`DocumentsSchema`).

`document.yaml` содержит 71 колонку, разбитую по группам источников. CLI использует секцию `export.columns`, поэтому порядок колонок в документации, схеме и выгрузке совпадает.【F:library/schemas/document_spec.py†L13-L118】【F:config/schema/document.yaml†L1-L200】 Для проверки можно выполнить `python -m library.schemas.document_spec` и вывести актуальный перечень.


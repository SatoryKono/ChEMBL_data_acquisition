# Data Schema

This document summarises the canonical input and output structures enforced by the pipelines. Output tables are validated with Pandera schemas located under `library/schemas`, ensuring the documentation mirrors the implementation.

## Input tables

| File | Purpose | Required columns |
|------|---------|------------------|
| `activity.csv` | Seed identifiers for `get_activity_data`. | `activity_chembl_id` |
| `assay.csv` | Seed identifiers for `get_assay_data`. | `assay_chembl_id` |
| `documents.csv` | Seed identifiers for `get_document_data` (`chembl`/`all`). | `document_chembl_id` |
| `targets.csv` | Seed identifiers for `get_target_data`. | `target_chembl_id` |
| `testitem.csv` | Seed identifiers for `get_testitem_data`. | `molecule_chembl_id` |

## Output tables

The following sections are generated directly from the Pandera schemas and therefore reflect the live expectations of the codebase.

### activity.csv (processed export)
- **Purpose:** Normalised activity measurements enriched with derived bounds, annotations and pipeline metadata.
- **Schema source:** `library/schemas/activities.py` (`ActivitiesSchema`).

| Column | Pandera dtype | Required | Nullable |
| --- | --- | --- | --- |
| `activity_id` | `None` | True | True |
| `molecule_chembl_id` | `str` | True | True |
| `assay_chembl_id` | `str` | True | True |
| `activity_comment` | `str` | False | True |
| `assay_description` | `str` | False | True |
| `assay_variant_accession` | `str` | False | True |
| `assay_variant_mutation` | `str` | False | True |
| `bao_format` | `str` | False | True |
| `bao_label` | `str` | False | True |
| `data_validity_comment` | `str` | False | True |
| `data_validity_description` | `str` | False | True |
| `document_chembl_id` | `str` | False | True |
| `pchembl_value` | `object` | False | True |
| `potential_duplicate` | `None` | False | True |
| `qudt_units` | `str` | False | True |
| `record_id` | `None` | False | True |
| `relation` | `object` | False | True |
| `src_assay_id` | `None` | False | True |
| `src_id` | `None` | False | True |
| `standard_relation` | `str` | False | True |
| `standard_units` | `str` | False | True |
| `type` | `str` | False | True |
| `units` | `str` | False | True |
| `value` | `object` | False | True |
| `standard_type` | `str` | False | True |
| `standard_value` | `float64` | True | True |
| `lower_value` | `float64` | False | True |
| `upper_value` | `float64` | False | True |
| `activity_properties` | `str` | False | True |
| `action_type` | `str` | False | True |
| `properties_hash` | `str` | False | True |
| `pipeline_version` | `str` | False | True |
| `timestamp_utc` | `str` | False | True |

> The CLI removes the intermediate `standard_lower_value`/`standard_upper_value` columns after bounds are resolved to keep the output header stable.【F:library/cli/entrypoints/activity.py†L1239-L1360】

### assay.csv (processed export)
- **Purpose:** Consolidated assay metadata with deterministic typing and timestamp bookkeeping.
- **Schema source:** `library/schemas/assays.py` (`AssaysSchema`).

| Column | Pandera dtype | Required | Nullable |
| --- | --- | --- | --- |
| `assay_chembl_id` | `str` | True | True |
| `accession` | `str` | False | True |
| `assay_cell_type` | `str` | False | True |
| `assay_subcellular_fraction` | `str` | False | True |
| `assay_group` | `str` | False | True |
| `assay_tissue` | `str` | False | True |
| `assay_strain` | `str` | False | True |
| `bao_format` | `str` | False | True |
| `description` | `str` | False | True |
| `document_chembl_id` | `str` | False | True |
| `isoform` | `None` | False | True |
| `mutation` | `str` | False | True |
| `target_chembl_id` | `str` | False | True |
| `year` | `None` | False | True |
| `pipeline_version` | `str` | False | True |
| `timestamp_utc` | `str` | False | True |

### target.csv (processed export)
- **Purpose:** Merged target attributes from ChEMBL, UniProt and IUPHAR with legacy column ordering required by downstream Power Query consumers.
- **Schema source:** `library/schemas/targets.py` (`TargetsSchema`).

| Column | Pandera dtype | Required | Nullable |
| --- | --- | --- | --- |
| `target_chembl_id` | `str` | True | True |
| `uniprot_id_primary` | `str` | False | True |
| `uniprot_ids_all` | `str` | False | True |
| `isoform_ids` | `str` | False | True |
| `isoform_names` | `str` | False | True |
| `isoform_synonyms` | `str` | False | True |
| `hgnc_id` | `str` | False | True |
| `gene_symbol` | `str` | False | True |
| `protein_name_canonical` | `str` | False | True |
| `protein_name_alt` | `str` | False | True |
| `organism` | `str` | False | True |
| `taxon_id` | `object` | False | True |
| `lineage_superkingdom` | `str` | False | True |
| `lineage_phylum` | `str` | False | True |
| `lineage_class` | `str` | False | True |
| `sequence_length` | `str` | False | True |
| `features_signal_peptide` | `object` | False | True |
| `features_transmembrane` | `object` | False | True |
| `features_topology` | `str` | False | True |
| `ptm_glycosylation` | `object` | False | True |
| `ptm_lipidation` | `object` | False | True |
| `ptm_disulfide_bond` | `object` | False | True |
| `ptm_modified_residue` | `object` | False | True |
| `xref_chembl` | `str` | False | True |
| `xref_uniprot` | `str` | False | True |
| `xref_ensembl` | `str` | False | True |
| `xref_iuphar` | `str` | False | True |
| `gtop_target_id` | `str` | False | True |
| `gtop_synonyms` | `str` | False | True |
| `gtop_natural_ligands_n` | `str` | False | True |
| `gtop_interactions_n` | `str` | False | True |
| `gtop_function_text_short` | `str` | False | True |
| `uniprot_last_update` | `str` | False | True |
| `uniprot_version` | `str` | False | True |
| `pipeline_version` | `str` | False | True |
| `timestamp_utc` | `str` | False | True |
| `pfam` | `str` | False | True |
| `interpro` | `str` | False | True |
| `xref_pdb` | `str` | False | True |
| `xref_alphafold` | `str` | False | True |
| `hgnc_name` | `str` | False | True |
| `uniProtkbId` | `str` | False | True |
| `secondaryAccessions` | `str` | False | True |
| `recommendedName` | `str` | False | True |
| `geneName` | `str` | False | True |
| `secondaryAccessionNames` | `str` | False | True |
| `molecular_function` | `str` | False | True |
| `cellular_component` | `str` | False | True |
| `subcellular_location` | `str` | False | True |
| `topology` | `str` | False | True |
| `transmembrane` | `object` | False | True |
| `intramembrane` | `object` | False | True |
| `glycosylation` | `object` | False | True |
| `lipidation` | `object` | False | True |
| `disulfide_bond` | `object` | False | True |
| `modified_residue` | `object` | False | True |
| `phosphorylation` | `object` | False | True |
| `acetylation` | `object` | False | True |
| `ubiquitination` | `object` | False | True |
| `signal_peptide` | `object` | False | True |
| `propeptide` | `object` | False | True |
| `GuidetoPHARMACOLOGY` | `str` | False | True |
| `family` | `str` | False | True |
| `SUPFAM` | `str` | False | True |
| `PROSITE` | `str` | False | True |
| `InterPro` | `str` | False | True |
| `Pfam` | `str` | False | True |
| `PRINTS` | `str` | False | True |
| `TCDB` | `str` | False | True |
| `pref_name` | `str` | False | True |
| `AddCellularitySmart ` | `str` | False | True |
| `tax_id` | `str` | False | True |
| `species_group_flag` | `str` | False | True |
| `target_components` | `str` | False | True |
| `protein_classifications` | `str` | False | True |
| `cross_references` | `str` | False | True |
| `gene_symbol_list` | `str` | False | True |
| `protein_synonym_list` | `str` | False | True |
| `reactions` | `str` | False | True |
| `reaction_ec_numbers` | `str` | False | True |
| `protein_class_pred_L1` | `str` | False | True |
| `protein_class_pred_L2` | `str` | False | True |
| `protein_class_pred_L3` | `str` | False | True |
| `protein_class_pred_rule_id` | `str` | False | True |
| `protein_class_pred_evidence` | `str` | False | True |
| `protein_class_pred_confidence` | `str` | False | True |
| `iuphar_target_id` | `str` | False | True |
| `iuphar_family_id` | `str` | False | True |
| `iuphar_type` | `str` | False | True |
| `iuphar_class` | `str` | False | True |
| `iuphar_subclass` | `str` | False | True |
| `iuphar_chain` | `str` | False | True |
| `iuphar_name` | `str` | False | True |
| `iuphar_full_id_path` | `str` | False | True |
| `iuphar_full_name_path` | `str` | False | True |

### testitems.csv (processed export)
- **Purpose:** Enriched molecule catalogue combining ChEMBL metadata, parent hierarchies and PubChem descriptors.
- **Schema source:** `library/schemas/testitems.py` (`TestitemsSchema`).

| Column | Pandera dtype | Required | Nullable |
| --- | --- | --- | --- |
| `molecule_chembl_id` | `str` | True | True |
| `parent_molecule_chembl_id` | `str` | False | True |
| `salt_chembl_id` | `str` | False | True |
| `natural_product` | `boolean` | False | True |
| `prodrug` | `boolean` | False | True |
| `polymer_flag` | `boolean` | False | True |
| `black_box_warning` | `None` | False | True |
| `first_approval` | `None` | False | True |
| `max_phase` | `str` | False | True |
| `canonical_smiles` | `str` | False | True |
| `standard_inchi` | `str` | False | True |
| `standard_inchi_key` | `str` | False | True |
| `molecule_type` | `str` | False | True |
| `oral` | `None` | False | True |
| `parenteral` | `None` | False | True |
| `pref_name` | `str` | False | True |
| `pubchem_canonical_smiles` | `str` | False | True |
| `pubchem_cid` | `object` | False | True |
| `pubchem_inchi` | `str` | False | True |
| `pubchem_inchikey` | `str` | False | True |
| `pubchem_isomeric_smiles` | `str` | False | True |
| `pubchem_iupac_name` | `str` | False | True |
| `pubchem_molecular_formula` | `str` | False | True |
| `structure_type` | `str` | False | True |
| `topical` | `None` | False | True |
| `pipeline_version` | `str` | False | True |
| `timestamp_utc` | `str` | False | True |

### documents.csv (processed export)
- **Purpose:** Consolidated bibliographic metadata across ChEMBL, PubMed, Semantic Scholar, OpenAlex and CrossRef.
- **Schema source:** Declarative YAML at `config/schema/document.yaml`, materialised through `DOCUMENT_DECLARATION` and exposed as `DocumentsSchema`.

`document.yaml` currently lists 71 ordered columns grouped by provider. The CLI exporters honour the `export.columns` projection defined in the same file, ensuring that documentation, schemas and runtime outputs remain aligned.【F:library/schemas/document_spec.py†L13-L118】【F:config/schema/document.yaml†L1-L200】 Run `python -m library.schemas.document_spec` to dump the effective column order if required for audits.


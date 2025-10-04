# Schema crosswalks

This appendix maps input columns to output columns and highlights key joins.

## Document pipeline

| Input column | Output column(s) | Notes |
|--------------|------------------|-------|
| `document_chembl_id` | `document_chembl_id` | Primary key propagated through all stages. |
| `PMID` | `pubmed_id`, `PubMed.PMID`, `scholar.PMID` | Used when `--mode` includes `pubmed`. |
| `DOI` (fallback file) | `doi_normalised`, `PubMed.DOI`, `crossref.*` | Fallback applied when external sources disagree. |

## Target pipeline

| Input column | Output column(s) | Notes |
|--------------|------------------|-------|
| `target_chembl_id` | `target_chembl_id` | Required for all modes. |
| `mapping_uniprot_id` | `uniprot_id_primary` | Default column for UniProt processing in `--mode all`. |
| `uniprot_id` | `uniprot_id_primary` (when `--column uniprot_id`) | Direct UniProt processing. |

### Joins

1. **ChEMBL stage** – groups by `target_chembl_id`, aggregates cross references
   and basic metadata.
2. **UniProt stage** – merges on `uniprot_id_primary` to enrich with sequences,
   GO terms and topology data.
3. **IUPHAR stage** – matches UniProt IDs to IUPHAR chains and families.

## Assay pipeline

| Input column | Output column | Notes |
|--------------|---------------|-------|
| `assay_chembl_id` | `assay_chembl_id` | Primary key. |
| `Target ID` | `target_chembl_id` | Mapped using ChEMBL metadata. |

## Activity pipeline

| Input column | Output column(s) | Notes |
|--------------|------------------|-------|
| `activity_chembl_id` | `activity_id` | Identifier is preserved. |
| `standard_type`, `standard_relation`, `standard_value`, `standard_units` | Same columns | Validated and normalised. |
| `action_type` | `action_type` | Recomputed when `activity_enrichment.action_type.enabled` is `true`. |
| `activity_properties` | `activity_properties`, `properties_hash` | Flattened summary + deterministic hash. |

## Test item pipeline

| Input column | Output column(s) | Notes |
|--------------|------------------|-------|
| `molecule_chembl_id` | `molecule_chembl_id` | Primary key. |
| `pref_name` | `pref_name` | Preserved; optional enrichment from PubChem. |
| `parent_molecule_chembl_id` (if provided) | `parent_molecule_chembl_id` | Overrides dictionary-derived parent. |

## Cross-pipeline relationships

- Activities link to assays via `assay_chembl_id`, targets via
  `assay_variant_accession`/`assay_variant_mutation` metadata, documents via
  `document_chembl_id`, and test items via `molecule_chembl_id`.
- The star schema is described in [`architecture/DATA_MODEL.md`](../architecture/DATA_MODEL.md).

When adding new columns ensure the crosswalk is updated to reflect how data
flows from inputs through to exports.

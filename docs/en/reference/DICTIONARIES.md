# Dictionary datasets

The pipelines rely on static resources bundled under `config/dictionary`. This
reference summarises each dataset, the consuming pipeline and important fields.

| Path | Used by | Purpose |
|------|---------|---------|
| `_target/_IUPHAR/_IUPHAR_target.csv` | Target (`iuphar`, `all`) | Mapping from UniProt accession to Guide to PHARMACOLOGY target metadata. |
| `_target/_IUPHAR/_IUPHAR_family.csv` | Target (`iuphar`, `all`) | Hierarchy of IUPHAR families (used to build `iuphar_full_*` columns). |
| `_target/_uniprot/*.json` | Target (`uniprot`, `all`) | Cached UniProt KB responses keyed by accession. Refresh with the UniProt downloader tools when schema changes. |
| `_target/targets_type.csv` | Target QA | Maps ChEMBL target types to canonical groupings for reporting. |
| `_testitem/molecule_catalog.csv` | Test item | Parent/child molecule relationships including `parent_molecule_chembl_id`. |
| `_testitem/molecule_hierarchy.csv` | Test item | Precomputed fallback hierarchy used when the API omits parent data. |
| `_document/fallback_doi_template.csv` | Document (optional) | Template structure for DOI override files consumed via `--fallback-doi-*`. |
| `_taxonomy/taxonomy.csv` | Assay | Maps `tax_id` values to canonical `assay_group` and `assay_strain` labels used during assay enrichment. |

### Maintaining dictionary data

- Update CSVs atomically and keep commits small; the QA suite compares hashes in
  smoke tests.
- When adding new columns record the meaning in this document and ensure the
  consuming pipeline handles nulls gracefully.
- For cached JSON responses store only the minimal subset required for tests to
  avoid large commits. The production pipeline should re-fetch data on demand.

### Parent molecule enrichment

The test item pipeline uses the dictionary files as the source of truth for
parent-child relationships and fallback lookups. Mandatory columns:

| File | Columns |
|------|---------|
| `_testitem/molecule_catalog.csv` | `molecule_chembl_id`, `parent_molecule_chembl_id`, optional PubChem fields |
| `_testitem/molecule_hierarchy.csv` | `molecule_chembl_id`, `parent_molecule_chembl_id`, `pref_name`, `level` |

Ensure new datasets include `parent_molecule_chembl_id`; missing values are
logged and, if `parent_fallback` is enabled, replaced using the hierarchy file.

### Target classification helpers

Post-processing does not introduce new dictionary files, but the modular
`library.postprocessing.targets` pipeline relies on existing columns to derive a
compact view for QA metrics. The helper reads `protein_classifications` from the
canonical export to populate the `target_class` and `protein_family` columns and
combines `pref_name`, component descriptions and alternative names into the
deterministic `synonyms` list.

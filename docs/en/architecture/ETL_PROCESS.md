# ETL process

The ETL flow is identical across languages; this document details the execution
within each pipeline stage.

## 1. Documents

1. Load the seed CSV (`--input`) and validate required columns (`document_chembl_id`,
   mode-specific identifier column).
2. Fetch ChEMBL metadata (`mode=chembl`): call the ChEMBL API in batches, merge
   results, normalise DOIs and authors, record fetch status.
3. Fetch PubMed / partner data (`mode=pubmed`): batch PMIDs, call PubMed,
   Semantic Scholar, OpenAlex, CrossRef respecting rate limits, merge JSON payloads.
4. Merge stages (`mode=all`): join ChEMBL and PubMed payloads, apply fallback DOI
   overrides, deduplicate by DOI/PMID.
5. Post-process: derive publication class, compute review scores, write CSV,
   metadata and quality artefacts.

## 2. Targets

1. Read identifier CSV (`target_chembl_id` and optional UniProt columns).
2. Stage **ChEMBL**: resolve target records, collect cross references, taxonomies
   and classification data.
3. Stage **UniProt**: for each primary accession load cached JSON or request the
   REST API, extract sequence features, GO annotations and PTMs.
4. Stage **IUPHAR**: map UniProt IDs to IUPHAR targets and families using local
   dictionaries.
5. Merge intermediate frames, normalise columns using `TARGETS_COLUMN_ORDER`,
   validate with Pandera and export.

## 3. Assays

1. Validate input columns (`assay_chembl_id`, `target_chembl_id`).
2. Fetch assay metadata from ChEMBL, including BAO annotations, correlation flags
   and document links.
3. Apply flexible type coercion for numeric/boolean fields, enrich with pipeline
   metadata and write outputs.

## 4. Test items

1. Validate input molecules (`molecule_chembl_id`).
2. Fetch ChEMBL molecule details (structure, administration flags, warnings).
3. Augment with PubChem using the configured resolution order (cache → SMILES →
   InChIKey → InChI → preferred name).
4. Reconcile parent molecules using dictionaries; emit warnings for missing
   parents and apply fallback hierarchy when enabled.
5. Export CSV with boolean coercion and metadata sidecars.

## 5. Activities

1. Validate identifiers (`activity_chembl_id`, `assay_chembl_id`,
   `molecule_chembl_id`).
2. Fetch ChEMBL activities with pagination and optional multithreading (`--workers`).
3. Normalise measurement columns (`standard_*`), derive bounds from relations,
   calculate action types via the enrichment rules, flatten activity properties.
4. Validate non-negative `standard_value`, allowed `action_type`, and write CSV,
   metadata and QC reports.

## Shared post-processing

- Deterministic row ordering using identifier columns and stable sorting.
- Column ordering via schema-specific constants (`DOCUMENT_SCHEMA_COLUMNS`,
  `TARGETS_COLUMN_ORDER`, etc.).
- Metadata sidecars capturing configuration origin, dependency versions and
  hashes.
- Table quality profiling controlled by `system.doc_quality`.

Refer to [`QC_DEPENDENCIES.md`](./QC_DEPENDENCIES.md) for QA dependencies and
[`DATA_MODEL.md`](./DATA_MODEL.md) for downstream analytical structure.

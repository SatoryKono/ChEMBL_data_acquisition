# Analytical Report on the ETL Workflow

## ETL Overview

The pipeline is orchestrated through the shared configuration file `config/config.yaml`, which defines base URLs for external APIs, retry parameters, RPS limits, logging settings, and directories for inputs and outputs. Command scripts in `scripts/` are organized by entity and rely on shared utilities for I/O, normalization, validation, and logging from the `library/` package, while auxiliary helpers live under `library.utils.cli_tools`, ensuring repeatable execution for documents, assays, activities, test items, and targets.

## Data Sources

| Source | Integration | Retrieved data |
| --- | --- | --- |
| **ChEMBL API** | The `ChemblClient` HTTP client manages headers and retries, while entity-specific functions `get_documents`, `get_assays`, `get_activities`, `get_testitem`, `get_targets` accept chunk-size/timeout parameters from the configuration and batch API calls. | Publication, assay, activity, test-item, and target metadata with key identifiers (`document_chembl_id`, `assay_chembl_id`, `target_chembl_id`, `molecule_chembl_id`, `pubmed_id`). |
| **PubMed / Semantic Scholar / OpenAlex / CrossRef** | The `fetch_pubmed_records` module and companion clients use delay, limit, and parallelism settings from the `document.pubmed` and `document.all` configuration sections, merging responses into a single dataframe. | PubMed XML metadata, publication types, DOIs, identifiers, and cross-links from Semantic Scholar, OpenAlex, and CrossRef. |
| **UniProt REST API** | `uniprot_library.process` retrieves JSON entries while observing RPS thresholds and reusing prepared accession CSVs; the `target all` mode wraps calls and caches results. | Full UniProt cards: names, taxonomy, sequences, secondary accessions, and related cross-references. |
| **IUPHAR and local CSV files** | `fetch_iuphar` reads local classifiers and, when needed, invokes REST services; `target_postprocessing` merges the reference data with ChEMBL and UniProt extracts. | IUPHAR classifications, families, HGNC identifiers, and predicted protein classes. |
| **PubChem PUG REST** | The `add_pubchem_data` function initializes a PubChem session, then issues rate-limited requests for unique SMILES via the `pubchem` configuration. | CID, IUPAC names, formulas, InChI, and SMILES used to enrich test items. |
| **Local dictionaries and prepared CSV/Excel files** | `io.read_ids` and entity modules load files from `dictionary/` and `data/` as input identifier lists or reference tables, reducing external calls. | Restrictive ID lists, target-type dictionaries, and initialization Excel workbooks. |

## Data Transformation

### Shared normalization and quality checks

All entities are processed through the unified `library.constants.normalization` layer, which replaces non-standard characters such as “μ” with “u”, aligns comparison operators (`<` → `<=`, `>` → `>=`), and trims identifiers across the dataframe, eliminating discrepancies before validation. After normalization each script adds the technical columns `pipeline_version` and `timestamp_utc` via `add_pipeline_metadata`, sourcing the version from `pyproject.toml` or the installed package to preserve export lineage.

Validation relies on `pandera`: when mismatches occur, the `SidecarErrors` helper aggregates row-level violations and persists them into a separate CSV with metadata, keeping the primary export intact. Export duties are handled by `io.write_csv`, which enforces deterministic ordering of rows and columns, creates directories, and generates a YAML sidecar containing the launch command, configuration snapshot, and file hash. After writing, each table is assessed by `analyze_table_quality`, delivering numeric profiles and type distributions.

### Document pipeline

1. **Collection and merge** — The `all` mode first extracts ChEMBL document metadata, normalizes DOI locally, then fetches PubMed identifiers where necessary and calls PubMed, Semantic Scholar, OpenAlex, and CrossRef, caching the DOI map for reuse. The `merge_with_chembl` function aligns `pubmed_id` types and joins the external metadata with the ChEMBL export while dropping duplicate fields.
2. **Source normalization** — `merge_metadata` harmonizes DOIs across sources, gathers unique publication types, computes weighted scores (review/experimental/unknown), and composes the final classification so terminology is unified before export.
3. **Post-processing** — `postprocess_documents` enforces consistent types on text and numeric columns, derives a `date_code`, sorts rows by date, adds a sequential index, extracts review flags from text, and establishes the final column order that mirrors the Power Query workflow.
4. **Export** — `_finalise_export` appends pipeline metadata, reshapes the dataframe to match the schema, casts non-numeric fields to strings, maps key columns to export names, and writes CSV alongside YAML metadata and a JSON quality report covering DOI coverage, publication classes, and source errors.

### Assay pipeline

1. **Extraction and aggregation** — The `get_assay_data` script reads identifiers from CSV, requests ChEMBL data, then runs `postprocess_assays` before validation to compute the `assay_with_same_target` counter per document-target pair for downstream analysis.
2. **Normalization and metadata** — After column alignment the dataframe is normalized (`normalize_assays`), enriched with technical fields, and validated with the `AssaysSchema`; any mismatches are routed to a sidecar file while clean rows proceed.
3. **Export** — Columns are ordered with schema fields first and additional ones alphabetically; CSV output includes metadata and a quality report, ensuring reproducible and self-documenting exports.

### Activity pipeline

1. **Data retrieval** — `get_activity_data` reads activity IDs, requests ChEMBL rows in batches, and applies `normalize_activities`, aligning operators and identifiers prior to validation.
2. **Validation and errors** — `validate_activities` with `return_result=True` returns both the cleaned dataframe and potential issues. Failed rows go into the sidecar with context, while logs capture warnings about missing optional columns.
3. **Export** — Columns are aligned to the schema, CSV is written with key sort fields, and YAML sidecar plus quality report are generated. The same template is used across entities, delivering a unified artifact format.

### Test-item pipeline

1. **ChEMBL + PubChem aggregation** — After loading ChEMBL data, `add_pubchem_data` iterates unique canonical SMILES, fetches CID and properties, then concatenates results; unmatched cases keep empty rows to preserve table structure.
2. **Normalization and validation** — `normalize_testitems` standardizes string fields, `add_pipeline_metadata` injects technical columns, and `validate_testitems` in lazy mode captures errors and persists them via `SidecarErrors` without halting the main flow.
3. **Export** — As with other pipelines, columns follow the schema-then-alphabet rule; CSV and YAML are written, SHA-256 is computed, and a quality report checks determinism and logs the outcome.

### Target pipeline

1. **Acquisition modes** — `get_target_data` supports standalone scenarios (`chembl`, `uniprot`, `iuphar`, `all`). In `chembl` mode the data is normalized (`normalize_targets`), augmented with metadata, and validated with a sidecar for issues; exports include YAML and quality analysis.
2. **Composite `all` mode** — Sequential functions `fetch_chembl`, `fetch_uniprot`, and `fetch_iuphar` gather data, join UniProt/IUPHAR with ChEMBL, then `merge_results` injects protein class predictions and passes the table to post- and final processing.
3. **Target post-processing** — `postprocess_targets` normalizes UniProt identifiers, uppercases gene synonyms, consolidates synonym and EC lists, fills optional fields with defaults, and restructures the column set using `TARGETS_COLUMN_ORDER` while preserving the origin of `target_chembl_id` and accession lists.
4. **Finalization** — `finalise_targets` removes duplicates and entries lacking UniProt IDs, invokes the taxonomy classifier to derive `type` and related flags from UniProt lineage plus ChEMBL species data, lowercases selected fields, and `validate_and_write` performs normalization, adds metadata, validates against `TargetsSchema`, replaces gaps with “-”, writes CSV with a fixed order, and emits a quality report.

## Data Export

Each script produces a bundle of artifacts: a primary CSV with deterministic row/column ordering, a YAML sidecar (`*.meta.yaml`) capturing configuration, launch command, row counts, and checksum, plus auxiliary reports—CSV with validation errors and JSON/text quality summaries (e.g., DOI coverage and publication class distribution). Target runs additionally create intermediate files (`*_chembl.csv`, `*_uniprot.csv`, `*_iuphar.csv`) in `all` mode to aid diagnostics for individual steps.

## Entity Relationships

* Assays reference documents (`document_chembl_id`) and targets (`target_chembl_id`), and receive a parallel-target counter during post-processing to highlight well-studied combinations.
* Activities bridge assays, molecules, and documents via the respective identifiers and, thanks to operator normalization, align consistently with downstream analytics.
* Test items connect ChEMBL and PubChem through SMILES and CID, adding chemical properties for future joins with external compound libraries.
* Targets aggregate UniProt, IUPHAR, and the built-in taxonomy classifier into a single table enriched with protein-class predictions and synonym lists, ready to join with activities and assays via `target_chembl_id`/`uniprot_id`.
* Documents include normalized DOIs, publication types, and review flags, enabling alignment of assay/activity results with source quality.

## Project Structure

* **`config/config.yaml`** — Central configuration for APIs, limits, paths, reference data, and parameters of each subsystem.
* **`scripts/`** — CLI wrappers for entity loading, quality reporting, and dictionary maintenance; every command covers input reading, client calls, normalization, validation, and export.
* **`library/`** — Core business logic: API clients, post-processing (documents, targets, assays), normalization, validation, logging, CSV operations, and sidecar handling.
* **`library/constants/`** — `pandera` schemas and normalization routines for every entity.
* **`dictionary/` and `data/`** — Local dictionaries, UniProt/IUPHAR caches, and input CSV/Excel files for launching pipelines.
* **`docs/`** — Documentation for configuration, execution, and outputs; this report extends it with an end-to-end ETL description.

The report covers the entire cycle—from data sources through normalization, post-processing, and export—highlighting quality controls and entity relationships. It equips newcomers to quickly understand, extend, or troubleshoot existing pipelines without compromising data integrity.

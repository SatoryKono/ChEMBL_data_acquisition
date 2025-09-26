# Configuration Reference

## Overview

* All command-line tools load their defaults from [`config.yaml`](../config.yaml) in the project root.
* Values are validated against [`config.schema.json`](../config.schema.json); unknown keys raise an error during start-up.
* Settings can be overridden via environment variables and CLI flags. Precedence is: `config.yaml` < environment variables < CLI arguments.

## Layout of `config.yaml`

| Section | Purpose |
| --- | --- |
| `sources` | Connectivity details for external services such as ChEMBL, UniProt, CrossRef, PubChem and PubMed. |
| `local` | Paths to local resources, CSV defaults and initialisation workbooks. |
| `system` | Logging, retry policy, global rate limiting and document classification weights. |

Sensitive values (API tokens, personal e-mails) should be injected via environment variables rather than committed to the repository.

## `sources.chembl`

### API client (`sources.chembl.api`)

| Key | Default | Description |
| --- | --- | --- |
| `chembl_base` | `https://www.ebi.ac.uk/chembl/api/data` | Base URL for the ChEMBL REST API. |
| `timeout_connect` | `5` | Connection timeout in seconds. |
| `timeout_read` | `30` | Read timeout in seconds. |
| `retries` | `3` | HTTP retry attempts handled by `requests`. |
| `backoff_factor` | `0.5` | Multiplier for exponential backoff between retries. |
| `rps` | `20` | Allowed requests per second for the rate limiter. |
| `burst` | `20` | Bucket size used by the token bucket limiter. |
| `user_agent` | `chembl-da/0.1 (mailto:contact@example.org)` | User-Agent header; replace the contact e-mail in production. |

### Response cache (`sources.chembl.cache`)

| Key | Default | Description |
| --- | --- | --- |
| `cache_ttl` | `3600` | Time-to-live for cached API responses in seconds. |
| `cache_maxsize` | `1024` | Maximum number of cached responses. |

<a id="sources-chembl-molecule-catalog"></a>
### Molecule catalogue (`sources.chembl.molecule_catalog`)

| Key | Default | Description |
| --- | --- | --- |
| `cache_path` | `data/cache/molecule_parent_catalog.json` | Location of the JSON cache storing molecule parent-child relationships reused by enrichment jobs. |
| `endpoint` | `molecule` | ChEMBL REST resource queried when the cache needs to be refreshed. |
| `child_field` | `molecule_chembl_id` | JSON field containing the child molecule identifier extracted from API responses. |
| `parent_field` | `parent_molecule_chembl_id` | JSON field containing the parent molecule identifier extracted from API responses. |
| `page_size` | `500` | Number of records requested per API page while rebuilding the catalogue. |

### Pipelines (`sources.chembl.pipelines`)

Each sub-section below defines defaults for the respective CLI utility. CLI arguments are merged back into the configuration before execution.

#### Activity pipeline (`activity`)

| Key | Default | Description |
| --- | --- | --- |
| `column` | `activity_chembl_id` | Input column containing activity identifiers. |
| `batch_size` | `50` | Batch size for API requests. |
| `timeout` | `30.0` | Request timeout in seconds. |
| `limit` | `null` | Optional cap on identifiers processed. |
| `dry_run` | `false` | Skip network calls and file generation when `true`. |

#### Activity bounds (`activity_bounds`)

| Key | Default | Description |
| --- | --- | --- |
| `enable_from_relation` | `true` | Derive bounds from `standard_value` and relation operators when explicit limits are absent. |
| `enable_from_uncertainty` | `false` | Parse `standard_text_value` expressions like `value ± delta` to generate bounds. |
| `rounding_digits` | `3` | Decimal digits used when rounding the derived limits. |
| `clamp_nonnegative` | `true` | Clamp negative bounds to zero for concentration-like metrics. |
| `log_unknown_relations` | `true` | Emit warnings when relation markers are not recognised by the pipeline. |

#### Activity enrichment (`activity_enrichment`)

The activity pipeline enriches raw ChEMBL payloads with canonical lower/upper bounds using the rules implemented by
`compute_activity_bounds` in `scripts/get_activity_data.py`. Configuration for the feature is stored in the `activity_bounds`
block and controls the following deterministic stages (executed in order for every row):【F:scripts/get_activity_data.py†L212-L353】【F:library/config.py†L371-L388】

1. Use `standard_lower_value`/`standard_upper_value` when both are populated.
2. Combine `standard_value` with the opposite explicit limit (for example `standard_upper_value`) and fill the missing bound.
3. Inspect `standard_relation` when `enable_from_relation` is `true`, mapping operators such as `=`, `≈`, `>=`, `<=`, `between` and `range` to appropriate bounds.
4. Parse `±` expressions from `standard_text_value` when `enable_from_uncertainty` is enabled.

Each completed step locks in previously derived numbers; disabling a stage simply skips it without mutating earlier results.

Toggles and precision defaults are configured via YAML or environment overrides:

```yaml
activity_bounds:
  enable_from_relation: false
  enable_from_uncertainty: true
  rounding_digits: 2
  clamp_nonnegative: true
```

```bash
export CHEMBL_DA__ACTIVITY_BOUNDS__ENABLE_FROM_RELATION=false
export CHEMBL_DA__ACTIVITY_BOUNDS__ROUNDING_DIGITS=2
```

The CLI only exposes high-level switches such as `--batch-size` or `--dry-run`; enrichment-specific options must be changed in the configuration file or via the corresponding `CHEMBL_DA__ACTIVITY_BOUNDS__*` variables. CLI values still win over file/env defaults for overlapping keys declared on the parser (column, batch size, limits).【F:scripts/get_activity_data.py†L536-L603】

#### Assay pipeline (`assay`)

| Key | Default | Description |
| --- | --- | --- |
| `column` | `assay_chembl_id` | Input column with assay identifiers. |
| `batch_size` | `50` | Batch size for API requests. |
| `timeout` | `30.0` | Request timeout in seconds. |
| `limit` | `null` | Optional cap on identifiers processed. |

#### Test item pipeline (`testitem`)

| Key | Default | Description |
| --- | --- | --- |
| `column` | `molecule_chembl_id` | Input column with compound identifiers. |
| `batch_size` | `50` | Batch size for API requests. |
| `timeout` | `30.0` | Request timeout in seconds. |
| `limit` | `null` | Optional cap on identifiers processed. |

#### Test item molecule enrichment (`testitem_molecule_enrichment`)

| Key | Default | Description |
| --- | --- | --- |
| `enable` | `true` | Master switch for the enrichment stage that derives salt identifiers and catalogue flags. |
| `sources.molecule_catalog_path` | `dictionary/molecule_catalog.csv` | CSV with `molecule_chembl_id` and the `natural_product`/`prodrug`/`polymer_flag` columns. |
| `sources.molecule_hierarchy_path` | `dictionary/molecule_hierarchy.csv` | CSV that maps derivatives to their parent molecule (`molecule_chembl_id`, `parent_molecule_chembl_id`). |
| `output.salt_as_null_when_absent` | `true` | Emit `null` (or `-` when set to `false`) when the compound is not a salt. |
| `flags.coerce_to_bool` | `true` | Normalise catalogue values such as `Y/N`, `1/0`, `yes/no` to pandas nullable booleans. |
| `flags.parent_fallback` | `true` | Reuse parent flag values when the child entry is missing. |
| `logging.warn_missing_molecule` | `true` | Log warnings when a molecule listed in the export is not present in the hierarchy or catalogue. |
| `logging.warn_inconsistent_flags` | `true` | Emit warnings when child and parent flags disagree before fallback. |

#### Document pipeline (`document`)

| Sub-section | Key | Default | Description |
| --- | --- | --- | --- |
| `pubmed` | `column` | `PMID` | Column with PubMed identifiers. |
|  | `sleep` | `5.0` | Delay between polling cycles in seconds. |
|  | `workers` | `1` | Number of worker threads. |
|  | `batch_size` | `5` | Number of IDs requested per batch. |
|  | `limit` | `null` | Optional cap on identifiers processed. |
| `chembl` | `column` | `document_chembl_id` | Column with document identifiers exported by ChEMBL. |
|  | `chunk_size` | `50` | Batch size for API requests. |
|  | `timeout` | `30.0` | Request timeout in seconds. |
|  | `limit` | `null` | Optional cap on identifiers processed. |
| `all` | `column` | `document_chembl_id` | Column with identifiers processed by the combined pipeline. |
|  | `chunk_size` | `5` | ChEMBL request size used when `run_all` calls `cl.get_documents`. |
|  | `sleep` | `5.0` | Delay between PubMed polling cycles during enrichment. |
|  | `workers` | `1` | Worker threads orchestrating ChEMBL fetch + PubMed enrichment. |
|  | `batch_size` | `5` | PubMed request size passed to `fetch_pubmed_records`. |
|  | `timeout` | `30.0` | Request timeout applied to both ChEMBL and PubMed calls. |
|  | `limit` | `null` | Optional cap on identifiers handled by the combined run. |

#### Target pipeline (`target`)

| Sub-section | Key | Default | Description |
| --- | --- | --- | --- |
| `uniprot` | `column` | `uniprot_id` | Column with UniProt identifiers. |
|  | `data_dir` | `dictionary/uniprot` | Directory holding cached UniProt JSON files. |
|  | `limit` | `null` | Optional cap on identifiers processed. |
| `chembl` | `column` | `target_chembl_id` | Column with ChEMBL target identifiers. |
|  | `chunk_size` | `5` | Batch size for API requests. |
|  | `timeout` | `30.0` | Request timeout in seconds. |
|  | `limit` | `null` | Optional cap on identifiers processed. |
| `iuphar` | `target_csv` | `dictionary/_IUPHAR/_IUPHAR_target.csv` | Lookup table with IUPHAR target metadata. |
|  | `family_csv` | `dictionary/_IUPHAR/_IUPHAR_family.csv` | Lookup table with IUPHAR family metadata. |
|  | `limit` | `null` | Optional cap on identifiers processed. |
| `all` | `data_dir` | `dictionary/uniprot` | Directory containing cached UniProt data. |
|  | `target_csv` | `dictionary/_IUPHAR/_IUPHAR_target.csv` | IUPHAR target reference data. |
|  | `family_csv` | `dictionary/_IUPHAR/_IUPHAR_family.csv` | IUPHAR family reference data. |
|  | `chunk_size` | `5` | Batch size when combining all sources. |
|  | `timeout` | `30.0` | Request timeout in seconds. |
|  | `organism_csv` | `dictionary/_Target/targets_type.csv` | Taxonomy and target type mapping. |
|  | `uniprot_column` | `uniprot_id` | Column used to join UniProt data. |
|  | `chembl_out` | `null` | Optional override for the combined ChEMBL output path. |
|  | `uniprot_out` | `null` | Optional override for the combined UniProt output path. |
|  | `iuphar_out` | `null` | Optional override for the combined IUPHAR output path. |
|  | `limit` | `null` | Optional cap on identifiers processed. |

## Other external sources (`sources.*`)

| Section | Default base URL | Key settings |
| --- | --- | --- |
| `openalex` | `https://api.openalex.org` | `timeout_connect=5`, `timeout_read=30`, `retries=3`, `rps=5`, `burst=5`, `mailto="contact@example.org"`. |
| `crossref` | `https://api.crossref.org` | Same structure as `openalex`; provide a personal `mailto`. |
| `uniprot.api` | `https://rest.uniprot.org` | `timeout_connect=5`, `timeout_read=30`, `rps=25`, `burst=25`, `delay=0.25` seconds. |
| `uniprot.mapping` | `https://rest.uniprot.org/idmapping` | `poll_interval=0.5` seconds, `timeout=300.0` seconds, `cache_ttl=null`. |
| `iuphar` | `https://www.guidetopharmacology.org/services` | `timeout_connect=5`, `timeout_read=30`, `rps=5`, `burst=5`. |
| `pubchem` | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | `timeout_connect=5`, `timeout_read=60`, `retries=3`, `rps=5`, `burst=5`, `delay=0.2`, `cache_ttl=3600`, `prefer_local_values=true`, `cache_ttl_hours=null`. |
| `pubmed` | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils` | `timeout_connect=5`, `timeout_read=10`, `retries=2`. |
| `semantic_scholar` | `https://api.semanticscholar.org/graph/v1` | `timeout_connect=5`, `timeout_read=10`, `retries=2`. |

All URLs must comply with the respective service usage policies, including rate limits and contact information requirements.

## Local resources (`local`)

### Reference data (`local.resources`)

| Key | Default | Description |
| --- | --- | --- |
| `dictionary_dir` | `dictionary` | Root directory with lookup tables. |
| `iuphar_target_csv` | `dictionary/_IUPHAR/_IUPHAR_target.csv` | IUPHAR target mapping table. |
| `iuphar_family_csv` | `dictionary/_IUPHAR/_IUPHAR_family.csv` | IUPHAR family mapping table. |
| `uniprot_data_dir` | `dictionary/uniprot` | Cached UniProt JSON responses. |
| `organism_csv` | `dictionary/_Target/targets_type.csv` | Organism and taxonomy mapping. |
| `targets_type_csv` | `dictionary/_Target/targets_type.csv` | Target type classification table. |

### I/O defaults (`local.io`)

| Key | Default | Description |
| --- | --- | --- |
| `output_dir` | `data/output` | Base directory for generated datasets. |
| `cache_dir` | `.cache` | Location of the HTTP cache. |
| `csv_sep` | `,` | Default delimiter when reading and writing CSV files. |
| `csv_encoding` | `utf-8-sig` | Default encoding for CSV exports. |
| `na_markers` | `["#N/A"]` | Extra values treated as missing identifiers when reading CSV files. |
| `exist_ok` | `true` | Create directories automatically when `true`. |

### Initialisation workbooks (`local.init`)

| Key | Default | Description |
| --- | --- | --- |
| `same_doc` | `data/input/ChEMBL/ChEMBL_same_document_20_05.xlsx` | Workbook with same-document pairs for initialisation. |
| `all_doc` | `data/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx` | Workbook with cross-document pairs for initialisation. |
| `output_dir` | `data/output/ChEMBL/processed` | Destination for pre-processed initialisation files. |

## System settings (`system`)

| Sub-section | Key | Default | Description |
| --- | --- | --- | --- |
| `log` | `level` | `INFO` | Default logging level. |
|  | `format` | `[%(asctime)s] %(levelname)s %(name)s: %(message)s` | Log message format. |
|  | `datefmt` | `%Y-%m-%d %H:%M:%S` | Timestamp format. |
| `rate` | `global_rps` | `100` | Global requests-per-second budget shared across clients. |
|  | `global_burst` | `100` | Global token bucket burst capacity. |
|  | `limiter_cache_maxsize` | `128` | Maximum cached limiter instances. |
|  | `limiter_cache_ttl` | `600` | TTL for cached limiters in seconds. |
| `retry` | `max_attempts` | `3` | Number of retry attempts for recoverable errors. |
|  | `backoff_factor` | `0.5` | Base multiplier for exponential backoff. |
|  | `status_forcelist` | `[429, 500, 502, 503, 504]` | HTTP status codes that trigger retries. |
| `doc_type` | `weights` | `{pubmed: 4, openalex: 3, scholar: 2}` | Weighting applied to document sources. |
|  | `thresholds` | `{review: 1, experimental: 1, unknown: 2}` | Minimum counts for document type classification. |
|  | `limit` | `null` | Optional limit on classified records. |

## Environment variable overrides

Environment variables follow the pattern `CHEMBL_DA__SECTION__SUBSECTION__KEY`. For example:

```bash
export CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT="project/1.0 (mailto:me@example.com)"
export CHEMBL_DA__LOCAL__IO__OUTPUT_DIR=/mnt/datasets
```

Common short aliases:

| Variable | Maps to |
| --- | --- |
| `CHEMBL_DA_BASE` | `sources.chembl.api.chembl_base` |
| `CHEMBL_DA_RPS` | `sources.chembl.api.rps` |
| `CHEMBL_DA_BURST` | `sources.chembl.api.burst` |
| `CHEMBL_DA_TIMEOUT_CONNECT` | `sources.chembl.api.timeout_connect` |
| `CHEMBL_DA_TIMEOUT_READ` | `sources.chembl.api.timeout_read` |
| `CHEMBL_DA_CACHE_TTL` | `sources.chembl.cache.cache_ttl` |
| `CHEMBL_DA_CACHE_MAXSIZE` | `sources.chembl.cache.cache_maxsize` |
| `CHEMBL_DA_OUTDIR` | `local.io.output_dir` |
| `CHEMBL_DA_CACHE_DIR` | `local.io.cache_dir` |
| `CHEMBL_DA_GLOBAL_RPS` | `system.rate.global_rps` |
| `CHEMBL_DA_GLOBAL_BURST` | `system.rate.global_burst` |
| `CHEMBL_DA_LIMITER_CACHE_MAXSIZE` | `system.rate.limiter_cache_maxsize` |
| `CHEMBL_DA_LIMITER_CACHE_TTL` | `system.rate.limiter_cache_ttl` |
| `CHEMBL_DA_LOG_LEVEL` | `system.log.level` |
| `CHEMBL_DA_LOG_FORMAT` | `system.log.format` |
| `CHEMBL_DA_LOG_DATEFMT` | `system.log.datefmt` |
| `CHEMBL_DA_RETRY_MAX_ATTEMPTS` | `system.retry.max_attempts` |
| `CHEMBL_DA_RETRY_BACKOFF_FACTOR` | `system.retry.backoff_factor` |
| `CHEMBL_DA_DICT_DIR` | `local.resources.dictionary_dir` |
| `CHEMBL_DA_UNIPROT_DATA_DIR` | `local.resources.uniprot_data_dir` |
| `CHEMBL_DA_IUPHAR_TARGET_CSV` | `local.resources.iuphar_target_csv` |
| `CHEMBL_DA_IUPHAR_FAMILY_CSV` | `local.resources.iuphar_family_csv` |
| `CHEMBL_DA_ORGANISM_CSV` | `local.resources.organism_csv` |
| `CHEMBL_DA_TARGETS_TYPE_CSV` | `local.resources.targets_type_csv` |
| `CHEMBL_DA_OPENALEX_BASE` | `sources.openalex.base` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_CONNECT` | `sources.openalex.timeout_connect` |
| `CHEMBL_DA_OPENALEX_TIMEOUT_READ` | `sources.openalex.timeout_read` |
| `CHEMBL_DA_OPENALEX_MAILTO` | `sources.openalex.mailto` |
| `CHEMBL_DA_CROSSREF_BASE` | `sources.crossref.base` |
| `CHEMBL_DA_CROSSREF_TIMEOUT_CONNECT` | `sources.crossref.timeout_connect` |
| `CHEMBL_DA_CROSSREF_TIMEOUT_READ` | `sources.crossref.timeout_read` |
| `CHEMBL_DA_CROSSREF_MAILTO` | `sources.crossref.mailto` |
| `CHEMBL_DA_UNIPROT_BASE` | `sources.uniprot.api.base` |
| `CHEMBL_DA_IUPHAR_BASE` | `sources.iuphar.base` |
| `CHEMBL_DA_PUBCHEM_BASE` | `sources.pubchem.base` |

Any other key can be targeted using the long `CHEMBL_DA__...` form.

## CLI overrides

* Supply `--config` to point at an alternative YAML file; defaults to `config.yaml`.
* Pass `--print-config` to print the effective configuration (after environment and CLI overrides) and exit.
* Any CLI argument mapped via `apply_config_overrides` updates the configuration. For example `--batch-size 25` sets `sources.chembl.pipelines.activity.batch_size` for the current run.
* Nested parameters are changed via `config.yaml` or environment variables such as `CHEMBL_DA__SOURCES__CHEMBL__API__RPS=10`. Flags like `--sources.…` are not defined by the parsers.

## Validation workflow

At start-up the configuration loader:

1. Reads `config.yaml`.
2. Applies environment overrides and CLI-derived overrides.
3. Rejects unknown keys according to the JSON schema.
4. Ensures `output_dir` and `cache_dir` exist (creating them when `local.io.exist_ok` is `true`).

Keep the schema and documentation in sync when adding new configuration options.

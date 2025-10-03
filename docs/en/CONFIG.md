# Configuration Reference

## Overview

* All command-line tools load their defaults from `config/config.yaml` in the project root.
* Values are validated by `library.config.load_config`, which calls `Config.model_validate` from Pydantic. `config.schema.json` documents the same structure for tooling but is not executed during start-up.
* Settings can be overridden via a `config.local.yaml` placed next to the primary configuration (including custom paths supplied via `--config`), environment variables, and CLI flags. Precedence is: `config/config.yaml` < `config.local.yaml` < environment variables < CLI arguments.

### Runtime compatibility

The validated configuration assumes the dependency ranges declared in `pyproject.toml`:

| Component | Supported range |
| --- | --- |
| Python | `==3.11.*` |
| numpy | `>=2.3.3,<3.0` |
| pandas | `>=2.3.3,<3.0` |
| requests | `>=2.32.5,<3.0` |
| PyYAML | `>=6.0.3,<7.0` |
| openpyxl | `>=3.1.5,<4.0` |
| pyarrow | `>=17.0.0,<18.0` |
| jsonschema | `>=4.25.1,<5.0` |
| pandera | `>=0.26.1,<0.27` |
| pydantic | `>=2.11.9,<3.0` |
| cachetools | `>=5.3.3,<6.0` |

## Layout of `config/config.yaml`

The configuration is structured into logical blocks. High-level keys group settings for external services, local file paths, and system-wide behaviours.

| Section | Purpose |
| --- | --- |
| `sources` | Connectivity details for external services such as ChEMBL, UniProt, CrossRef, PubChem, and PubMed. |
| `local` | Paths to local resources, CSV defaults, and initialisation workbooks. |
| `system` | Logging, retry policy, global rate limiting, table quality profiling, and document classification weights. |
| `activity_bounds` | Rules for deriving numerical bounds from activity data. |
| `activity_enrichment` | Settings for annotating activities with derived labels and properties. |
| `testitem_molecule_enrichment` | Configuration for enriching test items with parent molecule data. |

Sensitive values (API tokens, personal e-mails) should be injected via environment variables rather than committed to the repository.

Path settings may reference the `$CHEMBL_DA_BASE_PATH` placeholder. During runtime, it resolves to the CLI `--base-path` argument, the `CHEMBL_DA_BASE_PATH` environment variable, or by default, `~/.local/share/chembl-da`. User-home shortcuts (`~`) are expanded before relative paths are resolved against the configuration file.

## `sources.chembl`

### API client (`sources.chembl.api`)

| Key | Default | Description |
| --- | --- | --- |
| `chembl_base` | `https://www.ebi.ac.uk/chembl/api/data` | Base URL for the ChEMBL REST API. |
| `timeout_connect` | `5.0` | Connection timeout in seconds. |
| `timeout_read` | `30.0` | Read timeout in seconds. |
| `retries` | `3` | Maximum number of attempts for failed requests. |
| `backoff_factor` | `0.5` | Multiplier for exponential backoff between retries. |
| `rps` | `5` | Allowed requests per second for the rate limiter. |
| `burst` | `5` | Bucket size used by the token bucket limiter. |
| `user_agent` | `chembl-da/1.0 (mailto:chembl-data@ebi.ac.uk)` | User-Agent header. Replace the contact with your own address before production use; the validator rejects placeholder emails like `contact@example.org`. |

### Response cache (`sources.chembl.cache`)

| Key | Default | Description |
| --- | --- | --- |
| `cache_ttl` | `3600` | Time-to-live for cached API responses in seconds. |
| `cache_maxsize` | `1024` | Maximum number of cached responses. |

### Molecule catalogue (`sources.chembl.molecule_catalog`)

| Key | Default | Description |
| --- | --- | --- |
| `cache_path` | `.../cache/molecule_parent_catalog.json` | Location of the JSON cache for molecule parent-child relationships. |
| `sqlite_path` | `.../cache/molecule_parent_catalog.sqlite` | Location of the SQLite cache for parent-child lookups. |
| `endpoint` | `molecule` | ChEMBL REST resource for cache refresh. |
| `child_field` | `molecule_chembl_id` | Field for the child molecule identifier. |
| `parent_field` | `parent_molecule_chembl_id` | Field for the parent molecule identifier. |
| `hierarchy_lookup_path` | `dictionary/_testitem/molecule_hierarchy.csv` | Optional local CSV for seeding parent-child relationships. |
| `hierarchy_lookup_encoding` | `utf-8-sig` | Encoding for the hierarchy lookup CSV. |
| `hierarchy_lookup_delimiter` | `,` | Delimiter for the hierarchy lookup CSV. |
| `force_refresh_existing` | `false` | If `true`, rebuilds parent relationships even if they already exist in the data. |
| `fields` | `('molecule_chembl_id', 'parent_molecule_chembl_id')` | Fields requested from the ChEMBL API. |
| `filters` | `{'parent_molecule_chembl_id__isnull': 'false'}` | Filters for the API call to get only records with parents. |
| `page_size` | `500` | Records requested per API page. |
| `fallback_single_limit` | `null` | Cap on single-molecule fallback requests after a bulk fetch fails. |

### Pipelines (`sources.chembl.pipelines`)

Defaults for each pipeline CLI utility.

#### Activity pipeline (`activity`)

| Key | Default | Description |
| --- | --- | --- |
| `column` | `activity_id` | Input column containing activity identifiers. |
| `batch_size` | `5` | Batch size for API requests. |
| `workers` | `1` | Number of parallel worker processes. |
| `timeout` | `30.0` | Request timeout in seconds. |
| `limit` | `null` | Optional cap on identifiers processed. |
| `offset` | `0` | Number of identifiers to skip. |
| `dry_run` | `false` | If `true`, skip network calls and file generation. |

#### Assay pipeline (`assay`)

| Key | Default | Description |
| --- | --- | --- |
| `column` | `assay_chembl_id` | Input column with assay identifiers. |
| `batch_size` | `10` | Batch size for API requests. |
| `timeout` | `30.0` | Request timeout in seconds. |
| `limit` | `null` | Optional cap on identifiers processed. |

#### Test item pipeline (`testitem`)

| Key | Default | Description |
| --- | --- | --- |
| `column` | `molecule_chembl_id` | Input column with compound identifiers. |
| `batch_size` | `1000` | Batch size for API requests (max 1000). |
| `timeout` | `30.0` | Request timeout in seconds. |
| `limit` | `null` | Optional cap on identifiers processed. |
| `offset` | `0` | Starting position for paginated exports. |
| `fields` | `(...)` | Default list of ChEMBL/PubChem fields to request. |
| `request_limit` | `1000` | Hard cap on paginated requests per run. |
| `retries` | `null` | Max retries for failed API calls (uses global `system.retry` if `null`). |
| `backoff_factor` | `null` | Backoff multiplier (uses global `system.retry` if `null`). |
| `batch_retry` | (sub-model) | Settings for retrying with a smaller batch size. |
| `parent_watchdog_idle_minutes` | `0.0` | Timeout in minutes for parent catalog idle time. |
| `execution_budget_minutes` | `null` | Optional execution time limit for the pipeline. |

#### Document pipeline (`document`)

| Sub-section | Key | Default | Description |
| --- | --- | --- | --- |
| `pubmed` | `column` | `PMID` | Column with PubMed identifiers. |
| | `sleep` | `5.0` | Delay between polling cycles. |
| | `workers` | `1` | Number of worker threads. |
| | `batch_size` | `100` | Number of IDs per PubMed request. |
| | `limit` | `null` | Optional cap on identifiers. |
| `chembl` | `column` | `document_chembl_id` | Column with ChEMBL document identifiers. |
| | `chunk_size` | `50` | Batch size for API requests. |
| | `timeout` | `30.0` | Request timeout in seconds. |
| | `limit` | `null` | Optional cap on identifiers. |
| `all` | `column` | `document_chembl_id` | Column for the combined pipeline. |
| | `chunk_size` | `5` | ChEMBL request size. |
| | `sleep` | `5.0` | PubMed polling delay. |
| | `workers` | `1` | Worker threads for combined run. |
| | `batch_size` | `50` | PubMed request size. |
| | `timeout` | `30.0` | Request timeout. |
| | `limit` | `null` | Optional cap on identifiers. |

#### Target pipeline (`target`)

| Sub-section | Key | Default | Description |
| --- | --- | --- | --- |
| `uniprot` | `column` | `uniprot_id` | Column with UniProt identifiers. |
| | `data_dir` | `dictionary/_target/_uniprot` | Directory for cached UniProt JSON files. |
| | `limit` | `null` | Optional cap on identifiers. |
| `chembl` | `column` | `target_chembl_id` | Column with ChEMBL target identifiers. |
| | `chunk_size` | `5` | Batch size for API requests. |
| | `timeout` | `30.0` | Request timeout in seconds. |
| | `limit` | `null` | Optional cap on identifiers. |
| `iuphar` | `target_csv` | `.../_IUPHAR_target.csv` | Lookup table for IUPHAR target metadata. |
| | `family_csv` | `.../_IUPHAR_family.csv` | Lookup table for IUPHAR family metadata. |
| | `limit` | `null` | Optional cap on identifiers. |
| `all` | `uniprot_column` | `uniprot_id` | Column used to join UniProt data. |
| | `...` | | See `TargetAllCfg` model for other fields. |

## Other External Sources (`sources.*`)

| Section | Default `rps` | Default `burst` | Key settings |
| --- | --- | --- | --- |
| `openalex` | `4` | `5` | `base`, `timeout_*`, `retries`, `mailto` |
| `crossref` | `4` | `5` | `base`, `timeout_*`, `retries`, `mailto` |
| `uniprot.api` | `4` | `5` | `base`, `timeout_*`, `delay` |
| `uniprot.mapping` | `N/A` | `N/A` | `base`, `poll_interval`, `timeout`, `cache_ttl` |
| `iuphar` | `4` | `5` | `base`, `timeout_*` |
| `pubchem` | `3` | `5` | See `PubChemCfg` model for all settings. |
| `pubmed` | `null` | `null` | `base`, `timeout_*`, `retries` |
| `semantic_scholar` | `null` | `null` | `base`, `timeout_*`, `retries` |

### PubChem Lookups (`sources.pubchem`)

The `PubChemCfg` model contains numerous options for PubChem integration. Refer to the model in `library/config.py` for a complete and authoritative list.

## Local Resources (`local`)

### Reference Data (`local.resources`)

| Key | Default |
| --- | --- |
| `dictionary_dir` | `dictionary/` |
| `iuphar_target_csv` | `dictionary/_target/_IUPHAR/_IUPHAR_target.csv` |
| `iuphar_family_csv` | `dictionary/_target/_IUPHAR/_IUPHAR_family.csv` |
| `uniprot_data_dir` | `dictionary/_target/_uniprot` |
| `targets_type_csv` | `dictionary/_target/targets_type.csv` |

### I/O Defaults (`local.io`)

| Key | Default |
| --- | --- |
| `output_dir` | `$CHEMBL_DA_BASE_PATH/output` |
| `cache_dir` | `~/.cache/chembl-da` |
| `csv_sep` | `,` |
| `csv_encoding` | `utf-8-sig` |
| `csv_fallback_encodings`| `('utf-8', 'cp1252', ...)` |
| `csv_fallback_separators`| `('\t', ';')` |
| `na_markers` | `('#N/A',)` |
| `csv_chunksize` | `10000` |
| `exist_ok` | `true` |
| `keep_na_markers` | `false` |

## System Settings (`system`)

| Sub-section | Key | Default |
| --- | --- | --- |
| `log` | `level` | `INFO` |
| `rate` | `global_rps` | `8` |
| | `global_burst` | `8` |
| | `limiter_cache_maxsize` | `128` |
| | `limiter_cache_ttl` | `600` |
| `retry` | `max_attempts` | `3` |
| | `backoff_factor` | `0.5` |
| | `backoff_cap` | `null` |
| | `status_forcelist` | `[429, 500, 502, 503, 504]` |

## Environment Variable Overrides

Environment variables follow the pattern `CHEMBL_DA__SECTION__SUBSECTION__KEY`. For example: `CHEMBL_DA__SOURCES__CHEMBL__API__RPS=10`.
Many common options also have short aliases. Refer to the `_ALIAS_MAP` dictionary in `library/config.py` for the complete and authoritative list of all aliases.
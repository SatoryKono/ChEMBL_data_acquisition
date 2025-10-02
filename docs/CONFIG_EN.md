# Configuration Reference

## Overview

* All command-line tools load their defaults from [`config/config.yaml`](../config/config.yaml) in the project root.
* Values are validated by `library.config.load_config`, which calls `Config.model_validate` from Pydantic. [`config.schema.json`](../config.schema.json) documents the same structure for tooling but is not executed during start-up.
* Settings can be overridden via a `config.local.yaml` placed next to the primary configuration (including custom paths supplied via `--config`), environment variables and CLI flags. Precedence is: `config/config.yaml` < `config.local.yaml` < environment variables < CLI arguments.

## Layout of `config/config.yaml`

| Section | Purpose |
| --- | --- |
| `sources` | Connectivity details for external services such as ChEMBL, UniProt, CrossRef, PubChem and PubMed. |
| `local` | Paths to local resources, CSV defaults and initialisation workbooks. |
| `system` | Logging, retry policy, global rate limiting, table quality profiling and document classification weights. |

Sensitive values (API tokens, personal e-mails) should be injected via environment variables rather than committed to the repository.

Path settings may reference the ``$CHEMBL_DA_BASE_PATH`` placeholder. During
runtime it resolves to the CLI ``--base-path`` argument, the
``CHEMBL_DA_BASE_PATH`` environment variable or, by default,
``~/.local/share/chembl-da``. User-home shortcuts (``~``) are expanded before
relative paths are resolved against the configuration file.

## `sources.chembl`

### API client (`sources.chembl.api`)

| Key | Default | Description |
| --- | --- | --- |
| `chembl_base` | `https://www.ebi.ac.uk/chembl/api/data` | Base URL for the ChEMBL REST API. |
| `timeout_connect` | `5` | Connection timeout in seconds. |
| `timeout_read` | `30` | Read timeout in seconds. |
| `retries` | `3` | Maximum number of attempts performed by higher-level clients; the shared HTTP adapter does not retry automatically. |
| `backoff_factor` | `0.5` | Multiplier for exponential backoff between retries. |
| `rps` | `20` | Allowed requests per second for the rate limiter. |
| `burst` | `20` | Bucket size used by the token bucket limiter. |
| `user_agent` | `chembl-da/1.0 (mailto:chembl-data@ebi.ac.uk)` | User-Agent header. Replace the contact with your own address before production use; the validator still rejects the placeholder `contact@example.org`. Set via `CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT`. |


### Response cache (`sources.chembl.cache`)

| Key | Default | Description |
| --- | --- | --- |
| `cache_ttl` | `3600` | Time-to-live for cached API responses in seconds. |
| `cache_maxsize` | `1024` | Maximum number of cached responses. |

<a id="sources-chembl-molecule-catalog"></a>
### Molecule catalogue (`sources.chembl.molecule_catalog`)

| Key | Default | Description |
| --- | --- | --- |
| `cache_path` | `"$CHEMBL_DA_BASE_PATH/cache/molecule_parent_catalog.json"` | Location of the JSON cache storing molecule parent-child relationships reused by enrichment jobs; override via `CHEMBL_DA_MOLECULE_CATALOG_CACHE` (alias for `CHEMBL_DA__SOURCES__CHEMBL__MOLECULE_CATALOG__CACHE_PATH`). |
| `sqlite_path` | `"$CHEMBL_DA_BASE_PATH/cache/molecule_parent_catalog.sqlite"` | Location of the SQLite cache powering parent-child lookups; override via `CHEMBL_DA_SOURCES_CHEMBL_MOLECULE_CATALOG_SQLITE_PATH` (or the canonical `CHEMBL_DA__SOURCES__CHEMBL__MOLECULE_CATALOG__SQLITE_PATH`). |
| `endpoint` | `molecule` | ChEMBL REST resource queried when the cache needs to be refreshed. |
| `child_field` | `molecule_chembl_id` | JSON field containing the child molecule identifier extracted from API responses. |
| `parent_field` | `parent_molecule_chembl_id` | JSON field containing the parent molecule identifier extracted from API responses. |
| `force_refresh_existing` | `false` | When `true`, rebuilds parent relationships even if the incoming dataset already contains parent identifiers, ensuring the cache wins over source columns. |
| `fields` | `['molecule_chembl_id', 'parent_molecule_chembl_id']` | List of fields requested from the ChEMBL API when populating or refreshing the catalogue; extend to retrieve extra metadata alongside identifiers. |
| `filters` | `{'parent_molecule_chembl_id__isnull': 'false'}` | Query parameters appended to every API call; defaults keep only rows that already have parent assignments in ChEMBL. |
| `hierarchy_lookup_path` | `../dictionary/_testitem/molecule_hierarchy.csv` | Optional CSV used as an offline seed for parent-child relationships before querying ChEMBL; override when distributing a curated hierarchy snapshot or relocating the dictionary folder. |
| `hierarchy_lookup_encoding` | `utf-8-sig` | Text encoding applied when reading the hierarchy lookup CSV; change when the snapshot is saved with a different charset (for example Latin-1 from legacy exports). |
| `hierarchy_lookup_delimiter` | `,` | Delimiter expected by the hierarchy lookup loader; override for semicolon- or tab-separated snapshots produced by regional data teams. |
| `page_size` | `500` | Number of records requested per API page while rebuilding the catalogue. |
| `fallback_single_limit` | `null` | Caps the number of single-molecule fallback requests performed after bulk fetching fails; `null` keeps the fallback unlimited. |

### Pipelines (`sources.chembl.pipelines`)

Each sub-section below defines defaults for the respective CLI utility. CLI arguments are merged back into the configuration before execution.

#### Activity pipeline (`activity`)

| Key | Default | Description |
| --- | --- | --- |
| `column` | `activity_chembl_id` | Input column containing activity identifiers. |
| `batch_size` | `50` | Batch size for API requests. |
| `timeout` | `30.0` | Request timeout in seconds. |
| `workers` | `1` | Number of parallel worker processes handling API batches. |
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

This block controls derived annotations appended to every activity row. It works alongside the separate
`activity_bounds` configuration (documented below) that is responsible solely for calculating numerical limits.
`activity_enrichment` itself is divided into two independent sub-sections that can be toggled individually.

##### Action type labelling (`activity_enrichment.action_type`)

* `enabled` — master switch (`true` by default).
* `column` — name of the output column that stores the resolved label (`action_type`).
* Logging flags control diagnostics when the source data does not produce a value:
  * `log_missing` emits a warning when the label cannot be determined (`true`).
  * `log_distribution` prints summary statistics after enrichment (`true`).
* `metrics` maps activity measurements (IC50, EC50, etc.) to default labels such as `inhibition` or `activation`. Extend the
  mapping to support additional metrics.
* `triages`, `functionality` and `mechanism` are optional overrides for explicit text matches. The defaults keep `triages` and
  `mechanism` empty while normalising common functional annotations (agonist, antagonist, etc.).
* The helper field lists (`triage_fields`, `functionality_fields`, `mechanism_fields`) decide which columns are scanned for
  keywords or manual labels before applying metric defaults.
* `allowlist` enumerates labels allowed in the output. Values not present in the list fall back to `fallback` after logging the
  anomaly.
* `positive_label` and `negative_label` define human-readable aliases used when the data represents positive/negative modulators
  (`PAM`/`NAM`). The neutral fallback is `unknown`.

##### Activity properties flattening (`activity_enrichment.activity_properties`)

* `enabled` — feature switch (`true`).
* `column` — name of the raw JSON-like source column (`activity_properties`).
* `summary_column` — reserved for the future text renderer. The current implementation keeps the JSON payload in `column`,
  generates only the deterministic fingerprint in `hash_column` and does not emit a separate summary field.
* `name_field`, `value_field`, `units_field` identify keys within each property record (`type`, `value`, `units`).
* `separator` and `pair_separator` are legacy formatting knobs retained for backwards compatibility; the current
  JSON serialisation ignores them.

* `drop_source_column` removes the original structured column after summarisation (`true`).
* Logging flags default to `false`, muting missing/distribution reports unless troubleshooting is required.
* `allowlist` restricts which property groups are retained (measurement, assay, comments, effect_features, triage, mechanism,
  functionality).
* `hash_column` stores a deterministic fingerprint of the preserved properties (`properties_hash`), enabling change detection in
  downstream jobs.

##### Activity bounds (`activity_bounds`)

The activity pipeline enriches raw ChEMBL payloads with canonical lower/upper bounds using the rules implemented by
`compute_activity_bounds` in `library.processing.activity`. Configuration for this feature is stored in the
`activity_bounds` block (separate from `activity_enrichment`) and controls the following deterministic stages (executed
in order for every row):

1. Use `standard_lower_value`/`standard_upper_value` when both are populated.
2. Combine `standard_value` with the opposite explicit limit (for example `standard_upper_value`) and fill the missing bound.
3. Inspect `standard_relation` when `enable_from_relation` is `true`, mapping operators such as `=`, `≈`, `>=`, `<=`, `between`
   and `range` to appropriate bounds.
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

The CLI only exposes high-level switches such as `--batch-size` or `--dry-run`; enrichment-specific options must be changed in the configuration file or via the corresponding `CHEMBL_DA__ACTIVITY_BOUNDS__*` variables. CLI values still win over file/env defaults for overlapping keys declared on the parser (column, batch size, limits).

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
| `batch_size` | `1000` | Batch size for API requests. |
| `timeout` | `30.0` | Request timeout in seconds. |
| `limit` | `null` | Optional cap on identifiers processed. |
| `offset` | `0` | Starting position for paginated exports, enabling resumable runs. |
| `request_limit` | `1000` | Hard cap on the number of paginated requests performed in a single execution. |
| `retries` | `5` | Maximum number of retry attempts applied to failed API calls. |
| `backoff_factor` | `0.5` | Multiplier controlling exponential backoff between retry attempts. |
| `fields` | `['molecule_chembl_id', 'parent_molecule_chembl_id', 'pref_name', 'max_phase', 'molecule_type', 'first_approval', 'oral', 'parenteral', 'topical', 'black_box_warning', 'structure_type', 'molecule_structures.canonical_smiles', 'molecule_structures.standard_inchi', 'molecule_structures.standard_inchi_key', 'pubchem_cid', 'pubchem_iupac_name', 'pubchem_molecular_formula', 'pubchem_isomeric_smiles', 'pubchem_canonical_smiles', 'pubchem_inchi', 'pubchem_inchikey']` | List of ChEMBL and PubChem fields requested for each test item batch. |


Pagination now tracks both the starting `offset` and a `request_limit`, allowing operators to resume large exports and bound the
number of pages consumed per run. Network failures are handled via the shared retry policy (`retries`/`backoff_factor`) that
applies exponential backoff between attempts before surfacing the error.

#### Test item molecule enrichment (`testitem_molecule_enrichment`)

| Key | Default | Description |
| --- | --- | --- |
| `enable` | `true` | Master switch for the enrichment stage that derives salt identifiers and catalogue flags. |
| `sources.molecule_catalog_path` | `../dictionary/_testitem/molecule_catalog.csv` | CSV with `molecule_chembl_id` and the `natural_product`/`prodrug`/`polymer_flag` columns. |
| `sources.molecule_hierarchy_path` | `../dictionary/_testitem/molecule_hierarchy.csv` | CSV that maps derivatives to their parent molecule (`molecule_chembl_id`, `parent_molecule_chembl_id`). |
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
|  | `data_dir` | `../dictionary/_target/_uniprot` | Directory holding cached UniProt JSON files. |
|  | `limit` | `null` | Optional cap on identifiers processed. |
| `chembl` | `column` | `target_chembl_id` | Column with ChEMBL target identifiers. |
|  | `chunk_size` | `5` | Batch size for API requests. |
|  | `timeout` | `30.0` | Request timeout in seconds. |
|  | `limit` | `null` | Optional cap on identifiers processed. |
| `iuphar` | `target_csv` | `../dictionary/_target/_IUPHAR/_IUPHAR_target.csv` | Lookup table with IUPHAR target metadata. |
|  | `family_csv` | `../dictionary/_target/_IUPHAR/_IUPHAR_family.csv` | Lookup table with IUPHAR family metadata. |
|  | `limit` | `null` | Optional cap on identifiers processed. |
| `all` | `data_dir` | `../dictionary/_target/_uniprot` | Directory containing cached UniProt data. |
|  | `target_csv` | `../dictionary/_target/_IUPHAR/_IUPHAR_target.csv` | IUPHAR target reference data. |
|  | `family_csv` | `../dictionary/_target/_IUPHAR/_IUPHAR_family.csv` | IUPHAR family reference data. |
|  | `chunk_size` | `5` | Batch size when combining all sources. |
|  | `timeout` | `30.0` | Request timeout in seconds. |

|  | `uniprot_column` | `uniprot_id` | Column used to join UniProt data. |
|  | `chembl_out` | `null` | Optional override for the combined ChEMBL output path. |
|  | `uniprot_out` | `null` | Optional override for the combined UniProt output path. |
|  | `iuphar_out` | `null` | Optional override for the combined IUPHAR output path. |
|  | `limit` | `null` | Optional cap on identifiers processed. |

> Target taxonomy (`type` column and classifier flags) is computed by a built-in
> module that uses UniProt lineage fields (`genus`, `lineage_superkingdom`,
> `lineage_phylum`, `lineage_class`), the `taxon_id`, and ChEMBL's
> `species_group_flag`. No external organism lookup CSV is required.

## Other external sources (`sources.*`)

| Section | Default base URL | Key settings |
| --- | --- | --- |
| `openalex` | `https://api.openalex.org` | `timeout_connect=5`, `timeout_read=30`, `retries=3`, `rps=5`, `burst=5`, `mailto="chembl-data@ebi.ac.uk"` (replace with your own contact; placeholder domains such as `example.org` are rejected). |
| `crossref` | `https://api.crossref.org` | Same structure as `openalex`; provide a personal `mailto`. |
| `uniprot.api` | `https://rest.uniprot.org` | `timeout_connect=5`, `timeout_read=30`, `rps=25`, `burst=25`, `delay=0.25` seconds. |
| `uniprot.mapping` | `https://rest.uniprot.org/idmapping` | `poll_interval=0.5` seconds, `timeout=300.0` seconds, `cache_ttl=null`. |
| `iuphar` | `https://www.guidetopharmacology.org/services` | `timeout_connect=5`, `timeout_read=30`, `rps=5`, `burst=5`. |
| `pubchem` | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | See [detailed PubChem options](#pubchem-lookups-sourcespubchem). |
| `pubmed` | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils` | `timeout_connect=5`, `timeout_read=10`, `retries=2`, optional `rps`/`burst` overrides for document acquisition. |
| `semantic_scholar` | `https://api.semanticscholar.org/graph/v1` | `timeout_connect=5`, `timeout_read=10`, `retries=2`, optional `rps`/`burst` overrides for document acquisition. |

All URLs must comply with the respective service usage policies, including rate limits and contact information requirements.

### PubChem lookups (`sources.pubchem`)

PubChem augmentation is primarily used by the test item pipeline. Every key below mirrors [`config/config.yaml`](../config/config.yaml) and
can be overridden through environment variables (see [Environment variables](#environment-variables)). The table lists both the
auto-generated `CHEMBL_DA_SOURCES_PUBCHEM_*` aliases and the generic `CHEMBL_DA__SOURCES__PUBCHEM__*` form. The base URL also
supports the short alias documented in the [Environment variable aliases](#environment-variables) table.

| Key | Default | Description | Environment override(s) |
| --- | --- | --- | --- |
| `enable` | `true` | Master switch enabling PubChem enrichment for test item records. | `CHEMBL_DA_SOURCES_PUBCHEM_ENABLE`, `CHEMBL_DA__SOURCES__PUBCHEM__ENABLE` |
| `base` | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | PubChem PUG REST endpoint queried for metadata. | `CHEMBL_DA_PUBCHEM_BASE`, `CHEMBL_DA_SOURCES_PUBCHEM_BASE`, `CHEMBL_DA__SOURCES__PUBCHEM__BASE` |
| `timeout_connect` | `5` | Connection timeout (seconds) for establishing new HTTP sessions. | `CHEMBL_DA_SOURCES_PUBCHEM_TIMEOUT_CONNECT`, `CHEMBL_DA__SOURCES__PUBCHEM__TIMEOUT_CONNECT` |
| `timeout_read` | `60` | Read timeout (seconds) waiting for server responses. | `CHEMBL_DA_SOURCES_PUBCHEM_TIMEOUT_READ`, `CHEMBL_DA__SOURCES__PUBCHEM__TIMEOUT_READ` |
| `timeout_seconds` | `30.0` | Upper bound for a single CID resolution attempt, including retries. | `CHEMBL_DA_SOURCES_PUBCHEM_TIMEOUT_SECONDS`, `CHEMBL_DA__SOURCES__PUBCHEM__TIMEOUT_SECONDS` |
| `retries` | `3` | Number of attempts executed by the PubChem retry loop before giving up. | `CHEMBL_DA_SOURCES_PUBCHEM_RETRIES`, `CHEMBL_DA__SOURCES__PUBCHEM__RETRIES` |
| `rps` | `5` | Per-service request-per-second budget used by the rate limiter. | `CHEMBL_DA_SOURCES_PUBCHEM_RPS`, `CHEMBL_DA__SOURCES__PUBCHEM__RPS` |
| `burst` | `5` | Token bucket size paired with the `rps` limit. | `CHEMBL_DA_SOURCES_PUBCHEM_BURST`, `CHEMBL_DA__SOURCES__PUBCHEM__BURST` |
| `delay` | `0.2` | Fixed pause (seconds) inserted between retry attempts. | `CHEMBL_DA_SOURCES_PUBCHEM_DELAY`, `CHEMBL_DA__SOURCES__PUBCHEM__DELAY` |
| `backoff_initial_seconds` | `0.5` | Initial exponential backoff applied after 429/5xx responses. | `CHEMBL_DA_SOURCES_PUBCHEM_BACKOFF_INITIAL_SECONDS`, `CHEMBL_DA__SOURCES__PUBCHEM__BACKOFF_INITIAL_SECONDS` |
| `resolve_order` | `cache → smiles → inchikey → inchi → pref_name` | Order in which lookup strategies are attempted when resolving PubChem CIDs. | `CHEMBL_DA_SOURCES_PUBCHEM_RESOLVE_ORDER`, `CHEMBL_DA__SOURCES__PUBCHEM__RESOLVE_ORDER` |
| `cache_ttl` | `3600` | Lifespan (seconds) of the in-memory HTTP response cache. | `CHEMBL_DA_SOURCES_PUBCHEM_CACHE_TTL`, `CHEMBL_DA__SOURCES__PUBCHEM__CACHE_TTL` |
| `cache_maxsize` | `1024` | Maximum number of entries retained by the in-memory HTTP response cache. | `CHEMBL_DA_SOURCES_PUBCHEM_CACHE_MAXSIZE`, `CHEMBL_DA__SOURCES__PUBCHEM__CACHE_MAXSIZE` |
| `cache_ttl_hours` | `null` | Optional expiry (hours) for the persisted CID cache; `null` keeps entries indefinitely. | `CHEMBL_DA_SOURCES_PUBCHEM_CACHE_TTL_HOURS`, `CHEMBL_DA__SOURCES__PUBCHEM__CACHE_TTL_HOURS` |
| `cid_cache_path` | `"$CHEMBL_DA_BASE_PATH/cache/pubchem_cid_cache.json"` | Path to a JSON file storing resolved CIDs for re-use across runs. | `CHEMBL_DA_SOURCES_PUBCHEM_CID_CACHE_PATH`, `CHEMBL_DA__SOURCES__PUBCHEM__CID_CACHE_PATH` |
| `batch_size` | `50` | Number of rows processed per PubChem batch request; concurrency never exceeds `min(batch_size, rps)`. | `CHEMBL_DA_SOURCES_PUBCHEM_BATCH_SIZE`, `CHEMBL_DA__SOURCES__PUBCHEM__BATCH_SIZE` |
| `prefer_local_smiles` | `false` | Skip remote lookups when local SMILES/InChIKey columns are already populated. | `CHEMBL_DA_SOURCES_PUBCHEM_PREFER_LOCAL_SMILES`, `CHEMBL_DA__SOURCES__PUBCHEM__PREFER_LOCAL_SMILES` |
| `prefer_local_values` | `true` | Preserve existing `pubchem_*` columns when lookups return empty payloads. | `CHEMBL_DA_SOURCES_PUBCHEM_PREFER_LOCAL_VALUES`, `CHEMBL_DA__SOURCES__PUBCHEM__PREFER_LOCAL_VALUES` |
| `use_parent_for_salts` | `true` | Escalate to parent structures when salt-specific lookups fail. | `CHEMBL_DA_SOURCES_PUBCHEM_USE_PARENT_FOR_SALTS`, `CHEMBL_DA__SOURCES__PUBCHEM__USE_PARENT_FOR_SALTS` |
| `allow_polymer` | `false` | Permit lookups that resolve to polymer or mixture entries. | `CHEMBL_DA_SOURCES_PUBCHEM_ALLOW_POLYMER`, `CHEMBL_DA__SOURCES__PUBCHEM__ALLOW_POLYMER` |
| `write_not_found_literal` | `false` | Write the literal `Not Found` when PubChem fails to return a CID. | `CHEMBL_DA_SOURCES_PUBCHEM_WRITE_NOT_FOUND_LITERAL`, `CHEMBL_DA__SOURCES__PUBCHEM__WRITE_NOT_FOUND_LITERAL` |

> Tip: `resolve_order` accepts any combination of supported strategies—adjust the list to prioritise cached values or specific
> identifiers while keeping `cache` first to honour warm CID lookups.

## Local resources (`local`)

### Reference data (`local.resources`)

| Key | Default | Description |
| --- | --- | --- |
| `dictionary_dir` | `../dictionary` | Root directory with lookup tables. |
| `iuphar_target_csv` | `../dictionary/_target/_IUPHAR/_IUPHAR_target.csv` | IUPHAR target mapping table. |
| `iuphar_family_csv` | `../dictionary/_target/_IUPHAR/_IUPHAR_family.csv` | IUPHAR family mapping table. |
| `uniprot_data_dir` | `../dictionary/_target/_uniprot` | Cached UniProt JSON responses. |
| `targets_type_csv` | `../dictionary/_target/targets_type.csv` | Target type classification table. |


The `dictionary/_target` folder mirrors the current repository layout; all
IUPHAR and UniProt lookups are stored there by default.


### I/O defaults (`local.io`)

| Key | Default | Description |
| --- | --- | --- |
| `output_dir` | `"$CHEMBL_DA_BASE_PATH/output"` | Base directory for generated datasets. |
| `cache_dir` | `"~/.cache/chembl-da"` | Location of the HTTP cache. |
| `csv_sep` | `,` | Default delimiter when reading and writing CSV files. |
| `csv_fallback_separators` | `["\t", ";"]` | Additional delimiters tried when the primary separator does not expose the requested column. |
| `csv_encoding` | `utf-8-sig` | Default encoding for CSV exports. |
| `csv_chunksize` | `10000` | Rows processed per batch by deterministic CSV writers; see [`config/config.yaml`](../config/config.yaml). |
| `na_markers` | `["#N/A"]` | Extra values treated as missing identifiers when reading CSV files. |
| `keep_na_markers` | `false` | Preserve identifiers matching `na_markers` instead of filtering them out. |
| `exist_ok` | `true` | Create directories automatically when `true`. |

### Initialisation workbooks (`local.init`)

| Key | Default | Description |
| --- | --- | --- |
| `same_doc` | `"$CHEMBL_DA_BASE_PATH/input/ChEMBL/ChEMBL_same_document_20_05.xlsx"` | Workbook with same-document pairs for initialisation. |
| `all_doc` | `"$CHEMBL_DA_BASE_PATH/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx"` | Workbook with cross-document pairs for initialisation. |
| `output_dir` | `"$CHEMBL_DA_BASE_PATH/output/ChEMBL/processed"` | Destination for pre-processed initialisation files. |

Paths under `data/input/ChEMBL/*.xlsx` are placeholders included for local smoke tests. Replace them with the workbooks prepared by your organisation (or copy the manually supplied files into the desired location) before starting the initialisation routines. Refer to [docs/USAGE_EN.md](./USAGE_EN.md) for guidance on preparing the input templates.

## System settings (`system`)

| Sub-section | Key | Default | Description |
| --- | --- | --- | --- |
| `log` | `level` | `INFO` | Default logging level. Structured JSON output keeps a fixed schema handled by [`library/logging_setup.py`](../library/logging_setup.py), so message and timestamp formatting are not configurable here. |
| `rate` | `global_rps` | `100` | Global requests-per-second budget shared across clients. |
|  | `global_burst` | `100` | Global token bucket burst capacity. |
|  | `limiter_cache_maxsize` | `128` | Maximum cached limiter instances. |
|  | `limiter_cache_ttl` | `600` | TTL for cached limiters in seconds. |
| `retry` | `max_attempts` | `3` | Number of retry attempts for recoverable errors. |
|  | `backoff_factor` | `0.5` | Base multiplier for exponential backoff. |
|  | `status_forcelist` | `[429, 500, 502, 503, 504]` | HTTP status codes that trigger retries. |
| `doc_quality` | `enable` | `true` | Toggle generation of table quality reports. |
|  | `sample_rows` | `null` | Limit analysis to the first `N` rows; `null` processes the full dataset. |
|  | `include_columns` | `null` | Optional allow list of column names to profile. |
|  | `exclude_columns` | `null` | Optional deny list of column names to skip. |
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
| `CHEMBL_DA_PUBMED_RPS` | `sources.pubmed.rps` |
| `CHEMBL_DA_PUBMED_BURST` | `sources.pubmed.burst` |
| `CHEMBL_DA_SEMANTIC_SCHOLAR_RPS` | `sources.semantic_scholar.rps` |
| `CHEMBL_DA_SEMANTIC_SCHOLAR_BURST` | `sources.semantic_scholar.burst` |
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
| `CHEMBL_DA_RETRY_MAX_ATTEMPTS` | `system.retry.max_attempts` |
| `CHEMBL_DA_RETRY_BACKOFF_FACTOR` | `system.retry.backoff_factor` |
| `CHEMBL_DA_DICT_DIR` | `local.resources.dictionary_dir` |
| `CHEMBL_DA_UNIPROT_DATA_DIR` | `local.resources.uniprot_data_dir` |
| `CHEMBL_DA_IUPHAR_TARGET_CSV` | `local.resources.iuphar_target_csv` |
| `CHEMBL_DA_IUPHAR_FAMILY_CSV` | `local.resources.iuphar_family_csv` |
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
| `CHEMBL_DA_PUBCHEM_USER_AGENT` | `sources.pubchem.user_agent` |

Any other key can be targeted using the long `CHEMBL_DA__...` form.

> Structured logging details live in [`library/log.py`](../library/log.py) and [`library/logging_setup.py`](../library/logging_setup.py), explaining why the JSON layout (including date formatting) cannot be overridden via `config.yaml`.

## CLI overrides

* Supply `--config` to point at an alternative YAML file; defaults to `config/config.yaml`.
* Pass `--print-config` to print the effective configuration (after environment and CLI overrides) and exit.
* Any CLI argument mapped via `apply_config_overrides` updates the configuration. For example `--batch-size 25` sets `sources.chembl.pipelines.activity.batch_size` for the current run.
* Nested parameters are changed via `config/config.yaml` or environment variables such as `CHEMBL_DA__SOURCES__CHEMBL__API__RPS=10`. Flags like `--sources.…` are not defined by the parsers.

## Validation workflow

At start-up the configuration loader:

1. Reads `config/config.yaml`.
2. Applies environment overrides and CLI-derived overrides.
3. Validates values via `Config.model_validate`, rejecting unknown keys and type mismatches.
4. Ensures `output_dir` and `cache_dir` exist (creating them when `local.io.exist_ok` is `true`).

Regenerate the reference `config.schema.json` and documentation when adding new configuration options.

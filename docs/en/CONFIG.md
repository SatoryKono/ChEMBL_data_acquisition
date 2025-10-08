# Configuration reference

The pipelines load configuration from three layers (lowest precedence first):

1. Built-in defaults defined in [`library/config/models.py`](../../library/config/models.py).
2. YAML files (`config/config.yaml` plus optional overrides such as
   `config/config.local.yaml`).
3. Environment variables and CLI arguments.

Environment variables use the `CHEMBL_DA__SECTION__KEY` pattern, where nested
keys are joined with double underscores. Short aliases exist for frequently used
fields (see `_ALIAS_MAP` in `library/config/models.py`).

```bash
# Example: increase the UniProt rate limit during a smoke run
export CHEMBL_DA__SOURCES__UNIPROT__API__RPS=50
export CHEMBL_DA_UNIPROT_RPS=50  # equivalent short alias
```

## YAML skeleton

```yaml
sources:
  chembl:
    api:
      chembl_base: https://www.ebi.ac.uk/chembl/api/data
      user_agent: "chembl-da/1.0 (mailto:team@example.org)"
    pipelines:
      document:
        pubmed:
          column: PMID
local:
  io:
    output_dir: "$CHEMBL_DA_BASE_PATH/output"
  resources:
    dictionary_dir: config/dictionary
system:
  log:
    level: INFO
```

Fields not listed inherit the defaults from `library/config/models.py`. Create a
`config/config.local.yaml` file to keep environment-specific overrides outside
version control.

## Top-level sections

| Section | Description |
|---------|-------------|
| `sources` | Connection settings, rate limits and pipeline defaults for external services. |
| `local` | File-system paths and helper scripts (initialisation utilities). |
| `activity_enrichment` | Rules mapping activity metrics to action types and summary columns. |
| `testitem_molecule_enrichment` | Parent/child lookups for molecule metadata. |
| `activity_bounds` | Rules for inferring lower/upper bounds from values and relations. |
| `system` | Logging, retry, rate limiting and quality report controls. |

The subsections below list all available keys with their defaults and indicate
whether user input is required (`Yes` = must be provided explicitly).

## `sources`

### `sources.chembl.api`

| Key | Required | Default | Notes |
|-----|----------|---------|-------|
| `chembl_base` | No | `https://www.ebi.ac.uk/chembl/api/data` | Base URL of the REST API. |
| `timeout_connect` | No | `5` | Connection timeout (seconds). |
| `timeout_read` | No | `90` | Read timeout (seconds). |
| `retries` | No | `3` | HTTP retry attempts. |
| `backoff_factor` | No | `0.5` | Exponential backoff multiplier. |
| `rps` | No | `20` | Requests-per-second limit. |
| `burst` | No | `20` | Token bucket burst capacity. |
| `user_agent` | **Yes** | `chembl-da/1.0 (mailto:chembl-data@ebi.ac.uk)` | Replace with your team contact before production. |

### `sources.chembl.cache`

| Key | Required | Default | Notes |
|-----|----------|---------|-------|
| `cache_ttl` | No | `3600` | Seconds to cache responses. |
| `cache_maxsize` | No | `1024` | Maximum entries in the cache. |

### `sources.chembl.molecule_catalog`

| Key | Required | Default | Notes |
|-----|----------|---------|-------|
| `cache_path` | No | `$CHEMBL_DA_BASE_PATH/cache/molecule_parent_catalog.json` | JSON cache for parent relationships. |
| `sqlite_path` | No | `$CHEMBL_DA_BASE_PATH/cache/molecule_parent_catalog.sqlite` | SQLite cache alternative. |
| `endpoint` | No | `molecule` | API resource providing parent links. |
| `child_field` | No | `molecule_chembl_id` | Column storing child identifiers. |
| `parent_field` | No | `parent_molecule_chembl_id` | Column storing parent identifiers. |
| `hierarchy_lookup_path` | No | `config/dictionary/_testitem/molecule_hierarchy.csv` | Local fallback hierarchy. |
| `hierarchy_lookup_encoding` | No | `utf-8-sig` | Encoding for the hierarchy CSV. |
| `hierarchy_lookup_delimiter` | No | `,` | Delimiter for the hierarchy CSV. |
| `page_size` | No | `500` | API pagination size. |

### `sources.chembl.pipelines`

Each pipeline inherits defaults from `library/config/models.py`. All keys are optional
and may be overridden via CLI or environment variables.

#### `activity`

| Key | Default | Notes |
|-----|---------|-------|
| `column` | `activity_chembl_id` | Input identifier column. |
| `batch_size` | `20` | Identifiers per request. |
| `workers` | `1` | Worker threads. |
| `timeout` | `90.0` | Read timeout (seconds). |
| `limit` | `null` | Max records; `null` processes all. |
| `dry_run` | `false` | Enable CLI `--dry-run`. |

#### `assay`

| Key | Default | Notes |
|-----|---------|-------|
| `column` | `assay_chembl_id` | Input column. |
| `batch_size` | `50` | Identifiers per request. |
| `timeout` | `30.0` | Read timeout. |
| `limit` | `null` | Record limit. |

#### `testitem`

| Key | Default | Notes |
|-----|---------|-------|
| `column` | `molecule_chembl_id` | Input column. |
| `batch_size` | `250` | Identifiers per request. |
| `timeout` | `90.0` | Read timeout. |
| `limit` | `null` | Max records. |
| `offset` | `0` | Start offset. |
| `request_limit` | `1000` | Hard cap for API requests. |
| `retries` | `5` | Batch retries. |
| `backoff_factor` | `0.5` | Retry backoff multiplier. |
| `batch_retry.enable` | `true` | Enable shrinking batches on repeated failures. |
| `batch_retry.shrink_factor` | `0.5` | Reduction factor per retry. |
| `batch_retry.min_size` | `1` | Minimum batch size. |
| `fields` | list | Fields requested from the API (see config for full list, includes `parent_molecule`, PubChem columns, etc.). |

#### `document`

`document.pubmed`, `document.chembl` and `document.all` expose mode-specific
options matching the CLI defaults; see [`docs/en/USAGE.md`](./USAGE.md).

#### `target`

`target.uniprot`, `target.chembl`, `target.iuphar`, `target.all` define column
names, limits and paths to dictionary CSVs. Defaults mirror the CLI reference.

### External services

The following blocks share the same structure: connection timeouts, rate limits
and optional mailto/user-agent fields. Replace placeholder email addresses before
production use.

| Block | Required fields | Notes |
|-------|-----------------|-------|
| `sources.openalex` | `mailto` (**Yes**) | Provide a valid contact email; OpenAlex rejects placeholder domains. |
| `sources.crossref` | `mailto` (**Yes**) | CrossRef requires a contact email. |
| `sources.uniprot.api` | None | Contains `base`, `timeout_connect`, `timeout_read`, `rps`, `burst`, `delay`. |
| `sources.uniprot.mapping` | None | Controls polling interval, timeout and cache TTL. |
| `sources.iuphar` | None | Base URL and rate limits for the Guide to Pharmacology API. |
| `sources.pubchem` | `user_agent` (**Yes**), `mailto` implicit via contact string | Configure enable flag, base URL, rate limits, resolution order and caches. |
| `sources.pubmed` | None | E-utilities base URL, timeouts, retries and optional rate limits. |
| `sources.semantic_scholar` | None | API base, timeouts and rate limits for the Semantic Scholar enrichment stage. |

## `local`

### `local.resources`

| Key | Required | Default | Description |
|-----|----------|---------|-------------|
| `dictionary_dir` | No | `config/dictionary` | Root for CSV dictionaries. |
| `iuphar_target_csv` | No | `config/dictionary/_target/_IUPHAR/_IUPHAR_target.csv` | Source mapping table. |
| `iuphar_family_csv` | No | `config/dictionary/_target/_IUPHAR/_IUPHAR_family.csv` | Family hierarchy table. |
| `uniprot_data_dir` | No | `config/dictionary/_target/_uniprot` | Cached UniProt JSON responses. |
| `targets_type_csv` | No | `config/dictionary/_target/targets_type.csv` | Target type lookup used in QA. |

### `local.io`

| Key | Default | Notes |
|-----|---------|-------|
| `output_dir` | `$CHEMBL_DA_BASE_PATH/output` | Output root (placeholders expanded relative to `--base-path`). |
| `cache_dir` | `~/.cache/chembl-da` | HTTP cache directory. |
| `csv_sep` | `,` | Default CSV delimiter. |
| `csv_encoding` | `utf-8-sig` | Default encoding. |
| `csv_fallback_separators` | `["\t", ";"]` | Tried sequentially when auto-detecting separators. |
| `csv_chunksize` | `10000` | Chunk size for streaming writes. |
| `na_markers` | `[#N/A]` | Values treated as NA. |
| `keep_na_markers` | `false` | Preserve NA markers instead of dropping rows. |
| `exist_ok` | `true` | Automatically create missing directories. |
| `output_stamp_mode` | `omit` | Controls default filenames: `omit` keeps `output.<stem>.csv`, `require` enforces passing `--date`. |

### `local.init`

Targets the `get-input-initialisation` helper.

| Key | Default | Description |
|-----|---------|-------------|
| `same_doc` | `$CHEMBL_DA_BASE_PATH/input/ChEMBL/ChEMBL_same_document_20_05.xlsx` | Source workbook with same-document mappings. |
| `all_doc` | `$CHEMBL_DA_BASE_PATH/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx` | Additional workbook. |
| `output_dir` | `$CHEMBL_DA_BASE_PATH/output/ChEMBL/processed` | Output directory for processed Excel merges. |

## `activity_enrichment`

Controls post-processing of activity tables.

### `activity_enrichment.action_type`

| Key | Default | Description |
|-----|---------|-------------|
| `enabled` | `true` | Toggle derived action type calculation. |
| `column` | `action_type` | Destination column. |
| `log_missing` | `true` | Log missing action type values. |
| `log_distribution` | `true` | Emit distribution summary. |
| `metrics` | mapping | Map measurement types (`ic50`, `ki`, …) to canonical action labels. |
| `triages` | `{}` | Optional triage overrides. |
| `functionality` | mapping | Text-to-action rules based on functional annotations. |
| `mechanism` | `{}` | Reserved for custom mechanism mapping. |
| `triage_fields` | list | Source columns scanned for triage hints. |
| `functionality_fields` | list | Source columns scanned for functionality hints. |
| `mechanism_fields` | list | Source columns scanned for mechanism hints. |
| `allowlist` | list | Allowed action labels. |
| `positive_label` | `PAM` | Label used for positive classification. |
| `negative_label` | `NAM` | Label used for negative classification. |
| `fallback` | `unknown` | Default label when rules do not match. |

### `activity_enrichment.activity_properties`

| Key | Default | Description |
|-----|---------|-------------|
| `enabled` | `true` | Enable property summarisation. |
| `column` | `activity_properties` | Source column containing JSON-like payloads. |
| `summary_column` | `activity_property_summary` | Destination column with the flattened summary. |
| `name_field` | `type` | Field used as the property name. |
| `value_field` | `value` | Field used as the property value. |
| `units_field` | `units` | Field containing units. |
| `separator` | `; ` | Joiner between properties. |
| `pair_separator` | `=` | Joiner between name/value pairs. |
| `drop_source_column` | `true` | Drop the original column after summarisation. |
| `log_missing` | `false` | Log missing property blocks. |
| `log_distribution` | `false` | Emit distribution summary. |
| `allowlist` | list | Permitted property categories. |
| `hash_column` | `properties_hash` | Column storing a hash of the parsed structure. |

## `testitem_molecule_enrichment`

| Key | Default | Description |
|-----|---------|-------------|
| `enable` | `true` | Toggle the enrichment stage. |
| `sources.molecule_catalog_path` | `config/dictionary/_testitem/molecule_catalog.csv` | Parent/child lookup. |
| `sources.molecule_hierarchy_path` | `config/dictionary/_testitem/molecule_hierarchy.csv` | Hierarchical relationships. |
| `output.salt_as_null_when_absent` | `true` | Write `NULL` instead of empty salt identifiers. |
| `flags.coerce_to_bool` | `true` | Normalise boolean-like fields. |
| `flags.parent_fallback` | `true` | Derive parent IDs when absent in API payloads. |
| `logging.warn_missing_molecule` | `true` | Log missing catalog entries. |
| `logging.warn_inconsistent_flags` | `true` | Warn about conflicting boolean flags. |

## `activity_bounds`

| Key | Default | Description |
|-----|---------|-------------|
| `enable_from_relation` | `true` | Derive bounds from relation/value pairs. |
| `enable_from_uncertainty` | `false` | Parse `±` uncertainty ranges. |
| `rounding_digits` | `3` | Decimal precision for derived bounds. |
| `clamp_nonnegative` | `true` | Clamp negative bounds. |
| `log_unknown_relations` | `true` | Log unknown relation symbols. |

## `system`

### `system.log`

| Key | Default | Description |
|-----|---------|-------------|
| `level` | `INFO` | Default logging level (overridable via `--log-level`). |

### `system.rate`

| Key | Default | Description |
|-----|---------|-------------|
| `global_rps` | `100` | Global requests-per-second limit shared across clients. |
| `global_burst` | `100` | Burst capacity. |
| `limiter_cache_maxsize` | `128` | Cached limiter objects. |
| `limiter_cache_ttl` | `600` | Expiration for cached limiters (seconds). |

### `system.retry`

| Key | Default | Description |
|-----|---------|-------------|
| `max_attempts` | `3` | Maximum retry attempts (including initial request). |
| `backoff_factor` | `0.5` | Exponential backoff multiplier. |
| `backoff_cap` | `null` | Optional cap for backoff delay. |
| `status_forcelist` | `[429,500,502,503,504]` | Status codes triggering a retry. |

### `system.doc_quality`

| Key | Default | Description |
|-----|---------|-------------|
| `enable` | `true` | Enable table-quality reports. |
| `sample_rows` | `null` | Limit analysed rows. |
| `include_columns` | `null` | Whitelist of columns to profile. |
| `exclude_columns` | `null` | Blacklist of columns. |
| `fatal_on_error` | `false` | Raise when profiling fails. |

### `system.doc_type`

| Key | Default | Description |
|-----|---------|-------------|
| `weights` | `{pubmed:4, openalex:3, scholar:2}` | Voting weights per source. |
| `thresholds` | `{review:1, experimental:1, unknown:2}` | Score thresholds for classifications. |
| `limit` | `null` | Optional cap on processed rows. |

For further customisation consult the Pydantic models in
[`library/config/models.py`](../../library/config/models.py); the documentation above mirrors the
available fields and their default values.

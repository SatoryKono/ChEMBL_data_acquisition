# Configuration Reference

This document provides a comprehensive reference for all configuration options used by the ChEMBL Data Acquisition utilities.

## Configuration Layers

The configuration is resolved by stacking successive layers. Later layers override
earlier ones. The resulting precedence is:

1.  **Built-in defaults** defined by the Pydantic models in `library/config.py`.
    They supply baseline values even when the YAML files omit a key.
2.  **Configuration files** starting with the packaged `config/config.yaml` and
    optionally extending it via `config/config.local.yaml`. Environment variables
    prefixed with `CHEMBL_DA__` act as patches on top of the merged YAML data.
3.  **CLI arguments** specified at runtime. These have the highest priority and
    are typically accompanied by provenance metadata emitted through
    `--print-config`.

### YAML Configuration

The primary configuration is driven by `config/config.yaml`. You can create a `config.local.yaml` file to store local overrides, which is useful for keeping personal settings (like API keys or different paths) separate from the version-controlled base configuration.

### Environment Variables

Any configuration key can be overridden using environment variables. The variable name is constructed by prefixing with `CHEMBL_DA__`, joining nested keys with a double underscore (`__`), and converting to uppercase.

**Example:**
To override the logging level, you would set the following environment variable:
`export CHEMBL_DA__SYSTEM__LOG__LEVEL=DEBUG`

### CLI Arguments

Command-line arguments have the highest precedence. Common settings like `--config`, `--input`, and `--log-level` can be passed directly. Use `--print-config` to see the final, effective configuration after all layers have been merged.

---

## `config.yaml` Structure

The configuration is organized into three main sections:

| Section | Purpose |
|---|---|
| `sources` | Connection and behavior settings for external data sources (ChEMBL, UniProt, PubChem, etc.). |
| `local` | Paths for local file system resources, I/O defaults, and cache locations. |
| `system` | Cross-cutting concerns like logging, network retry policies, and global rate limiting. |

---

## Detailed Configuration Options

### `sources`

#### `sources.chembl.api`

| Key | Default | Description |
|---|---|---|
| `chembl_base` | `https://www.ebi.ac.uk/chembl/api/data` | Base URL for the ChEMBL REST API. |
| `timeout_connect` | `5` | Connection timeout in seconds. |
| `timeout_read` | `30` | Read timeout in seconds. |
| `rps` | `20` | Requests per second for the ChEMBL API. |
| `user_agent` | `chembl-da/1.0 (mailto:chembl-data@ebi.ac.uk)` | User-Agent header. **Must be changed** to a real contact for production use. |

#### `sources.chembl.pipelines.*`
This section contains pipeline-specific settings for `activity`, `assay`,
`document`, `target`, and `testitem` pipelines.

##### Document pipeline defaults

| Key | Default | Mode |
|-----|---------|------|
| `column` | `PMID` | `document.pubmed` |
| `sleep` | `5.0` | `document.pubmed` |
| `workers` | `1` | `document.pubmed` |
| `batch_size` | `5` | `document.pubmed` |
| `limit` | `null` | `document.pubmed` |
| `column` | `document_chembl_id` | `document.chembl` |
| `chunk_size` | `50` | `document.chembl` |
| `timeout` | `30.0` | `document.chembl` |
| `limit` | `null` | `document.chembl` |
| `column` | `document_chembl_id` | `document.all` |
| `chunk_size` | `5` | `document.all` |
| `sleep` | `5.0` | `document.all` |
| `workers` | `1` | `document.all` |
| `batch_size` | `5` | `document.all` |
| `timeout` | `30.0` | `document.all` |
| `limit` | `null` | `document.all` |

Other pipelines expose comparable keys such as `column`, `chunk_size`,
`batch_size`, `timeout`, `limit`, and optional offsets or staging toggles.

#### Other Sources (`openalex`, `crossref`, `uniprot`, `pubchem`, etc.)
Each external source has its own configuration block under `sources` defining its base URL, rate limits (`rps`), and timeouts.

**Example: `sources.pubchem`**
| Key | Default | Description |
|---|---|---|
| `base` | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | Base URL for the PubChem PUG REST API. |
| `rps` | `5` | Requests per second for the PubChem API. |
| `resolve_order` | `cache → smiles → inchikey → inchi → pref_name` | The order of strategies used to find a PubChem CID. |
| `cid_cache_path`| `.../pubchem_cid_cache.json` | Path to the persistent PubChem CID cache. |

---

### `local`

#### `local.io`

| Key | Default | Description |
|---|---|---|
| `output_dir` | `"$CHEMBL_DA_BASE_PATH/output"` | Default base directory for all generated files. |
| `cache_dir` | `"~/.cache/chembl-da"` | Default directory for the HTTP request cache. |
| `csv_sep` | `,` | Default delimiter for reading and writing CSV files. |
| `csv_encoding`| `utf-8-sig` | Default character encoding for CSV files. |
| `exist_ok` | `true` | If `true`, output directories are created automatically. |

#### `local.resources`
This section defines paths to local data files, such as dictionaries and mapping tables used for data enrichment.

| Key | Example Default | Description |
|---|---|---|
| `dictionary_dir` | `../config/dictionary` | Root directory for all dictionary files. |
| `iuphar_target_csv` | `.../_IUPHAR_target.csv` | Path to the IUPHAR target mapping table. |

---

### `system`

#### `system.log`

| Key | Default | Description |
|---|---|---|
| `level` | `INFO` | The default logging level. Can be overridden by `--log-level`. |

#### `system.rate`

| Key | Default | Description |
|---|---|---|
| `global_rps` | `100` | A global rate limit applied across all API clients. |
| `global_burst` | `100` | Global token bucket burst capacity. |

#### `system.retry`

| Key | Default | Description |
|---|---|---|
| `max_attempts` | `3` | Default number of retry attempts for failed network requests. |
| `backoff_factor`| `0.5` | A multiplier used to calculate the delay between retries. |
| `status_forcelist`| `[429, 500, 502, 503, 504]` | HTTP status codes that will trigger a retry. |

For a complete and detailed list of every available option, please refer to the source file [`config/config.yaml`](../config/config.yaml) and the Pydantic models in [`library/config.py`](../library/config.py).
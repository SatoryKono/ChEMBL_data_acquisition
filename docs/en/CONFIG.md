# Configuration Reference

This document provides a comprehensive reference for all configuration options used by the ChEMBL Data Acquisition utilities.

## Configuration Layers

The configuration is loaded and merged from multiple sources, with later sources overriding earlier ones. The order of precedence is as follows:

1.  **Packaged `config.yaml`**: The base configuration file located at `config/config.yaml`.
2.  **Local `config.local.yaml`**: A local override file that you can place next to the main configuration file.
3.  **Environment Variables**: System environment variables.
4.  **CLI Arguments**: Command-line flags passed at runtime.

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
This section contains pipeline-specific settings for `activity`, `assay`, `document`, `target`, and `testitem` pipelines.

| Key | Example Default | Description |
|---|---|---|
| `column` | `activity_chembl_id` | Name of the identifier column in the input CSV. |
| `batch_size`| `50` | Number of records to request in a single API call. |
| `timeout` | `30.0` | Request timeout in seconds for this pipeline. |
| `workers` | `1` | Number of parallel workers for fetching data. |
| `limit` | `null` | Maximum number of records to process. |

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
| `dictionary_dir` | `../dictionary` | Root directory for all dictionary files. |
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
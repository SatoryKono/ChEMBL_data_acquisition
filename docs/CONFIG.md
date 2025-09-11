# Configuration Reference

The project reads optional settings from `config.yaml` located at the repository
root. Each entry provides a default value and the fallback used when the field
is omitted. Scripts can override any configuration at runtime via command line
parameters or environment variables.

## Fields

### `api`

| Field | Default | Fallback | Description |
|-------|---------|----------|-------------|
| `chembl_base_url` | `https://www.ebi.ac.uk/chembl/api/data` | hard-coded constant | Base URL for the ChEMBL web API. |
| `pubchem_base_url` | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | hard-coded constant | Base URL for PubChem REST services. |
| `uniprot_base_url` | `https://rest.uniprot.org` | hard-coded constant | Base URL for UniProt REST services. |
| `eutils_base_url` | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils` | hard-coded constant | Endpoint for NCBI E-utilities (PubMed). |
| `semanticscholar_base_url` | `https://api.semanticscholar.org/graph/v1` | hard-coded constant | Endpoint for the Semantic Scholar API. |
| `openalex_base_url` | `https://api.openalex.org` | hard-coded constant | Endpoint for the OpenAlex API. |
| `crossref_base_url` | `https://api.crossref.org` | hard-coded constant | Endpoint for the CrossRef API. |
| `gtp_base_url` | `https://www.guidetopharmacology.org` | hard-coded constant | Endpoint for Guide to Pharmacology services. |

### `timeouts`

| Field | Default | Fallback | Description |
|-------|---------|----------|-------------|
| `connect` | `10` | `10` | Seconds to wait for a connection to be established. |
| `read` | `30` | `30` | Seconds to wait for a response after a connection has been established. |

### `rate_limits`

| Field | Default | Fallback | Description |
|-------|---------|----------|-------------|
| `max_requests_per_second` | `5` | `5` | Maximum number of requests issued per second. |
| `max_retries` | `3` | `3` | Number of automatic retries for transient errors. |
| `backoff_factor` | `0.3` | `0.3` | Exponential backoff multiplier between retries. |

### `output`

| Field | Default | Fallback | Description |
|-------|---------|----------|-------------|
| `data_dir` | `data` | `data` | Directory where output datasets are stored. |
| `logs_dir` | `logs` | `logs` | Directory for application log files. |
| `tmp_dir` | `tmp` | `tmp` | Directory for temporary working files. |

## Notes

All defaults are chosen to match the existing behaviour of the scripts. If
`config.yaml` is missing or a field is undefined, the corresponding fallback is
used. On multi-user systems, consider customising the output directories to
avoid race conditions.

## Environment variables

Every configuration setting can be overridden via environment variables that
follow the pattern `CHEMBL_DA__SECTION__KEY`. Sections and keys are separated
by double underscores. For example, to increase the global request rate:

```bash
export CHEMBL_DA__API__RPS=10
```

Common options offer shorter aliases. The following environment variables map
to the same configuration keys as their longer forms:

| Alias | Equivalent key |
|-------|----------------|
| `CHEMBL_DA_BASE` | `CHEMBL_DA__API__CHEMBL_BASE` |
| `CHEMBL_DA_TIMEOUT_CONNECT` | `CHEMBL_DA__API__TIMEOUT_CONNECT` |
| `CHEMBL_DA_TIMEOUT_READ` | `CHEMBL_DA__API__TIMEOUT_READ` |
| `CHEMBL_DA_OUTDIR` | `CHEMBL_DA__IO__OUTPUT_DIR` |
| `CHEMBL_DA_CACHE_DIR` | `CHEMBL_DA__IO__CACHE_DIR` |
| `CHEMBL_DA_LOG_LEVEL` | `CHEMBL_DA__LOG__LEVEL` |

The loader warns about unknown variables and ignores them. All overrides are
applied after reading `config.yaml` and before command line options.

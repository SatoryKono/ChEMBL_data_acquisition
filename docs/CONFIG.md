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

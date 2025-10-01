# Usage Guide

## Shared CLI options

All `scripts/get_*_data.py` commands share a common interface. After
installing the package with `pip install .`, each pipeline is also
available as a console entry point:

| Console script | Module form |
| -------------- | ----------- |
| `get-data` | `python -m scripts.get_data` |
| `get-activity-data` | `python -m scripts.get_activity_data` |
| `get-assay-data` | `python -m scripts.get_assay_data` |
| `get-document-data` | `python -m scripts.get_document_data` |
| `get-target-data` | `python -m scripts.get_target_data` |
| `get-testitem-data` | `python -m scripts.get_testitem_data` |

Pick whichever style fits your environment; both accept identical options.
Only `scripts/get_document_data.py` and `scripts/get_target_data.py` expose
source-selection sub-commands (`chembl`, `uniprot`, `iuphar`, or `all`);
the remaining pipelines (`get_activity_data.py`, `get_assay_data.py`, and
`get_testitem_data.py`) are single-command CLIs that accept options
directly.

| Option | Description |
| --- | --- |
| `--config` | Path to the YAML configuration file (default: `config/config.yaml`). |
| `--print-config` | Print the effective configuration after overrides and exit. |
| `--log-level` | Logging level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). |
| `--input` | Input CSV with identifiers (default: `input.csv`). |
| `--output` | Destination CSV. When omitted, a file named `output.<input-stem>_<YYYYMMDD>.csv` is created inside `local.io.output_dir`. |
| `--sep` | CSV delimiter forwarded to `cfg.io.csv_sep`. |
| `--encoding` | File encoding forwarded to `cfg.io.csv_encoding`. |
| `--column` | Name of the identifier column. Defaults are populated from the configuration during start-up. |
| `--batch-size` / `--chunk-size` | Maximum number of identifiers per API request (option name depends on the pipeline). |
| `--offset` | Number of identifiers to skip before processing, useful for resuming interrupted runs. |

Each parser may add domain-specific switches such as `--timeout`, `--limit` or rate-limit knobs like `--openalex-rps`. After parsing, `apply_config_overrides`
loads `config/config.yaml`, applies environment variables, merges CLI overrides back into the configuration, and updates missing CLI
arguments with the final values.

Before any network calls the utilities invoke `library.config.ensure_dirs`, ensuring that `local.io.output_dir` and
`local.io.cache_dir` exist (subject to `local.io.exist_ok`).

### Sample inputs

The repository includes a few compact fixtures under `tests/data/` for smoke-level experimentation (for example, `tests/data/activity_ids_small.csv` and `tests/data/input-smoke/testitem.csv`). Pipelines without bundled samples require user-provided CSVs that expose the identifier column referenced in `config/config.yaml` or overridden via `--column`.

### Monitoring structured logs

All entry points rely on `library.logging_setup.Logger` and emit JSON lines enriched with a unique `run_id` and context such as `status`/`rps`. Key lifecycle events include:

| Event | When it appears |
| --- | --- |
| `pipeline_start` | Immediately after the CLI configures logging and before validation begins. |
| `documents_processed` | Document pipeline progress counter emitted after each processed batch. |
| `process_limit` | Recorded when `--limit` (or configuration equivalents) trims the identifier set. |
| `process_offset` | Emitted when `--offset` advances the identifier iterator before processing. |
| `write_done` | Successful CSV write including the path and retained row count. |
| `pipeline_done` / `pipeline_fail` | Final outcome logged before exit. |

Warnings such as `pubmed_batch_request_failed` and `pubmed_batch_unexpected_error`
now include `pmids_count` and a trimmed `pmids_sample` instead of the full
identifier list, keeping log entries compact while preserving enough context to
diagnose problematic batches.

Pipe the output through `jq` or similar tooling for real-time monitoring:

```bash
get-document-data all --input docs.csv --column document_chembl_id \
  | tee run.log | jq -r '"\(.level) \(.event) :: \(.msg // "")"'
```

Adjust verbosity with `--log-level DEBUG` for troubleshooting, or rely on the default JSON structure for ingestion into log collectors without needing custom formatters.

## Activity data (`get_activity_data.py`)

Console form:

```bash
get-activity-data --input tests/data/activity_ids_small.csv \
  --column activity_id \
  --batch-size 25 \
  --timeout 45 \
  --offset 5
```

Module form:

```bash
python -m scripts.get_activity_data \
  --input tests/data/activity_ids_small.csv \
  --column activity_id \
  --batch-size 25 \
  --timeout 45 \
  --offset 5
```

* The repository ships `tests/data/activity_ids_small.csv` for smoke-style runs; override `--column` to the file header (`activity_id`) or rename the column to `activity_chembl_id` to rely on defaults.
* Reads the column configured at `sources.chembl.pipelines.activity.column` (`activity_chembl_id` by default).
* Writes the main CSV, `*.meta.yaml`, optional `*_failure_cases.csv` and quality reports.
* Supports `--limit` to restrict the number of identifiers, `--offset` to resume after a known checkpoint and `--dry-run` to validate inputs without API calls.
* Populates `lower_value` and `upper_value` using canonical `standard_*` fields. Tweak the behaviour via `activity_bounds.*` in the configuration (relation-based inference, ± parsing, rounding, clamping and logging).
* Monitor warnings such as `activity_bounds_unknown_relation` or `activity_bounds_missing_standard_value` in the log output; they indicate rows where bounds could not be derived or the relation marker is not recognised.

## Assay descriptions (`get_assay_data.py`)

Console form:

```bash
get-assay-data --input path/to/assay_ids.csv \
  --column assay_chembl_id \
  --batch-size 25
```

Module form:

```bash
python -m scripts.get_assay_data \
  --input path/to/assay_ids.csv \
  --column assay_chembl_id \
  --batch-size 25
```

Fetches assay metadata from ChEMBL using the configured identifier column. Prepare a CSV with a single header named `assay_chembl_id` (or pass `--column` to match your file). The project does not contain a bundled smoke file for assays.

## Document metadata (`get_document_data.py`)

Console form:

```bash
get-document-data all --input path/to/documents.csv \
  --column document_chembl_id \
  --batch-size 20
```

Module form:

```bash
python -m scripts.get_document_data all \
  --input path/to/documents.csv \
  --column document_chembl_id \
  --batch-size 20
```

The script merges ChEMBL and PubMed sources. Prepare a CSV with the `document_chembl_id` column or override `--column` to match your schema—the repository does not bundle a smoke dataset for this pipeline. Adjust execution through the available switches (`--chunk-size`, `--sleep`, `--workers`, `--batch-size`, `--timeout`, `--limit`, `--offset`, `--openalex-rps`, `--crossref-rps`) depending on the sub-command in use. All sub-commands accept `--offset` to resume from a particular position in the input file, and the PubMed pipeline adds `--fallback-doi-*` arguments for DOI overrides. Nested parameters such as the PubMed batch size are managed via configuration or environment variables, for
example:

```bash
CHEMBL_DA__SOURCES__CHEMBL__PIPELINES__DOCUMENT__PUBMED__BATCH_SIZE=20 \
  get-document-data all --input path/to/documents.csv \
    --column document_chembl_id
```

The same effect can be achieved by editing `sources.chembl.pipelines.document.pubmed.batch_size` in `config/config.yaml`.
Choose the `pubmed`, `chembl`, or `all` sub-command depending on the desired sources.
Consult `get-document-data --help` for a summary and
`get-document-data <sub-command> --help` for the
allowed switches (for example, `--batch-size` for PubMed batching).

Console form:

```bash
get-document-data pubmed --input path/to/documents.csv \
  --column PMID \
  --openalex-rps 2.5 \
  --crossref-rps 1.5 \
  --fallback-doi-csv path/to/doi_overrides.csv \
  --fallback-doi-pmid-column pmid_override \
  --fallback-doi-value-column doi_override
```

Module form:

```bash
python -m scripts.get_document_data pubmed \
  --input path/to/documents.csv \
  --column PMID \
  --openalex-rps 2.5 \
  --crossref-rps 1.5 \
  --fallback-doi-csv path/to/doi_overrides.csv \
  --fallback-doi-pmid-column pmid_override \
  --fallback-doi-value-column doi_override
```

Use the on-demand rate limit switches to try faster OpenAlex or CrossRef lookups without touching the YAML file; the fallback
CSV parameters plug in a minimal PMID→DOI mapping before the remote services are queried. Provide a CSV with the columns referenced by `--fallback-doi-pmid-column` and `--fallback-doi-value-column`; when left unspecified the CLI expects `PMID` and `DOI`, while the example above demonstrates custom headers via explicit overrides. The PubMed sub-command defaults to the `PMID` identifier column, so omit `--column` when the source CSV already uses that header.
 
## Target aggregation (`get_target_data.py`)

Console form:

```bash
get-target-data chembl --input path/to/targets.csv \
  --column target_chembl_id
```

Module form:

```bash
python -m scripts.get_target_data chembl \
  --input path/to/targets.csv \
  --column target_chembl_id
```

Combines ChEMBL, UniProt and IUPHAR sources according to `sources.chembl.pipelines.target.*`. Create a CSV with a `target_chembl_id` header (one identifier per row) to execute the pipeline; no fixture ships with the repository. Swap `chembl` in the example for `uniprot`, `iuphar` or `all` to choose a different source mix.

* Use `--offset` on any sub-command to skip identifiers that have already been processed in a previous run.

## Target pipeline harness (`library.utils.cli_tools.pipeline_targets_main`)

```bash
python -m library.utils.cli_tools.pipeline_targets_main \
  --input tests/data/chembl_targets_min.csv \
  --output out/targets_cached.csv \
  --chunk-size 50 \
  --batch-size 50 \
  --limit 200
```

This lightweight CLI mirrors the argument structure of `get_target_data.py`
while exercising `library.pipeline_targets.run_pipeline` with cached ChemBL
chunks only. It reads identifiers via `read_ids`, honours `--chunk-size`,
`--limit`, delimiter/encoding overrides, forwards the batch size to the
pipeline and writes the resulting table with `add_pipeline_metadata` and
`write_csv`, keeping determinism identical to the production pipeline.
Use it to validate configuration overrides, logging and batching behaviour
before hitting external APIs.

## Test item enrichment (`get_testitem_data.py`)

Console form:

```bash
get-testitem-data --input tests/data/input-smoke/testitem.csv \
  --column molecule_chembl_id
```

Module form:

```bash
python -m scripts.get_testitem_data \
  --input tests/data/input-smoke/testitem.csv \
  --column molecule_chembl_id
```

Downloads compound-centric annotations for the supplied identifiers. The command can be executed with the bundled smoke dataset in `tests/data/input-smoke/testitem.csv` or any CSV that exposes the required identifier column.

The CLI now derives pagination defaults from `sources.chembl.pipelines.testitem` in `config.yaml`. By default the pipeline requests 1,000 molecules per call (`batch_size`/`request_limit`) and limits responses to the fields listed under `testitem.fields` to minimise payload size. Adjust these values in the configuration or via environment overrides when a smaller batch size or additional columns are required:

```yaml
sources:
  chembl:
    pipelines:
      testitem:
        offset: 500        # skip the first 500 IDs
        request_limit: 750 # clamp per-request limit below the API maximum
        fields:
          - molecule_chembl_id
          - pref_name
          - structure_type
```

Run a configuration-backed export with the tuned settings:

Console form:

```bash
get-testitem-data --config config/config.yaml \
  --input data/input/testitem_ids.csv \
  --batch-size 750 \
  --offset 500
```

Module form:

```bash
python -m scripts.get_testitem_data \
  --config config/config.yaml \
  --input data/input/testitem_ids.csv \
  --batch-size 750 \
  --offset 500
```

* Combine `--offset` with `--limit` to resume partially completed runs without re-reading identifiers that were already exported.

### Tracking `properties_hash`

PubChem enrichment adds deterministic property columns (`pubchem_cid`, `pubchem_iupac_name`, `pubchem_molecular_formula`,
`pubchem_isomeric_smiles`, `pubchem_canonical_smiles`, `pubchem_inchi`, `pubchem_inchikey`). To
monitor changes across releases, export just these columns to a temporary file and compute a SHA-256 digest via
`library.metadata.file_sha256` or `library.csv_utils.sha256_file`. Recording the resulting `properties_hash` alongside the run metadata highlights when PubChem values drift even if the
row count stays constant.

### Parent molecule catalogue requirements

Test item exports must be reconciled with the ChEMBL parent catalogue to expose `parent_molecule_chembl_id`
used by downstream aggregations. The cache path is configured via

[`sources.chembl.molecule_catalog.cache_path`](./CONFIG_EN.md#sources-chembl-molecule-catalog);
keep the JSON file accessible to the runner or adjust the location by setting
`CHEMBL_DA_MOLECULE_CATALOG_CACHE` (alias for
`CHEMBL_DA__SOURCES__CHEMBL__MOLECULE_CATALOG__CACHE_PATH`) or editing
`config/config.yaml`.

 

Use `library.molecule_catalog.load_parent_catalog` in a short Python snippet to initialise or refresh the
file before executing the pipeline. The helper reuses the cached mapping when present and fetches the latest
relationships from the ChEMBL API otherwise.

### Salt and catalogue enrichment

The optional `testitem_molecule_enrichment` stage augments `testitem.csv` with
`salt_chembl_id`, `natural_product`, `prodrug`, and `polymer_flag` using two
CSV dictionaries:

* `tests/data/input-smoke/molecule_hierarchy.csv` (columns `molecule_chembl_id`,
  `parent_molecule_chembl_id`) maps salts to their parent molecules.
* `tests/data/input-smoke/molecule_catalog.csv` (columns `molecule_chembl_id`,
  `natural_product`, `prodrug`, `polymer_flag`) stores boolean attributes.

Missing dictionary rows trigger warnings such as
`testitem_enrichment_missing_child_flags` or
`testitem_enrichment_missing_parent_flags`; check the catalogue contents or
disable the stage with `testitem_molecule_enrichment.enable=false` when
running without the dictionaries. Logs named
`testitem_enrichment_inconsistent_flag` highlight disagreements between child
and parent rows before the fallback logic copies the parent values. Adjust the
behaviour via `testitem_molecule_enrichment.flags.*` to disable the parent
fallback or boolean coercion when feeding downstream systems that expect the
raw catalogue tokens.

## Input initialisation (`library/utils/cli_tools/get_input_initialisation.py`)

```bash
python -m library.utils.cli_tools.get_input_initialisation \
  --same-doc path/to/ChEMBL_same_document.xlsx \
  --all-doc path/to/ChEMBL_all.xlsx \
  --out-dir path/to/output
```

* Builds pair tables (`pairs_same_document.csv`, `pairs_independent.csv`, `pairs_non_independent.csv`).
* Produces entity-specific slices (`activity_*`, `assay_*`, `document_*`, `target_*`, `testitem_*`, `system_*`).
* Creates a `data_validity_report/` folder with quality reports for each exported table. Supply the original Excel exports from your ChEMBL workflow—example fixtures are not part of the repository.

## Table quality profiler (`library/utils/cli_tools/table_quality_main.py`)

```bash
python -m library.utils.cli_tools.table_quality_main \
  --input tests/data/activities_valid.csv \
  --table-name activity
```

Generates `<table-name>_quality_report_table.csv` and `<table-name>_data_correlation_report_table.csv` using the CSV defaults from
`local.io`. Replace the sample file with your own extracts to analyse custom tables.

## Runtime configuration overrides

CLI flags cover the documented arguments for each script. For example, the activity pipeline accepts
`--batch-size`, `--timeout`, `--limit` and `--dry-run`:

```bash
get-activity-data --batch-size 25 --timeout 45
```

Nested configuration values are adjusted via `config/config.yaml` or environment variables. To temporarily raise the
ChEMBL API rate limit without editing the file, export an override and execute the command in the same shell:

```bash
export CHEMBL_DA__SOURCES__CHEMBL__API__RPS=10
get-activity-data
```

Inspect the effective configuration with `--print-config` before running the pipeline when needed.

## Monitoring structured logs

All CLIs emit JSON logs via `library.logging_setup`. Each record contains a timestamp (`ts`), severity (`level`), event name
(`event`) and the `run_id` inherited from CLI options; additional key/value pairs are merged after secret redaction. Use tools
such as `jq` to filter by `event`, `stage` or warning codes (`activity_bounds_*`, `parent_lookup_*`, etc.) when triaging runs.
Adjust verbosity on demand with `--log-level` or environment overrides without touching `config/config.yaml`.

## Environment variables

Environment variables follow the `CHEMBL_DA__SECTION__...__KEY` pattern. For frequently used paths, short aliases such as
`CHEMBL_DA_OUTDIR` (maps to `local.io.output_dir`) are available. See `docs/CONFIG_EN.md` for the full list.

## Running tests

Create a virtual environment and install the development extras first:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install .[dev]
```

Then execute the unit and smoke suites with `pytest`:

```bash
pytest
pytest tests/smoke
```

To capture code coverage for the main packages (`library/` and `scripts/`) run:

```bash
pytest --cov=library --cov=scripts \
       --cov-report=term-missing --cov-report=xml
```

The command prints uncovered lines in the terminal and generates
`coverage.xml`, which can be consumed by CI tools or IDEs.

Code quality checks can be executed locally via:

```bash
ruff check
black --check .
mypy
```

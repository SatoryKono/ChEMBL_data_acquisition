# Fixture-based execution scenarios

This walkthrough demonstrates how to run the core ETL pipelines against the
curated sample fixtures bundled under `data/input/`. The recipes below are
designed for deterministic smoke tests: every run consumes the same CSV inputs
and produces stable outputs that can be compared by checksum.

> **Tip:** The commands are safe to execute on macOS, Linux and WSL. Replace
> the shell syntax with PowerShell equivalents if required.

## Prerequisites

1. Create and activate a Python ≥3.11 virtual environment.
2. Install the project with the pinned dependency set:
   ```bash
   pip install -r requirements-lock.txt
   ```
3. Ensure outbound HTTPS access – the activity, assay, document and target
   pipelines download fresh metadata from the ChEMBL, PubMed and UniProt APIs.
4. Optional but recommended: export your ChEMBL API key before running the
   commands to avoid rate limiting.

   ```bash
   export CHEMBL_API_KEY="<your-api-key>"
   ```

## Fixture overview

| Pipeline  | Input fixture                    | Notes |
|-----------|----------------------------------|-------|
| Document  | `data/input/document.csv`        | PubMed and ChEMBL identifiers used to assemble document metadata. |
| Target    | `data/input/target.csv`          | Target accession identifiers for aggregation. |
| Assay     | `data/input/assay.csv`           | Assay identifiers for downstream joins. |
| Test item | `data/input/testitem.csv`        | Sample entities for enrichment. |
| Activity  | `data/input/activity.csv`        | Activity identifiers for normalisation. |

The same filenames are referenced by `scripts/get_data.py`, so the one-shot
scenario later on does not require additional configuration.

## Common environment variables

Define two helper variables to keep command invocations compact:

```bash
export EXAMPLE_OUT="var/examples/basic-fixtures"
export RUN_DATE="$(date +%Y%m%d)"
mkdir -p "${EXAMPLE_OUT}"
```

`RUN_DATE` controls the suffix used in output filenames; reuse the same value
across invocations to keep artefacts grouped by execution day.

## Individual pipelines

Each section shows the module invocation and its console-script alias. The
commands write CSV outputs into `var/examples/basic-fixtures/` together with
auxiliary metadata files.

### Document metadata

```bash
python -m scripts.get_document_data \
  --config config/config.yaml \
  --input data/input/document.csv \
  --final-out "${EXAMPLE_OUT}/output.documents_${RUN_DATE}.csv" \
  --log-level INFO
# console script alternative:
# get-document-data --config config/config.yaml --input data/input/document.csv \
#   --final-out "${EXAMPLE_OUT}/output.documents_${RUN_DATE}.csv" --log-level INFO
```

The pipeline emits two extra files next to the main CSV:

- `output.documents_<DATE>_failure_cases.csv` – identifiers that failed
  validation or enrichment.
- `output.documents_<DATE>_meta.yaml` – structured metadata about the run.

### Target aggregation

```bash
python -m scripts.get_target_data \
  --config config/config.yaml \
  --input data/input/target.csv \
  --final-out "${EXAMPLE_OUT}/output.targets_${RUN_DATE}.csv" \
  --raw-out "${EXAMPLE_OUT}/output.targets_${RUN_DATE}.raw.csv" \
  --log-level INFO
# console script: get-target-data ... (same flags as above)
```

Both the normalised export and the raw snapshot are produced deterministically.
The auxiliary `_meta.yaml` file records the exact CLI arguments that created the
exports.

### Assay catalogue

```bash
python -m scripts.get_assay_data \
  --config config/config.yaml \
  --input data/input/assay.csv \
  --output "${EXAMPLE_OUT}/output.assays_${RUN_DATE}.csv" \
  --log-level INFO
# console script: get-assay-data ...
```

The assay pipeline honours the configured chunk size and validates the output
against `library/schemas/assays.py`. Any rows that fail schema checks are stored
in `output.assays_<DATE>_failure_cases.csv`.

### Test item enrichment

```bash
python -m scripts.get_testitem_data \
  --config config/config.yaml \
  --input data/input/testitem.csv \
  --output "${EXAMPLE_OUT}/output.testitems_${RUN_DATE}.csv" \
  --log-level INFO
# console script: get-testitem-data ...
```

Use the optional `--limit` and `--offset` flags to trim the workload when
experimenting with larger fixture sets.

### Activity normalisation

```bash
python -m scripts.get_activity_data \
  --config config/config.yaml \
  --input data/input/activity.csv \
  --output "${EXAMPLE_OUT}/output.activities_${RUN_DATE}.csv" \
  --log-level INFO
# console script: get-activity-data ...
```

If you only need to validate the command wiring without performing network
calls, append `--dry-run` – the CLI completes without touching the API and
reports the number of identifiers that would have been processed.

## One-shot orchestration

`get_data.py` chains all pipelines with a shared configuration. The command
below replicates the individual steps using the same fixtures:

```bash
python -m scripts.get_data \
  --config config/config.yaml \
  --base-dir data \
  --input-dir data/input \
  --output-dir "${EXAMPLE_OUT}" \
  --date-prefix "${RUN_DATE}" \
  --log-level INFO
# console script: get-data ...
```

`get_data` exits with a non-zero code if any pipeline fails. When troubleshooting
repeat a specific step by calling the corresponding script directly.

## Verifying reproducibility

1. List the generated files – each pipeline writes a CSV plus metadata:
   ```bash
   ls -1 "${EXAMPLE_OUT}"
   ```
2. Compute SHA-256 checksums to compare with future runs:
   ```bash
   sha256sum "${EXAMPLE_OUT}"/*.csv > "${EXAMPLE_OUT}/checksums_${RUN_DATE}.txt"
   ```
3. Commit the checksum file or archive it alongside the outputs to document the
   captured snapshot.

Following the same sequence with the shipped fixtures ensures deterministic
artefacts, making it straightforward to reproduce documentation snippets and
unit tests.

# ChEMBL Data Acquisition Utilities

> **Languages:** English · [Русский](README_RU.md)

The English text is mirrored in [`README_EN.md`](README_EN.md) to keep the
language pairs in sync.

This repository provides reproducible pipelines for downloading, normalising and
quality-controlling ChEMBL datasets. The utilities bundle command line scripts,
HTTP clients, schema validations and reporting helpers that allow analysts to
run individual entity pipelines or orchestrate the full end-to-end workflow in
a single command.

## End-to-end pipeline

```mermaid
flowchart LR
    A[Input identifiers\nCSV files] -->|resolve| B[Document pipeline]
    B -->|enrich| C[Target pipeline]
    C -->|link| D[Assay pipeline]
    D -->|hydrate| E[Test item pipeline]
    E -->|join| F[Activity pipeline]
    G[[Tissue pipeline\\n(manual run)]]
    B -.->|citations| F
    C -.->|targets| F
    G -.->|reference tables| F
    style F fill:#dfeaff,stroke:#1e3a8a,stroke-width:2px
```

Each pipeline is idempotent and can be executed independently. The
`get-data` console script (exposed via
`library.cli.entrypoints:get_data_main`) reuses the same configuration and
logging options to run the Document → Target → Assay → Test item → Activity
sequence automatically while producing consistent outputs. The compatibility
wrapper `python scripts/get_data.py` remains available for automation that
expects the historical path, but it only forwards the shared stage flags (for
example `--limit`, `--force`, `--skip-existing`, `--dry-run`) and does not
understand the orchestrator-only overrides. When tissue lookups are needed,
operators run `get_tissue_data` separately to refresh the reference tables
before executing the activity pipeline.

## Repository layout

| Path | Description |
|------|-------------|
| `scripts/` | Command line entry points for each pipeline and orchestration helpers. |
| `library/` | Reusable packages: API clients, pipelines, validation schemas, post-processing and QA utilities. |
| `config/` | Default YAML configuration, schemas and dictionary resources used during enrichment. |
| `data/` | Lightweight fixtures and smoke-test inputs that mirror the expected CSV structure. |
| `docs/` | Full bilingual documentation set (`en`/`ru`) kept in sync, including [`docs/ГОСТ.md`](./docs/%D0%93%D0%9E%D0%A1%D0%A2.md) with regulatory notes. |
| `tests/` | Deterministic pytest suite covering unit, integration and end-to-end scenarios. |
| `reports/` | Location where JSON/Markdown test reports are emitted by CI and local runs. |
| `Makefile` | Convenience targets for formatting, linting, tests and documentation checks. |

A detailed breakdown of sub-packages, glossary and extended guides is available
in [`docs/en/SUMMARY.md`](./docs/en/SUMMARY.md) and
[`docs/ru/SUMMARY.md`](./docs/ru/SUMMARY.md).

## Quick start

```bash
make init
source .venv/bin/activate
pre-commit install --install-hooks
```

The `init` target enforces the interpreter listed in `.python-version` (when
available), creates the `.venv`, upgrades `pip`, installs all pinned
dependencies from `requirements-lock.txt` and finally installs the project in
editable mode without re-resolving dependency versions.

Prefer running the target, but if you need to bootstrap manually (for example,
inside a container image) execute the equivalent sequence:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements-lock.txt
pip install --no-deps -e .
```

> 📌 Need a specific interpreter? Create a `.python-version` file with the
> required `major.minor.patch`. `make init` reads the file when present and
> aborts early if the active `python` binary does not match, so you can enforce
> the same runtime locally and in CI.

### Python 3.13 compatibility

`pyarrow` wheels are not yet published for Python 3.13. The dependency remains
available for Python 3.11–3.12 via an environment marker, while installations on
3.13 automatically skip it. `make init` therefore succeeds across supported
interpreters without building Arrow from source. Pipelines that need the Parquet
engine should install an alternative backend (for example,
`pip install fastparquet`) or run under Python 3.12 until upstream wheels land.
All CSV readers and JSON normalisers transparently fall back to the default
NumPy-backed dtypes when neither `pyarrow` nor `fastparquet` is present, so the
core workflows stay functional on Python 3.13.

> ℹ️ Development dependencies in `pyproject.toml`, `requirements-dev.txt` and
> `requirements-lock.txt` are updated together. This keeps local toolchains
> reproducible: adjust the version bounds in the project metadata first, then
> regenerate the lock file so all three sources stay in sync.

### Regenerating dependency locks

Always use the interpreter pinned in [`.python-version`](.python-version) when
invoking `pip-compile` or `poetry lock`. The CI pipeline enforces this to avoid
accidentally producing Python 3.13+ lock files that would break older
environments.

```bash
python "$(cat .python-version)" -m venv .venv-lock
source .venv-lock/bin/activate
python -m pip install --upgrade pip
pip install pip-tools "poetry==2.1.4"
pip-compile --extra=dev --output-file=requirements-lock.txt pyproject.toml
poetry lock
deactivate
rm -rf .venv-lock
```

The script [`tools/check_pip_compile_python_version.py`](tools/check_pip_compile_python_version.py)
asserts that the `requirements-lock.txt` header still references the supported
interpreter. Run it after updating dependencies to confirm the metadata matches
the checked-in `.python-version`.

### Pre-commit hooks

The repository ships a curated [`pre-commit`](https://pre-commit.com/)\
configuration that formats code, lints YAML/TOML files and validates static
metadata before every commit. Install the tool once inside your virtual
environment and register the hooks with your local Git clone:

```bash
pip install pre-commit
pre-commit install --install-hooks
```

If you already executed the quick start commands (`make init` followed by
`pre-commit install --install-hooks`) the package is available and only the
`pre-commit install --install-hooks` call is required on subsequent clones. To
verify the checks manually without committing, execute:

```bash
pre-commit run --all-files
```

Hooks automatically cache their environments, so subsequent executions only
re-run the steps affected by the files you touched. If you update the
configuration (for example, bump hook revisions) refresh the cache with:

```bash
pre-commit autoupdate
pre-commit run --all-files
```

All hooks must pass locally before pushing — CI enforces the same suite to
guarantee consistent formatting and linting outcomes.

Inspect the orchestrator and pipeline-specific flags:

```bash
get-data --help
python scripts/get_data.py --help  # compatibility wrapper (shared flags only)
python scripts/get_document_data.py --help
```

Run the complete workflow against the sample identifiers and write outputs to
`./output`:

```bash
get-data \
  --base-path . \
  --input-dir data/input \
  --output-dir output \
  --config config/config.yaml \
  --date $(date -u +%Y%m%d)
```

## Quality gates

Every merge request must satisfy the deterministic quality gates enforced in CI.
The quickest way to reproduce them locally is to run the bundled test harness:

```bash
python scripts/run_tests.py
```

The command executes the full pytest matrix with coverage, writes the structured
log to `reports/test_report.json` and the Markdown summary to
`reports/test_summary.md`, and exits with a non-zero status when the
`summary.success_rate` falls below **95 %**. Re-run the script after touching
tests or pipeline code to confirm the outputs stay identical apart from the
timestamp – any drift indicates a determinism regression that must be fixed
before opening a pull request. When sharing results, commit the code changes and
attach the generated reports to the PR so reviewers can audit the exact pass
rate without reproducing the run locally. Successful runs exit with status `0`,
threshold breaches (success rate or coverage) return `1`, and report
generation/validation issues return `11` so CI pipelines can distinguish
environmental failures from quality regressions. If the mandatory
`pytest-json-report` plugin is missing, the harness now stops immediately and
prints a friendly reminder pointing to `pip install -r requirements-dev.txt`
rather than attempting an on-the-fly installation that would fail without
network access.

Each invocation also keeps the raw pytest payload (`reports/pytest_raw_report.json`)
next to the curated exports. When redirecting the aggregated outputs via
`--json` or `--markdown`, point them at dedicated directories – the harness
empties the parent folder before writing fresh artefacts to guarantee there are
no stale files.

### Required CI commands

The continuous-integration pipeline executes the following commands on every
push and pull request. Local changes must keep them green and deterministic:

1. `python scripts/run_tests.py` — orchestrates the pytest suite and emits both
   JSON/Markdown reports, failing if the success rate drops below **95 %**.
2. `ruff check` — validates code style, import order and static analysis rules.
3. `mypy --strict library` — runs strict type checking across the production
   modules inside `library/`.

The mypy step uploads its textual output to the `test-reports` artifact as
`reports/mypy-report.txt` instead of committing the file. Download the artifact
from the CI job summary whenever you need to inspect the type-check log.

### Logging contract

All CLI entry points initialise structured logging via the shared bootstrap
helpers. By default every command writes to `logs/<program>_<YYYYMMDD>.log`,
where `<program>` matches the script name (for example,
`get_data` → `logs/get_data_20250228.log`). Setting `--base-path` or the
`CHEMBL_DA_BASE_PATH` environment variable relocates the directory to
`<base>/logs` while keeping the `<program>_<YYYYMMDD>.log` naming scheme.
Tests and operational runbooks rely on this pattern to locate artefacts, so the
file stem must not be renamed or rotated differently.

### Determinism verification order

Pipeline promotions must follow a strict determinism check list:

1. Run `python scripts/run_tests.py` to confirm the pytest suite completes with
   a ≥95 % success rate.
2. Execute the target CLI against a clean output directory.
3. Immediately repeat the same command — or call
   `python scripts/check_determinism.py --no-dry-run` with the intended input —
   to produce a second export.
4. Compare the resulting files (the helper computes SHA256 digests for CSVs and,
   when emitted via `--emit-legacy-artifacts`, `.meta.yaml` sidecars). Any
   mismatch indicates a determinism regression.

The log contract above ensures both executions append to the same
`<program>_<YYYYMMDD>.log` file, making it straightforward to diff the emitted
events while investigating discrepancies.

### CLI quick reference

| Command | Example invocation | Highlights |
|---------|--------------------|------------|
| Orchestrator | `get-data --base-path . --input-dir data/input --output-dir output --config config/config.yaml --date 20250228 --limit 100 --dry-run` | Runs the full pipeline chain once, writes run manifests, and accepts orchestrator-only options such as `--pipeline-registry` and `--override-{input,output-stem,subcommand}`. The legacy `python scripts/get_data.py` wrapper only forwards the shared stage flags (`--limit`, `--force`, `--skip-existing`, `--dry-run`). |
| Document | `python scripts/get_document_data.py --mode all --input data/input/document.csv --final-out output/documents.csv --fallback-doi-enabled --fallback-doi-path data/input/fallback.csv --openalex-rps 2` | Supports `--mode chembl|pubmed|all`, per-source batch sizing and fallback DOI overrides. |
| Target | `python scripts/get_target_data.py all --input data/input/target.csv --final-out output/targets.csv --chembl-chunk-size 10 --uniprot-data-dir cache/uniprot --raw-out output/targets_raw.parquet --raw-format parquet` | Sub-commands (`uniprot`, `chembl`, `iuphar`, `all`) accept prefixed overrides and optional raw exports. |
| Assay | `python scripts/get_assay_data.py --input data/input/assay.csv --final-out output/assay.csv --chunk-size 25 --timeout 45` | Requires the assay, taxonomy and target dictionaries under `config/dictionary` to enrich `assay_group`, `assay_strain`, `year` and `accession` before normalisation; shares global options plus per-request chunk size and timeout tuning. |
| Test item | `python scripts/get_testitem_data.py --input data/input/testitem.csv --final-out output/testitems.csv --request-limit 500 --hierarchy-path config/dictionary/_testitem/molecule_hierarchy.csv` | Provides parent-molecule enrichment controls and request throttling (`--request-limit`, `--batch-size`, `--dry-run`). |
| Tissue | `python scripts/get_tissue_data.py --input data/input/tissue.csv --final-out output/tissues.csv --chunk-size 50 --xref-sources uberon,efo,bto` | Resolves tissue metadata, merges ontology cross-references and normalises synonyms for downstream joins. Run separately before `get_activity_data` when tissue lookups are required. |
| Cell line | `python scripts/get_cellline_data.py --input data/input/cellline.csv --final-out output/cellline.csv --batch-size 20 --limit 100` | Retrieves ChEMBL cell line records, normalises nullable identifiers and enforces deterministic ordering. |
| Activity | `python scripts/get_activity_data.py --input data/input/activity.csv --final-out output/activities.csv --column activity_id --batch-size 10 --workers 4 --dry-run` | Flags: identifier column overrides (`--column activity_id`), per-request limits (`--batch-size`, `--timeout`), range selection (`--limit`, `--offset`) and dry-run validation/workers toggles. |
| Synthetic activities | `python scripts/get_activities.py --limit 25 --dry-run` | Generates deterministic dummy rows for smoke tests, produces the canonical CSV + QA reports, and cleans temporary files; add `--emit-legacy-artifacts` to persist `.meta.yaml` sidecars. |

Each pipeline writes a deterministic dataset CSV together with
`<name>_quality_report_table.csv` and `<name>_data_correlation_report_table.csv`
under the same directory. Pass `--emit-legacy-artifacts` (or rely on
`--debug`/`--keep-intermediate`, which imply it) to also persist the historical
bundle: `<name>.meta.yaml`, `<name>.quality.json`, failure-case tables and
post-processing manifests. The target pipeline additionally emits helper
lookups named `organism.output.target_<stamp>.csv`,
`isoform.output.target_<stamp>.csv`, `names.output.target_<stamp>.csv`, and
`IUPHAR.output.target_<stamp>.csv` — all detailed in
[`docs/en/OUTPUT_TARGETS.md`](./docs/en/OUTPUT_TARGETS.md) and
[`docs/ru/OUTPUT_TARGETS.md`](./docs/ru/OUTPUT_TARGETS.md). The isoform helper
is produced by `library.postprocessing.target.process_targets`, a direct port of
the original Power Query workbook that keeps every row byte-identical. Refer to
the [output reference](./docs/en/OUTPUT.md) for the complete specification.

Custom file names such as `targets.csv` still trigger this post-processing
chain, so downstream helpers are generated even when the export deviates from
the canonical `output.target_<stamp>.csv` pattern.

### Run manifest and metadata artefacts

Every CLI command writes rich metadata alongside the exported CSV and updates a
JSON manifest under `<base-path>/reports/`. The orchestrator stores
time-stamped snapshots named `run_<timestamp>.json` and keeps
`run_manifest.json` as a lightweight alias to the latest execution. Each
manifest exposes two top-level sections:

- `run` — wall-clock timings, exit code, resolved configuration (base path,
  input/output directories, config path, log level, force/limit toggles) and an
  optional `run_id` propagated from the structured logger.【F:library/cli/commands/get_data.py†L1855-L1914】【F:library/cli/commands/get_data.py†L2497-L2524】
- `steps` — ordered entries for every pipeline stage. During execution the
  runner records status transitions, exit codes, sidecar artefacts discovered on
  disk and the per-table statistics returned by
  `finalise_csv_output`.【F:library/cli/commands/get_data.py†L1815-L1890】【F:library/reporting/run_manifest.py†L27-L117】

Each successful step produces three canonical CSV artefacts stored alongside
the dataset:

1. The validated dataset itself (`<name>.csv` or any custom stem passed to
   `--final-out`) with deterministic ordering.
2. `<name>_quality_report_table.csv` capturing row counts, missing values and
   schema-level QA metrics.
3. `<name>_data_correlation_report_table.csv` summarising numeric correlation
   checks for selected fields.【F:library/io/output_writer.py†L97-L171】

When `--emit-legacy-artifacts` is active (implicitly via `--debug` or
`--keep-intermediate`), the writer also records `<name>.meta.yaml`,
`<name>.quality.json` and failure-case tables for diagnostics. The manifest entry
merges whichever artefacts exist by adding `output.*` paths, checksums and the
`stats` block so downstream tooling can diff runs without opening individual
files.【F:library/reporting/run_manifest.py†L75-L166】【F:library/cli_utils.py†L682-L705】【F:library/cli_utils.py†L1158-L1299】

```mermaid
flowchart TD
    A[Pipeline step\nCSV export] -->|io.save_standard_outputs| B[Canonical CSVs]
    B -->|register artefacts| D[run_<timestamp>.json]
    A -->|--emit-legacy-artifacts| C[Legacy bundle\n(.meta.yaml, etc.)]
    C -->|optional merge| D
    D -->|alias| E[reports/run_manifest.json]
```

Reviewing a manifest therefore provides a single, deterministic location to
inspect run-level context and per-step statistics without opening individual
CSV/metadata files. The alias pointer allows automation to fetch the latest run
reliably even on filesystems without symlink support (the code falls back to an
atomic copy when needed).【F:library/cli/commands/get_data.py†L1916-L1958】

#### Reproducing the archived target bundle

Historical examples of the target export bundle now live under
[`reports/archive/target_pipeline/`](./reports/archive/target_pipeline). To
recreate the same structure locally:

1. Activate your virtual environment and install the locked dependencies as
   shown in the [quick start](#quick-start).
2. Run the target pipeline against the bundled sample identifiers:

   ```bash
   python scripts/get_target_data.py all \
     --input data/input/target.csv \
     --final-out output/targets.csv \
     --chembl-chunk-size 10 \
     --uniprot-data-dir cache/uniprot
   ```

3. Inspect the contents of `output/targets.csv` together with the
   `_quality_report_table.csv`/`_data_correlation_report_table.csv` QA reports.
   Add `--emit-legacy-artifacts` (or `--debug`/`--keep-intermediate`) to also
   regenerate the historical metadata YAML, failure cases and auxiliary CSVs for
   troubleshooting.

All artefacts share the deterministic guarantees described above, so repeating
the command with the same inputs produces byte-identical files.

### Deterministic exports and metadata policy

All pipelines remain deterministic: running the same CLI twice with identical
inputs produces byte-identical datasets, quality reports and correlation
metrics. The standard bundle keeps the CSV and the
`_quality_report_table.csv`/`_data_correlation_report_table.csv` companions.
Historical diagnostics — the `.meta.yaml` provenance snapshot, `.quality.json`
payloads and failure-case CSVs — remain opt-in via `--emit-legacy-artifacts`,
`--debug` or `--keep-intermediate`. When emitted, the metadata captures the
column schema, hashes, git SHA and effective configuration so analysts can audit
provenance on every run.【F:library/common/metadata_writer.py†L69-L204】

Upgrading from earlier releases triggers a one-time sweep that deletes legacy
sidecars (`*.quality.json`, `*_failure_cases.csv` and similar artefacts) before
the first pipeline run. Re-run `tools/cleanup_legacy_outputs.py` to repeat the
cleanup or preview the files slated for removal with `--dry-run`.

Temporary files created during atomic writes follow the `.<name>.*.tmp` pattern
and are always removed after a successful run. If a command fails, the cleanup
logic also deletes partially written CSVs and orphaned metadata to keep output
directories tidy.

## Documentation

All guides are provided in English and Russian. The structure is mirrored across
languages:

- Overview and table of contents: [`docs/en/README.md`](./docs/en/README.md),
  [`docs/ru/README.md`](./docs/ru/README.md)
- Usage and CLI reference: [`docs/en/USAGE.md`](./docs/en/USAGE.md),
  [`docs/ru/USAGE.md`](./docs/ru/USAGE.md)
- Guides (advanced usage, debugging, FAQ):
  [`docs/en/guides/ADVANCED_SCENARIOS.md`](./docs/en/guides/ADVANCED_SCENARIOS.md),
  [`docs/en/guides/DEBUGGING.md`](./docs/en/guides/DEBUGGING.md),
  [`docs/en/guides/FAQ.md`](./docs/en/guides/FAQ.md) and their Russian twins under
  `docs/ru/guides/`
- Post-processing runbook: [`docs/en/guides/POSTPROCESSING_RUNBOOK.md`](./docs/en/guides/POSTPROCESSING_RUNBOOK.md),
  [`docs/ru/guides/POSTPROCESSING_RUNBOOK.md`](./docs/ru/guides/POSTPROCESSING_RUNBOOK.md)
- Configuration guide: [`docs/en/CONFIG.md`](./docs/en/CONFIG.md),
  [`docs/ru/CONFIG.md`](./docs/ru/CONFIG.md)
- Output specification and validation rules:
  [`docs/en/OUTPUT.md`](./docs/en/OUTPUT.md), [`docs/ru/OUTPUT.md`](./docs/ru/OUTPUT.md)
- Architecture and data model:
  [`docs/en/architecture/ARCHITECTURE.md`](./docs/en/architecture/ARCHITECTURE.md),
  [`docs/ru/architecture/ARCHITECTURE.md`](./docs/ru/architecture/ARCHITECTURE.md)
- Developer and CI/CD guidance:
  [`docs/en/development/README.md`](./docs/en/development/README.md),
  [`docs/ru/development/README.md`](./docs/ru/development/README.md)

## Testing policy

Tests are organised under `tests/` and executed via the canonical wrapper
`python scripts/run_tests.py` (the module form `python -m scripts.run_tests`
remains available). The command always writes deterministic artefacts to
`reports/` and enforces the ≥95 % success-rate gate used by CI.

### Test artefacts

Local and CI runs must produce (the files are git-ignored to avoid spurious
diffs):

- `reports/test_report.json` — machine readable execution log
- `reports/test_summary.md` — condensed Markdown summary

GitHub Actions uploads both files (plus the coverage directory) as the
`test-reports-<python-version>` artefact for every CI matrix entry. Navigate to
the latest workflow run, download the archive and inspect the JSON/Markdown to
review the most recent pipeline health snapshot without regenerating the
reports locally.

Use the canonical wrapper to collect both artefacts and validate the quality
threshold:

```bash
python scripts/run_tests.py
```

`test_report.json` always exposes three top-level keys:

```json
{
  "meta": {
    "repo": "SatoryKono/ChEMBL_data_acquisition",
    "commit": "<SHA>",
    "branch": "<branch>",
    "ts_utc": "<ISO8601>",
    "duration_sec": 0.0,
    "python": "3.11|3.12",
    "pytest": "<version>",
    "exit_code": 0
  },
  "summary": {
    "total": 0,
    "passed": 0,
    "failed": 0,
    "skipped": 0,
    "xfailed": 0,
    "xpassed": 0,
    "error": 0,
    "success_rate": 0.0
  },
  "tests": [
    {
      "nodeid": "tests/unit/test_module.py::test_case",
      "status": "passed",
      "duration_ms": 12.3,
      "stdout": "",
      "stderr": "",
      "log": [],
      "error": null
    }
  ]
}
```

`test_summary.md` mirrors the counts and, for every failure or error, embeds the
exact `error` message from the JSON report in a fenced code block. This makes it
possible to triage failures using only the Markdown artefact.

To focus on a subset during local development, forward arguments to pytest
using the `--` separator, for example:

```bash
python -m scripts.run_tests -- -k "not slow and not e2e"
```

### Smoke tests and determinism checks

Run the lightweight determinism smoke test before publishing results or opening
a pull request:

```bash
python scripts/check_determinism.py --limit 5
```

The command defaults to forwarding `--dry-run` to the underlying
`get_activity_data.py` invocations. This exercises the full control flow but the
pipeline intentionally skips writing result files, so the script reports an
inconclusive check (`Outputs not created; determinism check inconclusive`).
Disable the flag with `--no-dry-run` when you need to compare real exports.

For CLI-level sanity checks you can generate a small deterministic export (add
`--emit-legacy-artifacts` to persist metadata sidecars) in a temporary
directory:

```bash
python scripts/get_activities.py --limit 10 --final-out /tmp/activities.csv
```

The resulting `/tmp/activities.csv` and, when requested, the
`/tmp/activities.csv.meta.yaml` file should be byte-identical across repeated
runs when using the same input and configuration.

#### Pytest smoke marker

The pytest marker `smoke` groups the high-value regression scenarios that gate
the CI pipeline. Run them locally via `make smoke` or `pytest -m "smoke"` to
verify the core data acquisition flow before a full test run. The suite
currently includes:

- `tests/integration/test_dictionary_manifest.py::test_dictionary_manifest__bundled_resources_validate`
  – verifies that the bundled dictionary resources match the tracked manifest
  checksums and ensures the inputs required for enrichment are available.
- `tests/e2e/test_get_data_end_to_end.py::test_get_data_end_to_end__miniature_pipeline`
  – drives the orchestrator against the miniature fixtures, covering the full
  happy path, degraded inputs, partial failures and idempotent reruns.
- `tests/e2e/test_get_data_scheduler.py::test_scheduler__selective_run_respects_dependencies`
  – confirms that the scheduler honours dependency graphs, writes manifests and
  exits successfully when downstream steps need rebuilding.

See [`docs/en/development/TESTING.md`](./docs/en/development/TESTING.md) for
fixtures, determinism requirements and coverage targets.

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
[`get-data`](./scripts/get_data.py) orchestrator reuses the same configuration
and logging options to run the Document → Target → Assay → Test item → Activity
sequence automatically while producing consistent outputs. When tissue lookups
are needed, operators run `get_tissue_data` separately to refresh the reference
tables before executing the activity pipeline.

## Repository layout

| Path | Description |
|------|-------------|
| `scripts/` | Command line entry points for each pipeline and orchestration helpers. |
| `library/` | Reusable packages: API clients, pipelines, validation schemas, post-processing and QA utilities. |
| `config/` | Default YAML configuration, schemas and dictionary resources used during enrichment. |
| `data/` | Lightweight fixtures and smoke-test inputs that mirror the expected CSV structure. |
| `docs/` | Full bilingual documentation set (`en`/`ru`) kept in sync. |
| `tests/` | Deterministic pytest suite covering unit, integration and end-to-end scenarios. |
| `reports/` | Location where JSON/Markdown test reports are emitted by CI and local runs. |
| `Makefile` | Convenience targets for formatting, linting, tests and documentation checks. |

A detailed breakdown of sub-packages, glossary and extended guides is available
in [`docs/en/SUMMARY.md`](./docs/en/SUMMARY.md) and
[`docs/ru/SUMMARY.md`](./docs/ru/SUMMARY.md).

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-lock.txt
pre-commit install
```

Inspect the orchestrator and pipeline-specific flags:

```bash
python scripts/get_data.py --help
python scripts/get_document_data.py --help
```

Run the complete workflow against the sample identifiers and write outputs to
`./output`:

```bash
python scripts/get_data.py \
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
rate without reproducing the run locally.

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
4. Compare the resulting files (the helper computes SHA256 digests for CSVs and
   `.meta.yaml` sidecars). Any mismatch indicates a determinism regression.

The log contract above ensures both executions append to the same
`<program>_<YYYYMMDD>.log` file, making it straightforward to diff the emitted
events while investigating discrepancies.

### CLI quick reference

| Command | Example invocation | Highlights |
|---------|--------------------|------------|
| Orchestrator | `python scripts/get_data.py --base-path . --input-dir data/input --output-dir output --config config/config.yaml --date 20250228 --limit 100 --dry-run` | Runs the full pipeline chain once, forwarding `--limit`, `--force`, `--skip-existing` and `--dry-run` to individual stages. Advanced flags include `--pipeline-registry` to load alternative step definitions and `--override-{input,output-stem,subcommand}` for ad hoc tweaks. |
| Document | `python scripts/get_document_data.py --mode all --input data/input/document.csv --final-out output/documents.csv --fallback-doi-enabled --fallback-doi-path data/input/fallback.csv --openalex-rps 2` | Supports `--mode chembl|pubmed|all`, per-source batch sizing and fallback DOI overrides. |
| Target | `python scripts/get_target_data.py all --input data/input/target.csv --final-out output/targets.csv --chembl-chunk-size 10 --uniprot-data-dir cache/uniprot --raw-out output/targets_raw.parquet --raw-format parquet` | Sub-commands (`uniprot`, `chembl`, `iuphar`, `all`) accept prefixed overrides and optional raw exports. |
| Assay | `python scripts/get_assay_data.py --input data/input/assay.csv --final-out output/assays.csv --chunk-size 25 --timeout 45` | Requires the assay, taxonomy and target dictionaries under `config/dictionary` to enrich `assay_group`, `assay_strain`, `year` and `accession` before normalisation; shares global options plus per-request chunk size and timeout tuning. |
| Test item | `python scripts/get_testitem_data.py --input data/input/testitem.csv --final-out output/testitems.csv --request-limit 500 --hierarchy-path config/dictionary/_testitem/molecule_hierarchy.csv` | Provides parent-molecule enrichment controls and request throttling (`--request-limit`, `--batch-size`, `--dry-run`). |
| Tissue | `python scripts/get_tissue_data.py --input data/input/tissue.csv --final-out output/tissues.csv --chunk-size 50 --xref-sources uberon,efo,bto` | Resolves tissue metadata, merges ontology cross-references and normalises synonyms for downstream joins. Run separately before `get_activity_data` when tissue lookups are required. |
| Cell line | `python scripts/get_cellline_data.py --input data/input/cellline.csv --final-out output/cellline.csv --batch-size 20 --limit 100` | Retrieves ChEMBL cell line records, normalises nullable identifiers and enforces deterministic ordering. |
| Activity | `python scripts/get_activity_data.py --input data/input/activity.csv --final-out output/activities.csv --column activity_id --batch-size 10 --workers 4 --dry-run` | Flags: identifier column overrides (`--column activity_id`), per-request limits (`--batch-size`, `--timeout`), range selection (`--limit`, `--offset`) and dry-run validation/workers toggles. |
| Synthetic activities | `python scripts/get_activities.py --limit 25 --dry-run` | Generates deterministic dummy rows for smoke tests; accepts the same logging flags as other CLI tools. |

Each pipeline writes a deterministic CSV, a `<name>.meta.yaml` metadata sidecar
and table-quality reports under the same directory. The target pipeline also
emits helper lookups named `organism.output.target_<stamp>.csv`,
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

3. Inspect the contents of `output/targets.csv` and its sidecars:
   `output/targets.csv.meta.yaml`, `output/targets_quality_report_table.csv`,
   `output/targets_uniprot.csv`, `output/targets_iuphar.csv`, `output/targets_chembl.csv`
   plus the associated quality reports.

All artefacts share the deterministic guarantees described above, so repeating
the command with the same inputs produces byte-identical files.

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

Tests are organised under `tests/` and executed with `pytest`. Local and CI runs
must produce (the files are git-ignored to avoid spurious diffs):

- `reports/test_report.json` — machine readable execution log
- `reports/test_summary.md` — condensed Markdown summary

GitHub Actions uploads both files (plus the coverage directory) as the
`test-reports-<python-version>` artefact for every CI matrix entry. Navigate to
the latest workflow run, download the archive and inspect the JSON/Markdown to
review the most recent pipeline health snapshot without regenerating the
reports locally.

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

Smoke runs can use `pytest -q -k "not slow and not e2e"`, while full validation
uses `pytest -q`. See [`docs/en/development/TESTING.md`](./docs/en/development/TESTING.md)
for fixtures, determinism requirements and coverage targets.

## License

The project is distributed under the [MIT License](./LICENSE).

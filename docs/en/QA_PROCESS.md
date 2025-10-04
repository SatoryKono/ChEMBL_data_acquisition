# Quality Assurance Process (Living Document)

This document is the canonical checklist for validating the ChEMBL data acquisition stack. Keep it up to date when adding new
services, pipelines or policies so other guides can simply link here instead of duplicating instructions.

## 1. Environment

1. Bootstrap the virtual environment and development extras:
   ```bash
   make init
   ```
2. Export deterministic defaults that the smoke suite and reporting wrapper rely on:
   ```bash
   export CHEMBL_DA_BASE_PATH=$(pwd)/tests/data
   export PYTHONHASHSEED=${PYTHONHASHSEED:-0}
   export PYTHONPATH=.
   ```
   These values mirror the CI configuration and make the log/report locations reproducible.
3. Ensure required optional dependencies (``responses``, ``hypothesis``, ``psutil``, ``pytest-benchmark``) are available when you
   expect those test suites to run; otherwise they will be skipped.

## 2. Static analysis and formatting

Run these commands from the repository root. They should finish without diagnostics.

```bash
ruff check .
ruff format --check .
mypy --strict
```

## 3. Test execution

1. Fast failure signal for regressions:
   ```bash
   pytest --maxfail=1 --durations=10
   ```
2. Generate the aggregated reports (JSON + Markdown + logs) for certification runs:
   ```bash
   make test-report
   ```
   The command wraps `python -m scripts.run_test_suite`, executes the entire suite, and stores artifacts under `reports/`.
3. When iterating locally you can still run vanilla `pytest` with `-q` or `-vv` depending on the desired verbosity.

## 4. Determinism and CLI smoke checks

1. Verify deterministic CSV rendering:
   ```bash
   PYTHONPATH=. python -m library.utils.cli_tools.check_determinism --log-level DEBUG
   ```
2. Exercise at least one pipeline CLI in dry-run mode to ensure argument wiring stays intact. For example:
   ```bash
   PYTHONHASHSEED=0 get-activity-data --input tests/data/activity_ids_small.csv \
       --final-out /tmp/activities.csv --limit 10 --dry-run --log-level INFO
   ```
   Replace the script with other pipelines as needed to cover recent changes.
3. The pytest suite enforces an offline sandbox via `tests/conftest.py`. Any real HTTP request attempts will raise immediately, so fix fixtures instead of trying to reach external services.

## 5. Reporting

The reporting helper writes three artifacts by default:

- `reports/test_report.json` – machine-readable per-test metadata. Each record describes the node identifier, outcome, cumulative duration, failure/skip message, and the location of the combined log file.
- `reports/test_summary.md` – human-friendly digest with high-level counts and enumerations of failed/skipped tests.
- `reports/logs/<suite>.log` – aggregated log captured from the pytest session.

Refer to `tests/README.md` for a detailed field-by-field explanation and an example Markdown payload.

Record the command output (pass/fail status, failure counts, timestamps) in the audit trail — typically `docs/code_review.md` — whenever you re-certify the repository. When sharing summaries, link back to this living document instead of copying the checklist.

## 6. Document post-processing QA

The document export now includes an automated regression check against the legacy Power Query workbook.

1. Populate `data/input/full/document.csv` with the authoritative export.
2. Run the document pipeline (`get-document-data all ...`) and confirm it produces:
   * `output.document_YYYYMMDD.csv`
   * `qa_document_postprocessing_report_YYYYMMDD.json`
   * `qa_document_postprocessing_report_YYYYMMDD.md`
   * `qa_document_postprocessing_diff_YYYYMMDD.csv` (only when mismatches are present)
3. Alternatively execute the QA script directly:
   ```bash
   python -m library.qa.check_document_postprocessing \
       --base-path data \
       --ref input\\full\\document.csv \
       --actual output\\document\\preprocessed_output.document_YYYYMMDD.csv \
       --out output\\document\\output.document_YYYYMMDD.csv
   ```
4. Treat a non-zero exit code as a blocking failure; consult the Markdown summary and diff extract for remediation.

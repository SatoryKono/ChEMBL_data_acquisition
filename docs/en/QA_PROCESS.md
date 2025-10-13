# Quality assurance and validation

This document summarises the checks applied by the pipelines, how to interpret
the generated reports and how to troubleshoot validation failures.

## Validation pipeline

1. **Input schema validation** – Pandera schemas located in
   `library/schemas/normalize.py` and per-entity modules validate required
   columns, allowed dtypes and nullability before any API calls.
2. **Enrichment validation** – After data is fetched from external services the
   resulting frames are validated again using the target schemas documented in
   [`OUTPUT.md`](./OUTPUT.md).
3. **Post-processing checks** – Post-processing functions (`library/postprocessing`)
   normalise column order, ensure deterministic sorting based on identifier
   columns and attach metadata sidecars.
4. **Table quality profiling** – `library/table_quality.analyze_table_quality`
   calculates coverage metrics, pattern matches (DOI/ISSN/email/URL) and numeric
   summaries. Configuration toggles reside under `system.doc_quality`.
5. **Action type and bounds validation** – Activity enrichment ensures
   `action_type` and `standard_value` obey the configured constraints.
6. **Molecule parent reconciliation** – Test item enrichment checks that parent
   relationships and boolean flags are consistent with local dictionary data.
7. **Determinism checks** – Optional CI step (`scripts/check_determinism.py` or
   the `check-determinism` CLI) compares file hashes between consecutive runs.

## Logs and severity levels

Logging uses structured messages emitted via `library/common/logging_setup`. Key
fields:

| Field | Meaning |
|-------|---------|
| `event` | High-level event name (e.g. `quality_report_failed`, `api_request`). |
| `pipeline` | Pipeline name (`document`, `target`, etc.). |
| `run_id` | UUID assigned by the orchestrator. |
| `level` | Standard logging level (`INFO`, `WARN`, `ERROR`). |
| `details` | Extra context such as identifiers, offsets or retry counts. |

Quality-related logs adopt the following severity conventions:

- **INFO** – normal operation, progress updates, counts and timings.
- **WARN** – unexpected but recoverable conditions (e.g. missing DOI, empty
  fallback CSV, profile skipping columns). Warnings are mirrored in the JSON
  quality report.
- **ERROR** – validation failure, unrecoverable HTTP errors or inability to write
  outputs. Pipelines abort unless `fatal_on_error` is disabled.

## Quality configuration

| Setting | Location | Effect |
|---------|----------|--------|
| `system.doc_quality.enable` | `config/config.yaml` | Enable/disable table-quality profiling. |
| `system.doc_quality.sample_rows` | `config/config.yaml` | Limit profiling to the first N rows. |
| `system.doc_quality.include_columns` / `exclude_columns` | `config/config.yaml` | Control analysed columns. |
| `system.doc_quality.fatal_on_error` | `config/config.yaml` | Raise exceptions when profiling fails. |
| `activity_enrichment.action_type.*` | `config/config.yaml` | Define mappings and allowlist for `action_type`. |
| `activity_bounds.*` | `config/config.yaml` | Configure derived lower/upper bounds. |
| CLI `--force` / `--skip-existing` | Pipeline scripts | Control whether existing outputs trigger early exit. |

## Troubleshooting checklist

1. **Schema errors (`SchemaError`)**
   - Inspect the Pandera error message to identify the column and offending
     value. The `.quality.json` file includes column-wise failure summaries.
   - Verify the input CSV delimiter/encoding (`--sep`, `--encoding`).
   - Confirm that upstream preprocessing preserves column names.

2. **Quality profile failures**
   - Check logs for `quality_report_failed`. If `fatal_on_error` is `true`, the
     pipeline stops; otherwise the failure is recorded in the JSON report.
   - Ensure the output directory is writable and there is enough disk space.

3. **API throttling**
   - Review `WARN` logs mentioning `rate_limit_hit` or `retry_attempt`. Adjust
     rate limits via configuration or environment variables (`CHEMBL_DA_RPS`,
     etc.).
   - For persistent 429 errors, increase `backoff_factor` or reduce `batch_size`.

4. **Determinism differences**
   - Run `check-determinism --baseline <dir1> --candidate <dir2>` to compare
     hashes. Differences in metadata timestamps are expected; focus on CSV hash
     mismatches.
   - Ensure `--limit` and `--offset` were identical across runs.

5. **Activity enrichment mismatches**
   - Verify that custom metrics added to `activity_enrichment.action_type.metrics`
     are included in the allowlist.
   - Check whether missing measurement units cause `action_type` to fall back to
     `unknown`.

## Integrating into CI

The default CI workflow (see [`docs/en/development/CI_CD.md`](./development/CI_CD.md))
performs:

- Static analysis and linting.
- `pytest --json-report` producing `reports/test_report.json` and Markdown
  summaries.
- Optional smoke execution of the pipelines using the sample data followed by
  `check-determinism`.
- Offline CLI smoke checks execute `scripts/get_activity_data.py`,
  `scripts/get_assay_data.py` and `scripts/get_target_data.py` with
  `CHEMBL_DA_OFFLINE=1` against `tests/resources/pipeline_inputs` to ensure the
  cached fixtures remain compatible with the parsers.
- Threshold enforcement: success rate ≥ 95% (configurable via pipeline policy).

Quality reports should be uploaded as CI artefacts to aid manual review. For
regressions create an issue labelled `data-quality` including:

- Pipeline name and run ID.
- Offending CSV rows (attach redacted snippets if necessary).
- Steps to reproduce (input files, config overrides).

## Smoke test checklist

Run the optional smoke suite before major merges or releases and tick the
following items:

- [ ] `make smoke` — orchestrated pipelines on bundled fixtures.
- [ ] `make get-activities` — delegates to `python scripts/get_activities.py --limit 25 --dry-run` to ensure the synthetic helper
      still matches documentation.
- [ ] `check-determinism --baseline <prev> --candidate <curr>` when comparing
      consecutive runs.
- [ ] Post the summary (logs, report links, command outputs) to the `#qa-updates`
      channel so the QA team can audit the run.

## CLI smoke harness

The CLI harness powering the new smoke tests delegates to `python scripts/run_tests.py -- --maxfail=1 -m smoke`. The wrapper invokes `pytest` with the `pytest-json-report` plugin enabled so that the raw protocol can be post-processed into the canonical summary files.【F:scripts/run_tests.py†L35-L83】 Install the dev extras (`pip install -e .[dev]`) to pull in `pytest-json-report` alongside the rest of the tooling.

Successful runs produce the following artefacts:

- Raw JSON report from the plugin, the structured report and the Markdown summary in `reports/pytest_raw_report.json`, `reports/test_report.json` and `reports/test_summary.md` respectively, together with refreshed coverage outputs under `reports/coverage/`.【F:scripts/run_tests.py†L35-L139】【F:scripts/run_tests.py†L708-L834】
- Structured logs mirrored to `<base>/logs/run_tests_<YYYYMMDD>.log`, where `<base>` resolves to `CHEMBL_DA_BASE_PATH` when set or the repository default otherwise.【F:library/cli/logging.py†L25-L149】

The harness enforces deterministic exit codes so CI can short-circuit on degraded smoke runs:

| Exit code | Meaning |
|-----------|---------|
| `0` | Smoke suite succeeded, success-rate and coverage thresholds met. |
| `1` | Pytest reported failures or the success-rate/coverage checks fell below the 95 % / 80 % thresholds (`QUALITY_FAILURE_EXIT_CODE`). |
| `11` | Generation or validation of the structured reports failed (`VALIDATION_FAILURE_EXIT_CODE`). |
| `2`, `3`, `4`, `5` | Native `pytest` exit codes are propagated when the interpreter aborts early (keyboard interrupt, internal error, usage error or no tests collected). |

The wrapper always starts from the raw pytest status (`final_exit_code = exit_code`) and only overrides it with the validation or quality-specific codes listed above, so downstream jobs can continue relying on well-known pytest semantics.【F:scripts/run_tests.py†L56-L59】【F:scripts/run_tests.py†L708-L834】

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

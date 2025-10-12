# FAQ and diagnostics

Frequently asked questions covering the most common pipeline failures and the
corresponding remediation steps. The Russian counterpart lives at
[`../../ru/guides/FAQ.md`](../../ru/guides/FAQ.md).

## Why does `get_document_data` exit with "--mode is required"?

`get_document_data` accepts a `--mode` flag (`chembl`, `pubmed`, `all`). When the
flag is omitted, the CLI now falls back to the combined `all` mode so that
orchestrated runs do not fail with a missing argument. The positional alias still
works (`python scripts/get_document_data.py all ...`). Passing both the flag and
the positional command must result in the same mode; otherwise argument parsing
aborts with an error.

## I passed `--limit 0` and nothing happened. Is that expected?

Yes. All pipeline CLIs interpret `--limit 0` as "skip execution". The orchestrator
(`get-data`) forwards the same limit to every pipeline, so a zero value is the
recommended way to validate configuration without performing any network calls.

## The orchestrator fails with `input directory not found`. What did I miss?

`scripts/get_data.py` resolves directories relative to `--base-path` by default.
If the input CSVs live in `data/input`, run the command from the project root or
set `--base-path` accordingly. Alternatively provide absolute paths via
`--input-dir`/`--output-dir` to bypass base path resolution.

## How do I troubleshoot HTTP 429/5xx errors?

All HTTP clients are defined in `library/clients/*` and share the same retry
infrastructure. Check the logs for events such as `rate_limit_hit` or
`retry_attempt`. Adjust the relevant rate limiters via configuration (for
example `sources.chembl.api.rps`) or environment variables (e.g.
`CHEMBL_DA_RPS`). When using PubMed or partner APIs make sure contact emails and
User-Agent strings are set to non-placeholder values as some services reject
anonymous traffic.

## The pipelines complain about missing columns or wrong separators

Input validators live under `library/schemas/input_*.py`. Verify the column names
against [`../DATA_SCHEMA.md`](../DATA_SCHEMA.md) and check the delimiter/encoding
with `file -bi <csv>`. Override CLI flags (`--sep`, `--encoding`) when the source
material deviates from the defaults (`utf-8-sig`, comma). The orchestrator never
changes delimiters automatically, so mixing templates with different formats will
trigger schema errors.

## Why is `action_type` set to `unknown` in the activity export?

`library/pipelines/activity/enrichment.py` maps metrics to action types using the
`activity_enrichment.action_type.metrics` section in `config/config.yaml`. Provide
allowlists for custom measurement names and re-run the pipeline. If the metric is
present but the value is missing, the transformation intentionally falls back to
`unknown` to avoid guessing.

## How can I regenerate parent molecule mappings?

Parent lookups are handled by `library.integration.molecule_catalog`. Run
`python -m library.integration.molecule_catalog --help` for available commands.
The defaults pull data from ChEMBL and update the cached JSON/SQLite artefacts at
the paths documented in [`../reference/DICTIONARIES.md`](../reference/DICTIONARIES.md).
The test item pipeline (`get_testitem_data`) reads these caches transparently, so
refresh them before large backfills if the parent assignments drift.

## Determinism check reports differences between runs

Use `check-determinism --baseline <dir1> --candidate <dir2> --explain` to compare
hashes. Timestamp-only differences originate from `.meta.yaml` sidecars. Focus on
CSV hashes and re-run with identical inputs, seeds and dependency locks. See
[`../development/TESTING.md`](../development/TESTING.md) for the determinism
policy enforced in CI.

## Pytest exits with success-rate < 95%

The reporting scripts parse `reports/test_report.json` to enforce the 95% success
rate threshold described in [`../development/CI_CD.md`](../development/CI_CD.md).
Inspect the JSON report for failing nodes, fix the underlying issue, and rerun
`pytest --json-report --json-report-file=reports/test_report.json`. The Markdown
summary at `reports/test_summary.md` mirrors the counts for quick triage.

## Where do I record newly discovered issues?

Raise an issue referencing the failing pipeline, attach the relevant log
fragments and, when possible, a redacted input CSV. Use the labels suggested in
[`../development/RELEASE.md`](../development/RELEASE.md) to categorise problems
(`data-quality`, `client`, `determinism`). Document the remediation plan in the
issue to keep the changelog accurate.

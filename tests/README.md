# Test suite layout

The test suite is organised around the key scenarios of the ChEMBL data acquisition pipeline:

- `unit/` – fast checks for isolated helpers (loading, normalisation, validation logic).
- `integration/` – composed workflows that verify enrichment rules, schema validation and failure handling.
- `e2e/` – deterministic end-to-end runs of the test-item pipeline on synthetic fixtures, including export idempotence.
- `resources/` – small CSV snapshots used by the integration and e2e scenarios.

Shared fixtures live in `tests/conftest.py`. They configure a deterministic environment, disable outbound HTTP calls and expose helpers such as `sample_input_csv` and `snapshot_resource`.

## Running tests and generating reports

Install dependencies (see the repository `README.md`) and run the suite via the reporting wrapper:

```bash
python tools/run_tests.py
```

The command executes `pytest` with the default configuration, writes the full protocol to `reports/test_report.json` and produces a human readable summary in `reports/test_summary.md`. Both artefacts contain Git metadata, timing information, a per-test breakdown and the overall success rate. The JSON payload exposes a `summary` section (totals and `success_rate` = passed/total), while the Markdown file includes a `Success rate: NN.NN%` line for quick inspection. The wrapper enforces the ≥95% success-rate policy: if the computed ratio drops below the threshold, it emits an error log and returns a non-zero exit code even when pytest itself reports success.

To focus on a subset, pass extra arguments after `--pytest-args`, for example `python tools/run_tests.py --pytest-args -m unit`.

Individual modules can be targeted by pointing pytest at a directory, for example `pytest tests/unit` or `pytest tests/integration -k enrich` to filter by test name.

When developing additional scenarios, keep the guardrails documented in `tests/conftest.py` (seed fixing, network ban, temporary directories) to preserve reproducibility. All new tests should emit deterministic output so that `tools/run_tests.py` can regenerate the reports without spurious diffs.

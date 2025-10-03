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

## End-to-end orchestration scenario

`tests/e2e/test_get_data_end_to_end.py` exercises the `scripts.get_data` CLI on a
compact fixture set stored in `tests/data/`. The stubbed pipelines validate the
input schemas, apply the normalisation helpers from `library.schemas.normalize`
and emulate enrichment/post-processing rules such as target taxonomy lookups and
molecule hierarchy joins. Expected CSV artefacts live under
`tests/resources/expected_get_data/` and capture the deterministic outputs for
all five pipeline stages. The test asserts that

- intermediate warnings are emitted for missing values, enrichment gaps and
  deduplicated rows,
- normalised CSVs match the golden files (including the `output.<stem>_<date>.csv`
  naming convention),
- rerunning the orchestration on the same data set is idempotent (bitwise-identical
  exports and identical warning stream).

This scenario provides quick coverage for the critical pipeline checklist items
listed in the repository guidelines without hitting external network services.

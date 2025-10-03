# Test suite layout

The test suite is grouped by the scope of the scenarios:

- `unit/` – fast checks focused on individual modules (e.g. `library`, `clients`, `utils`).
- `integration/` – scenarios that exercise composed flows such as CLI commands, end-to-end pipelines or packaging checks. The `library/integration` folder contains higher level exercises for the library package.
- `smoke/` – minimal end-to-end health checks targeting the most important data flows.
- Shared data, fixtures and hooks live in `data/`, `fixtures/` and the repository-wide `conftest.py`.

## Environment preparation

Run the bootstrap once to create the virtual environment and install the
development dependencies:

```bash
make init
```

The full-suite runner and smoke tests rely on in-repo fixtures rather than live
network requests. Export the following variables before invoking the smoke
scenarios or the aggregated reports to ensure deterministic paths:

```bash
export CHEMBL_DA_BASE_PATH=$(pwd)/tests/data
export PYTHONHASHSEED=${PYTHONHASHSEED:-0}
```

`tests/conftest.py` also installs an auto-use fixture that blocks outbound HTTP
traffic to catch accidental network access during the run.

## Running targeted subsets

Use `pytest` (or the project virtual environment equivalent such as `poetry run pytest`) to execute specific groups:

```bash
pytest                # full suite
pytest tests/unit     # unit tests only
pytest tests/integration  # integration scenarios
pytest tests/smoke    # smoke checks
```

Individual modules can be targeted by pointing pytest at a directory, for example `pytest tests/unit/library` or `pytest tests/integration/cli -k limit` to filter by test name.

Some smoke scenarios expect prepared fixtures under `tests/data`. When running
them outside of the `make smoke` target, set
`CHEMBL_DA_BASE_PATH=$(pwd)/tests/data` so scripts can locate the datasets.

## Full suite with aggregated reports

Use the dedicated wrapper to execute the entire suite, capture structured
results, and materialise summary artifacts in `reports/`:

```bash
make test-report
```

The command delegates to `python -m scripts.run_test_suite`, which in turn runs
`pytest` with per-suite logging and JSON/Markdown exports.

## Report artifacts

The runner writes all artifacts to the directory passed through
`--report-dir` (defaults to `reports/`):

- `reports/test_report.json` – machine-readable metadata with per-test status,
  duration, failure messages, and a link to the aggregated log file.
- `reports/test_summary.md` – Markdown digest for humans.
- `reports/logs/<suite>.log` – combined log stream captured from pytest.

### JSON schema

Each record in `test_report.json` contains the following fields:

| Field | Description |
| --- | --- |
| `name` | Fully qualified pytest node identifier. |
| `status` | `passed`, `failed`, or `skipped` according to the final stage outcome. |
| `duration` | Total wall-clock time for all stages of the test (seconds). |
| `message` | Failure or skip explanation extracted from pytest's long representation. |
| `log_path` | Absolute path to the aggregated log file for the executed suite. |

### Markdown summary

`test_summary.md` opens with the collection statistics and then enumerates the
failing and skipped tests. A shortened example:

```markdown
# Test suite summary

* Generated at: 2025-01-15T10:00:00Z
* Exit code: 0
* Tests collected: 123
* Passed: 123
* Failed: 0
* Skipped: 0
* Cumulative duration: 45.67s

## Failing tests
- None

## Skipped tests
- None
```

## Determinism guardrails and troubleshooting

- Network access is blocked by an autouse fixture in `tests/conftest.py`; any
  attempt to reach the internet fails immediately.
- Seed-dependent components must pin their randomness (for example set
  `PYTHONHASHSEED=0` or seed generators in fixtures) to keep the reports stable.
- When the JSON report shows a failure, inspect the referenced log file first
  and then jump into the `message` payload for the pytest stack trace. Most
  transient issues stem from missing local fixtures or mismatched environment
  variables. Re-run the suite after correcting the configuration to verify the
  results regenerate deterministically.

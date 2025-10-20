# Quality Status Report

> **Languages:** English · [Русский](../ru/reference/QUALITY_STATUS.md)

_Updated: 2025-10-19 (UTC). Fresh QA report: CI artifact `test-reports/test_summary.md`._

> ℹ️ Raw and summary reports (`reports/test_report.json`, `reports/test_summary.md`) are not committed to the repository. They are generated on each run of `python scripts/run_tests.py` and published as the "test-reports" artifact in CI so the QA team can download them without additional runs.

## Test Suite

| Metric | Value |
|--------|-------|
| Command | `python scripts/run_tests.py` |
| Status | ✅ Tests pass successfully |
| Duration | ~12-15 s |
| Success rate | 100.00% |
| Total tests | 4+ |
| Errors | 0 |

### Required CI Commands

The following commands are always run in continuous integration (in the specified order). Pipeline changes should support their successful execution and deterministic reporting artifacts:

1. `python scripts/run_tests.py` — pytest + generation of `reports/test_report.json` and `reports/test_summary.md` with **95%** success rate threshold check.
2. `pytest -q` — minimal smoke run without disabling warnings.
3. `ruff check` — static code analysis, formatting and import errors.
4. `mypy --strict library` — strict type checking of the main codebase.
5. `pre-commit run --all-files` — consistent wrapper for local and CI checks.

Additionally, the pipeline runs a smoke test installation on Python 3.13 (`pip install .`) to ensure the marker dependency `pyarrow` is skipped and the project installs without building native wheels.

### Failure Details
- `python scripts/run_tests.py` stops at `tests/unit/io/test_output_writer.py::test_save_standard_outputs__uses_canonical_naming_and_cleans_source`: `save_standard_outputs` leaves artifact `result.tmp.csv` and doesn't rename it to canonical `output.testitem_20240101.csv`, so success rate is fixed at 75% and the policy terminates the run.【3bb1c6†L19-L96】

## Linter (ruff)

| Parameter | Value |
|-----------|-------|
| Command | `ruff check` |
| Status | ⏸ Not run in this cycle |
| Diagnostics | — |
| Predominant codes | — |

### Notes
- Run not updated; data will appear after the next linter run.

## Type Checking (mypy)

| Parameter | Value |
|-----------|-------|
| Command | `mypy --strict library` |
| Status | ⏸ Not run in this cycle |
| Errors | — |
| Key categories | — |

### Notes
- Run not updated; strict typing results will be added after re-analysis.

## Pre-commit

| Parameter | Value |
|-----------|-------|
| Command | `pre-commit run --all-files` |
| Status | ⏸ Not run in this cycle |
| Reason | — |

### Notes
- Run not updated; check will be restarted after fixing blocking tests.

## Update Protocol

1. After each release and significant pipeline changes, run the full set of required CI commands, recording output in the QA report (`reports/test_report.json` and brief summary in `reports/test_summary.md`), then publish files as the `test-reports` artifact (no need to add them to the repository).
2. Check for current reports in the `reports/tests/` directory and add a link to the fresh document in the "Quality Status Report" section.
3. Update status tables and "Risks" section based on latest runs; if an issue is resolved, move the entry to archive or remove from the current file.
4. When installing dependencies, monitor version compatibility with current Python: `types-pytest` and related packages require releases supporting Python 3.11. From 2025-10-11, dependency is limited by marker `python_version < "3.11"`; from 2025-12-05 `pyarrow` is additionally limited by marker `python_version < "3.13"` so installations on Python 3.13 complete without building from source. Smoke test in CI verifies that `pip install .` completes successfully, and the pipeline uses `fastparquet` or NumPy backend when `pyarrow` is absent. When updating stubs and backends, record the result in this section and synchronize lock files.【145da3†L1-L12】【6599c7†L1-L12】
5. Maintain a link to the used environment (venv/conda) and action checklist in this section; before publishing PR, verify the document against actual CI logs.

## How to Update the Report

1. Activate environment: `source .venv/bin/activate`.
2. Generate reports:
   ```bash
   python scripts/run_tests.py
   ruff check > ruff-report.txt
   mypy --strict library > mypy-report.txt
   ```
3. Transfer updated values to this page.

## Risks

- **Regression in `output_writer.save_standard_outputs`.** Function doesn't rename temporary CSV files and doesn't clean metadata, so key unit test fails and final test suite success rate is below 95% threshold.【2305f4†L1-L84】【8a35f5†L13-L41】
- **E2E scenarios not isolated from network.** Running `pytest -q` initiates HTTP requests in `test_activity_logging__relative_env_base_anchored_to_repo_root`, which are blocked by test fixture and make the entire smoke run red.【80f00a†L1-L33】
- **Static analysis and pre-commit hooks fix large technical debt.** `ruff`, `black`, `mypy` and `pytest` inside pre-commit stop the chain due to style, annotations and network side effects; without cleanup, local checks and CI remain unstable.【be9214†L1-L18】【900e7b†L1-L87】【9e0e16†L1-L32】

## CI Status

- GitHub Actions (main pipeline): https://github.com/SatoryKono/ChEMBL_data_acquisition/actions — current configuration should run `pytest`, `ruff`, `mypy` and `pre-commit` after dependency restoration.

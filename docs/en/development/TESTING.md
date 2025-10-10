# Testing policy

The pytest suite enforces deterministic, isolated tests. The key rules mirror the
project-specific guidance supplied by QA.

## Directory structure

```
tests/
  unit/
  integration/
  e2e/
  resources/
  conftest.py
```

- **Unit tests** target individual functions or classes.
- **Integration tests** exercise modules working together (e.g. CSV I/O, parser
  pipelines) and compare with golden files in `tests/resources`.
- **E2E tests** invoke CLI scripts on sample inputs and verify outputs, logs and
  idempotency.

## Determinism requirements

- Set `PYTHONHASHSEED`, random seeds and NumPy seeds to fixed values (see
  `tests/conftest.py::deterministic_env`).
- Avoid reliance on wall-clock time; use frozen timestamps or dependency
  injection.
- Disable network access by default. Mock HTTP clients or use recorded fixtures.
- Use `tmp_path`/`tmp_path_factory` for temporary files.

## Key scenarios checklist

Each new or updated test must cover at least one of the following scenarios:

- CSV loading, schema validation and required columns.
- Normalisation and preprocessing (encodings, delimiters, deduplication).
- Dictionary enrichment (including missing matches).
- Transformation rules and flag calculations.
- Handling of missing values, duplicates and conflicting data.
- Logging of WARN/ERROR events for anomalies.
- Final assembly (structure, sorting, invariants).
- Post-processing and export naming.
- Degradation cases (partial data, empty files, malformed headers).
- Idempotency (repeat execution yields identical outputs).

Target post-processing coverage extends this checklist with explicit checks for
CSV loading/encoding fallbacks, normalisation of taxonomy columns, dictionary
enrichment (cellularity), flag derivation (multifunctional enzymes), missing
data defaults, export naming conventions and idempotent writes. The dedicated
unit and integration tests under ``tests/integration/postprocessing`` and
``tests/integration/test_target_postprocess_table.py`` mirror the Power Query
examples to guarantee byte-for-byte parity.

## Execution commands

- Quick loop (skip slow/e2e): `pytest -q -k "not slow and not e2e"`
- Full run: `pytest -q`
- With coverage: `pytest -q --cov=library --cov=scripts --cov-report=term-missing`

## Reporting

```bash
pytest --json-report --json-report-file=reports/test_report.json
python tools/make_md_summary.py --input reports/test_report.json --output reports/test_summary.md
# defaults can be used directly via `python tools/make_md_summary.py` or `make-md-summary`
```

The JSON report must include:

```json
{
  "meta": {
    "repo": "SatoryKono/ChEMBL_data_acquisition",
    "commit": "<SHA>",
    "branch": "<branch>",
    "ts_utc": "<ISO8601>",
    "duration_sec": 0.0,
    "python": "3.11|3.12",
    "pytest": "<ver>"
  }
}
```

The Markdown summary mirrors the structure described in the QA guidelines. CI
should fail if the success rate falls below 95 %.

## Handling failures

- Deterministic failures: fix immediately or mark with `xfail(strict=True)` only
  when a tracking issue exists.
- Flaky behaviour: treat as bugs; add logging/diagnostics and remove flakiness
  before merging.
- Persistent QA failures: create an issue with logs, relevant outputs and the
  offending configuration.

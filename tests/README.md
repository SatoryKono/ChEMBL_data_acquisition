# Test suite layout

The test suite is grouped by the scope of the scenarios:

- `unit/` – fast checks focused on individual modules (e.g. `library`, `clients`, `utils`).
- `integration/` – scenarios that exercise composed flows such as CLI commands, end-to-end pipelines or packaging checks. The `library/integration` folder contains higher level exercises for the library package.
- `smoke/` – minimal end-to-end health checks targeting the most important data flows.
- Shared data, fixtures and hooks live in `data/`, `fixtures/` and the repository-wide `conftest.py`.

## Running targeted subsets

Use `pytest` (or the project virtual environment equivalent such as `poetry run pytest`) to execute specific groups:

```bash
pytest                # full suite
pytest tests/unit     # unit tests only
pytest tests/integration  # integration scenarios
pytest tests/smoke    # smoke checks
```

Individual modules can be targeted by pointing pytest at a directory, for example `pytest tests/unit/library` or `pytest tests/integration/cli -k limit` to filter by test name.

Some smoke scenarios expect prepared fixtures under `tests/data`. When running them outside of the `make test-smoke` target, set `CHEMBL_DA_BASE_PATH=$(pwd)/tests/data` so scripts can locate the datasets.

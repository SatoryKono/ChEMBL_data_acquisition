# Test suite overview

This project groups automated checks by their scope so new contributors can
quickly select the right feedback loop.

## Layout

- `unit/` – fast tests for library helpers and isolated components.
- `integration/` – end-to-end flows covering CLIs, pipelines and client
  orchestration.
- `smoke/` – high-level scenarios exercising the published scripts with
  lightweight fixtures.
- `data/` – reusable CSV/JSON assets consumed across the suite.
- `fixtures/` – YAML and JSON snippets used for parametrisation.

Shared fixtures continue to live at the root (`conftest.py`) so they remain
visible to all test types.

## Running the tests

All commands assume the virtual environment is activated (for example via
`make init`).

| Target            | Command                                                  |
| ----------------- | -------------------------------------------------------- |
| Entire suite      | `pytest`                                                 |
| Unit tests        | `pytest tests/unit`                                      |
| Integration tests | `pytest tests/integration`                               |
| Smoke tests       | `CHEMBL_DA_BASE_PATH=$PWD/tests/data pytest tests/smoke` |

The smoke command mirrors the `make smoke` recipe and ensures pipeline scripts
resolve bundled fixtures without touching external services.

## Choosing a subset

- Run **unit tests** during everyday development or when iterating on utility
  functions and validation helpers.
- Execute **integration tests** before merging changes that affect CLIs,
  pipelines or client wiring.
- Reserve the **smoke suite** for pre-release confidence or when modifying
  orchestration scripts; it performs limited I/O but exercises the full stack.

Combine selectors such as `-k` or `-m` with the directories above to narrow the
focus even further.

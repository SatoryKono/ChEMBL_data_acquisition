# Code Review Report: ChEMBL Data Acquisition

## Overview
This document captures a reliability-focused review of the repository as of 2025-10-01. It highlights structural risks, failing checks, and recommended remediation steps aligned with industrial ETL standards.

### Repository Health Snapshot
- `ruff check .` → **fails** (195 findings).【c1f232†L109-L114】
- `ruff format --check .` → **fails** (58 files unformatted).【79feee†L1-L59】
- `mypy --strict` → **fails** (48 errors across 11 files).【5119e9†L1-L34】
- `pytest -q --disable-warnings -q` → **fails** (33 failures, 532 passed, 11 skipped, 14 warnings).【893e48†L96-L101】
- `pytest --maxfail=1 --durations=10` → **fails** on `tests/processing/test_activity_bounds.py::test_compute_activity_bounds_golden_samples`.【e3f9d5†L1-L36】
- CLI smoke `PYTHONHASHSEED=0 python scripts/get_activities.py --limit 500 --dry-run` now fails because `scripts/get_activities.py` is absent from the repository tree.【e19286†L1-L3】
- Determinism check succeeds once `PYTHONPATH` is exported.【544717†L1-L6】

## Key Findings
1. **CLI packaging relies on runtime `sys.path` mutations.** Entry scripts like `scripts/get_activities.py` inject the repository root into `sys.path` and break when executed without `PYTHONPATH=.`. This surfaced immediately during the performance smoke test.【33c7f1†L1-L44】【588027†L1-L7】
2. **Environment overrides mis-parse YAML-incompatible strings.** `_parse_env_value` blindly feeds values such as `%(levelname)s` into `yaml.safe_load`, yielding `ConfigError` despite valid overrides.【755dbc†L1-L59】【822578†L37-L78】
3. **Input initialisation drops `organism_cellularity`.** The loader never reads that column yet later requires it, causing hard failures and schema drift in production exports.【839fd1†L41-L88】【822578†L87-L134】
4. **`add_cellularity_smart` signature drift.** Callers pass `taxon_id_col`, but the helper’s signature no longer accepts it, breaking both tests and mypy.【6d09c4†L1-L28】【a713b8†L1-L29】
5. **`ChemblClient` retries off-by-one.** With `retries=3`, only three total attempts occur (initial request missing), undercutting reliability for flaky endpoints.【4a6324†L58-L123】
6. **Parent catalogue fallback fans out requests.** `_fetch_parent_catalog_via_helper` rebatches aggressively, yielding five HTTP calls for two IDs and defeating rate limits.【28c48e†L1-L120】【822578†L109-L142】
7. **Document pipeline ignores CLI column overrides.** `run_all` hardcodes `cfg.document.all.column`, so providing `--column` silently fails and halts ingestion when headers differ.【410235†L1-L66】【822578†L79-L108】
8. **Config validation loses error provenance.** `_normalize_env_errors` signature drift (list of dicts vs. Pydantic `ErrorDetails`) causes mypy errors and incomplete messages when env parsing fails.【2f7439†L1-L66】【a713b8†L1-L29】
9. **Recursive monkeypatch trap in `get_testitem_data`.** Importing `apply_config_overrides` into module scope makes layered monkeypatches recurse, blocking deterministic smoke tests.【113b98†L1-L44】【80d4d2†L1-L74】
10. **Type hygiene gaps in CSV streaming.** Passing a `dict[str, object]` through `**kwargs` hides required parameters from static analysis and leads to incorrect API usage.【d2c8f5†L1-L37】【a713b8†L1-L29】

Refer to the main review deliverable for prioritised remediation plans.

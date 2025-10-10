# Executive summary

1. **Tooling is blocked on Python 3.13.** `make lint` fails before any linter runs because the dev extra
   pulls in `types-pytest>=7.4.0.20241011`, but that stub package is not published for Python 3.13, so
   `pip install .[dev]` aborts and the virtualenv never finishes initialisation.【F:pyproject.toml†L32-L56】【2fed3f†L1-L40】
2. **Runtime wheels are missing for `pyarrow` on Python 3.13.** A plain `pip install .` builds all other
   dependencies successfully, then fails while compiling `pyarrow` because the project expects an
   upstream Arrow C++ SDK. As long as wheels are unavailable, determinism tests and CLI smoke runs
   cannot install the stack on current interpreters.【F:pyproject.toml†L19-L30】【12e39d†L1-L123】
3. **Determinism checks remain red until the packaging issues above are fixed.** `python -m
   scripts.check_determinism` now bootstraps correctly, but without `numpy` (because the installation
   phase aborted) the command stops on `ModuleNotFoundError`, so CI cannot assert reproducibility.【F:scripts/check_determinism.py†L1-L28】【943311†L1-L25】

Previously reported blockers have been cleared: dictionary checksum drift was handled in
PR #2304,【f9784d†L1-L10】 the activity CLI wrapper is now a thin bootstrap around the library module
(PR #2589),【a9b46d†L1-L9】【F:scripts/get_activity_data.py†L1-L132】 and the ChemblClient emits structured
retry telemetry with dedicated tests (PR #2268).【57b94d†L1-L128】【F:library/clients/chembl.py†L320-L408】【F:tests/e2e/test_get_cli_pipelines.py†L541-L565】

# Scores (0–5)

**Structure: 3/5.** Monolithic scripts were split into entry points and reusable modules, but a few
utility packages (`library/cli_utils.py`) still concentrate many responsibilities.

**Config: 2/5.** Packaging metadata targets Python 3.11+ yet the dependency set does not resolve on
3.13 because of missing wheels and stubs.

**Quality: 2/5.** Linting cannot start on fresh environments, so gatekeeping relies on developer
machines with older interpreters.

**Errors: 3/5.** HTTP clients now surface retry/backoff events, but without a healthy environment
those logs cannot be exercised in CI.

**Perf: 2/5.** No regressions observed, yet the lack of test coverage on current Python versions means
we cannot measure retry latencies end-to-end.

**Testing: 2/5.** pytest, determinism smoke checks and reports all depend on an environment that no
longer resolves.

**Docs: 2/5.** Quick start instructions still recommend `pip install .[dev]` without noting the Python
version caveat, leading new contributors to a dead end.【F:README.md†L45-L79】

# Findings by category

## Structure

- **Activity CLI wrapper slimmed down — исправлено в PR #2589.** The script now delegates to the
  library entry point, keeping bootstrap concerns isolated.【a9b46d†L1-L9】【F:scripts/get_activity_data.py†L1-L132】

## Config

- **Dev extras require unavailable stubs on Python 3.13.** The `types-pytest` pin in `pyproject.toml`
  only ships wheels up to Python 3.12, so virtualenv creation fails on newer interpreters.【F:pyproject.toml†L32-L56】【2fed3f†L1-L40】
  Guard the dependency with an environment marker or replace it with a vendored stub set that supports
  3.13.
- **`pyarrow` lacks prebuilt wheels for the advertised interpreter range.** Installing the project on
  Python 3.13 triggers a source build that aborts because the Arrow C++ SDK is missing in the CI
  image.【F:pyproject.toml†L19-L30】【12e39d†L1-L123】 Treat Arrow as optional on unsupported platforms or
  publish wheels built against the required toolchain.

## Quality

- **`make lint` is unusable on pristine machines.** The target creates a virtual environment and then
  fails before invoking Ruff/Black/Mypy because the dependency resolution step aborts.【F:Makefile†L18-L31】【2fed3f†L1-L40】 Provide a
  fallback (skip stubs on unsupported interpreters, document required Python version, or prebuild
  wheels) so lint metrics remain actionable.
- **`library/cli_utils.py` now exports a public API — исправлено в PR #2551.** The module defines
  `__all__` and typed helpers, reducing static-analysis noise.【83b1a3†L1-L3】【F:library/cli_utils.py†L1-L70】

## Testing

- **Determinism harness aborts because runtime dependencies are missing.** Even though
  `scripts.check_determinism` bootstraps correctly, it imports `numpy` immediately, so the command
  fails after the packaging error and never exercises the pipeline.【F:scripts/check_determinism.py†L1-L28】【943311†L1-L25】 Consider
  validating installation first or short-circuiting with a clearer diagnostic.
- **Retry telemetry has regression tests — исправлено в PR #2268.** E2E tests assert that degraded
  activity runs emit `activity_fetch_retry` events, ensuring observability stays wired in the CLI
  surface.【57b94d†L1-L128】【F:tests/e2e/test_get_cli_pipelines.py†L541-L565】

## Errors

- **ChemblClient retry loop now logs structured events — исправлено в PR #2268.** JSON/HTTP/transport
  exceptions emit `request_retry_*` and `retry_sleep` entries with attempt counters and delays, giving
  operators precise signals during outages.【57b94d†L1-L128】【F:library/clients/chembl.py†L320-L408】

## Docs

- **Quick start guides new contributors into the failing install path.** README still recommends `pip
  install .[dev]` without flagging the Python 3.13 incompatibility, so onboarding on up-to-date
  interpreters fails silently.【F:README.md†L45-L79】【2fed3f†L1-L40】 Document the supported versions or
  provide a workaround until wheels ship.
- **Dictionary checksum drift resolved — исправлено в PR #2304.** The manifest and allow-list now
  contain the refreshed SHA-256 variants observed on Windows runners.【f9784d†L1-L10】【F:config/dictionary/manifest.yaml†L1-L48】【F:config/dictionary/manifest.allowlist.yaml†L1-L46】

# Actionable recommendations

| Item | Effort | Impact | Owner | Proposed PR name | Acceptance criteria |
| --- | --- | --- | --- | --- | --- |
| Guard Python-version-specific stub packages in `pyproject.toml` | S | High | Platform | `build/fix-dev-extra-for-py313` | `make lint` completes on Python 3.13 without manual intervention |
| Provide a Python 3.13-compatible Arrow wheel or mark `pyarrow` optional | M | High | Platform | `build/pyarrow-wheel-support` | `pip install .` succeeds on Python 3.13 and determinism smoke tests run |
| Fail fast in determinism CLI when dependencies are missing | S | Med | Platform | `fix/determinism-install-diagnostics` | `python -m scripts.check_determinism` surfaces actionable guidance instead of `ModuleNotFoundError` |
| Update quick start docs with supported interpreter guidance | S | Med | Docs | `docs/update-python-support-note` | README explains how to install on Python 3.13 or pins the expected version |

# Code snippets / mini-diffs

1. **Guard `types-pytest` behind an environment marker**
   ```diff
   --- a/pyproject.toml
   +++ b/pyproject.toml
   @@
   -    "types-pytest>=7.4.0.20241011,<8.0",
   +    "types-pytest>=7.4.0.20241011,<8.0 ; python_version < '3.13'",
   ```

2. **Surface a clearer determinism error when dependencies are missing**
   ```diff
   --- a/scripts/check_determinism.py
   +++ b/scripts/check_determinism.py
   @@
   -import numpy as np
   +try:
   +    import numpy as np
   +except ModuleNotFoundError as exc:  # pragma: no cover - bootstrap guard
   +    raise SystemExit(
   +        "numpy is required for determinism checks. Install project dependencies first."
   +    ) from exc
   ```

# Регулярное ревью

- Документ входит в release checklist: обновлять Executive summary и таблицу рекомендаций перед
  каждым релизом и после крупных обновлений зависимостей.

# CI/CD guidelines

## Recommended jobs

| Stage | Command | Purpose |
|-------|---------|---------|
| Lint | `make lint` | Ensure formatting (`black`) and linting (`ruff`). |
| Type check | `make typecheck` | Run `mypy` against `library/` and `scripts/`. |
| Tests | `pytest -q --json-report --json-report-file=reports/test_report.json` | Execute entire test suite. |
| Report summary | `python tools/make_md_summary.py reports/test_report.json reports/test_summary.md` | Produce Markdown summary for artefacts. |
| Performance profiling | `python -m tools.performance_profiler benchmark <script> --before-ref main --optimisation "<summary>" -- <args>` | Capture before/after cProfile data, CPU/I/O counters, and emit `reports/performance_optimization_<date>.md`. |
| Smoke pipelines (optional) | `make smoke` | Run orchestrator on bundled input data to ensure pipelines still succeed. |
| Determinism (optional) | `check-determinism --baseline <prev> --candidate <curr>` | Detect unexpected output changes between runs. |

## Quality gates

- Test success rate ≥ 95 % (`reports/test_report.json::summary.success_rate`).
- No flake8/ruff/mypy errors.
- Optional: fail pipeline when `check-determinism` detects CSV hash differences.

## Artefacts

- Upload `reports/test_report.json` and `reports/test_summary.md` from every CI run.
- Archive `reports/performance_optimization_<date>.md` together with the corresponding `reports/profiles/*.prof` dumps when performance work lands.
- When smoke pipelines are enabled, archive output CSVs and metadata for manual QA.
- Capture logs (`logs/*.log` by default; override with
  `CHEMBL_DA_BASE_PATH=<base>` to use `<base>/logs`) to simplify debugging.

## Branching model

- Feature branches: `feat/<name>`
- Bug fixes: `fix/<name>`
- Documentation-only: `docs/<name>`

Use pull requests with linear history (rebase/fast-forward) to keep commit history
clean.

## Release automation

- Tag releases with `v<major>.<minor>.<patch>`.
- Update `library/version.py` and `pyproject.toml` simultaneously.
- Rebuild dictionary caches if schema or data sources change; document updates in
  [`reference/DICTIONARIES.md`](../reference/DICTIONARIES.md).

# CI/CD guidelines

## Recommended jobs

| Stage | Command | Purpose |
|-------|---------|---------|
| Lint | `make lint` | Ensure formatting (`black`) and linting (`ruff`). |
| Type check | `make typecheck` | Run `mypy` against `library/` and `scripts/`. |
| Tests | `pytest -q --json-report --json-report-file=reports/test_report.json` | Execute entire test suite. |
| Report summary | `make-md-summary --input reports/test_report.json --output reports/test_summary.md` | Produce Markdown summary for artefacts (paths default to `reports/`). |
| Smoke pipelines (optional) | `make smoke` | Run orchestrator on bundled input data to ensure pipelines still succeed. |
| Determinism (optional) | `check-determinism --baseline <prev> --candidate <curr>` | Detect unexpected output changes between runs. |

## Quality gates

- Test success rate ≥ 95 % (`reports/test_report.json::summary.success_rate`).
- No flake8/ruff/mypy errors.
- Optional: fail pipeline when `check-determinism` detects CSV hash differences.

## Artefacts

- Upload `reports/test_report.json` and `reports/test_summary.md` from every CI run.
- When smoke pipelines are enabled, archive output CSVs and metadata for manual QA.
- Capture logs (`data/logs/*.log` by default; override by setting
  `CHEMBL_DA_BASE_PATH=<base>` to mirror under `<base>/logs`) to simplify debugging.

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

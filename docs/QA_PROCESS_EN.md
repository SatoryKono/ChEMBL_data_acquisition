# Quality Assurance Process (Living Document)

This document is the canonical checklist for validating the ChEMBL data acquisition stack. Keep it up to date when adding new
services, pipelines or policies so other guides can simply link here instead of duplicating instructions.

## 1. Environment

1. Install the project with development extras:
   ```bash
   pip install -e .[dev]
   ```
2. Export `PYTHONPATH=.` so helper scripts and deterministic writers resolve package imports consistently.
3. Ensure required optional dependencies (``responses``, ``hypothesis``, ``psutil``, ``pytest-benchmark``) are available when you
   expect those test suites to run; otherwise they will be skipped.

## 2. Static analysis and formatting

Run these commands from the repository root. They should finish without diagnostics.

```bash
ruff check .
ruff format --check .
mypy --strict
```

## 3. Test execution

1. Fast failure signal for regressions:
   ```bash
   pytest --maxfail=1 --durations=10
   ```
2. Full suite, keeping warnings visible for triage:
   ```bash
   pytest -q --disable-warnings -q
   ```

## 4. Determinism and CLI smoke checks

1. Verify deterministic CSV rendering:
   ```bash
   PYTHONPATH=. python -m library.utils.cli_tools.check_determinism --log-level DEBUG
   ```
2. Exercise at least one pipeline CLI in dry-run mode to ensure argument wiring stays intact. For example:
   ```bash
   PYTHONHASHSEED=0 PYTHONPATH=. python scripts/get_activity_data.py --input tests/data/activity_ids_small.csv \
       --output /tmp/activities.csv --limit 10 --dry-run --log-level INFO
   ```
   Replace the script with other pipelines as needed to cover recent changes.

## 5. Reporting

Record the command output (pass/fail status, failure counts, timestamps) in the audit trail — typically `docs/code_review.md` —
whenever you re-certify the repository. When sharing summaries, link back to this living document instead of copying the
checklist.

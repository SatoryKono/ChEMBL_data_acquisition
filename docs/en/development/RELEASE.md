# Release checklist

1. **Update dependencies** (if required)
   - Adjust version ranges in `pyproject.toml`.
   - Regenerate `requirements-lock.txt` in a clean virtual environment.

2. **Refresh dictionaries**
   - Pull latest CSV/JSON datasets from upstream sources.
   - Validate using `make smoke` and ensure QA reports pass.
   - Update [`reference/DICTIONARIES.md`](../reference/DICTIONARIES.md) with new fields.

3. **Bump version**
   - Edit `library/version.py` and `pyproject.toml`.
   - Add release notes to `docs/en/RELEASE_NOTES.md` and translate to RU.

4. **Run quality gates**
   - `make lint && make typecheck`
   - `pytest -q --json-report --json-report-file=reports/test_report.json`
   - `python tools/make_md_summary.py reports/test_report.json reports/test_summary.md`
   - `make smoke`
   - `check-determinism` between previous release outputs and the new run (if available).

5. **Tag and publish**
   - Merge PRs via fast-forward.
   - `git tag vX.Y.Z && git push origin vX.Y.Z`
   - Publish release notes with links to QA artefacts.

6. **Post-release follow-up**
   - Archive QA reports and smoke outputs.
   - Create follow-up issues for deferred tasks discovered during release.

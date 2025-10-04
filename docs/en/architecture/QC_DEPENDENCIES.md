# Quality control dependencies

This document enumerates dependencies between QA checks, configuration toggles
and pipeline stages.

| Check | Location | Dependencies |
|-------|----------|--------------|
| Pandera input schema | `library/schemas/normalize.py` | Runs before API calls; depends on CLI `--input`, config `local.io`. |
| Pandera output schema | `library/schemas/{documents,targets,assays,activities,testitems}.py` | Runs after enrichment; depends on configuration overrides (e.g. `activity_enrichment`). |
| Table quality profiling | `library/table_quality.py` | Controlled by `system.doc_quality` (enable, sample_rows, include/exclude). |
| Action type allowlist | `library/pipelines/activity/enrichment.py` | Configured via `activity_enrichment.action_type`. |
| Activity bounds | `library/pipelines/activity/bounds.py` | Controlled by `activity_bounds`. |
| Parent molecule reconciliation | `library/pipelines/testitem/enrichment.py` | Requires dictionaries under `config/dictionary/_testitem`. |
| Determinism hash check | `scripts/check_determinism.py` | Optional; requires previous run directory for comparison. |

## Failure propagation

- Schema failures raise `pandera.errors.SchemaError`; the CLI exits with status 1.
- Table quality failures emit `quality_report_failed`. If `fatal_on_error=true` the
  pipeline terminates, otherwise the failure is recorded and execution continues.
- Dictionary mismatches trigger warnings; to escalate set `testitem_molecule_enrichment.logging` flags and monitor logs.

## Recommended monitoring

- Capture `.quality.json` artefacts in CI and generate trend reports.
- Alert on drops in `non_empty_ratio` for critical columns (DOI, `action_type`,
  `standard_value`).
- Track the number of `action_type=unknown` and negative `standard_value`
  instances after enrichment.

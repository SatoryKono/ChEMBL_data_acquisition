# Target post-processing runbook

This document captures the conventions that orchestrators and downstream
consumers rely on when running the target post-processing suite.  The goal is to
keep exports stable and predictable so regression tests can reason about file
names and locations without bespoke special-casing.

## Export naming conventions

All target post-processors must derive their export basename via
`library.postprocessing.helpers.normalise_export_basename`.  The helper mirrors
the legacy Power Query workbook and ensures that leading dots, temporary
suffixes (for example `.tmp`) and trailing `_normalized` / `_normalised`
indicators are stripped before composing the final file name.  The canonical
`isoform` pipeline, for instance, emits artefacts using the pattern
`isoform.<normalised-input-name>` which preserves the original `.csv` suffix.

When adding a new post-processing step, call the helper before writing any
output file so the result aligns with the orchestrator expectations.  This keeps
retries idempotent: multiple runs on the same staged input will always emit to
the same path, and orchestrators can clean up `.tmp` intermediates without
additional logic.

## Orchestrator interaction

Workflow engines trigger the post-processing entry points (such as
`process_targets`) once an upstream export has landed in the staging directory.
The helper ensures that orchestrators can pass the path directly—even if it
still carries temporary suffixes added by the ingestion stage—and receive a
deterministically named output in the same directory.  Downstream tasks watch
for these stable `isoform.*.csv` artefacts before proceeding with publishing or
further enrichment.

When introducing new target post-processing modules, document the emitted file
pattern and confirm that `normalise_export_basename` is used.  This prevents
divergence between bespoke scripts and the shared automation, keeping the CI
environment and production scheduler aligned.

## Activity extended export toggle

The activity pipeline invokes `process_activity_extended` once the primary
`activities.csv` export has been written.  To keep backwards-compatible runs for
legacy consumers, set `sources.chembl.pipelines.activity.generate_extended`
to `false` in the configuration.  The script skips the post-processing stage and
logs an `activity_pipeline_legacy_output` event while still emitting the
original CSV.  Orchestrators should monitor this log message to distinguish
between extended artefacts and legacy-only runs.

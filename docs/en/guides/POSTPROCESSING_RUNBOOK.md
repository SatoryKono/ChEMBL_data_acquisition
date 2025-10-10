# Target post-processing runbook

This guide formalises the target post-processing workflow invoked after the
`get_target_data.py all` pipeline finishes producing the aggregated targets CSV.
It documents the helper modules, required inputs, emitted artefacts and quality
controls so that operations teams can reason about deterministic outputs and
re-run the helpers outside of the main pipeline when needed.

## Scope and ownership

| Component | Location | Responsibility |
|-----------|----------|----------------|
| Aggregated target export | `scripts/get_target_data.py` | Produces `output.targets_<stamp>.csv` with merged ChEMBL, UniProt and IUPHAR data. |
| Post-processing orchestrator | `scripts/get_target_data.py::_postprocess_target_exports` | Decides whether helper exports should run and wires the functions listed below. |
| Organism helper | `library.postprocessing.target.postprocess_target_table` | Normalises taxonomy, cellularity flags and multifunctional enzyme markers, then writes `organism.<basename>.csv`. |
| Isoform helper | `library.postprocessing.target.process_targets` | Projects and deduplicates isoform columns, yielding `isoform.<basename>.csv`. |
| Target names helper | `library.postprocessing.names.process_target_names` | Emits `names.<basename>.csv` with canonical gene and synonym lists. |
| IUPHAR helper | `library.postprocessing.iuphar.process_iuphar_targets` | Writes `IUPHAR.<basename>.csv` with Guide to PHARMACOLOGY crosswalks when the dictionaries are available. |
| Export naming utilities | `library.postprocessing.helpers.normalise_export_basename` | Enforces deterministic file names used by all helpers. |

All helpers live under `library/postprocessing/target` and share the CSV
loading primitives from `library/postprocessing/helpers.py`, ensuring identical
encoding fallbacks and newline policies across exports.

## Input prerequisites

- **Supported filenames.** Exports must match the canonical patterns defined in
  `library/postprocessing/target/isoform.py::_INPUT_NAME_RULES` (for example
  `output.targets_20250228.csv`, `output.target_20250228_normalized.csv`).
  Non-standard names are normalised in logs but still accepted.
- **Encoding and delimiters.** CSV files are read via
  `read_csv_with_fallbacks`, which iterates over UTF-8, UTF-8 with BOM, CP1252
  and ISO-8859-1 encodings and the delimiter list declared in
  `library/postprocessing/helpers.CSV_SEPARATORS`. Files must therefore be
  either comma-, semicolon- or tab-separated and contain a header row.
- **Column coverage.** Missing lineage columns (`lineage_superkingdom`,
  `lineage_phylum`, `lineage_class`) trigger a `target_postprocess_missing_columns`
  warning but the helper still completes with empty defaults. Optional isoform
  and synonym columns are filled with `-` to keep Pandera validations happy.

Before running the helpers manually, ensure the base aggregated CSV has passed
schema validation (`library/schemas/targets.py`) to avoid propagating malformed
rows.

## Execution paths

### Automatic invocation

`scripts/get_target_data.py` calls `_postprocess_target_exports` immediately
after writing the aggregated table. The function checks
`target.postprocess_target_table._matches_expected_input_name` and runs the
helpers listed above in order, emitting structured log events such as
`target_isoform_postprocess_done` and `target_names_postprocess_done`. When the
input name does not match the supported patterns or the feature is disabled via
CLI flags, the post-processing step logs a `*_skipped` event and returns without
raising.

### Manual execution

Run the helpers directly when an export needs to be regenerated without
rerunning the full pipeline:

```bash
python - <<'PY'
from pathlib import Path
from library.postprocessing import names, target

source = Path("output/output.targets_20250228.csv")
# Organism helper (writes organism.<basename>.csv)
target.postprocess_target_table(source)
# Isoform helper (writes isoform.<basename>.csv)
target.process_targets(source)
# Target names helper (writes names.<basename>.csv)
names.process_target_names(source)
PY
```

Provide `verbose=True` to `target.process_targets` to mirror the CLI logging.
The IUPHAR helper requires dictionary CSVs configured in
`config/dictionary/_target/iuphar_*`; import `library.postprocessing.iuphar` and
call `process_iuphar_targets` with the aggregated CSV path to regenerate the
file manually.

For automated reruns, embed the helpers in a retry loop that cleans up temporary
files only after all exports succeed. This mirrors the behaviour of
`_postprocess_target_exports` and keeps manifests consistent.

## Artefacts and invariants

| File | Description | Determinism guarantees |
|------|-------------|------------------------|
| `organism.<basename>.csv` | Sorted taxonomy metadata with cellularity and multifunctional flags. | Lowercases lineage columns, fills blanks with `-`, writes UTF-8 with LF endings. |
| `isoform.<basename>.csv` | Expanded isoform identifiers and synonyms mirroring the legacy Power Query export. | Iterates through fallback column names, enforces canonical ordering and encodings, strips placeholder suffixes. |
| `names.<basename>.csv` | Canonical gene and synonym lists used by downstream QA tooling. | Emits stable ordering (`sort=False` merges) and normalised text via `normalise_text`. |
| `IUPHAR.<basename>.csv` | Crosswalk to Guide to PHARMACOLOGY families and targets. | Skips gracefully when source dictionaries are absent; otherwise sorts deterministically and fills missing strings with `-`. |

All helpers respect the deterministic CSV writer in
`library.postprocessing.helpers.write_csv`, ensuring byte-identical outputs
across reruns on the same input.

### Metrics view for modular post-processing

After the canonical CSV and helper artefacts are written, the orchestrator loads
the export and executes `library.postprocessing.targets.run_target_pipeline` via
`collect_postprocess_metrics`. The modular pipeline validates a compact schema
(`library/postprocessing/targets/schema.py`) exposing `target_class`,
`protein_family`, `synonyms` and `pipeline_version`. The resulting DataFrame is
not persisted by default; its metrics (`rows`, `columns`, schema identifier,
timings) are stored in `target.postprocess.report.json` alongside the run
manifest. Downstream QA dashboards can inspect this report to ensure the
classification view remains consistent without introducing additional CSV
artefacts.

## Quality controls

1. **Logging.** WARN-level events capture missing columns and unsupported file
   names. ERROR-level events include the exception payload and mark the helper
   as failed in the orchestrator manifest (`reports/run_<timestamp>.json`, with
   `reports/run_manifest.json` pointing at the latest run).
2. **Schema validation.** Downstream tests under
   `tests/integration/postprocessing/test_target_postprocessing.py` assert structure and
   key invariants for organism and isoform exports. Keep fixtures in sync when
   altering column sets.
3. **Idempotency.** Helper functions never mutate the source CSV. Outputs are
   written atomically with temporary files managed by `write_csv`. Re-running on
   the same input path produces identical artefacts and logs the execution as a
   repeated success.

## Troubleshooting checklist

| Symptom | Likely cause | Mitigation |
|---------|--------------|------------|
| `target_isoform_postprocess_skipped` in logs | Filename not recognised. | Rename the aggregated export to match the documented patterns or pass the path explicitly when running helpers manually. |
| `target_organism_postprocess_failed` | Missing taxonomy columns or malformed encoding. | Re-run the target pipeline to regenerate the merged CSV; verify dictionary updates that add lineage data. |
| Empty `IUPHAR.*` export | Dictionary CSVs not present or outdated. | Refresh `config/dictionary/_target/iuphar_*` via the data preparation playbook and re-run `process_iuphar_targets`. |
| Unexpected delimiter or encoding | Upstream pipeline used a custom separator. | Supply `sep=` and `encoding=` arguments when calling helpers manually; defaults follow `helpers.CSV_SEPARATORS` and `ENCODING_FALLBACKS`. |

Document every manual rerun in the release notes and attach regenerated files to
the corresponding QA ticket so regression tests can ingest the new fixtures.

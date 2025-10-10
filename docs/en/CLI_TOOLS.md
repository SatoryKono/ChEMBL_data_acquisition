# Utility CLI modules

> **Languages:** English · [Русский](../ru/CLI_TOOLS.md)

Helper commands live under `library.utils.cli_tools` and can be executed with
`python -m`. This keeps the distribution tidy while preserving module execution
for debugging.

| Module | Example command | Purpose |
|--------|-----------------|---------|
| `library.utils.cli_tools.check_determinism` | `python -m library.utils.cli_tools.check_determinism --baseline out1 --candidate out2` | Compare CSV hashes between runs. |
| `library.utils.cli_tools.chunk_io_main` | `python -m library.utils.cli_tools.chunk_io_main --input data.csv --final-out copy.csv` | Stream CSV files in deterministic order with Unix newlines. |
| `library.utils.cli_tools.csv_utils_main` | `python -m library.utils.cli_tools.csv_utils_main --input data.csv --final-out clean.csv` | Re-serialise CSVs with consistent ordering. |
| `library.utils.cli_tools.dtype_inspector_main` | `python -m library.utils.cli_tools.dtype_inspector_main --log-level INFO` | Inspect pandas dtypes emitted by pipelines. |
| `library.utils.cli_tools.get_activities` | `python -m library.utils.cli_tools.get_activities --limit 10 --final-out output/activities_smoke.csv` | Generate synthetic activity rows and emit deterministic CSV + `.meta.yaml` artefacts for smoke runs. |
| `library.utils.cli_tools.get_document_type` | `python -m library.utils.cli_tools.get_document_type --input docs.csv` | Classify document rows with bundled heuristics. |
| `library.utils.cli_tools.get_input_initialisation` | `python -m library.utils.cli_tools.get_input_initialisation --same-doc a.xlsx --all-doc b.xlsx` | Merge Excel workbooks into canonical CSV. |
| `library.utils.cli_tools.mapper_batch_main` | `python -m library.utils.cli_tools.mapper_batch_main --input ids.csv --final-out mapped.csv` | Batch mapping from ChEMBL IDs to UniProt accessions. |
| `library.utils.cli_tools.mapper_main` | `python -m library.utils.cli_tools.mapper_main --input ids.csv --final-out mapped.csv` | Interactive mapper for diagnostics. |
| `library.utils.cli_tools.pipeline_targets_main` | `python -m library.utils.cli_tools.pipeline_targets_main --input targets.csv` | Re-run target pipeline harness with cached data. |
| `library.utils.cli_tools.table_quality_main` | `python -m library.utils.cli_tools.table_quality_main --input data.csv --table-name data` | Produce column-level quality reports. |

Both mapping CLIs honour `io.na_markers` and `io.keep_na_markers` from
[`CONFIG.md`](./CONFIG.md). Each module exposes a `main(argv)` function for reuse in tests.

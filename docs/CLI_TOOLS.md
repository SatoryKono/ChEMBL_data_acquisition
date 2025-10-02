# Utility CLI modules

> **Languages:** [English](./CLI_TOOLS.md) · [Русский](./CLI_TOOLS_RU.md)

The lightweight helper commands previously exposed as loose modules under
`scripts/` now live in the `library.utils.cli_tools` package so that they can be
invoked with `python -m <module>`. The relocation keeps the distribution tidy
while preserving direct module execution for ad-hoc debugging.

| Module | Typical command | Purpose |
|--------|-----------------|---------|
| `library.utils.cli_tools.check_determinism` | `python -m library.utils.cli_tools.check_determinism` | Validate that deterministic CSV writing still produces matching hashes. |
| `library.utils.cli_tools.chunk_io_main` | `python -m library.utils.cli_tools.chunk_io_main --input data.csv --output copy.csv` | Stream input CSV files in chunks while preserving deterministic ordering. |
| `library.utils.cli_tools.csv_utils_main` | `python -m library.utils.cli_tools.csv_utils_main --input data.csv` | Re-serialise arbitrary CSV files with deterministic ordering. |
| `library.utils.cli_tools.dtype_inspector_main` | `python -m library.utils.cli_tools.dtype_inspector_main --log-level INFO` | Inspect pandas dtypes emitted by the core `get_*_data` pipelines. |
| `library.utils.cli_tools.get_activities` | `python -m library.utils.cli_tools.get_activities --limit 10` | Emit synthetic activity rows to verify logging and CLI wiring. |
| `library.utils.cli_tools.get_document_type` | `python -m library.utils.cli_tools.get_document_type --input docs.csv` | Classify document rows with the bundled heuristics for unit tests. |
| `library.utils.cli_tools.get_input_initialisation` | `python -m library.utils.cli_tools.get_input_initialisation --same-doc path.xlsx --all-doc path.xlsx` | Combine Excel workbooks that initialise input pairs for downstream QA. |
| `library.utils.cli_tools.mapper_batch_main` | `python -m library.utils.cli_tools.mapper_batch_main --input ids.csv --output mapped.csv` | Map ChEMBL identifiers to UniProt accessions using batch configuration. |
| `library.utils.cli_tools.mapper_main` | `python -m library.utils.cli_tools.mapper_main --input ids.csv --output mapped.csv` | Lightweight interactive mapper for quick lookups and diagnostics. |
| `library.utils.cli_tools.pipeline_targets_main` | `python -m library.utils.cli_tools.pipeline_targets_main --input targets.csv` | Run the cached target pipeline harness to refresh stored artefacts and exercise staging flags (`--raw-out`, `--raw-format`, `--id-cols`, `--no-reindex-raw`, `--normalize-at-export` / `--no-normalize-at-export`). |
| `library.utils.cli_tools.table_quality_main` | `python -m library.utils.cli_tools.table_quality_main --input data.csv --table-name data` | Generate column-level quality reports for arbitrary CSV datasets with optional sampling and column filters. |

Both mapping CLIs honour the [`io.na_markers`](CONFIG_EN.md#io) list when filtering
placeholder identifiers and use [`io.keep_na_markers`](CONFIG_EN.md#io) to decide
whether to keep those placeholders in the mapping input.

All modules continue to expose a `main` function so they can still be wired into
`pyproject.toml` entry points. When invoking them programmatically, import the
module from `library.utils.cli_tools` and call `main(argv)` to reuse the command
line interfaces in tests.

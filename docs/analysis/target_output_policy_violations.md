# get_target_data.py output anomalies

## Expected behaviour

Per the Output Standardization Policy every `get_*_data.py` command must emit exactly three CSV artefacts via `save_standard_outputs`: the consolidated dataset, a correlation report and a QC report. No other files should remain after a standard run.

## Why extra artefacts are produced

### 1. Legacy outputs are force-enabled by the CLI harness

`scripts/get_data.py` forcibly sets `args.keep_intermediate = True` before delegating to the underlying target CLI. This implicitly turns on `emit_legacy_artifacts` inside `run_cli_command`, because the helper maps either `--debug` or `--keep-intermediate` to the legacy output bundle. As a result, even when the operator does not pass `--emit-legacy-artifacts`, every orchestrated run keeps the raw dumps, failure CSVs and sidecars. 【F:scripts/get_data.py†L63-L82】【F:library/cli_utils.py†L289-L309】

### 2. Raw-mode pipelines call `run_pipeline` without disabling legacy writers

When the ChEMBL stage is executed in "raw dump" mode (`normalize_at_export=False`), `_run_pipeline_with_meta` delegates to `run_pipeline` with its default parameters. The latter assumes `emit_legacy_artifacts=True` and therefore writes both the canonical bundle and the historical artefacts (raw CSV, `_failure_cases.csv`, `.meta.yaml`). That reproduces files such as `target_raw.csv` and `target_raw_failure_cases.csv` alongside the standard outputs. 【F:library/cli/commands/get_target_data.py†L2184-L2268】【F:library/cli_utils.py†L524-L583】

### 3. The validation/export step still persists raw and metadata files when legacy mode is active

`validate_and_write` keeps the backwards-compatible branches for `emit_legacy_artifacts`. Once that flag is enabled (see #1), the function writes the raw frame to `raw_out`, stores metadata sidecars and leaves the validation failure CSV in place instead of pruning it. These branches are what generate files such as `target_raw.csv`, `target_raw.meta.yaml` and `target_failure_cases.csv`. 【F:library/cli/commands/get_target_data.py†L3636-L3777】

### 4. Optional post-processing creates yet another export

If `--postprocess` is supplied, `run_target_postprocess_if_requested` materialises `output_postprocessed.targets.csv` next to the canonical outputs. Because the standard output bundle is already written, this extra CSV is effectively a duplicate unless post-processing is explicitly required. 【F:library/cli/commands/get_target_data.py†L429-L505】

## Remediation suggestions

* Remove the forced `keep_intermediate` toggle in `scripts/get_data.py` or gate it behind an explicit debug flag so regular runs keep legacy artefacts disabled.
* Pass `emit_standard_outputs=True, emit_legacy_artifacts=False` when invoking `_run_pipeline_with_meta` for raw-mode fetches; the raw CSV can still be streamed by the custom writer without the pipeline producing additional files.
* Simplify `validate_and_write` by dropping the legacy save branches (or guarding them behind an explicit flag) so that only `save_standard_outputs` is executed in normal runs.
* Consider writing the post-processed dataset to an in-memory buffer and reusing the canonical writer instead of leaving a fourth CSV on disk.

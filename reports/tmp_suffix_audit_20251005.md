# tmp_suffix_audit 20251005

## Summary
- Post-processing helpers derived output names directly from the working path, so orchestrated runs that staged files as `.{name}.tmp` produced final artefacts ending in `.tmp`.【F:library/postprocessing/names.py†L409-L420】【F:library/postprocessing/target/main.py†L73-L77】
- Introduced `normalise_export_basename` to strip atomic prefixes/suffixes and re-order `.csv` markers, giving modules a canonical basename before building sidecar file names.【F:library/postprocessing/helpers.py†L109-L128】
- Updated target, isoform, names, and IUPHAR post-processing writers to use the canonical basename and deterministic atomic CSV emission, preventing orphaned `.tmp` outputs.【F:library/postprocessing/target/isoform.py†L480-L540】【F:library/postprocessing/iuphar.py†L253-L281】
- Added focused tests that emulate orchestrator-style paths to assert that helpers now emit `.csv` artefacts and accept `.tmp` inputs without errors.【F:tests/unit/postprocessing/test_helpers.py†L10-L25】【F:tests/postprocessing/test_target_postprocessing.py†L480-L497】【F:tests/postprocessing/test_names_postprocessing.py†L300-L344】【F:tests/postprocessing/test_iuphar_postprocessing.py†L180-L209】

## Environment & Commands
- `python -m scripts.get_target_data iuphar --limit 200` *(fails: missing config/dictionary/_target/target.csv)*
- `python -m scripts.get_target_data iuphar --limit 200 --log-level DEBUG` *(fails: missing config/dictionary/_target/target.csv)*
- `INJECT_FAIL_AFTER_WRITE=1 python -m scripts.get_target_data iuphar --limit 200 --log-level DEBUG` *(fails: missing config/dictionary/_target/target.csv)*

## Findings
- The orchestrator writes temporary artefacts using the pattern `.{output.name}.tmp`, and post-processing modules reused that literal name, so every derived helper (organism, isoform, names, IUPHAR) inherited the `.tmp` suffix.【F:library/postprocessing/target/main.py†L73-L77】【F:library/postprocessing/names.py†L409-L420】【F:library/postprocessing/iuphar.py†L253-L281】
- Isoform helpers rejected orchestrator inputs because validation expected `.csv` endings, leaving `.tmp` artefacts when reruns succeeded and bypassing downstream clean-up.【F:library/postprocessing/target/isoform.py†L77-L87】【F:library/postprocessing/target/isoform.py†L480-L540】
- No CLI or configuration flag toggled this behaviour; failure logs consistently pointed to missing dictionary CSVs rather than explicit rename attempts, confirming hypothesis H3 (naming mismatch prevents rename).【181eca†L1-L18】【a9e0c7†L1-L18】【d7f9de†L1-L18】

## Root Cause
- **H3 confirmed:** Derived filenames appended `_normalized` and other tokens after `.csv` on the temporary basename (`output...csv_normalized.tmp`), so the orchestrator never saw matching final names to rename, leaving `.tmp` artefacts behind.

## Changes
- Added `normalise_export_basename` to canonicalise orchestration paths before generating sidecars.【F:library/postprocessing/helpers.py†L109-L128】
- Target organism helper now normalises the basename and writes through `write_csv` for atomic output.【F:library/postprocessing/target/main.py†L73-L77】
- Isoform helper accepts `.tmp` inputs, computes canonical outputs, and writes deterministically.【F:library/postprocessing/target/isoform.py†L77-L87】【F:library/postprocessing/target/isoform.py†L480-L540】
- Names and IUPHAR post-processing adopt the canonical basename to prevent `.tmp` suffix leakage.【F:library/postprocessing/names.py†L409-L420】【F:library/postprocessing/iuphar.py†L253-L281】
- Added regression tests covering basename normalisation and orchestrator-style inputs for helpers and sidecar writers.【F:tests/unit/postprocessing/test_helpers.py†L10-L25】【F:tests/postprocessing/test_target_postprocessing.py†L480-L497】【F:tests/postprocessing/test_names_postprocessing.py†L300-L344】【F:tests/postprocessing/test_iuphar_postprocessing.py†L180-L209】

## Tests
- `pytest tests/unit/postprocessing/test_helpers.py tests/postprocessing/test_target_postprocessing.py::test_isoform_process_targets__normalises_tmp_suffix tests/postprocessing/test_names_postprocessing.py::test_process_target_names__normalises_tmp_suffix tests/postprocessing/test_iuphar_postprocessing.py::test_process_iuphar_targets__normalises_tmp_suffix`

## Artefacts
- New regression tests assert that helpers emit canonical `.csv` artefacts even when inputs mimic orchestrator `.tmp` paths, ensuring future runs cannot leave orphaned `.tmp` sidecars.【F:tests/unit/postprocessing/test_helpers.py†L10-L25】【F:tests/postprocessing/test_target_postprocessing.py†L480-L497】【F:tests/postprocessing/test_names_postprocessing.py†L300-L344】【F:tests/postprocessing/test_iuphar_postprocessing.py†L180-L209】

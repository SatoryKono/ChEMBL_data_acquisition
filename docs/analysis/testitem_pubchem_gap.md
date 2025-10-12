# Test Item PubChem Enrichment Gap

## Summary
- Running the dedicated CLI `scripts/get_testitem_data.py` with a proper `--base-path`
  causes `library/cli/parser.prepare_io_paths` to resolve the data root and the
  subsequent call to `apply_config_overrides` passes that path into the
  configuration loader. This keeps placeholders such as
  `$CHEMBL_DA_BASE_PATH/cache/pubchem_cid_cache.json` pointing at the repository
  cache directory and PubChem augmentation works as expected.
- The orchestration command implemented in
  `library/cli/commands/get_data.py` launches `get_testitem_data.py` as a
  subprocess, but the helper `PipelineStep.build_arguments` that assembles the
  forwarded CLI parameters omits `--base-path`, `--input-dir`, and
  `--output-dir`. As a consequence the subprocess falls back to the default
  placeholder base (`~/.local/share/chembl-da`) and looks for cache and
  dictionary resources outside the repository checkout.
- With the base path redirected to the home directory the PubChem CID cache and
  other lookup tables referenced via `$CHEMBL_DA_BASE_PATH` are missing, so the
  enrichment routine `library/pipelines/testitem/pubchem.add_pubchem_data` does
  not find any CID hints and all `pubchem_*` columns remain empty when the
  orchestrator is used.

## Evidence

| Context | Relevant code |
|---------|----------------|
| Base path resolution in the standalone CLI | `prepare_io_paths` normalises `--base-path` and feeds it into the config loader, so direct executions use the repository data directory.【F:library/cli/parser.py†L838-L897】【F:library/cli/parser.py†L635-L684】 |
| Placeholder fallback | `load_config` resolves `$CHEMBL_DA_BASE_PATH` against the provided base path; if none is supplied the helper defaults to `~/.local/share/chembl-da`.【F:library/config/loader.py†L403-L445】【F:library/config/env.py†L31-L44】 |
| PubChem paths depend on `$CHEMBL_DA_BASE_PATH` | The default configuration stores the CID cache and related artefacts under `$CHEMBL_DA_BASE_PATH/cache`, so the wrong base path leaves the subprocess without these resources.【F:config/config.yaml†L180-L213】 |
| Orchestrator subprocess arguments | `_run_testitem_subprocess` relies on `PipelineStep.build_arguments`, which only forwards `--config`, `--input`, and `--final-out` – none of the directory options reach the child process.【F:library/cli/commands/get_data.py†L1552-L1620】【F:library/pipelines/registry.py†L59-L120】 |
| Working implementation in the legacy script | The compatibility wrapper in `scripts/get_data.py` still adds `--base-path`, `--input-dir`, and `--output-dir`, which explains why running the legacy script preserves PubChem enrichment.【F:scripts/get_data.py†L240-L316】 |

## Root cause
The orchestrator in `library/cli/commands/get_data.py` does not forward the
resolved base path (nor the derived input/output directories) to the
`get_testitem_data.py` subprocess. Without these arguments the subprocess loads
`config/config.yaml` with `base_path=None`, so all `$CHEMBL_DA_BASE_PATH` placeholders resolve to the fallback directory under the user’s home. The PubChem
pipeline expects the CID cache and cached lookups in the repository `data/cache`
folder and therefore cannot populate any `pubchem_*` columns.

## Recommended fixes
1. Extend `_run_testitem_subprocess` to forward `--base-path`, `--input-dir`, and
   `--output-dir` (mirroring the behaviour in `scripts/get_data.py`).
2. Alternatively export `CHEMBL_DA_BASE_PATH` in the orchestrator environment
   before spawning the subprocess, keeping placeholders intact for all child
   commands.
3. Add an integration test that runs the orchestrator against the fixture input
   and asserts that at least one `pubchem_*` column contains non-null values, so
   regressions in argument forwarding are caught early.

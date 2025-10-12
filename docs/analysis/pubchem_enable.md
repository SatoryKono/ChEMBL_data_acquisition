# PubChem Enable Flag Audit

## Configuration defaults

- `config/config.yaml` ships with `sources.pubchem.enable` set to `true`, keeping PubChem enrichment active unless callers explicitly opt out through the documented `CHEMBL_DA_PUBCHEM_ENABLE` variable.【F:config/config.yaml†L184-L208】
- The sample override file `config/example.yaml` reiterates that the flag should stay enabled and still defaults it to `true` for local development scenarios.【F:config/example.yaml†L14-L23】
- The Pydantic model for `PubChemCfg` defines `enable` with a default of `True`, so programmatic instantiations of `Config` inherit the same behaviour even without loading YAML files.【F:library/config/models.py†L505-L543】

## CLI behaviour

- The `get_testitem_data` command maps CLI options such as timeout and batch size but does not expose `sources.pubchem.enable`, so calling the pipeline via the official entry point preserves the configuration defaults.【F:library/cli/commands/get_testitem_data.py†L909-L958】
- Within the pipeline, the `pubchem_enabled` flag gates augmentation using the value from configuration, ensuring enrichment only stops when the config explicitly disables it.【F:library/pipelines/testitem/cli.py†L756-L768】

## Environment and CI

- The shared CI workflow exports a minimal set of environment variables and does not override `CHEMBL_DA_PUBCHEM_ENABLE`, preventing accidental deactivation in automated runs.【F:.github/workflows/ci.yml†L1-L98】
- The provided `.env.example` file also refrains from defining the enable flag, helping local developers keep PubChem augmentation active by default.【F:.env.example†L1-L15】

## Regression safeguards

- `tests/unit/test_config_pubchem_enable.py` asserts that the plain `Config` model, the shipped YAML configuration, and the default CLI invocation all resolve to `sources.pubchem.enable = True`, adding coverage that will fail if a future change flips the toggle off by default.【F:tests/unit/test_config_pubchem_enable.py†L1-L31】

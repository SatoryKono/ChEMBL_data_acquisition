# Release Notes

> **Languages:** [English](./RELEASE_NOTES.md) · [Русский](../ru/RELEASE_NOTES.md)

## Unreleased

- Clarified the minimum supported Python version as 3.11 across documentation and runtime checks.
- Documented the requirement to configure `api.user_agent` with a real contact, including validation behaviour and the
  `CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT` environment override.
- Downgraded ChEMBL and PubChem cache hit/miss logging to DEBUG to reduce noise in default INFO logs.
- Added UniProt isoform fallbacks, richer cross-reference extraction (AlphaFold, PDB, Ensembl), and additional logging to
  trace identifier coverage through the target merge pipeline.
- Normalised CLI path overrides so `--base-path`, `--input-dir`, `--output-dir`, and `--cache-dir` propagate into
  `local.io` configuration fields.


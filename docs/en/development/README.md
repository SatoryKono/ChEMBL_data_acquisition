# Developer handbook

This folder complements the public documentation with engineering-specific
information.

- [`LOCAL_SETUP.md`](./LOCAL_SETUP.md) — local environment setup, tooling and
  conventions.
- [`TESTING.md`](./TESTING.md) — pytest structure, determinism policy, reporting.
- [`CI_CD.md`](./CI_CD.md) — recommended CI jobs and quality gates.
- [`RELEASE.md`](./RELEASE.md) — checklist for releasing new versions and updating
  dictionary data.

## CLI Inconsistencies

There is a known inconsistency in the command-line interface design across different scripts. Some scripts, like `get_target_data.py`, use subcommands to select functionality (e.g., `uniprot`, `chembl`), while others, such as `get_document_data.py`, use a `--mode` flag.

This is a point of architectural drift that may be addressed in a future refactoring. For now, developers should be aware of this difference when working with the various ETL scripts.

Contributions should follow the coding standards in `pyproject.toml` and respect
the deterministic testing policy described below.

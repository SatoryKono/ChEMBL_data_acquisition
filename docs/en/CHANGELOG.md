# Changelog

> **Languages:** [English](./CHANGELOG.md) · [Русский](../ru/CHANGELOG.md)
## [0.1.10] - 2025-10-14
- Fixed the duplicate registration of `--log-level` and `--input-dir` in
  `scripts/get_target_data.py`, restoring CLI execution on Windows
  orchestrations that rely on those flags.
## [0.1.9] - 2025-10-21
- Extended `scripts/get_target_data.py` with CLI options for `--log-level` and
  `--input-dir`, aligning the metadata output with orchestrator expectations and
  normalising UTC defaults.
## [0.1.7] - 2025-10-19
- Documentation and protocol refresh covering PR-0 through PR-8: updated pipeline diagrams, CLI examples with `--emit-legacy-artifacts`, refreshed configuration notes and regenerated the bilingual protocol DOCX.
- Added a dedicated `document_schema` module exposing `DocumentSchema` and `validate_document_frame` for enriched CrossRef/OpenAlex metadata validation with strict typing.
## [0.1.9] - 2025-10-20
- Added legacy-compatible options to `scripts/get_target_data.py`, accepting
  positional commands and `--log-level`/`--input-dir` flags so orchestration
  runners continue to work without wrapper scripts.

## [0.1.8] - 2025-10-20
- Restored CLI bootstrapping for `scripts/get_target_data.py` so direct Windows
  invocations resolve internal modules without manual `PYTHONPATH` tweaks.


All notable changes to this project will be documented in this file.


## [0.1.7] - 2025-10-19
- Documentation and protocol refresh covering PR-0 through PR-8: updated pipeline diagrams, CLI examples with `--emit-legacy-artifacts`, refreshed configuration notes and regenerated the bilingual protocol DOCX.


## [0.1.7] - 2025-10-19
- Added retention-aware project cleanup with ``--check``/``--days`` flags that track
  planned deletions and skip fresh files by default.
- Aligned the orchestration summary with the canonical ``data/output`` directory and
  logged per-date_tag CSV inventories to guarantee the expected 15 artefacts.
- Standardised the target postprocessing pipeline with a thin CLI that writes
  deterministic CSV artefacts and metadata sidecars.
- Added a dedicated Pandera schema together with UniProt and GtoPdb enrichment
  using resilient HTTP clients and structured logging.
- Extended automated coverage with merge, integration, and CLI tests to keep
  the target export reproducible.
- Added a focused assay postprocessing module deriving UTC timestamps and year
  attributes directly from ChEMBL metadata.
- Introduced a dedicated Pandera schema for the streamlined assay export and
  regression tests covering timestamp handling and dictionary enrichment.

## [0.1.6] - 2025-10-18
  Added structured logging utilities with run identifiers and stage duration
  tracking, plus resilient HTTP retry wrappers adopted by the Chembl, PubChem,
  UniProt, and Crossref clients.
- Enforced Pandera-based validation for activity, assay, document, target, and
  test item postprocessing pipelines with structured logging on failures.
- Restored the `list_output_files` helper used by the orchestrator summary to
  avoid runtime failures after successful pipeline executions.

## [0.1.7] - 2025-10-19
- Introduced a modern document enrichment step that retrieves ChEMBL records,
  normalises DOIs, augments metadata from CrossRef/OpenAlex, and emits QC/
  correlation reports validated via the new `fetch_normalize_document` helper.


## [0.1.5] - 2025-10-18
- Hardened the orchestrator by exiting with a detailed error when CSV artefacts are missing, listing the discovered files for faster triage.
- Normalised byte-valued dictionary resource paths before ``Path`` conversion so type checking accepts manifest references provided as bytes.

## [0.1.4] - 2025-10-17
- Regenerated the bilingual documentation set after the October audit: README pairs
  now stay in lockstep, the protocol DOCX export is reproducible via
  `scripts/convert_md_to_docx.py`, and the change log reflects the approved
  CHEMBL-DM01/2.1 baseline.

## [0.1.3] - 2025-10-16
- Refined hierarchy resolution to favour curated relationships and introduce safe fallbacks when source data is incomplete.
- Extended structured INFO-level logging to cover fallback decisions so hierarchy handling remains auditable.

## [0.1.2] - 2025-10-09
- Added consistent INFO-level logging across `get_document_data`, `get_assay_data`, `get_testitem_data`, and `get_activity_data` to make ETL progress traceable.

## [0.1.1] - 2025-10-02
- Moved all top-level Markdown documentation into the `docs/` directory and updated internal links.
- Bumped the package version to reflect the documentation reorganisation.

## [0.1.0] - 2025-09-28
- Initial public release.

# Changelog

> **Languages:** [English](./CHANGELOG.md) · [Русский](../ru/CHANGELOG.md)

All notable changes to this project will be documented in this file.


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

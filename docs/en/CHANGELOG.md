# Changelog

> **Languages:** [English](./CHANGELOG.md) · [Русский](../ru/CHANGELOG.md)

All notable changes to this project will be documented in this file.

## [0.1.4] - 2025-10-21
- Harmonised the pytest layout by relocating post-processing regression suites under `tests/integration/postprocessing` and
  updating module markers.
- Refreshed the documentation set in both languages to reference the new structure, clarify reporting outputs and remove stale
  references to the deprecated `tests/postprocessing` path.
- Synced the testing handbook and runbook with the current QA checklist and added status notes to completed improvement
  proposals.

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

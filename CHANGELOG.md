# Changelog / Журнал изменений

All notable changes to this project are documented here following
[Semantic Versioning](https://semver.org/) and a bilingual format. Each
entry mirrors the repository documentation in English and Russian.

## [0.2.0] - 2025-10-02
### Added / Добавлено
- Bilingual version banners across the README and detailed guides to keep
  Russian and English content synchronised.
- Top-level changelog consolidating release history and pointing to the
  QA checklist and reference manuals.

### Changed / Изменено
- Raised the declared project version to `0.2.0` in packaging metadata
  and documentation.
- Reworked README_EN.md and README_RU.md to present a consistent
  structure (overview, installation, pipelines, configuration,
  datasets, QA/testing, support) with updated command examples.
- Updated all reference manuals in `docs/` to share common terminology,
  explicitly reference the canonical CLI modules, and surface the current
  version at the top of each file.
- Converted `docs/RELEASE_NOTES.md` into a bilingual summary that links
  back to this changelog.

### Fixed / Исправлено
- Removed outdated statements about unreleased features and clarified the
  scope of `--raw-out` and related flags across the documentation.
- Ensured every cross-reference points to existing modules, directories,
  or configuration files.
- Restored the exported `MAX_DIFF_KEY_EXPORT` constant in the QA helper to keep
  regression tests functional.

## [0.1.0] - 2024-??-??
### Added / Добавлено
- Initial public release with CLI pipelines, configuration management,
  validation schemas, and the QA suite.
- Baseline English/Russian documentation set and reproducible packaging
  artefacts.


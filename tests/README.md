# Test suite layout

The test suite is organised around the key scenarios of the ChEMBL data acquisition pipeline:

- `unit/` – fast checks for isolated helpers (loading, normalisation, validation logic).
- `integration/` – composed workflows that verify enrichment rules, schema validation and failure handling. Snapshot-backed
  regression suites that previously lived in `tests/postprocessing/` now reside under
  `tests/integration/postprocessing/` to keep all multi-module scenarios together.
- `e2e/` – deterministic end-to-end runs of the test-item pipeline on synthetic fixtures, including export idempotence.
- `resources/` – small CSV snapshots used by the integration and e2e scenarios.

Standalone smoke checks that previously lived alongside the legacy `tests/run_tests.py`
wrapper were moved into the directories above:

- activity column filtering helpers now reside in `tests/unit/test_activity_output_columns.py`;
- assay output pruning checks live in `tests/unit/test_assay_output_columns.py`;
- CLI logging scenarios for `get_activity_data` were relocated to `tests/e2e/test_activity_logging.py` and
  `tests/e2e/test_logging_get_activity_data.py`.

When adding new coverage, select the directory that matches the scope (unit, integration, postprocessing or e2e) instead of
placing the file at the repository root. This keeps the suite hierarchy stable and avoids accidental collection gaps.

## Naming conventions

- Test modules follow `test_<module>.py` (for example, `test_normalize_rules.py`).
- Individual tests use `test_<unit_of_work>__<case>()` to encode both the subject and the scenario variant.
- Apply `@pytest.mark.unit`, `@pytest.mark.integration`, `@pytest.mark.e2e`, `@pytest.mark.slow` and
  `@pytest.mark.network` consistently so that subsets can be selected via `-m` filters.

Shared fixtures live in `tests/conftest.py`. They configure a deterministic environment, disable outbound HTTP calls and expose helpers such as `sample_input_csv` and `snapshot_resource`.

## Post-processing helper coverage

`tests/integration/postprocessing/test_target_postprocessing.py` exercises the helper modules (`helpers.py`, `cellularity.py`, `multifunctional.py`) that reproduce the Power Query logic for the target lookup. The tests assert deterministic CSV loading, taxonomy normalisation, cellularity labels, multifunctional flag derivation and byte-identical exports using the snapshots under `tests/resources/target_postprocess_power_query_*.csv`.

`tests/integration/test_target_postprocess_table.py` feeds the same fixtures through `helpers.postprocess_target_table_file` to validate file I/O, export naming and idempotency.
`tests/integration/test_pipeline_quality_matrix.py` exercises a miniature CSV pipeline backed by golden fixtures in `tests/resources/pipeline_quality/`. The scenario validates schema enforcement, preprocessing, enrichment fallbacks, warning emission, export invariants and idempotence in line with the QA checklist below.

## Integration enrichment checklist

- [x] PubChem augmentation: cache hits/misses, polymer handling and TTL expiry (`tests/integration/test_pubchem_augmentation.py`).

## CLI pipeline orchestration

Declarative pipeline scenarios are now covered in both integration and end-to-end suites. The integration module `tests/integration/test_pipeline_config_loading.py` exercises `load_pipeline_config`, asserting that YAML-defined steps, runners and metrics definitions are resolved into deterministic orchestrator payloads. The complementary end-to-end suite `tests/e2e/test_get_cli_pipelines.py` drives every `scripts/get_*` entrypoint with those declarative configurations, patches the network-heavy stages with deterministic stubs and validates that the orchestrator:

- loads miniature CSV fixtures and normalises core columns,
- emits structured WARN/ERROR events for missing or duplicated records,
- records pipeline metadata such as `pipeline_version`, per-step metric payloads, runner identifiers and a deterministic
  `generated_at` timestamp derived from the CLI invocation,
- writes the derived tables with deterministic ordering and derived fields, and
- honours the edge-case flags `--skip-existing`, `--limit` and `--dry-run` without invoking unintended side effects.

The scenarios exercise success and failure paths for `get_testitem_data`, `get_document_data`, `get_target_data`, `get_assay_data`, `get_tissue_data` and `get_activity_data`, bringing both the declarative configuration loading and CLI surface under the deterministic test umbrella.

## Running tests and generating reports

Install dependencies (see the repository `README.md`) and run the suite via the canonical reporting wrapper:

```bash
python scripts/run_tests.py
```

(`python -m scripts.run_tests` remains available for backwards compatibility.)

The command executes `pytest` with the default configuration, writes the full protocol to `reports/test_report.json` and produces a human readable summary in `reports/test_summary.md`. Both artefacts contain Git metadata, timing information, a per-test breakdown and the overall success rate. The JSON payload exposes a `summary` section (totals and a `success_rate` ratio computed as `(passed + xfailed) / max(1, total - skipped)`, ranging from 0.0 to 1.0), while the Markdown file includes a `Success rate: NN.NN%` bullet for quick inspection. The wrapper enforces the ≥95% success-rate policy: if the computed ratio drops below the threshold, it emits an error log and returns a non-zero exit code even when pytest itself reports success. All invocations also configure structured logging – log events are mirrored to `data/logs/run_tests_<YYYYMMDD>.log` (or the directory defined by `CHEMBL_DA_BASE_PATH`). Pass `--verbose` to lift the logger to DEBUG and forward the same verbosity to pytest’s log capture.

The generated JSON/Markdown files are git-ignored; CI publishes them together with the coverage directory as the `test-reports-<python-version>` artefact so that the latest results can be downloaded directly from GitHub Actions.

To guard against stalled executions, the wrapper accepts a `--pytest-timeout <seconds>` flag (or the `CHEMBL_DA_PYTEST_TIMEOUT` environment variable) that aborts the pytest run once the limit is exceeded. When triggered, the subprocess is terminated, the partial output is logged and the wrapper returns a non-zero exit code.

To focus on a subset, forward additional arguments to pytest after the `--` separator, for example `python -m scripts.run_tests -- -m unit`. Combine `--verbose` with the forwarding flag to observe detailed DEBUG events in both the console and the generated log file.

Individual modules can be targeted by pointing pytest at a directory, for example `pytest tests/unit` or `pytest tests/integration -k enrich` to filter by test name.

When developing additional scenarios, keep the guardrails documented in `tests/conftest.py` (seed fixing, network ban, temporary directories) to preserve reproducibility. All new tests should emit deterministic output so that `scripts/run_tests.py` can regenerate the reports without spurious diffs. The legacy `tests/run_tests.py` entry point remains available for now but issues a `DeprecationWarning` and will be retired once downstream jobs migrate to the canonical wrapper.

## End-to-end scenario checklist

`tests/e2e/test_get_data_end_to_end.py` drives the `scripts.get_data` orchestrator against the miniature fixtures in `tests/resources/pipeline_inputs`. The stubbed pipelines validate inputs, perform deterministic normalisation/post-processing and emit the canonical filenames in a temporary directory. The scenario also re-runs the workflow after intentionally corrupting `document.csv` (column removed) to verify schema validation, error logging and cleanup behaviour. The assertions cover the full QA checklist:

- [x] Загрузка входных CSV и валидация схемы/типов/обязательных колонок
- [x] Нормализация и предобработка (включая кодировки, разделители)
- [x] Обогащение данными из словарей/справочников (в т.ч. отсутствие соответствий)
- [x] Правила трансформации и расчета флагов/категорий
- [x] Обработка пропусков, дублей, конфликтных значений
- [x] Детализация логирования и уровни WARN/ERROR при аномалиях
- [x] Итоговая сборка результатов: структура, сортировка, инварианты
- [x] Постобработка и экспорт: корректность форматов/имён/путей
- [x] Деградационные кейсы: частичные данные, пустые файлы, неверный заголовок
- [x] Идемпотентность: повторный запуск на тех же входах даёт тот же результат
- [x] Метаданные пайплайна: `pipeline_version`, метрики по шагам и идентификаторы раннеров сохраняются и проверяются
- [x] Пограничные флаги CLI: `--skip-existing`, `--limit`, `--dry-run` не нарушают детерминированность и корректность

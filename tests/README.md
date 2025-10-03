# Test suite layout

The test suite is organised around the key scenarios of the ChEMBL data acquisition pipeline:

- `unit/` – fast checks for isolated helpers (loading, normalisation, validation logic).
- `integration/` – composed workflows that verify enrichment rules, schema validation and failure handling.
- `e2e/` – deterministic end-to-end runs of the test-item pipeline on synthetic fixtures, including export idempotence.
- `resources/` – small CSV snapshots used by the integration and e2e scenarios.

Shared fixtures live in `tests/conftest.py`. They configure a deterministic environment, disable outbound HTTP calls and expose helpers such as `sample_input_csv` and `snapshot_resource`.

## Integration enrichment checklist

- [x] PubChem augmentation: cache hits/misses, polymer handling and TTL expiry (`tests/integration/test_pubchem_augmentation.py`).

## Running tests and generating reports

Install dependencies (see the repository `README.md`) and run the suite via the reporting wrapper:

```bash
python tools/run_tests.py
```

The command executes `pytest` with the default configuration, writes the full protocol to `reports/test_report.json` and produces a human readable summary in `reports/test_summary.md`. Both artefacts contain Git metadata, timing information, a per-test breakdown and the overall success rate. The JSON payload exposes a `summary` section (totals and `success_rate = passed/total`), while the Markdown file includes a `Success rate: NN.NN%` bullet for quick inspection. The wrapper enforces the ≥75% success-rate policy: if the computed ratio drops below the threshold, it emits an error log and returns a non-zero exit code even when pytest itself reports success.

To focus on a subset, pass extra arguments after `--pytest-args`, for example `python tools/run_tests.py --pytest-args -m unit`.

Individual modules can be targeted by pointing pytest at a directory, for example `pytest tests/unit` or `pytest tests/integration -k enrich` to filter by test name.

When developing additional scenarios, keep the guardrails documented in `tests/conftest.py` (seed fixing, network ban, temporary directories) to preserve reproducibility. All new tests should emit deterministic output so that `tools/run_tests.py` can regenerate the reports without spurious diffs.

## End-to-end scenario checklist

`tests/e2e/test_get_data_end_to_end.py` drives the `scripts.get_data` orchestrator against the miniature fixtures in `tests/data`. The stubbed pipelines validate inputs, perform deterministic normalisation/post-processing and emit the canonical filenames in a temporary directory. The scenario also re-runs the workflow after intentionally corrupting `document.csv` (column removed) to verify schema validation, error logging and cleanup behaviour. The assertions cover the full QA checklist:

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

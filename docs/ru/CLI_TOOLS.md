# Вспомогательные CLI

> **Языки:** [English](../en/CLI_TOOLS.md) · Русский

Мини-утилиты расположены в `library.utils.cli_tools` и запускаются через
`python -m`.

| Модуль | Пример | Назначение |
|--------|--------|------------|
| `library.utils.cli_tools.check_determinism` | `python -m library.utils.cli_tools.check_determinism --baseline out1 --candidate out2` | Сравнение хэшей CSV между запусками. |
| `library.utils.cli_tools.chunk_io_main` | `python -m library.utils.cli_tools.chunk_io_main --input data.csv --final-out copy.csv` | Потоковая обработка CSV с сохранением порядка. |
| `library.utils.cli_tools.csv_utils_main` | `python -m library.utils.cli_tools.csv_utils_main --input data.csv --final-out clean.csv` | Перезапись CSV с детерминированным порядком. |
| `library.utils.cli_tools.dtype_inspector_main` | `python -m library.utils.cli_tools.dtype_inspector_main --log-level INFO` | Просмотр типов pandas. |
| `library.utils.cli_tools.get_activities` | `python -m library.utils.cli_tools.get_activities --limit 10 --final-out output/activities_smoke.csv` | Синтетические активности с записью детерминированного CSV и `.meta.yaml` для smoke-запусков. |
| `library.utils.cli_tools.get_document_type` | `python -m library.utils.cli_tools.get_document_type --input docs.csv` | Классификация документов. |
| `library.utils.cli_tools.get_input_initialisation` | `python -m library.utils.cli_tools.get_input_initialisation --same-doc a.xlsx --all-doc b.xlsx` | Объединение Excel-шаблонов. |
| `library.utils.cli_tools.mapper_batch_main` | `python -m library.utils.cli_tools.mapper_batch_main --input ids.csv --final-out mapped.csv` | Пакетное сопоставление ChEMBL↔UniProt. |
| `library.utils.cli_tools.mapper_main` | `python -m library.utils.cli_tools.mapper_main --input ids.csv --final-out mapped.csv` | Интерактивный маппер. |
| `library.utils.cli_tools.pipeline_targets_main` | `python -m library.utils.cli_tools.pipeline_targets_main --input targets.csv` | Повторный запуск таргет-харнеса на кэше. |
| `library.utils.cli_tools.table_quality_main` | `python -m library.utils.cli_tools.table_quality_main --input data.csv --table-name data` | Создание профилей качества. |

Мапперы учитывают `io.na_markers` и `io.keep_na_markers` из [`CONFIG.md`](./CONFIG.md).
Каждый модуль предоставляет функцию `main(argv)` для повторного использования в тестах.

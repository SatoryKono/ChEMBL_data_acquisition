# Вспомогательные CLI-модули

> **Версия проекта:** 0.2.0 (2025-10-02)
>
> **Языки:** [English](CLI_TOOLS.md) · [Русский](CLI_TOOLS_RU.md)

Легковесные вспомогательные команды, которые ранее поставлялись как отдельные модули в каталоге `scripts/`, теперь находятся в пакете `library.utils.cli_tools`, чтобы их можно было запускать через `python -m <module>`. Этот перенос упорядочивает дистрибутив и при этом сохраняет возможность прямого запуска модулей для разовой отладки.

| Модуль | Типовая команда | Назначение |
|--------|-----------------|------------|
| `library.utils.cli_tools.check_determinism` | `python -m library.utils.cli_tools.check_determinism` | Проверяет, что детерминированная запись CSV по-прежнему выдаёт совпадающие хэши. |
| `library.utils.cli_tools.chunk_io_main` | `python -m library.utils.cli_tools.chunk_io_main --input data.csv --output copy.csv` | Потоково обрабатывает входные CSV по чанкам с сохранением детерминированного порядка. |
| `library.utils.cli_tools.csv_utils_main` | `python -m library.utils.cli_tools.csv_utils_main --input data.csv` | Пересериализует произвольные CSV с детерминированным порядком строк. |
| `library.utils.cli_tools.dtype_inspector_main` | `python -m library.utils.cli_tools.dtype_inspector_main --log-level INFO` | Анализирует типы pandas, которые возвращают основные пайплайны `get_*_data`. |
| `library.utils.cli_tools.get_activities` | `python -m library.utils.cli_tools.get_activities --limit 10` | Генерирует синтетические строки активностей для проверки логирования и CLI-обвязки. |
| `library.utils.cli_tools.get_document_type` | `python -m library.utils.cli_tools.get_document_type --input docs.csv` | Классифицирует строки документов с использованием встроенных эвристик для модульных тестов. |
| `library.utils.cli_tools.get_input_initialisation` | `python -m library.utils.cli_tools.get_input_initialisation --same-doc path.xlsx --all-doc path.xlsx` | Объединяет Excel-книги, которые готовят пары входных данных для последующего QA. |
| `library.utils.cli_tools.mapper_batch_main` | `python -m library.utils.cli_tools.mapper_batch_main --input ids.csv --output mapped.csv` | Сопоставляет идентификаторы ChEMBL с акцессиями UniProt по пакетной конфигурации. |
| `library.utils.cli_tools.mapper_main` | `python -m library.utils.cli_tools.mapper_main --input ids.csv --output mapped.csv` | Лёгкий интерактивный маппер для быстрых запросов и диагностики. |
| `library.utils.cli_tools.pipeline_targets_main` | `python -m library.utils.cli_tools.pipeline_targets_main --input targets.csv` | Запускает закешированную обвязку таргет-пайплайна и позволяет проверить стадийные флаги (`--raw-out`, `--raw-format`, `--id-cols`, `--no-reindex-raw`, `--normalize-at-export` / `--no-normalize-at-export`). |
| `library.utils.cli_tools.table_quality_main` | `python -m library.utils.cli_tools.table_quality_main --input data.csv --table-name data` | Формирует отчёты о качестве по столбцам для произвольных CSV-наборов. |

Обе утилиты маппинга учитывают список [`io.na_markers`](CONFIG_RU.md#io) при фильтрации placeholder-значений и опцию [`io.keep_na_markers`](CONFIG_RU.md#io), которая определяет, сохранять ли такие идентификаторы во входных данных.

Все модули по-прежнему экспортируют функцию `main`, поэтому их можно подключать к точкам входа в `pyproject.toml`. При программном вызове импортируйте модуль из `library.utils.cli_tools` и передавайте `main(argv)`, чтобы переиспользовать интерфейс командной строки в тестах.

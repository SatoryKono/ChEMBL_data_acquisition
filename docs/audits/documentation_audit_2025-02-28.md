# Отчёт по аудиту документации (28.02.2025)

## Таблица расхождений

| № | Документ | Раздел / файл | Содержание в документации | Реализация в коде | Тип расхождения | Рекомендация |
|---|----------|----------------|---------------------------|-------------------|-----------------|--------------|
| 1 | `docs/ChEMBL Data Acquisition Protocol v0.2.1.md` (устаревший) | §§4.2–4.6 (CLI параметры) | Использовались ключи `--input-csv`, `--output-csv`, отсутствовали `--base-path`, `--override-*`, не упоминались пайплайны tissue/cell line. | Все CLI строятся на `library/cli/parser.add_common_arguments` и ожидают `--input`, `--final-out`, `--base-path`; оркестратор поддерживает `--pipeline-registry` и `--override-*`; существуют `scripts/get_tissue_data.py` и `scripts/get_cellline_data.py`. | Устаревшее описание CLI | Заменить документ новой версией `v2.2` с актуальными параметрами и последовательностью шагов.【F:library/cli/parser.py†L126-L204】【F:library/pipelines/registry.py†L80-L128】 |
| 2 | `docs/en/CLI_TOOLS.md`, `docs/ru/CLI_TOOLS.md`, `docs/en/guides/USAGE.md`, `docs/ru/guides/USAGE.md` | Разделы про `get_activities` | Пример вызова содержал флаг `--output-csv`. | Скрипт использует общие аргументы `--final-out`/`--output` из `library/cli.parser`. | Несоответствие CLI аргументов | Обновить примеры на `--final-out`, чтобы совпадало с `argparse`-интерфейсом.【F:library/utils/cli_tools/get_activities.py†L21-L111】 |
| 3 | `docs/en/DATA_SCHEMA.md`, `docs/ru/DATA_SCHEMA.md` | Описания выходных таблиц | Таблицы колонок не совпадали со схемами: указывались `standard_lower_value`, `ASSAY_ID`, устаревшие поля документов и отсутствовали Pandera-типы. | Актуальные схемы определены в `library/schemas/*.py` и `config/schema/document.yaml`; финальный экспорт удаляет `standard_lower_value`/`standard_upper_value`. | Несоответствие схем данных | Синхронизировать документацию с Pandera-схемами и сослаться на YAML-декларацию документов.【F:library/schemas/activities.py†L31-L83】【F:library/cli/entrypoints/activity.py†L1239-L1360】 |

## Итог по направлениям

### Извлечение и оркестрация
- Протокол обновлён до версии 2.2, перечислены флаги `--base-path`, `--pipeline-registry`, `--override-*` и добавлены упоминания вспомогательных пайплайнов тканевых и клеточных справочников.【F:docs/ChEMBL Data Acquisition Protocol v2.2.md†L33-L118】

### Постобработка и схемы
- Документация по схемам (EN/RU) теперь отображает Pandera-типы и nullable-статусы для `activities`, `assays`, `targets`, `testitems`; для `documents` даны ссылки на `config/schema/document.yaml` и `DOCUMENT_DECLARATION`.【F:docs/en/DATA_SCHEMA.md†L1-L111】【F:docs/ru/DATA_SCHEMA.md†L1-L118】
- Протокол фиксирует последовательность шагов postprocessing и ссылку на `library/postprocessing/common/run_steps`.【F:docs/ChEMBL Data Acquisition Protocol v2.2.md†L141-L172】

### CLI и руководства
- Все ссылки на `get_activities` используют корректный аргумент `--final-out` как в английских, так и в русских гайдах.【F:docs/en/CLI_TOOLS.md†L9-L15】【F:docs/en/guides/USAGE.md†L250-L259】

### Диаграммы и ER-модель
- Протокол содержит актуальные Mermaid-диаграммы для ER-связей и ETL-потока, отражающие зависимость activity от документ/target/assay/testitem измерений.【F:docs/ChEMBL Data Acquisition Protocol v2.2.md†L21-L70】

## Обновлённые и добавленные документы

- `docs/ChEMBL Data Acquisition Protocol v2.2.md`
- `docs/en/CLI_TOOLS.md`, `docs/ru/CLI_TOOLS.md`
- `docs/en/guides/USAGE.md`, `docs/ru/guides/USAGE.md`
- `docs/en/DATA_SCHEMA.md`, `docs/ru/DATA_SCHEMA.md`
- `docs/audits/file_hashes_20250228.txt`
- Настоящий отчёт `docs/audits/documentation_audit_2025-02-28.md`

## Приложения

### A. Контрольные хэши файлов
Файл `docs/audits/file_hashes_20250228.txt` содержит SHA-256 для проверенных артефактов (CLI-скрипты, схемы и обновлённые документы).【F:docs/audits/file_hashes_20250228.txt†L1-L16】

### B. Версионность документации
- Протокол (`v2.2`) включает таблицу версий в §8, фиксируя дату, автора и ключевые изменения.【F:docs/ChEMBL Data Acquisition Protocol v2.2.md†L173-L186】


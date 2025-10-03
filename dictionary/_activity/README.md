# `_activity` reference

> Changelog
> - 2024-05-08 — добавлено описание каталога.

Каталог содержит эталонный экспорт `activity.csv`, который формируется пайплайном `scripts/get_activity_data.py` из ChEMBL API и используется в тестах/QA как эталонная таблица нормализованных активностей. Структура колонок и описание полей задокументированы в [docs/en/reference/DATA_SCHEMA.md](../../docs/en/reference/DATA_SCHEMA.md#activitycsv-processed-export).

## Файлы
| Имя | Назначение |
| --- | --- |
| `activity.csv` | Нормализованные активности с расширенными границами, JSON-представлением свойства `activity_properties`, флагами действий и метаданными пайплайна. |

## Источник и обновление
1. Сформируйте входной список идентификаторов активностей (см. раздел «Input Tables» в [DATA_SCHEMA.md](../../docs/en/reference/DATA_SCHEMA.md#input-tables)).
2. Запустите `python -m scripts.get_activity_data chembl --input <path>/activity.csv --output <path>/activity.csv` с параметрами окружения, описанными в [docs/en/guides/README.md](../../docs/en/guides/README.md#scriptsget_activity_datapy).
3. Проверьте качество через `python -m library.utils.cli_tools.table_quality_main --input <path>/activity.csv --table-name chembl_activity` и убедитесь, что схема совпадает с документацией.
4. Обновите `activity.csv` и отразите изменения в этом README (дата загрузки, источники, отклонения от схемы).

## Проверка качества
- Сравнивайте количество строк и хэши `properties_hash` между версиями, чтобы отследить неожиданные изменения.
- Перед коммитом прогоняйте unit-тесты, использующие экспорт активностей, чтобы убедиться в обратной совместимости.

# Зависимости контроля качества

Таблица показывает связь между проверками, конфигурацией и этапами пайплайна.

| Проверка | Где реализована | Зависимости |
|----------|-----------------|-------------|
| Входная схема Pandera | `library/schemas/normalize.py` | Выполняется до обращений к API; зависит от CLI (`--input`) и `local.io`. |
| Выходная схема Pandera | `library/schemas/{documents,targets,assays,activities,testitems}.py` | Запускается после обогащения; учитывает конфигурацию (например, `activity_enrichment`). |
| Профили качества | `library/table_quality.py` | Управляются `system.doc_quality` (enable, sample_rows, include/exclude). |
| Allowlist action type | `library/pipelines/activity/enrichment.py` | Настраивается через `activity_enrichment.action_type`. |
| Границы активностей | `library/pipelines/activity/bounds.py` | Зависит от блока `activity_bounds`. |
| Родительские молекулы | `library/pipelines/testitem/enrichment.py` | Требуются словари из `config/dictionary/_testitem`. |
| Проверка детерминизма | `scripts/check_determinism.py` | Опциональна; сравнивает с предыдущей директорией. |

## Поведение при сбоях

- Нарушение схемы вызывает `pandera.errors.SchemaError` и завершение со статусом 1.
- Сбой профиля качества логирует `quality_report_failed`; при `fatal_on_error=true`
  конвейер останавливается, иначе ошибка фиксируется и выполнение продолжается.
- Несоответствия в словарях формируют предупреждения; при необходимости усиливайте
  логирование параметрами `testitem_molecule_enrichment.logging`.

## Рекомендации по мониторингу

- Собирайте `.quality.json` в CI и строите отчёты по трендам.
- Настраивайте алерты на падение `non_empty_ratio` для ключевых колонок (DOI,
  `action_type`, `standard_value`).
- Отслеживайте количество `action_type=unknown` и отрицательных `standard_value`
  после обогащения.

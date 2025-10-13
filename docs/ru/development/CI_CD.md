# Рекомендации по CI/CD

## Рекомендуемые джобы

| Этап | Команда | Назначение |
|------|---------|------------|
| Lint | `make lint` | Проверка формата (`black`) и стиля (`ruff`). |
| Type check | `make typecheck` | `mypy` для `library/` и `scripts/`. |
| Tests | `pytest -q --json-report --json-report-file=reports/test_report.json` | Полный прогон тестов. |
| Report summary | `make-md-summary --input reports/test_report.json --output reports/test_summary.md` | Генерация Markdown-отчёта (пути по умолчанию — `reports/`). |
| Smoke (опция) | `make smoke` | Запуск оркестратора на тестовых данных. |
| Determinism (опция) | `check-determinism --baseline <prev> --candidate <curr>` | Проверка изменений CSV между прогоном. |

## Контроль качества

- Успешность тестов ≥ 95 % (`summary.success_rate` из JSON).
- Отсутствие ошибок lint/typecheck.
- При необходимости падение пайплайна при несоответствии `check-determinism`.

## Артефакты

- Загружайте `reports/test_report.json` и `reports/test_summary.md`.
- При smoke-запуске архивируйте выходные CSV и метаданные.
- Сохраняйте логи `data/logs/*.log` (по умолчанию) или `<base>/logs`, если задана
  `CHEMBL_DA_BASE_PATH`, для диагностики.

## Ветвление

- Фичи: `feat/<name>`
- Исправления: `fix/<name>`
- Документация: `docs/<name>`

Рекомендуется linear history (rebase/FF).

## Релизы

- Тэг `v<major>.<minor>.<patch>`.
- Обновляйте `library/version.py` и `pyproject.toml` одновременно.
- При изменении схем словарей пересобирайте кэши и обновляйте
  [`reference/DICTIONARIES.md`](../reference/DICTIONARIES.md).

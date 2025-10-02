# Release Notes / Список релизов

> **Project version / Версия проекта:** 0.2.0 (2025-10-02)
>
> Detailed history lives in [`CHANGELOG.md`](../CHANGELOG.md). The sections
> below highlight the key points for each tagged release.

## 0.2.0 — 2025-10-02

### Highlights / Основные изменения
- Синхронизированы английская и русская версии руководств, добавлены версии в
  заголовки и обновлены примеры запуска CLI.
- README переработан в двуязычный портал с картой документации и блоком по
  управлению релизами.
- Уточнено, что набор ключей `--raw-out` доступен в полной мере только для
  пайплайна таргетов; остальные команды игнорируют их до расширения функционала.
- Проект повышен до версии `0.2.0`, changelog переведён в формат Keep a Changelog
  с двумя языками.

### Testing focus / Тестовый контур
- `pre-commit run --all-files`
- `pytest`
- `python -m library.utils.cli_tools.check_determinism --limit 10`

## 0.1.0 — (date not recorded / дата не зафиксирована)

### Highlights / Основные изменения
- Initial public release with activity, assay, document, target and test item
  pipelines, configuration framework and QA harness.

### Testing focus / Тестовый контур
- Baseline pytest suite and deterministic IO smoke tests.

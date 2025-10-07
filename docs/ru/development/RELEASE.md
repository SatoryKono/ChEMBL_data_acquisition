# Чек-лист релиза

1. **Зависимости** (при необходимости)
   - Обновите диапазоны версий в `pyproject.toml`.
   - Пересоберите `requirements-lock.txt` в чистом окружении.

2. **Словари**
   - Скачайте свежие CSV/JSON.
   - Запустите `make smoke`, убедитесь, что QC проходит.
   - Обновите [`reference/DICTIONARIES.md`](../reference/DICTIONARIES.md).

3. **Версия**
   - Измените `library/version.py` и `pyproject.toml`.
   - Добавьте заметки в `docs/en/RELEASE_NOTES.md` и переведите в RU.

4. **Качество**
   - `make lint && make typecheck`
   - `pytest -q --json-report --json-report-file=reports/test_report.json`
   - `python tools/make_md_summary.py --input reports/test_report.json --output reports/test_summary.md`
     (аргументы можно опустить; консольный скрипт: `make-md-summary`)
   - `make smoke`
   - `check-determinism` против предыдущих выгрузок (если доступны).

5. **Тег и публикация**
   - Мёрдж через fast-forward.
   - `git tag vX.Y.Z && git push origin vX.Y.Z`
   - Опубликуйте релиз с ссылками на QA-артефакты.

6. **После релиза**
   - Архивируйте отчёты QA и smoke-выгрузки.
   - Создайте задачи по найденным проблемам.

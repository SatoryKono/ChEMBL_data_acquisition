# Процесс обеспечения качества (живой документ)

Этот документ — каноничный чек-лист для валидации стека выгрузки данных ChEMBL. Обновляйте его при добавлении новых сервисов,
конвейеров или политик, чтобы другим гайдам не приходилось дублировать инструкции, а достаточно было сослаться сюда.

## 1. Среда

1. Установите зафиксированные зависимости:
   ```bash
   pip install -r requirements-lock.txt
   ```
2. При необходимости добавьте пакет в editable-режиме для локальных entry points:
   ```bash
   pip install -e .
   ```
3. Экспортируйте `PYTHONPATH=.` — так вспомогательные скрипты и детерминированные писатели будут одинаково находить пакетные импорты.
4. Убедитесь, что опциональные зависимости (``responses``, ``hypothesis``, ``psutil``, ``pytest-benchmark``) установлены там, где вы
   планируете запускать соответствующие тестовые наборы; иначе они будут пропущены.

## 2. Статический анализ и форматирование

Запускайте следующие команды из корня репозитория. Они должны завершаться без диагностик.

```bash
ruff check .
ruff format --check .
mypy --strict
```

## 3. Запуск тестов

1. Быстрый сигнал о регрессиях:
   ```bash
   pytest --maxfail=1 --durations=10
   ```
2. Полный набор с видимыми предупреждениями для триажа.
   В повседневной сертификации используйте тихий режим (`-q`).
   Переключайтесь в подробный (`-vv`) вывод, когда нужна детальная информация об ошибках:
   ```bash
   pytest -q --disable-warnings
   ```

## 4. Проверка детерминизма и smoke-тесты CLI

1. Подтвердите детерминированный рендеринг CSV:
   ```bash
   PYTHONPATH=. python -m library.utils.cli_tools.check_determinism --log-level DEBUG
   ```
2. Прогоните хотя бы один конвейерный CLI в режиме dry-run, чтобы убедиться в целостности аргументов. Например:
   ```bash
   PYTHONHASHSEED=0 PYTHONPATH=. python scripts/get_activity_data.py --input tests/data/activity_ids_small.csv \
       --final-out /tmp/activities.csv --limit 10 --dry-run --log-level INFO
   ```
   При необходимости замените сценарий на другие конвейеры, чтобы покрыть недавние изменения.

## 5. Отчётность

Фиксируйте вывод команд (статус выполнения, число ошибок, временные метки) в журнале аудита — обычно `docs/code_review.md` — при
каждой повторной сертификации репозитория. Делитесь итогами через ссылки на этот живой документ, а не копированием чек-листа.

## 6. Контроль качества пост-обработки документов

Экспорт документов теперь сопровождается автоматической регрессионной проверкой относительно легаси-книги Power Query.

1. Заполните `data/input/full/document.csv` эталонной выгрузкой.
2. Запустите конвейер документов (`python -m scripts.get_document_data ...`) и убедитесь, что он создаёт:
   * `output.document_YYYYMMDD.csv`
   * `qa_document_postprocessing_report_YYYYMMDD.json`
   * `qa_document_postprocessing_report_YYYYMMDD.md`
   * `qa_document_postprocessing_diff_YYYYMMDD.csv` (только при наличии расхождений)
3. Либо выполните QA-скрипт напрямую:
   ```bash
   python -m qa.check_document_postprocessing \
       --base-path data \

       --ref input\\full\\document.csv \
       --actual output\\document\\preprocessed_output.document_YYYYMMDD.csv

       --out output\\document\\output.document_YYYYMMDD.csv

   ```
4. Рассматривайте ненулевой код возврата как блокирующую ошибку; для устранения смотрите Markdown-резюме и CSV с диффами.

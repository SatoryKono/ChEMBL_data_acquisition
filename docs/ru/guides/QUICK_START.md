# Быстрый старт

Краткая инструкция по локальному запуску пайплайна на тестовых данных из репозитория.

## 1. Подготовка окружения

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-lock.txt
pre-commit install
```

Дополнительно:

- `pip install -r requirements-dev.txt` — для линтеров и mypy.
- `export PYTHONPATH=$(pwd)` — при запуске модулей напрямую.

## 2. Проверка конфигурации

Основной файл — `config/config.yaml`. Для локальных переопределений создайте
`config/config.local.yaml`:

```yaml
sources:
  chembl:
    api:
      user_agent: "ChEMBL-ETL/2.1 (mailto:team@example.org)"
```

Проверить итоговые значения можно командой:

```bash
python scripts/get_document_data.py --mode chembl --print-config | less
```

## 3. Smoke-запуск

В каталоге `data/input` лежат минимальные CSV для всех сущностей. Запускайте
конвейеры по отдельности или через оркестратор:

```bash
# Отдельный конвейер
python scripts/get_document_data.py --mode all \
  --input data/input/document.csv \
  --final-out output/documents.csv

# Полная цепочка
poetry run get-data \
  --base-path . \
  --input-dir data/input \
  --output-dir output \
  --config config/config.yaml \
  --date $(date -u +%Y%m%d)
```

Готовые файлы окажутся в `output/`. Проверьте `<name>.meta.yaml` для сведений о
версии и параметрах запуска.

## 4. Запуск тестов

```bash
pytest -q --disable-warnings
pytest -q --json-report --json-report-file=reports/test_report.json
python tools/make_md_summary.py --input reports/test_report.json --output reports/test_summary.md
# аргументы можно опустить при использовании путей по умолчанию:
# python tools/make_md_summary.py
# make-md-summary
```

В `reports/test_summary.md` ожидается ≥95 % успешных тестов. Артефакты из `reports/`
нужно прикладывать к CI/issue.

## 5. Проверка детерминизма (опционально)

```bash
poetry run get-data --output-dir output/run1
poetry run get-data --output-dir output/run2
check-determinism --baseline output/run1 --candidate output/run2
```

Несовпадения хэшей сигнализируют о недетерминированном поведении или изменениях в
конфигурации. Для диагностики воспользуйтесь [`QA_PROCESS.md`](../QA_PROCESS.md).

# ChEMBL Data Acquisition Utilities

Набор ETL-инструментов на Python для получения, нормализации и экспорта
данных ChEMBL и сопутствующих источников (PubMed, UniProt, OpenAlex,
CrossRef, IUPHAR и др.). Скрипты ориентированы на воспроизводимый сбор
биологической информации: от выгрузки идентификаторов до подготовки
аналитических витрин в формате CSV/Parquet с метаданными и контролем
качества.

## Цели проекта

* Централизованный доступ к ChEMBL API и вспомогательным сервисам с
  уважением ограничений RPS и политик идентификации.
* Единая конфигурация для всех пайплайнов и CLI-утилит через `config.yaml`
  и переменные окружения.
* Проверяемые схемы Pandera и pydantic для детерминированной структуры
  выгрузок и sidecar-файлов.
* Построение витрин активности, биологических тестов, документов,
  таргетов и тест-айтемов с дополнительными справочниками.

## Структура репозитория

| Директория | Назначение |
|------------|------------|
| `scripts/` | CLI-команды `get_*_data.py`, утилиты проверки качества, конвертации, буферизации чанков. Каждый скрипт использует единый каркас `library/cli` и конфигурацию `config.yaml`. |
| `library/` | Повторно используемые модули: клиенты API (`chembl_client.py`, `pubmed_library.py`, `uniprot_library.py`), обёртки для ввода-вывода (`io.py`, `chunk_io.py`), нормализация и постобработка (`activity_extraction.py`, `document_postprocessing.py`, `target_postprocessing.py`), вспомогательные структуры и схемы. |
| `dictionary/` | Локальные справочники и заготовки (IUPHAR, классификаторы таргетов, словари документов, кэш UniProt). Служат источником обогащения и валидации данных. |
| `schemas/` | Pandera-схемы для основных таблиц (`activities.py`, `assays.py`, `documents.py`, `targets.py`, `testitems.py`) и модели метаданных (`meta.py`). |
| `tests/data/` | Минимальные тестовые наборы идентификаторов, CSV и конфигураций для локальных прогонов и unit-тестов. |

## Установка и запуск

1. **Подготовьте окружение**

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
   python -m pip install --upgrade pip setuptools wheel
   ```

2. **Установите зависимости проекта**

   ```bash
   pip install -e .[dev]
   ```

   Пакет включает `pandas`, `requests`, `PyYAML`, `pandera`, а также инструменты
   разработчика (`black`, `ruff`, `mypy`, `pytest`).

3. **Сконфигурируйте доступ**

   * Скопируйте `config.yaml` и адаптируйте под свои пути/лимиты (см.
     [CONFIG.md](CONFIG.md)).
   * Чувствительные параметры (токены, логины, ключи) храните в `.env` или
     переменных окружения `CHEMBL_DA__SECTION__KEY`.

4. **Запустите нужный скрипт**

   ```bash
   python -m scripts.get_activity_data \
     --input tests/data/activity_ids_small.csv \
     --output data/output/activities.csv \
     --config config.yaml \
     --log-level INFO
   ```

   Аналогично вызываются `get_assay_data.py`, `get_document_data.py`,
   `get_target_data.py`, `get_testitem_data.py`. Любой скрипт принимает
   флаги `--input`, `--output`, `--config`, `--limit`, `--log-level`, а
   также специфичные параметры (например, `--document-source` для
   публикаций).

5. **Проверьте результаты**

   Скрипты записывают основную таблицу, YAML-sidecar с метаданными и, при
   наличии ошибок валидации, CSV с описанием отказов (`*_failure_cases.csv`).

## Тестирование и покрытие

Для проверки регресса запускайте тесты из корня проекта. Pytest настроен через
`pyproject.toml` и автоматически подхватывает каталоги `tests/` и `library/`:

```bash
pytest
```

Чтобы собрать отчёт по покрытию библиотечных модулей и CLI-скриптов, выполните
последовательность команд:

```bash
coverage run -m pytest
coverage report
```

Полученный отчёт позволяет отследить узкие места. На текущей ревизии
наиболее низкое покрытие отмечено у `library/chembl_target.py`,
`library/iuphar_library.py` и `scripts/get_document_data.py`, поэтому при
доработках этих частей рекомендуется добавлять дополнительные тесты.

## Общая схема ETL

```
Идентификаторы → Сырые API-ответы → Нормализация/слияния → Валидация → Выходные таблицы + sidecar
```

1. **Источники данных**
   * Входные CSV/XLSX с идентификаторами (`data/input`, `tests/data`).
   * Внешние API: ChEMBL, PubMed, UniProt (REST и ID Mapping), IUPHAR,
     OpenAlex, CrossRef, Semantic Scholar, PubChem.

2. **Преобразования**
   * Загрузка чанками с управлением rate limit (token bucket, глобальные
     ограничения `rate.*`).
   * Обогащение дополнительными справочниками (`dictionary/`).
   * Нормализация и расчёт derived-полей (`library/normalization.py`,
     `document_postprocessing.py`, `target_postprocessing.py`).
   * Валидация Pandera и анализ качества (`library/table_quality.py`).

3. **Выход**
   * Таблицы активности, ассев, документов, таргетов, тест-айтемов.
   * YAML sidecar с метаданными генерации (`schemas/meta.py`).
   * Отчёты контроля качества, файлы ошибок и журналы выполнения.

## Дополнительно

* Для локальных проверок используйте `pytest` и `pre-commit run --all-files`.
* Описание конфигурации — в [CONFIG.md](CONFIG.md).
* Структура выходных таблиц — в [OUTPUT.md](OUTPUT.md).
* Пошаговая инструкция запуска — в [USAGE.md](USAGE.md).

# Практическое использование

Документ описывает настройку `config.yaml`, запуск основных скриптов и работу с
тестовыми данными. Предполагается, что зависимости уже установлены (см.
[README.md](README.md)).

## Шаг 1. Подготовка `config.yaml`

1. Скопируйте базовый файл:
   ```bash
   cp config.yaml config.local.yaml
   ```
2. Отредактируйте ключевые блоки:
   * `api.user_agent` — укажите e-mail для идентификации в ChEMBL/OpenAlex/CrossRef.
   * `io.output_dir` и `io.cache_dir` — директории для результатов и кэша.
   * `resources.*` — пути к локальным словарям, если структура каталога изменилась.
   * При необходимости обновите таймауты и лимиты (`api.rps`, `rate.global_rps`).
3. Секреты (токены, ключи) загрузите в окружение:
   ```bash
   export CHEMBL_DA__API__USER_AGENT="chembl-da/1.0 (mailto:you@org.org)"
   ```
4. Перед запуском убедитесь, что все каталоги существуют или включен флаг
   `io.exist_ok: true` (по умолчанию так).

## Шаг 2. Запуск скриптов

Общий шаблон команд:

```bash
python -m scripts.<module> \
  --input <путь-к-входу> \
  --output <путь-к-выходу> \
  --config config.local.yaml \
  --log-level INFO \
  [дополнительные флаги]
```

### Activity

*Получение активностей по списку `activity_chembl_id`.*

```bash
python -m scripts.get_activity_data \
  --input tests/data/activity_ids_small.csv \
  --output data/output/activities.csv \
  --config config.local.yaml \
  --limit 50 \
  --log-level INFO
```

Результаты: `activities.csv`, `activities.meta.yaml`, при ошибках —
`activities_failure_cases.csv`.

### Assay

*Загрузка ассеев по `assay_chembl_id`.*

```bash
python -m scripts.get_assay_data \
  --input data/input/assays.csv \
  --output data/output/assays.csv \
  --config config.local.yaml \
  --log-level DEBUG
```

Дополнительно можно запустить анализ качества:

```bash
python -m scripts.table_quality_main \
  --input data/output/assays.csv \
  --output data/output/assays_quality.csv
```

### Document

*Агрегация публикаций (ChEMBL + внешние источники).* Укажите источник через
флаг `--document-source` (`chembl`, `pubmed`, `all`). Пример полного конвейера:

```bash
python -m scripts.get_document_data \
  --document-source all \
  --input tests/data/documents_postprocess.csv \
  --output data/output/documents.csv \
  --config config.local.yaml \
  --log-level INFO
```

Выход: `documents.csv`, sidecar и файлы с ошибками/классификацией.

### Target

*Комбинированная витрина таргетов.* При запуске укажите режим (`chembl`,
`uniprot`, `iuphar`, `all`). Пример для объединённого режима:

```bash
python -m scripts.get_target_data \
  --target-source all \
  --input tests/data/chembl_targets_min.csv \
  --output data/output/targets.csv \
  --config config.local.yaml \
  --log-level INFO
```

Скрипт может сформировать дополнительные файлы (`targets_chembl.csv` и т.п.)
если соответствующие пути заданы в `config.target.all`.

### Test item

*Обогащение молекул данными ChEMBL/PubChem.*

```bash
python -m scripts.get_testitem_data \
  --input tests/data/iuphar_targets_min.csv \
  --output data/output/testitems.csv \
  --config config.local.yaml \
  --log-level INFO
```

## Где искать результаты

* Основные CSV — в `io.output_dir` (по умолчанию `data/output`).
* YAML sidecar и отчёты валидации создаются в той же директории, что и основная
  таблица.
* Журналы (`*.log`) можно перенаправить через флаг `--log-file` (см. `library/cli`).

## Тестовые данные

Каталог `tests/data/` содержит минимальные примеры для локальной отладки:

| Файл | Назначение |
|------|-----------|
| `activity_ids_small.csv` | Список активностей для smoke-теста `get_activity_data`. |
| `activities_valid.csv` / `activities_invalid.json` | Примеры корректных и ошибочных структур для модулей валидации. |
| `chembl_targets_min.csv`, `uniprot_targets_min.csv`, `iuphar_targets_min.csv` | Срезы справочников таргетов. |
| `documents_postprocess.csv` | Заготовка публикаций для отработки постобработки. |
| `pmids.csv` | Минимальный набор PubMed ID. |
| `sample_config.yaml` | Пример конфигурации для модульных тестов. |
| `io_types.csv`, `csv_utils_input.csv`, `meta_input.csv` | Таблицы для тестов ввода-вывода и генерации метаданных. |

Используйте эти файлы для быстрой проверки инфраструктуры без обращения к
производственным API.

## Типичные сценарии

1. **Быстрый smoke-тест:**
   ```bash
   python -m scripts.get_activity_data --input tests/data/activity_ids_small.csv --limit 5 --dry-run --log-level DEBUG
   ```
   Проверяет чтение входа и конфигурацию без сетевых запросов.

2. **Обновление справочника таргетов:**
   ```bash
   python -m scripts.get_target_data --target-source all --input tests/data/chembl_targets_min.csv --output data/output/targets.csv
   ```
   После выполнения проверьте `targets.meta.yaml` и запустите `check_determinism.py`.

3. **Обновление публикаций:**
   ```bash
   python -m scripts.get_document_data --document-source all --input tests/data/documents_postprocess.csv --output data/output/documents.csv
   python -m scripts.table_quality_main --input data/output/documents.csv --output data/output/documents_quality.csv
   ```

4. **Генерация тест-айтемов и отчёт по качеству:**
   ```bash
   python -m scripts.get_testitem_data --input tests/data/iuphar_targets_min.csv --output data/output/testitems.csv
   python -m scripts.table_quality_main --input data/output/testitems.csv --output data/output/testitems_quality.csv
   ```

## Отладка и контроль качества

* Флаг `--limit` помогает сократить время тестового прогона.
* При проблемах с сетью увеличьте `timeout` в соответствующем блоке конфигурации.
* Команда `python -m scripts.check_determinism --log-level DEBUG` проверяет
  совпадение свежих выгрузок с эталонами в `tests/data/golden/`.
* Используйте `pytest` для проверки модулей: `pytest -k activity` и т.д.

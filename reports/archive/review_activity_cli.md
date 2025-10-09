# Отчёт код-ревью: ChEMBL activity pipeline

## 1. Таблица замечаний

| № | Компонент | Проблема / Риск | Причина | Рекомендация / Исправление | Критичность |
|---|-----------|-----------------|---------|----------------------------|-------------|
| 1 | `library/cli/entrypoints/activity.py` | Бесконечное ожидание при кешировании `pref_name` блокирует пайплайн | Цикл ожидания после отметки `_CACHE_IN_PROGRESS` крутится с `sleep(0)` без таймаута; при зависшем воркере или исключении значение остаётся `IN_PROGRESS`, CPU нагружается и поток не выходит | Заменить busy-wait на `Condition`/`Event` с таймаутом и форсировать сброс флага после ошибок; добавить журналирование и фолбек на `None` после истечения ожидания | Высокая【F:library/cli/entrypoints/activity.py†L760-L819】 |
| 2 | `library/pipelines/activity/runner.py` | Запуск мутирует конфигурацию `cfg.activity` и влияет на последующие шаги | При клэмпе batch size/timeout значения записываются обратно в объект конфига, который разделяется между пайплайнами | Работать с копией (`cfg.model_copy(deep=True)`), возвращать новые значения через `ActivityCommandOptions` и не менять глобальное состояние | Высокая【F:library/pipelines/activity/runner.py†L150-L184】 |
| 3 | `library/cli/parser.py` | Нестабильный `run_id` ломает детерминизм логов и отчётов | `create_logger_config` генерирует UUID на каждый запуск без возможности переопределить | Принимать `seed`/`run_id` из CLI/ENV, по умолчанию детерминировать (например, хэш входного файла) | Высокая【F:library/cli/parser.py†L28-L42】 |
| 4 | `library/common/fetch_retry.py` | Статистика ошибок хранит полный список ID и раздувает sidecar | `ChunkFailureTracker.stats` собирает и сортирует все идентификаторы без ограничений | Ограничить объём (например, первые 100 ID + счётчики), остальное писать в отдельный файл/sidecar | Средняя【F:library/common/fetch_retry.py†L22-L58】 |
| 5 | `scripts/get_activity_data.py`, `library/cli/entrypoints/activity.py` | Полное дублирование логики CLI → риск рассинхронизации | Скрипт и entrypoint содержат одинаковые ~1600 строк | Удалить дубликат, оставить один источник истины (скрипт вызывает модуль или наоборот) | Средняя【F:scripts/get_activity_data.py†L946-L1195】【F:library/cli/entrypoints/activity.py†L942-L1195】 |
| 6 | `library/cli/entrypoints/activity.py` | Итоговое сообщение без имени события осложняет мониторинг | `_emit_completion_message` логирует человекочитаемую строку без структурированных полей | Логировать `logger.info("pipeline_done", rows=..., duration=..., output=...)`, текст при необходимости дублировать в stdout | Средняя【F:library/cli/entrypoints/activity.py†L468-L499】 |
| 7 | `library/cli/entrypoints/activity.py` | Потеря значений `pref_name=""` и `0` при обогащении | Фильтр `if value and ...` отбрасывает валидные «ложные» значения | Явно проверять `value is not None` и приводить к строке, не теряя пустые строки | Средняя【F:library/cli/entrypoints/activity.py†L803-L808】 |

## 2. План исправлений

### Структура и поддерживаемость
- **Файлы:** `scripts/get_activity_data.py`, `library/cli/entrypoints/activity.py`
- **Действие:** Рефакторинг — удалить дублирование, выделить общую реализацию в библиотечный модуль.
- **Объём:** 1 PR.
- **Пример commit message:** `refactor: de-duplicate activity CLI entrypoint`.

### Надёжность/ETL логика
- **Файлы:** `library/cli/entrypoints/activity.py`
- **Действие:**
  - заменить busy-wait при ожидании кеша на синхронизацию с таймаутом;
  - корректно обрабатывать falsy `pref_name`.
- **Объём:** 1 PR.
- **Пример commit message:** `fix: stabilise molecule pref_name enrichment concurrency`.

### Конфигурация и детерминизм
- **Файлы:** `library/pipelines/activity/runner.py`, `library/cli/parser.py`.
- **Действие:**
  - перестать мутировать `cfg.activity`, использовать копии;
  - добавить параметризированный `run_id`, сделать детерминированный дефолт.
- **Объём:** 1 PR.
- **Пример commit message:** `feat: make activity pipeline configuration immutable and deterministic`.

### Наблюдаемость и отчётность
- **Файлы:** `library/cli/entrypoints/activity.py`, `library/common/fetch_retry.py`.
- **Действие:**
  - структурировать completion-лог;
  - ограничить размер статистики chunk failures.
- **Объём:** 1 PR.
- **Пример commit message:** `chore: harden pipeline observability and failure reporting`.

## 3. Черновой план PR-пакетов

1. **PR-01_refactor_cli** — унификация CLI и логов (`scripts/get_activity_data.py`, `library/cli/entrypoints/activity.py`).
2. **PR-02_activity_reliability** — фиксы кеша `pref_name`, корректная обработка пустых значений.
3. **PR-03_config_determinism** — иммутабельная конфигурация и детерминированный `run_id`.
4. **PR-04_observability** — структурный лог, ограничение `chunk_fetch_failure_ids`.

## 4. Контрольные проверки перед merge

- `pytest` — все тесты проходят.
- `mypy` — без ошибок.
- `ruff check` — без предупреждений.
- `python -m scripts.get_activity_data --help` — CLI парсинг и help работают.
- Сравнение выходных CSV с golden-файлами — детерминированный результат.

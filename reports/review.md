Executive summary
- [Errors] CLI skip-existing ветка логирует событие через глобальный logger из `library.common.log`, поэтому заглушка в e2e тестах не видит `pipeline_skip_existing` и сценарий падает. → В бою не сработает контроль idempotency при повторных запусках. → Прокинуть в `run_activity_pipeline` экземпляр logger из `scripts.get_activity_data` или пробрасывать `LoggerConfig` до уровня модуля и логировать через него.【F:library/cli/commands/get_activity_data.py†L181-L194】【F:tests/e2e/test_get_cli_pipelines.py†L1212-L1234】
- [Perf] `scripts/check_determinism.py` последовательно запускает pipeline и получает разные SHA256 даже на одинаковых входах. → CSV-выгрузки недетерминированы, что ломает reproducibility и проверку регрессий. → Устранить источники случайности (фиксированный run_id, сортировка выходов, замена `datetime.now` в sidecar) и покрыть регресс-тестом.【F:scripts/check_determinism.py†L114-L137】【99cce3†L1-L3】
- [Quality] Lint/typing слой фактически не работает: `ruff check` находит 890 нарушений, `mypy --strict` — 403 ошибки. → Без статанализа накапливаются дефекты и деградирует поддерживаемость. → Почистить импорты, добавить аннотации, завести pre-commit с `ruff`/`mypy` и сделать прогон обязательным в CI.【791a9b†L1-L118】【18648e†L1-L38】
- [Config] В `run_activity_pipeline` значения таймаута и batch size клэмпятся через прямую мутацию `cfg.activity`, из-за чего одиночный запуск меняет глобальную конфигурацию для следующих шагов. → Соседние pipeline получают неожиданно изменённые лимиты. → Работать с `cfg.model_copy(deep=True)` и возвращать изменённую копию вместо модификации оригинала.【F:library/cli/commands/get_activity_data.py†L125-L155】
- [Config] `create_logger_config` генерирует случайный `run_id` на каждый запуск CLI. → Run manifest и логи становятся недетерминированными, что блокирует сравнение результатов и ломает отчёты. → Разрешить передавать фиксированный `run_id` через конфиг/ENV, дефолт — детерминированный (например, от input hash).【F:library/cli/parser.py†L29-L44】
- [Quality] `_emit_completion_message` пишет человекочитаемую строку без имени события, тогда как structured-логер ожидал `event`. → Метрики и алерты не ловят завершение пайплайна, появляются ложные срабатывания. → Поменять на `logger.info("pipeline_done", rows=..., duration=..., output=...)` и оставить человекочитаемый текст в stdout при необходимости.【F:scripts/get_activity_data.py†L445-L466】
- [Structure] Скрипт `scripts/get_activity_data.py` и модуль `library/cli/commands/get_activity_data.py` содержат дублирующийся код подготовки `Namespace`, клампов и финализации. → Любое изменение нужно вносить дважды, легко рассинхронизировать поведение CLI и библиотечного API. → Свести реализацию к одному источнику (скрипт вызывает модуль или наоборот) и удалить дубли.【F:scripts/get_activity_data.py†L1560-L1660】【F:library/cli/commands/get_activity_data.py†L90-L194】
- [Errors] `ChunkFailureTracker.stats` сериализует полный список проваленных ID, что при крупных батчах приводит к огромным sidecar-ам и OOM при переработке. → Контроль качества становится неработоспособным на больших данных. → Хранить только первые N идентификаторов и счётчики, полный список писать в sidecar-файл на диск.【F:library/common/fetch_retry.py†L43-L58】
- [Docs/Observability] Metadata sidecar записывает `generated_at` текущим временем, но нигде не фиксируется старт пайплайна и версия конфига, выводится исходный словарь без маскировки вложенных секретов. → Трудно сопоставлять артефакты и есть риск утечки токенов из вложенных структур. → Добавить поля `run_id`/`started_at`, прогонять `_mask_secrets` по всему дереву и документировать формат `.meta.yaml`.【F:library/common/metadata.py†L137-L200】
- [Testing] Покрытие критичных сценариев ограничено: e2e падает, а property-based проверки schema/pandera отсутствуют. → Легко пропустить регрессии в трансформациях и нарушить схемы. → Добавить pandera-валидаторы для основных DataFrame и property-based тесты для мапперов (activities/assays/documents/targets/testitems).【F:tests/e2e/test_get_cli_pipelines.py†L1212-L1234】【F:library/postprocessing/pipeline/activities/schema.py†L1-L120】

Scores (0–5)
Structure: 2/5 — Дублируется бизнес-логика между `scripts/*` и `library/cli/commands/*`, что усложняет сопровождение и вносит рассинхрон (см. раздел Structure).【F:scripts/get_activity_data.py†L1560-L1660】
Config: 2/5 — Конфигурация мутируется на лету (клэмпы таймаутов) и генерирует случайные идентификаторы, что противоречит политике детерминизма.【F:library/cli/commands/get_activity_data.py†L125-L155】【F:library/cli/parser.py†L29-L44】
Quality: 1/5 — Массовые нарушения lint/typing показывают отсутствие работающего статанализа.【791a9b†L1-L118】【18648e†L1-L38】
Errors: 2/5 — Ключевые события (`pipeline_skip_existing`) теряются, sidecar-ы могут переполняться списками ID, обработка ошибок нестабильна.【F:library/cli/commands/get_activity_data.py†L181-L194】【F:library/common/fetch_retry.py†L43-L58】
Perf: 2/5 — Невозможность гарантировать детерминизм CSV и потенциально взрывающиеся sidecar-ы негативно влияют на SLA и повторные запуски.【F:scripts/check_determinism.py†L114-L137】【99cce3†L1-L3】
Testing: 2/5 — E2E сценарий падает, отсутствуют property-based тесты и детальные схемы pandera для основных пайплайнов.【F:tests/e2e/test_get_cli_pipelines.py†L1212-L1234】
Docs: 3/5 — Базовая документация есть, но .meta.yaml не описан и не гарантирует защищённость секретов.【F:library/common/metadata.py†L137-L200】

Findings by category

Structure
- Проблема: Дублирование логики между `scripts/get_activity_data.py` и `library/cli/commands/get_activity_data.py`.
  - Пример: подготовка `Namespace`, клэмпы batch/timeout и финализация реализованы в двух местах.【F:scripts/get_activity_data.py†L1560-L1660】【F:library/cli/commands/get_activity_data.py†L90-L194】
  - Почему важно: расхождения в этих путях приводят к разному поведению CLI и библиотечных сценариев, что усложняет тестирование и сопровождение.
  - Исправление: вынести общую реализацию в один модуль (например, оставить `library/cli/commands/...` единственным исполнителем, а скрипту оставить только `main()`).
- Проблема: `ChunkFailureTracker.stats` возвращает весь список идентификаторов ошибок.
  - Пример: метод строит `unique_ids = sorted({...})` без ограничений по размеру.【F:library/common/fetch_retry.py†L43-L58】
  - Почему важно: большие выгрузки создают гигантские sidecar-файлы и потребление памяти при сериализации.
  - Исправление: ограничить список первыми N значениями и остальное писать в отдельный файл.

Config
- Проблема: Мутация конфигурации в `run_activity_pipeline` через прямое изменение `cfg.activity.*`.
  - Пример: `cfg.activity.batch_size = MAX_ACTIVITY_CHUNK_SIZE` и `cfg.activity.timeout = minimum_timeout` внутри функции.【F:library/cli/commands/get_activity_data.py†L125-L155】
  - Почему важно: изменение глобального `Config` затрагивает последующие шаги пайплайна, нарушая принцип иммутабельности конфигураций.
  - Исправление: работать с глубокой копией `cfg.model_copy(deep=True)` и возвращать обновлённую структуру, а не менять оригинал.
- Проблема: `LoggerConfig` генерирует случайный `run_id`.
  - Пример: `return LoggerConfig(run_id=uuid.uuid4().hex, level=level)`.【F:library/cli/parser.py†L29-L44】
  - Почему важно: идентификатор запуска попадает в метаданные и логи, делая их недетерминированными и ломая сравнение хэшей.
  - Исправление: позволить переопределять `run_id` через ENV/CLI и использовать детерминированный дефолт (например, хэш входного файла).

Quality
- Проблема: Массовые lint-нарушения (`ruff` сообщает 890 ошибок).
  - Пример: несортированные импорты и неиспользуемые переменные в `library/cli/utils.py` и других файлах.【791a9b†L1-L118】
  - Почему важно: без автоформата и lint-проверок сложно поддерживать единый стиль, растёт риск скрытых багов.
  - Исправление: настроить `ruff --fix`, `ruff format` и включить их в CI/pre-commit.
- Проблема: `mypy --strict` выдаёт 403 ошибки, включая неверные `TypedDict` ключи и возвращаемые типы.
  - Пример: `scripts/get_target_data.py` передаёт несуществующие ключи в `TypedDict` и аргументы неправильных типов.【18648e†L1-L38】
  - Почему важно: отсутствие статической типизации скрывает ошибки при рефакторинге.
  - Исправление: пройтись по ошибкам, добавить аннотации, обновить `__all__`, скорректировать сигнатуры.

Errors
- Проблема: `pipeline_skip_existing` логируется через глобальный logger, не попадая в заглушку e2e теста.
  - Пример: `logger.info("pipeline_skip_existing", ...)` вызывается из модуля, который не переопределяется в тесте.【F:library/cli/commands/get_activity_data.py†L181-L194】
  - Почему важно: мониторинг и алерты не видят, что пайплайн пропущен по `--skip-existing`, тест падает.【F:tests/e2e/test_get_cli_pipelines.py†L1212-L1234】
  - Исправление: передавать экземпляр logger из скрипта в функцию или экспортировать structured logger через CLI-конфигурацию.
- Проблема: `_emit_completion_message` выводит только человекочитаемую строку.
  - Пример: `logger.info(message)` без имени события и структурированных полей.【F:scripts/get_activity_data.py†L445-L466】
  - Почему важно: системы наблюдаемости не могут построить метрики, а тексты сложнее парсить.
  - Исправление: логировать событие `pipeline_done` с полями и, при необходимости, дублировать текст в stdout.

Perf
- Проблема: `scripts/check_determinism.py` фиксирует несовпадающие SHA256.
  - Пример: при сравнении двух запусков выводится `Mismatch detected` с разными хэшами.【F:scripts/check_determinism.py†L114-L137】【99cce3†L1-L3】
  - Почему важно: отсутствие детерминизма делает невозможной верификацию выпусков и сравнение результатов между окружениями.
  - Исправление: стабилизировать run-id, порядок строк/столбцов и любые случайные элементы перед записью.
- Проблема: Нет доступного бинаря `/usr/bin/time`, поэтому perf-smoke выполняется без измерения времени.
  - Пример: попытка `/usr/bin/time` возвращает «No such file or directory».【a3b0b6†L1-L3】
  - Почему важно: без измерения трудно контролировать бюджеты исполнения.
  - Исправление: использовать встроенный `time` shell или Python-измерение (`perf_counter`) в собственном обёртке.

Testing
- Проблема: E2E тест на `--skip-existing` падает из-за отсутствующего события.
  - Пример: утверждение `pipeline_skip_existing in events` не выполняется.【F:tests/e2e/test_get_cli_pipelines.py†L1212-L1234】
  - Почему важно: критический сценарий не покрыт и не проходит.
  - Исправление: после фикса логгера обновить тесты и добавить проверки structured логов.
- Проблема: Нет property-based тестов для схем pandera и валидации DataFrame.
  - Пример: схемы в `library/postprocessing/pipeline/activities/schema.py` не тестируются на граничных значениях.【F:library/postprocessing/pipeline/activities/schema.py†L1-L120】
  - Почему важно: данные из внешних API часто содержат крайние случаи; без property-based тестов легко пропустить ошибки.
  - Исправление: добавить Hypothesis-тесты для парсеров/мапперов и валидаторов pandera.

Docs
- Проблема: `.meta.yaml` описан не полностью, секреты маскируются только по верхним ключам.
  - Пример: `_mask_secrets` проверяет только имена ключей на верхнем уровне, вложенные словари остаются с токенами.【F:library/common/metadata.py†L142-L200】
  - Почему важно: риск утечки ключей API при публикации артефактов.
  - Исправление: расширить маскирование на вложенные структуры и задокументировать формат `.meta.yaml`.

Actionable recommendations
Item | Effort (S/M/L) | Impact (Low/Med/High) | Owner | Proposed PR name | Acceptance criteria
--- | --- | --- | --- | --- | ---
Прокинуть logger из CLI в `run_activity_pipeline` | M | High | Data Platform | chore/fix-activity-logger-plumbing | `pipeline_skip_existing` фиксируется в тесте и логах; e2e проходит.
Детерминировать run_id и выводы | M | High | Data Platform | feat/deterministic-run-metadata | `check_determinism.py` возвращает одинаковый хэш дважды; meta содержит фиксированные таймстемпы.
Включить lint/typing в CI | S | Medium | Platform Eng | chore/enable-ruff-mypy-ci | `ruff check` и `mypy --strict` проходят без ошибок; обновлён `pyproject` и pre-commit.
Ограничить размер `chunk_fetch_failure_ids` | S | Medium | Data Platform | fix/chunk-failure-stats-cap | Metadata содержит не более N ID; sidecar сохраняет полный список отдельно.
Добавить pandera/property-based тесты | M | Medium | QA/Data | test/add-pandera-hypothesis-suites | ≥5 новых тестов по ключевым сценариям, проходят локально и в CI.
Документировать и маскировать `.meta.yaml` | S | Medium | Technical Writers | docs/meta-yaml-spec | README/meta guide обновлён; секреты маскируются на всех уровнях.

Code snippets / mini-diffs
1. Logger прокидываем через аргумент:
```diff
-    if args.skip_existing and output_path.exists() and not args.force:
-        logger.info("pipeline_skip_existing", output=str(output_path))
+    event_logger = logger if run_logger is None else run_logger
+    if args.skip_existing and output_path.exists() and not args.force:
+        event_logger.info("pipeline_skip_existing", output=str(output_path))
```
2. Фиксируем run_id от входного файла:
```diff
-def create_logger_config(level: str) -> LoggerConfig:
-    return LoggerConfig(run_id=uuid.uuid4().hex, level=level)
+def create_logger_config(level: str, *, seed: str | None = None) -> LoggerConfig:
+    run_id = seed or uuid.uuid5(uuid.NAMESPACE_URL, level).hex
+    return LoggerConfig(run_id=run_id, level=level)
```
3. Ограничение списка ID в `ChunkFailureTracker`:
```diff
-        unique_ids = sorted({identifier for failure in self._failures for identifier in failure.chunk_ids})
-        return {
-            "chunk_fetch_failure_chunks": len(self._failures),
-            "chunk_fetch_failure_ids": unique_ids,
-        }
+        unique_ids = sorted({identifier for failure in self._failures for identifier in failure.chunk_ids})
+        sample = unique_ids[:50]
+        return {
+            "chunk_fetch_failure_chunks": len(self._failures),
+            "chunk_fetch_failure_ids": sample,
+            "chunk_fetch_failure_overflow": max(0, len(unique_ids) - len(sample)),
+        }
```
4. Работать с копией конфигурации:
```diff
-    batch_size = getattr(cfg.activity, "batch_size", None)
+    cfg_mut = cfg.model_copy(deep=True)
+    batch_size = getattr(cfg_mut.activity, "batch_size", None)
```
5. Структурированный completion лог:
```diff
-    logger.info(message)
+    logger.info(
+        "pipeline_done",
+        rows=resolved_rows,
+        null_fraction=null_fraction_value,
+        duration=duration_s,
+        output=str(output_path) if output_path else "n/a",
+        mode=mode,
+    )
+    logger.info(message)
```

Risk & rollback
- Исправление логгера `run_activity_pipeline`: риск — пропустить существующие интеграции, ожидающие глобальный logger. Мониторим structured-логи и e2e. Rollback — вернуть старый вызов, если обнаружены отсутствующие события у сторонних потребителей.
- Детерминизация run_id и метаданных: риск — сторонние системы завязаны на случайный run_id. Наблюдаем по CI отчётам и внутренним дашбордам. Rollback — временно включить флаг «legacy_run_id», возвращающий прежнее поведение.
- Ограничение `chunk_fetch_failure_ids`: риск — аналитики потеряют полный список ID в JSON. Мониторинг — сравниваем sidecar размеры и отчёты QA. Rollback — поднять лимит или вернуть полный список, используя отдельный файл.

PR plan
1. **Fix activity logger plumbing** — Протянуть logger/LoggerConfig в `run_activity_pipeline`, обновить e2e тесты. Тест-план: `pytest -q --disable-warnings` + e2e `test_get_activity_run_skip_existing`. Успех: тесты проходят, structured лог содержит событие.
2. **Deterministic run metadata** — Ввести детерминированный run_id, стабилизировать `write_meta_yaml`, добавить регрессионный тест `scripts/check_determinism.py`. Тест-план: локально `python scripts/check_determinism.py`. Успех: одинаковый SHA, скрипт завершаетcя с 0.
3. **Static analysis enablement** — Исправить импорты/аннотации в трогаемых модулях, включить `ruff format` и `mypy` в CI. Тест-план: `ruff check .`, `mypy --strict --ignore-missing-imports .`. Успех: команды проходят, CI зелёный.
4. **Failure stats compaction** — Ограничить списки ID, добавить unit-тест на сериализацию. Тест-план: `pytest tests/unit/test_fetch_retry.py`. Успех: sidecar логика сохраняет не более 50 ID.
5. **Schema validation tests** — Добавить pandera/property-based тесты, покрывающие ключевые DataFrame. Тест-план: `pytest -q`. Успех: ≥5 новых тестов, покрытия ключевых сценариев достигнуто.

Acceptance criteria
- Все обязательные проверки запущены; где не получилось (отсутствует `/usr/bin/time`) дано объяснение.【a3b0b6†L1-L3】
- По каждому приоритетному issue предложен конкретный фикс или дифф (см. Code snippets).
- Политика детерминизма CSV описана: источник проблемы и ожидаемый фикс (см. Executive summary, Perf).
- Ошибки конфиг-валидации продемонстрированы и план фикса указан (mutation run_id).【F:library/cli/commands/get_activity_data.py†L125-L155】
- Запланировано ≥5 тестов для ключевых потоков (см. Actionable recommendations / PR plan).

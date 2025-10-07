Generated: 2025-01-01T00:00:00Z

Executive summary

1. [Errors] Модуль декларативного постпроцессинга требует функцию `resolve_dotted_path`, но она отсутствует в `library/postprocess/common/import_utils.py`, что ломает импорт `library.postprocess` во всех CLI и тестах → запуск скриптов и пайплайнов завершается `ImportError`, сборка падает на старте → добавить совместимую реализацию `resolve_dotted_path` (обёртку над `import_by_path`) и покрыть загрузку конфигурации регрессионным тестом.【F:library/postprocess/common/config.py†L14-L176】【F:library/postprocess/common/import_utils.py†L12-L84】
2. [Structure] Скрипты `scripts/get_target_data.py` и др. вынуждены динамически патчить legacy-модули (`setattr`, `_override_cli_meta_writer`) из-за несогласованности между `library/postprocess` и `library/postprocessing` → жёсткие зависимости от внутренних деталей, высокая хрупкость при обновлениях → перенести слой совместимости в пакет `library.postprocess` и отказаться от monkey-patching в CLI.【F:scripts/get_target_data.py†L95-L140】
3. [Quality] В `library/postprocessing/names.py` одновременно определены две разные версии `process_target_names`, возвращающие строку и словарь → противоречивый публичный API, нарушенные сигнатуры и предупреждения mypy → оставить одну функцию с чёткой сигнатурой, выделить legacy-обёртку с явным названием.【F:library/postprocessing/names.py†L390-L425】【F:library/postprocessing/names.py†L592-L619】
4. [Testing] `pytest` останавливается ещё на импорте `tests/e2e/test_activity_logging.py` по той же причине (`resolve_dotted_path`), поэтому критические сценарии не покрываются вообще → восстановить импорт и добавить e2e-smoke, проверяющий запуск CLI.【aae7a8†L1-L18】
5. [Testing] Скрипт `scripts/check_determinism.py` падает сразу после старта из-за импортной ошибки, и детерминизм CSV не проверяется → починить загрузку постпроцессинга и добавить фикстуру, подтверждающую равенство хэшей в CI.【8b9edb†L1-L15】
6. [Performance] Производительный smoke `PYTHONHASHSEED=0 python scripts/get_activity_data.py --limit 500 --dry-run` не выполняется по той же ошибке, а штатная команда `/usr/bin/time` отсутствует → после восстановления импорта задействовать встроенный `time.perf_counter()` или `python -m timeit` и зафиксировать бюджет выполнения в тесте-профиле.【26dcd5†L1-L3】【16d81a†L1-L13】
7. [Quality] `ChunkFailureTracker.stats()` возвращает ссылку на общий `_EMPTY_STATS` (`{}`), поэтому любое постобработка (например, `stats_extra.setdefault`) мутирует глобальное состояние и ломает последующие вызовы → возвращать новый словарь при отсутствии ошибок или `copy()` результата.【F:library/common/fetch_retry.py†L22-L70】
8. [Structure] `run_pipeline` вытягивает `output_path`/`failure_path` через `locals()` из-за отражательных вызовов и патчинга сигнатуры → код трудно читать, а оптимизаторы/стат-анализаторы видят «магические» побочные эффекты → сделать параметры позиционными и удалить обходные манёвры; для бэкендов оставить тонкую адаптацию в месте вызова.【F:library/cli_utils.py†L269-L288】
9. [Quality] `ruff check`, `ruff format` и `mypy --strict` завершаются сотнями ошибок, что подтверждает отсутствие линтинга и статпроверок в CI → зафиксировать набор правил в pre-commit, довести код до «зелёного» состояния и включить проверки в pipeline.【130b7c†L1-L120】【ad8e81†L1-L181】【aaca19†L1-L120】
10. [Config] Из-за отсутствия резолвера шагов конфигурации (`resolve_dotted_path`) YAML-файлы постпроцессинга невалидируемы и фактически неиспользуемы; опечатка в callable не будет поймана на этапе загрузки → после восстановления резолвера добавить тест, грузящий `config/pipeline/*.yaml`, и схему на основе pandera/pydantic.【F:library/postprocess/common/config.py†L60-L176】

Scores (0–5)

Structure: 2/5 — два параллельных стека (`library/postprocess` и `library/postprocessing`) плюс обилие monkey-patching в CLI создают неявные связи и ломают инкапсуляцию, что видно по костылям в `scripts/get_target_data.py` и `run_pipeline`.
Config: 2/5 — конфиги богаты комментариями, но ключевой механизм (декларативные шаги) не работает из-за отсутствующего резолвера; валидация YAML и env-override не покрыта тестами.【F:library/postprocess/common/config.py†L14-L176】
Quality: 1/5 — дубли в `library/postprocessing/names.py`, возврат глобальных мутабельных объектов и сотни lint-ошибок демонстрируют низкую поддерживаемость.【F:library/postprocessing/names.py†L390-L619】【F:library/common/fetch_retry.py†L22-L70】
Errors: 1/5 — базовые команды (`pytest`, `get_activity_data`, `check_determinism`) падают на импорте, что указывает на отсутствие fail-fast тестов и интеграций.【aae7a8†L1-L18】【16d81a†L1-L13】
Perf: 2/5 — есть rate-limiter и кеш в клиентах, но perf-smoke не запускается, а профилировщик `tools/profile_activity_pipeline.py` сам не проходит линтеры и mypy.【16d81a†L1-L13】【aaca19†L1-L120】
Testing: 1/5 — структура каталога тестов присутствует, но фактически тесты не выполняются; отчётность JSON/MD реализована, однако не используется из-за падения на импортных ошибках.【aae7a8†L1-L18】【b093dc†L1-L160】
Docs: 3/5 — README и комментарии в `config/config.yaml` подробны, но отсутствует актуальный гайд по переходу с legacy-пакетов `library.postprocessing` на новые модули.

Findings by category

### Structure
- Проблема: Два конкурирующих пакета постпроцессинга и прокси-CLI заставляют скрипты патчить модули во время исполнения.【F:scripts/get_target_data.py†L95-L140】  
  Почему важно: любое изменение внутреннего API ломает рабочие пайплайны, потому что нет единого места совместимости.  
  Исправление: перенести адаптеры в `library.postprocess` (alias-функции, re-export) и запретить monkey-patching в CLI; добавить интеграционный тест на вызов `library.postprocess.activities.run_activity_pipeline`.
- Проблема: `run_pipeline` зависит от `locals()` и неявных побочных эффектов для обязательных аргументов.【F:library/cli_utils.py†L269-L288】  
  Почему важно: сложно анализировать типами, трудно рефакторить, риск сбоев при оптимизациях.  
  Исправление: изменить сигнатуру `run_pipeline` на явные позиционные параметры (`output_path: Path`, `failure_path: Path`) и адаптировать вызовы; добавить mypy-тест.

### Config
- Проблема: `load_pipeline_config` вызывает несуществующий `resolve_dotted_path`, из-за чего YAML шага нельзя загрузить даже при корректной конфигурации.【F:library/postprocess/common/config.py†L14-L176】  
  Почему важно: конфигурации с ошибками не валидируются заранее, пайплайн падает только в рантайме.  
  Исправление: реализовать `resolve_dotted_path` и добавить тест, проходящий по всем файлам `config/pipeline/*.yaml`.
- Проблема: env-overrides не покрыты тестами — `_apply_env_overrides` парсит любую строку через YAML и тихо глотает ошибки, что может приводить к неожиданным типам (например, `"false"` → `False`).【F:library/config/env.py†L239-L258】  
  Почему важно: опечатка в переменной окружения может менять тип параметра без диагностики.  
  Исправление: добавить property-based тесты с Hypothesis и строгий контроль допустимых типов (строка/число/булево).

### Quality
- Проблема: Две разные реализации `process_target_names` в одном модуле с разными типами возвращаемых значений.【F:library/postprocessing/names.py†L390-L619】  
  Почему важно: статический анализ и IDE не понимают API, высокий риск вызвать «не ту» версию.  
  Исправление: оставить современную реализацию, legacy-вариант перенести в отдельную функцию `_legacy_process_target_names` с явным предупреждением.
- Проблема: `ChunkFailureTracker.stats` возвращает мутабельный singleton `_EMPTY_STATS`, что приводит к неожиданным побочным эффектам при модификации результата.【F:library/common/fetch_retry.py†L22-L70】  
  Почему важно: дополнительные метрики (`stats_extra`) могут навсегда изменить глобальное состояние, ломая последующие вызовы.  
  Исправление: возвращать новый словарь при каждом вызове (`return {}`) и добавить тест, проверяющий неизменяемость исходника.

### Errors
- Проблема: Импорт `library.postprocess` падает при загрузке, вслед за ним ломаются CLI и тесты.【F:library/postprocess/common/import_utils.py†L12-L84】【aae7a8†L1-L18】  
  Почему важно: ключевые команды недоступны, пайплайн не стартует.  
  Исправление: реализовать отсутствующий резолвер и покрыть smoke-тестом `python -m library.postprocess.activities`.
- Проблема: Скрипт `check_determinism.py` не выполняется, поэтому никто не замечает изменения в CSV-диалекте.【8b9edb†L1-L15】  
  Почему важно: без ежедневной проверки теряется гарантия детерминизма, что критично для downstream-систем.  
  Исправление: после починки импорта добавить сравнение хэшей в CI и опубликовать результат в `reports/test_summary.md`.

### Performance
- Проблема: perf-smoke нельзя запустить из-за отсутствия системной `time` и импортной ошибки.【26dcd5†L1-L3】【16d81a†L1-L13】  
  Почему важно: нет контроля SLA на сборку 500 записей, нельзя ловить регрессии.  
  Исправление: перейти на `time.perf_counter()` внутри Python-скрипта и добавить pytest-benchmark для мини-набора.

### Testing
- Проблема: pytest и линтеры падают до выполнения тестов, поэтому требования по 95% успеха и отчётности недостижимы.【aae7a8†L1-L18】【130b7c†L1-L120】  
  Почему важно: основная цель тестового контура (стабильная проверка пайплайна) не выполняется.  
  Исправление: починить импорты, включить линтеры в CI и добавить smoke для каждого CLI.

### Docs
- Проблема: Нет документации по миграции от `library.postprocessing.*` к `library.postprocess.*`, из-за чего разработчики продолжают использовать legacy-модули и добавлять костыли.【F:scripts/get_target_data.py†L95-L140】  
  Почему важно: копирование кода и патчинг вместо использования нового API.  
  Исправление: добавить раздел в README и docstring-и в новые пакеты, описывающие рекомендуемый путь.

Actionable recommendations

| Item | Effort | Impact | Owner | Proposed PR name | Acceptance criteria |
|---|---|---|---|---|---|
| Реализовать `resolve_dotted_path` и smoke-тесты загрузки pipeline YAML | S | High | Core backend | "postprocess-resolver-fix" | `library.postprocess` импортируется, pytest запускается до конца, determinism-check проходит |
| Удалить monkey-patching в CLI, вынести совместимость в `library.postprocess` | M | High | Core backend | "cli-compat-layer" | Все CLI вызывают новые API без `setattr`, smoke e2e тесты зелёные |
| Дедуплицировать `process_target_names` и стабилизировать API | M | Med | Data ops | "target-names-api-cleanup" | Есть одна функция с docstring, unit-тесты покрывают возвращаемый тип |
| Починить `ChunkFailureTracker.stats` и добавить тест | S | Med | Data ops | "chunk-failure-stats-fix" | Новые тесты подтверждают отсутствие мутаций, mypy/ruff зелёные |
| Включить ruff/black/mypy/pytest в pre-commit и CI | M | High | Platform | "ci-lint-and-type-gates" | Все команды (`ruff`, `mypy`, `pytest`) проходят локально и в CI, отчёты JSON/MD прикрепляются |

Code snippets / mini-diffs

1. Реализация отсутствующего резолвера:
```diff
--- a/library/postprocess/common/import_utils.py
+++ b/library/postprocess/common/import_utils.py
@@
-__all__ = ["ImportResolutionError", "import_by_path"]
+def resolve_dotted_path(path: str, expected_type: type | tuple[type, ...] | object = _DEFAULT_SENTINEL) -> Any:
+    """Backward-compatible alias used by legacy YAML configs."""
+
+    return import_by_path(path, expected_type=expected_type)
+
+
+__all__ = ["ImportResolutionError", "import_by_path", "resolve_dotted_path"]
```
2. Возврат нового словаря в `ChunkFailureTracker.stats`:
```diff
@@
-        if not self._failures:
-            return _EMPTY_STATS
+        if not self._failures:
+            return {}
```
3. Дедупликация API `process_target_names`:
```diff
--- a/library/postprocessing/names.py
+++ b/library/postprocessing/names.py
@@
-def process_target_names(...):
-    # legacy body
-    return str(output_path)
+def process_target_names(...):
+    # современная реализация, возвращающая словарь с путём и сводкой
+    ...
+
+def process_target_names_legacy(...):
+    warnings.warn("process_target_names_legacy is deprecated", DeprecationWarning)
+    return str(process_target_names(...)["path"])
```
4. Удаление monkey-patching в CLI:
```diff
-    with _override_cli_meta_writer():
-        return run_pipeline(...)
+    definition = dataclasses.replace(definition, stats_callback=_capture_stats)
+    return run_pipeline(definition=definition, ...)
```
5. Добавление smoke-теста конфигурации:
```diff
+def test_pipeline_configs_load() -> None:
+    for path in (CONFIG_PIPELINE_DIR).glob("*.yaml"):
+        cfg = load_pipeline_config(path.stem)
+        assert cfg.steps, path
```

Risk & rollback

- `postprocess-resolver-fix`: Риск — ошибочное разрешение путей вызовет runtime при загрузке YAML. Мониторинг — запуск pytest и smoke `check_determinism`. Откат — вернуть старый модуль и временно зафиксировать зависимость в requirements-lock.
- `cli-compat-layer`: Риск — новые импорты могут не покрыть все legacy-сценарии. Мониторинг — e2e тесты CLI, метрики логирования WARN. Откат — временно вернуть monkey-patching и открыть issue с перечнем несовместимых модулей.
- `target-names-api-cleanup`: Риск — внешние потребители могли зависеть от строкового возвращаемого значения. Мониторинг — обновление документации и поиск прямых импортов через `rg`. Откат — восстановить legacy-обёртку и добавить DeprecationWarning.

PR plan

1. `postprocess-resolver-fix` — область: модуль `library.postprocess.common`. Тест-план: unit-тест `test_pipeline_configs_load`, smoke `pytest -q`. Метрики успеха: `pytest`, `scripts/check_determinism.py` завершаются успешно.
2. `cli-compat-layer` — область: CLI-скрипты и `library.postprocess` алиасы. Тест-план: e2e тесты `tests/e2e/*`. Метрики: успешный запуск `get_activity_data --dry-run` и отсутствие monkey-patching.
3. `target-names-api-cleanup` — область: `library/postprocessing/names.py`. Тест-план: unit/integration тесты на генерацию имен. Метрики: mypy без ошибок, вызовы возвращают ожидаемый словарь.
4. `ci-lint-and-type-gates` — область: tooling (`pyproject.toml`, workflows). Тест-план: `ruff`, `mypy`, `pytest`, `scripts/check_determinism.py`. Метрики: все проверки зелёные, отчёты JSON/MD публикуются.

Acceptance criteria

- Все обязательные проверки из задания выполнены либо обоснованно отмечены (см. вывод команд ниже).
- Для каждого приоритетного issue предложен конкретный фикс (см. мини-диффы 1–5).
- Политика детерминизма: после фикса `check_determinism.py` сравнивает SHA256, скрипт выполняется дважды и подтверждает совпадение.
- Ошибки валидации конфигов: после появления `resolve_dotted_path` загрузка YAML покрыта тестом `test_pipeline_configs_load`.
- Планируется ≥5 тестов: загрузка конфигов, smoke CLI, determinism hash, unit для `ChunkFailureTracker`, unit для `process_target_names`.

Обязательные проверки (фактический вывод)
- `ruff check .` — ошибки (см. лог).【130b7c†L1-L120】
- `ruff format --check .` — ошибки форматирования.【ad8e81†L1-L181】
- `mypy --strict --ignore-missing-imports .` — 308 ошибок.【aaca19†L1-L120】
- `pytest -q --disable-warnings` — ImportError на `resolve_dotted_path`.【e9bf0a†L1-L18】
- `pytest --maxfail=1 --durations=10` — тот же ImportError.【aae7a8†L1-L18】
- `python -m pip list` — список окружения (numpy/pandas/pandera и пр.).【64112f†L1-L42】
- `python -c "import sys,platform; print(sys.version); print(platform.platform())"` — Python 3.11.12, Linux 6.12.13.【113ba4†L1-L4】
- `PYTHONHASHSEED=0 time python scripts/get_activities.py --limit 500 --dry-run` — `time` отсутствует; см. ниже замену.【26dcd5†L1-L3】
- `PYTHONHASHSEED=0 python scripts/get_activity_data.py --limit 500 --dry-run` — ImportError (`resolve_dotted_path`).【16d81a†L1-L13】
- `python scripts/check_determinism.py` — ImportError, determinism не подтверждён.【8b9edb†L1-L15】

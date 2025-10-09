# Executive summary

1. [Config] Манифест словарей `dictionary_root` не содержит актуальный SHA256 `3d2b7a7d…`, из-за чего любой запуск CLI/тестов немедленно падает при валидации ресурсов. → Полный стоп пайплайна. → Добавить новый хэш в `config/dictionary/manifest.yaml` и зеркально в `manifest.allowlist.yaml`, пересобрав артефакт и дополнив документацию по обновлению словарей.
2. [Testing] Проверка детерминизма вызывает `scripts/get_activity_data.py` напрямую, но скрипт лишён bootstrap-кода, поэтому подпроцесс не видит пакет `library` и падает. → Авто-тесты детерминизма не запускаются. → Вернуть bootstrap-блок в `scripts/get_activity_data.py` и запускать CLI через `python -m scripts.get_activity_data` либо прокладывать `PYTHONPATH` внутри `check_determinism.py`.
3. [Errors] Все e2e тесты и smoke-команды падают на той же проверке checksum, что зафиксировано в выводах pytest/CLI. → Невозможность проверить пайплайн перед релизом. → Исправить словарные checksum'ы и добавить health-check в `scripts/run_tests.py`, который валидирует манифест до запуска pytest.
4. [Structure] `scripts/get_activity_data.py` (~1600+ строк) смешивает CLI, сетевые вызовы, постобработку и генерацию отчётов в одном модуле без bootstrap; повторное использование почти невозможно. → Трудоёмкая поддержка и высокий риск регрессий. → Разбить модуль на подмодули (`cli`, `pipelines`, `postprocess`) и опубликовать явное API в `library/cli/commands`.
5. [Quality] Модуль `library/cli_utils.py` создаёт глобалы без аннотаций и публичного API (`required_cols`, `schema_columns_dict` и др.). → mypy/ruff фиксации падают, реальное API остаётся неочевидным. → Добавить `__all__`, типы и вынести подготовку схем в отдельные dataclass/TypedDict структуры.
6. [Config] Автоматическое расширение конфигурации не обрабатывает случай, когда `--config` отсутствует (Path(None) из argparse не ловится). → Падение в рантайме с `TypeError` вместо дружелюбного сообщения. → В `library/cli/parser.apply_config_overrides` перед вызовом `load_config` валидировать `config_path` и выводить actionable error.
7. [Errors] `ChemblClient` не логирует первопричину таймаутов/HTTP 429 и не делает структурированные события об истощении rate limiter. → Сложно отлавливать деградации API. → Добавить структурированные WARN/ERROR с полями `event=retry`, `backoff_s`, `status`, `retry_count`.
8. [Performance] `library/utils/cli_tools/get_activities._frame_from_records` материализует генератор целиком в память. → При `limit` > десятков тысяч CLI падает по памяти. → Писать чанками через `pd.DataFrame.from_records(..., columns=...)` или стримить в CSV без полного materialize.
9. [Testing] `scripts/run_tests.py` не проверяет success-rate ≥95% если pytest упал раньше (как сейчас), а JSON/MD отчёты не создаются. → CI не выдаёт консистентный отчёт. → Обернуть генерацию отчётов в `try/finally` и при ошибке создавать отчёты с признаком `error`.
10. [Docs] Нет описанной процедуры обновления словарей/allowlist и критериев детерминизма рядом с данными. → Сбои вроде текущего остаются незадокументированными. → Добавить `config/dictionary/README.md` раздел «Как обновить manifest/allowlist» и ссылку из основного README.

# Scores (0–5)

**Structure: 2/5.** Разделение между CLI, библиотекой и данными формальное: `scripts/get_activity_data.py` монолитен, клиенты/мапперы перемешаны, отсутствуют явные boundary-слои (`library/pipelines` напрямую импортируют CLI utils).

**Config: 2/5.** Есть `config.yaml` и алиасы ENV, но валидация ресурсов ломает пайплайн, нет graceful-degradation и bootstrap для CLI, allowlist не синхронизирован.

**Quality: 2/5.** Ruff/black/mypy не проходят (1000+ ошибок), в коде много неаннотированных глобалов и дублирования; публичное API модулей неочевидно.

**Errors: 1/5.** Таймауты/ретраи есть, но отсутствуют понятные лог-события, падения по checksum/Path(None) не перехватываются, CLI не fail-fast с подсказками.

**Perf: 2/5.** Внутренние утилиты не используют чанки/streaming, `get_activities` грузит всё в память, нет метрик по rate limiter.

**Testing: 1/5.** pytest падает до запуска, determinism check ломается, обязательные отчёты не генерируются при ошибке, часть тестов всё ещё бьётся о реальные ресурсы.

**Docs: 2/5.** README описывает запуск, но нет свежей инструкции по обновлению словарей, нет описания политик детерминизма/allowlist.

# Findings by category

## Structure

- **Манифест словарей и allowlist не синхронизированы.**
  - Примеры: `config/dictionary/manifest.yaml:1-107`, `config/dictionary/manifest.allowlist.yaml:1-65`.
  - Почему важно: любые CLI/тесты, импортирующие словари, падают на старте (см. pytest и CLI логи), нет обходного пути.
  - Исправление: пересобрать словари, добавить новый checksum в оба файла, задокументировать процедуру.
- **`scripts/get_activity_data.py` нарушает границы слоёв.**
  - Пример: `scripts/get_activity_data.py:10-120` импортирует почти всё приложение, отсутствует bootstrap.
  - Почему важно: модуль гигантский, сложно тестировать/переиспользовать; любые изменения требуют массовых правок.
  - Исправление: вынести CLI оболочку в `library/cli/commands/get_activity_data.py`, оставить в `scripts/` только thin wrapper.
- **`library/cli_utils.py` выступает «свалкой» общего кода.**
  - Примеры: `library/cli_utils.py:303-360`, `library/cli_utils.py:660-719`.
  - Почему важно: глобальные переменные без типов, функции не экспортируются через `__all__`, сложно понять контракт.
  - Исправление: выделить модули `library/cli/schema.py`, `library/cli/stats.py`, добавить типы и публичные экспортируемые объекты.

## Config

- **Checksum ресурса `dictionary_root` устарел.**
  - Примеры: `config/dictionary/manifest.yaml:6-107`, лог ошибки pytest `2b0f5f†L5-L27`.
  - Почему важно: ни один сценарий, требующий словари, не запускается.
  - Исправление: добавить хэш `3d2b7a7da5380896972b4ccac5ceaad1ccdaf19e2e2f7da995e70770ab75579a`, актуализировать allowlist и снапшоты.
- **`apply_config_overrides` допускает `config_path=None`.**
  - Пример: `library/cli/parser.py:486-571`.
  - Почему важно: при вызове CLI без `--config` (например, из API) получаем `TypeError: expected str, bytes or os.PathLike`.
  - Исправление: перед вызовом `load_config` проверить `config_path`, при отсутствии — показать сообщение с путём до `DEFAULT_CONFIG_PATH`.
- **ENV-алиасы есть, но нет валидации на неожиданные ключи.**
  - Пример: `library/config/env.py:86-137` не проверяет, что путь существует в модели.
  - Почему важно: опечатка в переменной окружения silently игнорируется.
  - Исправление: добавить сбор предупреждений для неизвестных путей и включить их в отчёт загрузки конфигурации.

## Quality

- **Ruff/black не проходят.**
  - Примеры: `ruff check` вывод `63e634†L1-L112`, `ruff format --check` вывод `c4f73c†L1-L111`.
  - Почему важно: кодстайл разъезжается, сложно ревьюить.
  - Исправление: завести pre-commit, отформатировать код и включить проверку в CI.
- **Mypy выдаёт 374 ошибки, включая отсутствующие stubs и неверные типы.**
  - Пример: `mypy --strict` вывод `44ff25†L1-L115`.
  - Почему важно: невозможно полагаться на статическую проверку.
  - Исправление: установить `types-requests`, `types-PyYAML`, добавить аннотации в CLI/utils, включить скрипты в область mypy.
- **`library/cli_utils.py` отсутствует `__all__`, что ломает `from library.cli_utils import ...`.**
  - Пример: `scripts/get_activity_data.py:45-56` (mypy жалуется на attr-defined).
  - Почему важно: статика и IDE не видят публичные функции.
  - Исправление: определить `__all__` и/или разнести функциональность по специализированным модулям.

## Errors

- **Determinism check падает из-за отсутствия bootstrap.**
  - Примеры: `scripts/check_determinism.py:95-112`, лог `25ac04†L1-L8`.
  - Почему важно: ключевая гарантия детерминизма не проверяется.
  - Исправление: запускать CLI через `-m` и добавить `bootstrap_cli` в `scripts/get_activity_data.py`.
- **Таймауты/ретраи ChemblClient не логируются.**
  - Пример: `library/clients/chembl.py:175-317`.
  - Почему важно: при деградациях (429/504) нет структурированных событий.
  - Исправление: в `_request_with_retry` писать WARN/ERROR с полями `retry`, `backoff`, `status`, `elapsed`.
- **Нет fail-fast для отсутствующего словаря.**
  - Пример: `library/resources/dictionaries.py:500-505` выбрасывает исключение без подсказки.
  - Почему важно: пользователи не понимают, как исправить несоответствие.
  - Исправление: расширить исключение советом «запустите tools/build_dictionary_resources.py».

## Performance

- **Materialize DataFrame в памяти.**
  - Пример: `library/utils/cli_tools/get_activities.py:79-100`.
  - Почему важно: рост лимита приводит к всплеску памяти.
  - Исправление: стримить записи через writer, не создавая полный список.
- **Отсутствует троттлинг логов и метрик по rate limiter.**
  - Пример: `library/common/rate_limiter.py` (нет счётчиков/логов).
  - Почему важно: сложно понять, когда лимит исчерпан.
  - Исправление: добавить счётчики, интегрировать с `logging`.

## Testing

- **pytest падает на импорте словарей.**
  - Пример: `2b0f5f†L5-L27`.
  - Почему важно: тестовый контур не запускается.
  - Исправление: обновить словари/manifest и добавить smoke-тест manifest'а.
- **`scripts/run_tests.py` не выпускает отчёт при падении.**
  - Пример: `scripts/run_tests.py:124-200` — отчёты формируются только после успешного pytest.
  - Почему важно: CI остаётся без JSON/MD отчётов.
  - Исправление: окружить генерацию `try/finally`, выдавать отчёт с `success_rate=0` при ошибке.

## Docs

- **Нет инструкции по обновлению словарей и allowlist.**
  - Пример: `config/dictionary/README.md` отсутствует.
  - Почему важно: разработчики не знают, как синхронизировать checksum.
  - Исправление: добавить раздел в README.
- **README не упоминает детерминизм/политику отчётности.**
  - Пример: `README.md` — нет раздела про `reports/test_report.json` и determinism check.
  - Почему важно: новые разработчики не запускают обязательные проверки.
  - Исправление: дописать раздел «Quality Gates».

# Actionable recommendations

| Item | Effort | Impact | Owner | Proposed PR name | Acceptance criteria |
| --- | --- | --- | --- | --- | --- |
| Обновить manifest/allowlist словарей и добавить новую контрольную сумму | S | High | Data Eng | `fix/dictionary-manifest-2025-10` | `dictionary_root` принимает новый SHA, pytest/CLI проходят манифест-чек |
| Вернуть bootstrap в `scripts/get_activity_data.py` и поправить `check_determinism.py` | S | High | Platform | `fix/cli-bootstrap-determinism` | `python scripts/check_determinism.py` завершаетcя 0, smoke-CLI работает |
| Добавить pre-commit с ruff/black/mypy и привести код | M | Med | Platform | `chore/code-quality-gates` | `ruff check`, `ruff format --check`, `mypy --strict` успешны |
| Расщепить `scripts/get_activity_data.py` на слой CLI и библиотеку | M | High | ETL | `refactor/activity-cli-modularisation` | Новый модуль `library.cli.commands.get_activity_data` покрыт unit/integration тестами |
| Добавить логирование retry/timeout в `ChemblClient` и метрики rate limiter | M | Med | Platform | `feat/chembl-client-observability` | WARN/ERROR события при повторных попытках, интеграционный тест фиксирует лог |
| Улучшить отчётность тестов при падении pytest | S | Med | QA | `fix/test-report-failures` | Даже при провале pytest создаётся JSON/MD отчёт с success_rate<95% |

# Code snippets / mini-diffs

1. **Актуализация SHA словарей**
```diff
--- a/config/dictionary/manifest.yaml
+++ b/config/dictionary/manifest.yaml
@@
-      - "ac5176986b0fd769a190182d91c69a2ab5e62606608ccf7d9704413fb39ef55b"
+      - "ac5176986b0fd769a190182d91c69a2ab5e62606608ccf7d9704413fb39ef55b"
+      - "3d2b7a7da5380896972b4ccac5ceaad1ccdaf19e2e2f7da995e70770ab75579a"
```

2. **Зеркальное обновление allowlist**
```diff
--- a/config/dictionary/manifest.allowlist.yaml
+++ b/config/dictionary/manifest.allowlist.yaml
@@
-  dictionary_root:
+  dictionary_root:
     - "9f0497f849122a4e625722b23b02b9aadc422ddbfc7cabe17ee252951e1e4a15"
     - "bb98601cdc63ee4aeab49dac849f545e516b2a0a9b720174444af8975115a0b2"
     - "bccf4cfc745addb3966efe9db8c3cd0f537ef3f5025d059d9cdaa412b2867092"
     - "db25392613353b15acb21c88c057f6422d8cd32aea1a3fc710e5a0c4d060b91b"
     - "564f3b40ddde94f6ec9c5b8124e494c2116cdb686be130eb0c1a151e7ddd246f"
     - "387d8a4b45d8960e5f899b85199a1013d3029258b8b75f42c6a0365f402023db"
     - "ac5176986b0fd769a190182d91c69a2ab5e62606608ccf7d9704413fb39ef55b"
+    - "3d2b7a7da5380896972b4ccac5ceaad1ccdaf19e2e2f7da995e70770ab75579a"
```

3. **Возврат bootstrap в CLI**
```diff
--- a/scripts/get_activity_data.py
+++ b/scripts/get_activity_data.py
@@
-from __future__ import annotations
-
-# Bootstrap code removed - not needed
+from __future__ import annotations
+
+if __package__ in {None, ""}:
+    from _bootstrap import bootstrap_cli
+else:  # pragma: no cover
+    from ._bootstrap import bootstrap_cli
+
+bootstrap_cli(__package__, __file__)
```

4. **Исправление determinism check**
```diff
--- a/scripts/check_determinism.py
+++ b/scripts/check_determinism.py
@@
-    cmd = [
-        sys.executable,
-        str(Path(__file__).resolve().parents[0] / "get_activity_data.py"),
+    cmd = [
+        sys.executable,
+        "-m",
+        "scripts.get_activity_data",
         "--limit",
         str(limit),
@@
-    return subprocess.run(cmd, text=True, capture_output=True, env=env)
+    env.setdefault("PYTHONPATH", str(Path(__file__).resolve().parents[1]))
+    return subprocess.run(cmd, text=True, capture_output=True, env=env)
```

5. **Fail-fast для пустого `config_path`**
```diff
--- a/library/cli/parser.py
+++ b/library/cli/parser.py
@@
-    try:
-        base_path_arg = getattr(args, "base_path", None)
+    if config_path in (None, argparse.SUPPRESS):
+        raise ConfigError("Configuration path is missing; pass --config or install defaults")
+
+    try:
+        base_path_arg = getattr(args, "base_path", None)
```

# Risk & rollback

- **Обновление manifest/allowlist.** Риск: разработчики с локальными модификациями словарей снова получат mismatch. Мониторинг: проверять `scripts/check_determinism.py` и `get_activity_data --dry-run`. Откат: вернуть предыдущие хэши и заблокировать обновление словарей.
- **Bootstrap в CLI.** Риск: если пакет устанавливается из PyPI, двойной bootstrap может менять `sys.path`. Мониторинг: smoke `python -m scripts.get_activity_data --dry-run`. Откат: условно активировать bootstrap только в режиме `__package__ in {None, ""}`.
- **Fail-fast config.** Риск: сторонние интеграции, передававшие `None`, начнут падать раньше. Мониторинг: отслеживать новые ошибки в логах `config_error`. Откат: ослабить проверку и вернуть прежнее поведение с предупреждением.

# PR plan

1. **`fix/dictionary-manifest-2025-10`** — синхронизация manifest/allowlist, smoke-тест dictionary hash. Метрика: pytest import проходит.
2. **`fix/cli-bootstrap-determinism`** — bootstrap и determinism check; тест-план: `python scripts/check_determinism.py`, `PYTHONHASHSEED=0 python scripts/get_activities.py --dry-run`.
3. **`chore/code-quality-gates`** — включить pre-commit (ruff/black/mypy), привести код. Метрика: все статические проверки зелёные.
4. **`refactor/activity-cli-modularisation`** — вынести CLI слой. Тест-план: unit для `library.cli.commands.get_activity_data`, e2e `tests/e2e/test_get_activity_data_cli.py`.
5. **`feat/chembl-client-observability`** — добавить логирование retry/лимитеров. Тест-план: мок `requests` с 429 и проверка логов.

# Acceptance criteria

- Все обязательные проверки (`ruff check`, `ruff format --check`, `mypy --strict`, `pytest`, runtime sanity, determinism) выполняются или зафиксирован обоснованный ским.
- По каждому приоритетному issue указан конкретный фикс (см. mini-diffs выше).
- Политика детерминизма описана: determinism check после фикса запускается успешно.
- Пример валидации конфигов (`apply_config_overrides`) демонстрирует ошибку и её исправление.
- Добавить ≥5 новых тестов (unit для bootstrap/config, integration для determinism, e2e для словарей) в рамках предложенных PR.

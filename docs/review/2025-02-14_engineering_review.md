# Инженерное ревью ChEMBL data acquisition (2025-02-14)

## Итоговое состояние

| Категория | Балл | Ключевые наблюдения |
|-----------|------|---------------------|
| Structure | 3/5 | Чёткая конфигурация и слои клиентов/пост-обработки, но ключевой модуль `get_data.py` разросся до 2200+ строк и смешивает CLI, оркестрацию и QA, что усложняет сопровождение. 【F:library/cli/commands/get_data.py†L1-L2223】 |
| Config    | 3/5 | Pydantic-модели и `config.yaml` богато документированы, однако вспомогательные конвертеры, такие как `_as_path`, некорректно работают с `bytes`, что ломает конфиг-алиасы и манифесты. 【F:library/resources/dictionaries.py†L671-L699】 |
| Quality   | 2/5 | Линтер фиксирует десятки ошибок (пересечения импортов, дублирование `document_steps`), а mypy выдаёт критичные несоответствия типов в HTTP-клиентах и CLI-пайплайнах, сигнализируя о деградации статики. 【9eb0cd†L1-L178】【b86cc7†L1-L8】 |
| Errors    | 2/5 | Логгирование шагов построено на строковых шаблонах; парсер логов падает на апострофах и валит E2E тест `get_data`, что говорит об отсутствии fail-fast на уровне QA. 【F:tests/helpers/logs.py†L28-L67】【f09520†L1-L83】 |
| Perf      | 3/5 | Есть батчинг и лимитеры в HTTP-клиентах, но отсутствие pyarrow делает векторизацию хрупкой (пакет обязателен), а fallback'ы в pandas_utils скрывают дорогостоящие даунгрейды. 【F:pyproject.toml†L19-L55】【5f7b24†L12-L74】 |
| Testing   | 2/5 | E2E сценарий `test_get_data_end_to_end` падает, pytest суммарно не проходит, а pre-commit/mypy останавливаются на первых ошибках — нет зелёного baseline. 【f09520†L1-L83】【ca04ed†L1-L27】 |
| Docs      | 4/5 | Руководства и README детализированы (в т.ч. русские версии) и соответствуют структуре тестов, но нет актуальной сводки по техническому долгу/CI. |

## Категорийный разбор

### Structure (3/5)
*Плюсы*: модульная конфигурация, вынесенные HTTP-клиенты, post-processing пайплайны реализованы отдельно от CLI.
*Минусы*: `library/cli/commands/get_data.py` содержит и CLI, и orchestration, и QA-хуки, что затрудняет локализацию дефектов и рефакторинг. Дубли экспорта модулей (например, `document_steps`) усложняют импортный граф. 【F:library/cli/commands/get_data.py†L1-L2223】【9eb0cd†L13-L104】

### Config (3/5)
*Плюсы*: `config/config.yaml` документирован, ENV-алиасы описаны в `models.py`.
*Минусы*: `_as_path` пытается построить `Path` из `bytes`, что в рантайме даёт `TypeError` и ломает разрешение путей в manifest-провайдерах. Нужно явное декодирование. 【F:library/resources/dictionaries.py†L671-L699】

### Quality (2/5)
*Минусы*: Ruff показывает 73 ошибки, включая переопределение импортов и style-задачи; mypy фиксирует отсутствие модулей `_bootstrap`, ошибки в HTTP-клиентах и CLI, что указывает на отсутствие контроля типобезопасности. 【9eb0cd†L1-L178】【b86cc7†L1-L8】

### Errors & Reliability (2/5)
*Плюсы*: Есть `ChunkFailureTracker`, `compute_backoff_delay`, retry-конфигурации.
*Минусы*: Логи шагов сериализуются строками с `%s`, парсер логов (`shlex.split`) не переносит экранированные кавычки и рушит end-to-end тест, таким образом QA-пайплайн не выдерживает собственные логи. 【F:tests/helpers/logs.py†L28-L67】【b25592†L1-L84】

### Performance (3/5)
*Плюсы*: Клиенты ограничены rate-limiter'ами и кэшами (`TTLCache`).
*Минусы*: `pyarrow` — обязательная зависимость, и при установке на Python 3.13 идёт сборка из исходников и падение (нет бинарей), что делает деплой нестабильным и лишает pandas `dtype_backend`. 【F:pyproject.toml†L19-L55】【8ba756†L1-L122】

### Testing (2/5)
*Минусы*: Pytest падает на `test_get_data_end_to_end`, pre-commit не проходит из-за mypy, покрытие ключевых сценариев не подтверждено. Нет зелёного baseline для CI. 【f09520†L1-L83】【ca04ed†L1-L27】

### Documentation (4/5)
*Плюсы*: README/DEVELOPER/QA-руководства подробны.
*Минусы*: нет отдельной страницы о статусе CI и открытых долговых задачах.

## Ключевые проблемы и рекомендации

| № | Файл / строки | Проблема | Почему важно | Рекомендация |
|---|----------------|----------|--------------|--------------|
| 1 | `tests/helpers/logs.py` 28-67 | `shlex.split` роняет парсер на кавычках в сообщениях (`log.exception`), из-за чего e2e-тест падает. | QA теряет fail-fast, CI не зелёный. 【F:tests/helpers/logs.py†L28-L67】【f09520†L1-L83】 | Либо переключить `shlex.split` в `posix=False` и ловить `ValueError`, либо заменять на JSON-структурирование логов и обновить парсер. Добавить негативный тест на строку с вложенными апострофами. |
| 2 | `pyproject.toml` 19-29 | `pyarrow` объявлен обязательной зависимостью; на Python 3.13 нет wheel — сборка падает. | Деплой/CI не устанавливаются с новыми Python, ломая reproducible builds. 【F:pyproject.toml†L19-L29】【8ba756†L1-L122】 | Вынести `pyarrow` в optional extra (`pip install chembl-data-acquisition[pyarrow]`) либо добавить маркеры `python_version < "3.13"`. Зафиксировать бинарную версию в lock-файле и документировать fallback на pandas-native types. |
| 3 | `library/resources/dictionaries.py` 671-699 | `_as_path` создаёт `Path` из `bytes` → `TypeError`. | ENV-алиасы и тестовые ресурсы могут быть байтовыми (например, через `os.environb`), что ломает конфиг. 【F:library/resources/dictionaries.py†L671-L699】【efbb88†L1-L10】 | Явно декодировать байтовые пути (`os.fsdecode`) перед `Path`, добавить тест с байтовым путём и линтерный guard. |
| 4 | `library/cli/commands/get_document_data.py` 89-95 | Двойной импорт `document_steps`, нарушает ruff и маскирует реальное определение. | Сигнализирует об отсутствующем контроле качества и усложняет ревью. 【9eb0cd†L13-L104】 | Удалить дублирующий импорт, включить ruff в CI с `--fix` и завести задачу на разбор backlog-а (73 ошибки). |
| 5 | `library/clients/chembl.py` 22-30 | Mypy фиксирует обращение к `urllib3.exceptions.NameResolutionError`, которого нет в свежих версиях. | При апдейте urllib3 логика обработки DNS-ошибок перестанет работать. 【b86cc7†L1-L8】 | Проверить текущие классы исключений (`urllib3.exceptions.MaxRetryError`) и переписать маппинг с явной проверкой наличия атрибута. |
| 6 | `library/common/pandas_utils.py` 45-68 | Fallback при `ImportError` silently убирает `dtype_backend`, что ломает векторизацию и маскирует отсутствие pyarrow. | Без pyarrow pipeline работает медленнее и неочевидно деградирует. 【5f7b24†L12-L74】 | При невозможности импортировать pyarrow логировать WARN и требовать явного флага `--disable-pyarrow`. |
| 7 | `scripts/get_data.py` 1-48 | Обёртки через `_bootstrap` ломают mypy/ruff (нет stubs) и усложняют CLI entrypoints. | Нарушена проверяемость: static analysis падает, CLI сложно тестировать. 【0be492†L1-L46】【b86cc7†L1-L8】 | Вынести bootstrap-логику в единый helper и добавить `py.typed`/stubs для `_bootstrap`. Настроить mypy на исключение этих файлов или переписать на `importlib.resources`. |
| 8 | `tests/e2e/test_get_data_end_to_end.py` | Тест падает из-за парсера логов, а fallback не обрабатывает `SchemaValidationError`. | Потеря проверки ключевых сценариев. 【f09520†L1-L83】【b25592†L1-L84】 | После фикса парсера добавить ассерт на WARN/ERROR c кавычками и убедиться, что повторный запуск идемпотентен. |
| 9 | `pyproject.toml` dev-зависимости | Нет фиксированной версии `pytest-json-report` и генерация отчётов не интегрирована в CI. | Тестовые отчёты могут ломаться при апдейте плагина. 【F:pyproject.toml†L31-L55】 | Пиновать плагины, обновить `scripts/run_tests.py` для graceful fallback, добавить smoke на генерацию `test_report.json`. |
|10 | `library/qa/reporting.py` / `library/cli/pipeline_definition.py` | Mypy ошибки (`partial[...]` возвращает DataFrame вместо Callable). | QA-хуки могут работать не так, как задумывалось, типы не контролируются. 【22618d†L1-L34】 | Переписать API `build_table_quality_hook`, явные типы для коллбеков, покрыть unit-тестом. |

## План действий и PR-пакеты

1. **PR A – Логи и тесты**
   - Исправить `parse_log_lines` (posix=False + fallback) и добавить e2e-тест с кавычками.
   - Подкорректировать логгер (`structlog`/JSON) для `log.exception`, чтобы сообщения не дублировались. 【F:tests/helpers/logs.py†L28-L67】【b25592†L1-L84】

2. **PR B – Зависимости и установка**
   - Вынести `pyarrow` в optional extra или добавить маркеры совместимости; обновить документацию по установке.
   - Добавить явный WARN в `pandas_utils` при откате к pandas-native dtypes. 【F:pyproject.toml†L19-L55】【5f7b24†L12-L74】

3. **PR C – Config / ресурсы**
   - Исправить `_as_path`, покрыть тестами байтовые пути и ENV-алиасы.
   - Обновить `config.models` для корректной типизации optional полей (устранить mypy warnings). 【F:library/resources/dictionaries.py†L671-L699】

4. **PR D – Статика и bootstrap**
   - Удалить дублирующие импорты, включить `ruff --fix` в CI.
   - Обновить `_bootstrap` и CLI-скрипты так, чтобы mypy не падал, добавить stubs. 【9eb0cd†L1-L178】【b86cc7†L1-L8】

5. **PR E – QA- хуки**
   - Починить `build_table_quality_hook` и связанные partial'ы; добавить unit-тест.
   - Обновить документацию QA по новым API. 【22618d†L1-L34】

6. **PR F – HTTP клиенты**
   - Перепроверить обработку DNS/timeout-ошибок (`NameResolutionError`) и добавить интеграционный тест с responses/pytest-httpserver.
   - Задокументировать политику таймаутов в README. 【b86cc7†L1-L8】

7. **PR G – get_data модуль**
   - Разбить `get_data.py` на слой orchestration и слой CLI, выделить pipeline builder. Добавить smoke-тест на новый API.

8. **PR H – CI отчётность**
   - Зафиксировать версии pytest-плагинов, автоматизировать генерацию `reports/test_report.json` и проверку success rate ≥95%.

9. **PR I – Performance guardrails**
   - Добавить health-check на наличие pyarrow в runtime и явное предупреждение при деградации производительности.

10. **PR J – Документация по техническому долгу**
   - Создать страницу `docs/QUALITY_STATUS.md` с текущими метриками (pytest/ruff/mypy) и процессом устранения долга.

## Статус обязательных проверок

- `pytest -q --disable-warnings` — **FAIL** (падает на `test_get_data_end_to_end__miniature_pipeline`). 【f09520†L1-L83】
- `ruff check .` — **FAIL** (73 ошибки). 【9eb0cd†L1-L178】
- `mypy library/ scripts/` — **FAIL** (`_bootstrap`, HTTP-клиенты, CLI). 【b86cc7†L1-L8】
- `pre-commit run --all-files` — **FAIL** (остановлено из-за mypy и других хуков, прервано после многократных ошибок). 【ca04ed†L1-L27】

## Дополнительно

- Установка `pip install -e .[dev]` прерывается сборкой `pyarrow` на Python 3.13. Нужно предусмотреть бинарные зависимости или минимальную версию Python. 【8ba756†L1-L122】


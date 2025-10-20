# Отчёт по ревью ChEMBL_data_acquisition

## Резюме
- Блокирующие дефекты CLI устранены: `scripts/check_determinism.py` собирается без
  синтаксических ошибок и корректно обрабатывает dry-run/real-run режимы, что
  подтверждено текущей реализацией скрипта.【F:scripts/check_determinism.py†L133-L215】
- Smoke-команды требуют полного dev-набора: без `pytest-cov` `python scripts/run_tests.py`
  завершается с кодом 4, а `pytest -q --disable-warnings` падает на импортах
  `freezegun`/`hypothesis`/`responses`, поэтому перед запуском нужен `pip install -e .[dev]`.【2abddb†L1-L13】【0c4395†L1-L31】
- Линтер `ruff` фиксирует 35 нарушений (I001/F401/UP017) — в основном несортированные
  импорты и неиспользуемые хелперы; SyntaxError/RuntimeError не наблюдаются.【82bf11†L1-L9】

## Проверки
- `python scripts/run_tests.py` — **ошибка конфигурации**: `pytest` не распознаёт
  флаги покрытия без плагина `pytest-cov` и завершает работу с кодом 4 до старта тестов.【2abddb†L1-L13】
- `pytest -q --disable-warnings` — **ошибка окружения**: отсутствуют `freezegun`,
  `hypothesis`, `responses`, коллекция тестов обрывается на импорте фикстур.【0c4395†L1-L31】
- `mypy --strict library` — **522 ошибки**: нет stubs для `pandas`/`yaml`/`requests`,
  множество CLI-функций возвращают `Any` или неправильные типы.【2d2413†L1-L4】
- `ruff check .` — **35 ошибок**: несортированные импорты и неиспользуемые
  хелперы в CLI и тестах (I001/F401/UP017).【82bf11†L1-L9】
- `pre-commit run --all-files` — **цепочка падает**: `ruff`/`black` переформатируют
  файлы, `mypy` и `pytest` повторяют ошибки отсутствующих зависимостей, запуск
  завершается статусом 1.【a28a08†L1-L38】【6dce54†L1-L35】【9a9fe8†L1-L120】【77ee15†L1-L33】

## Оценка категорий
| Категория | Балл | Комментарий |
|-----------|------|-------------|
| Structure | 4/5 | Основные CLI и отчётные утилиты работают без патчей `sys.path`; типизированные модели `Config` покрывают alias-ENV схему. Единичный профилировщик всё ещё мутирует `sys.path`, что усложняет запуск вне репозитория.【F:library/config/models.py†L1-L140】【F:scripts/check_determinism.py†L133-L215】【F:tools/make_md_summary.py†L1-L73】【F:tools/profile_activity_pipeline.py†L1-L40】 |
| Config    | 4/5 | `config/config.yaml` остаётся единой точкой настройки с подробными комментариями и управлением чанкингом/лимитами; отсутствует артефакт версии `.meta`, что осложняет отслеживание релизов.【F:config/config.yaml†L1-L122】 |
| Quality   | 3/5 | `ruff` сигнализирует о 35 проблемах (несортированные импорты, неиспользуемые `scripts`-хелперы, устаревшие `pandas`-аннотации); критических синтаксических ошибок нет.【82bf11†L1-L9】 |
| Errors    | 3/5 | CLI по‑прежнему корректно логируют, однако `scripts/run_tests.py` без `pytest-cov` падает до старта тестов, оставляя отчётность пустой.【F:scripts/check_determinism.py†L133-L215】【2abddb†L1-L13】 |
| Perf      | 3/5 | Пайплайны используют чанкинг и rate limiter из конфигурации, однако профилировщик по-прежнему тянет полный `pandas` и правит `sys.path`, что мешает лёгким экспериментам.【F:config/config.yaml†L10-L104】【F:tools/profile_activity_pipeline.py†L1-L40】 |
| Testing   | 2/5 | Структура tests/unit|integration|e2e и фикстура `deterministic_env` сохранены, добавлен модуль `tests/unit/scripts/test_run_tests_script_renamed.py` для контроля артефактов, но без dev-зависимостей (`freezegun`, `hypothesis`, `responses`) smoke-прогон `pytest` не выполняется.【F:tests/conftest.py†L1-L86】【F:tests/unit/scripts/test_run_tests_script_renamed.py†L1-L120】【0c4395†L1-L31】 |
| Docs      | 4/5 | README и локализованные руководства содержат QA-регламент и инструкции запуска тестов, но единый статус CI ещё не встроен в документацию.【F:README_RU.md†L1-L166】 |

## Детали по категориям
### Structure
- Типизированные модели `Config` и обновлённые CLI (`scripts/check_determinism.py`,
  `tools/make_md_summary.py`) не требуют ручного управления путями и используют
  штатные импорты.【F:library/config/models.py†L1-L140】【F:scripts/check_determinism.py†L133-L215】【F:tools/make_md_summary.py†L1-L73】
- Профилировщик `tools/profile_activity_pipeline.py` до сих пор вставляет корень
  репозитория в `sys.path`, что мешает упаковке и изолированному профилированию.【F:tools/profile_activity_pipeline.py†L1-L40】

### Config
- Централизованный `config/config.yaml` документирует лимиты, чанкинг и алиасы
  окружения, обеспечивая воспроизводимость пайплайнов.【F:config/config.yaml†L1-L122】
- Отсутствует обязательный артефакт версии пайплайна (`.meta`), поэтому
  восстановление конкретного релиза требует сверки git-истории вручную.

### Quality
- `ruff` фиксирует 35 нарушений — в основном несортированные импорты (I001)
  и неиспользуемые/устаревшие символы (`F401`, `UP017`).【82bf11†L1-L9】
- Несогласованные импорты и строковые аннотации сохранились в CLI и тестовых
  утилитах (`library/cli/commands/get_document_data.py`, `library/pipelines/testitem/cli.py`).【F:library/cli/commands/get_document_data.py†L30-L70】【F:library/pipelines/testitem/cli.py†L640-L660】

### Errors
- `scripts/check_determinism.py` корректно логирует и сигнализирует ошибки сухого
  и полного прогонов, удаляя временные файлы по завершении.【F:scripts/check_determinism.py†L133-L215】
- `python scripts/run_tests.py` без установленного `pytest-cov` падает на этапе
  парсинга аргументов `--cov*`, поэтому smoke-инструкции должны явно
  предусматривать `pip install -e .[dev]`.【2abddb†L1-L13】

### Perf
- Конфигурация по-прежнему управляет чанкингом и backoff'ом, что делает пайплайны
  масштабируемыми.【F:config/config.yaml†L10-L104】
- Профилировщик использует полноценный `pandas` и ручные патчи путей, из-за чего
  быстрый анализ узких мест невозможен без тяжёлых зависимостей.【F:tools/profile_activity_pipeline.py†L1-L40】

### Testing
- Тестовые каталоги и фикстура `deterministic_env` поддерживают детерминированную
  структуру запусков; для контроля артефактов добавлен модуль
  `tests/unit/scripts/test_run_tests_script_renamed.py`.【F:tests/conftest.py†L1-L86】【F:tests/unit/scripts/test_run_tests_script_renamed.py†L1-L120】
- Без `freezegun`, `hypothesis`, `responses` smoke-запуск `pytest -q --disable-warnings`
  обрывается на импортах, поэтому dev-набор обязателен перед CI-прогоном.【0c4395†L1-L31】

### Docs
- README описывает контроль качества, запуск `scripts/run_tests.py` и требования
  к success rate, что помогает синхронизировать локальные и CI-проверки.【F:README_RU.md†L1-L166】
- Для прозрачности стоит добавить в документацию ссылку на актуальный статус CI.

## Findings by Category
| № | Файл / строки | Проблема | Причина | Исправление |
|---|---------------|----------|---------|-------------|
| 1 | `scripts/run_tests.py` L125-L140 | Запуск падает с кодом 4: `pytest` не понимает флаги `--cov*` без установленного `pytest-cov`, отчёты не генерируются.| Скрипт без проверки выставляет аргументы покрытия и не подтверждает наличие плагина.| Перед запуском проверять наличие `pytest-cov` (или устанавливать его через dev-экстры) и фиксировать инструкцию `pip install -e .[dev]` в CI.【F:scripts/run_tests.py†L125-L140】【2abddb†L1-L13】 |
| 2 | `tests/e2e/test_get_tissue_data_cli.py` L1-L40; `tests/unit/test_fetch_retry.py` L1-L20; `tests/unit/test_openalex_client.py` L1-L20 | Smoke-прогон `pytest -q --disable-warnings` падает на импортах `freezegun`/`hypothesis`/`responses`, ключевые сценарии не покрываются.| Dev-зависимости не устанавливаются перед запуском тестов.| Обновить README/CI-скрипты: обязательный шаг `pip install -e .[dev]`, либо маркировать тесты как `skip` без зависимостей.【F:tests/e2e/test_get_tissue_data_cli.py†L1-L40】【F:tests/unit/test_fetch_retry.py†L1-L20】【F:tests/unit/test_openalex_client.py†L1-L20】【0c4395†L1-L31】 |
| 3 | `library/cli/commands/get_document_data.py` L30-L70 | `ruff` выдаёт I001/UP017: несортированные импорты и строковые аннотации `pandas`.| Миграция CLI оставила смешанный блок импортов и устаревшие типы.| Отсортировать импорты (`ruff --fix`) и заменить строковые типы на `from __future__ import annotations` + явные импорты stubs.【F:library/cli/commands/get_document_data.py†L30-L70】【82bf11†L1-L5】 |
| 4 | `library/pipelines/testitem/cli.py` L640-L660 | `ruff` фиксирует F821: используется `chain`, но не импортирован, и реэкспортирует неиспользуемые хелперы.| Обновление CLI не синхронизировало `__all__` и набор импортов.| Добавить `from itertools import chain`, пересобрать `__all__` и удалить лишние импорты, затем покрыть тестами re-export.【F:library/pipelines/testitem/cli.py†L640-L660】【a28a08†L33-L48】 |
| 5 | `library/schemas/activities.py` L40-L70 | `mypy --strict` сообщает о неизвестных `pa.Check`/`pa.DataFrameSchema`.| Аннотации завязаны на `pandera`, но типы не импортированы и не защищены `TYPE_CHECKING`.| Добавить явные импорты под `TYPE_CHECKING` и/или подключить `pandera` stubs в dev-зависимости.【F:library/schemas/activities.py†L1-L80】【9a9fe8†L1-L40】 |

## План исправлений и PR-пакетов
### Минимально-инвазивные (PR1)
1. Зафиксировать установку dev-экстры (`pytest-cov`, `freezegun`, `hypothesis`, `responses`) в smoke-скриптах и README, добавить проверку наличия `pytest-cov` в `scripts/run_tests.py` перед запуском.
2. Привести в порядок импорты и реэкспорт `library/cli/commands/get_document_data.py` и `library/pipelines/testitem/cli.py`, устранив замечания `ruff` (I001/F821) и покрыв изменения unit-тестами.
3. Обновить `tests/unit/scripts/test_run_tests_script_renamed.py` и сопутствующие smoke-инструкции так, чтобы отчёты удалялись даже при раннем выходе pytest.

### Среднесрочные (PR2)
4. Добавить `TYPE_CHECKING`-импорты для `pandera` в схемах (`library/schemas/activities.py` и связанные модули), подключить stubs и сократить ключевые ошибки `mypy --strict`.
5. Настроить автоматический прогон `ruff --fix`/`mypy` в CI после подготовки dev-зависимостей и публиковать отчёты в `reports/`.

### Долгосрочные (PR3)
6. Пересмотреть зависимость профилировщика от ручных правок `sys.path`, вынести лёгкий модуль профилирования без тяжёлых пакетов.
7. Внедрить артефакт версионирования `.meta` и отобразить статус CI в документации, чтобы разработчики видели актуальные результаты smoke-проверок.

## Регламент обновления отчёта
- **Ответственный:** QA-инженер проекта.
- **Частота:** после каждого релиза пайплайна или не реже одного раза в месяц.
- **Процедура:** выполняются команды из раздела «Проверки» в свежем окружении с `pip install -e .[dev]`; результаты и новые риски фиксируются в документе. При отсутствии релизов в течение месяца QA-инженер подтверждает актуальность, оставляя запись в git-журнале.

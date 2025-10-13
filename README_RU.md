# Утилиты ChEMBL Data Acquisition

> **Языки:** [English](README_EN.md) · Русский

Текст на русском полностью синхронизирован с английской версией и поддерживается
в актуальном состоянии вместе с `README_EN.md`.

## Сквозной конвейер

```mermaid
flowchart LR
    A[Входные идентификаторы\nCSV-файлы] -->|resolve| B[Document pipeline]
    B -->|enrich| C[Target pipeline]
    C -->|link| D[Assay pipeline]
    D -->|hydrate| E[Test item pipeline]
    E -->|join| F[Activity pipeline]
    G[[Tissue pipeline\\n(ручной запуск)]]
    B -.->|citations| F
    C -.->|targets| F
    G -.->|справочники| F
    style F fill:#dfeaff,stroke:#1e3a8a,stroke-width:2px
```

Каждый пайплайн идемпотентен и может запускаться независимо. Консольная команда
`get-data` (entry point `library.cli.entrypoints:get_data_main`) использует
единую конфигурацию и настройки логирования, чтобы выполнить цепочку
«документы → таргеты → ассайи → тестовые объекты → активности» автоматически и
воспроизводимо. Совместимая обёртка `python scripts/get_data.py` сохранена для
автоматизации, которая полагается на исторический путь, но она прокидывает
только общие флаги этапов (например, `--limit`, `--force`, `--skip-existing`,
`--dry-run`) и не понимает оркестратор-специфичные переопределения. Когда нужны
связи по тканям, `get_tissue_data` запускают отдельно, чтобы обновить
справочные таблицы перед запуском пайплайна активностей.

## Структура репозитория

| Путь | Описание |
|------|----------|
| `scripts/` | Точки входа CLI для каждого пайплайна и вспомогательные оркестраторы. |
| `library/` | Переиспользуемые модули: API-клиенты, пайплайны, схемы валидации, пост-обработка и QA-утилиты. |
| `config/` | Базовый YAML-конфиг, схемы и словари для обогащения. |
| `data/` | Небольшие фикстуры и тестовые входные данные, повторяющие структуру CSV. |
| `docs/` | Полный комплект документации на английском и русском языках, включая [`docs/ГОСТ.md`](./docs/%D0%93%D0%9E%D0%A1%D0%A2.md) с заметками по соответствию нормативам. |
| `tests/` | Детерминированный набор тестов pytest (unit, integration, e2e). |
| `reports/` | Каталог для JSON/Markdown отчётов о запуске тестов. |
| `Makefile` | Удобные команды для форматирования, линтинга, тестирования и проверки документации. |

Детальная структура пакетов, глоссарий и расширенные руководства описаны в
[`docs/ru/SUMMARY.md`](./docs/ru/SUMMARY.md) и
[`docs/en/SUMMARY.md`](./docs/en/SUMMARY.md).

## Быстрый старт

```bash
make init
source .venv/bin/activate
pre-commit install --install-hooks
```

Цель `init` проверяет версию интерпретатора из `.python-version` (если файл
существует), создаёт виртуальное окружение `.venv`, обновляет `pip`,
устанавливает все зависимости из `requirements-lock.txt`, а затем ставит проект
в editable-режим без повторного пересчёта зависимостей.

Предпочтителен запуск цели `make init`, но при необходимости можно повторить
последовательность вручную (например, внутри контейнера):

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements-lock.txt
pip install --no-deps -e .
```

> 📌 Нужна конкретная версия интерпретатора? Создайте файл `.python-version` с
> требуемым `major.minor.patch`. `make init` считает его и завершится с ошибкой,
> если активный `python` не совпадает, что помогает выровнять окружения локально
> и в CI.

### Совместимость с Python 3.13

Сборки `pyarrow` для Python 3.13 ещё не опубликованы. Зависимость остаётся
доступной для Python 3.11–3.12 благодаря маркеру окружения, а установки на 3.13
пропускают её автоматически. Поэтому `make init` успешно выполняется на всех
поддерживаемых интерпретаторах без попыток собрать Arrow из исходников. Если
нужен паркетный движок, установите альтернативу (например,
`pip install fastparquet`) или временно используйте Python 3.12, пока не выйдут
официальные колёса. Все CSV-ридеры и нормализаторы JSON прозрачно откатываются
к NumPy-бэкенду при отсутствии `pyarrow`/`fastparquet`, так что базовые
пайплайны работают на Python 3.13.

> ℹ️ Dev-зависимости в `pyproject.toml`, `requirements-dev.txt` и
> `requirements-lock.txt` обновляются одновременно. Сначала скорректируйте
> диапазоны версий в метаданных проекта, затем пересоздайте lock-файл, чтобы все
> три источника оставались синхронизированными.

### Pre-commit хуки

Репозиторий поставляется с настроенным
[`pre-commit`](https://pre-commit.com/) конфигом, который форматирует код,
проверяет YAML/TOML и валидирует статические метаданные перед каждым коммитом.
Установите инструмент внутри виртуального окружения и зарегистрируйте хуки в
локальном клоне:

```bash
pip install pre-commit
pre-commit install --install-hooks
```

Если вы уже выполнили команды из раздела «Быстрый старт» (`make init`, затем
`pre-commit install --install-hooks`), пакет установлен и на последующих клонах
потребуется лишь повторный вызов `pre-commit install --install-hooks`. Для
ручной проверки без коммита выполните:

```bash
pre-commit run --all-files
```

Хуки кэшируют свои окружения, поэтому повторные запуски анализируют только
затронутые файлы. После обновления конфигурации (например, смены версий хуков)
очистите кэш командами:

```bash
pre-commit autoupdate
pre-commit run --all-files
```

Все проверки должны успешно проходить локально перед push — CI выполняет тот же
набор, чтобы гарантировать единообразное форматирование и линтинг.

Изучите флаги оркестратора и отдельных пайплайнов:

```bash
get-data --help
python scripts/get_data.py --help  # обёртка совместимости (только общие флаги)
python scripts/get_document_data.py --help
```

Полный прогон на образцах идентификаторов с выводом в `./output`:

```bash
get-data \
  --base-path . \
  --input-dir data/input \
  --output-dir output \
  --config config/config.yaml \
  --date $(date -u +%Y%m%d)
```

> 💡 Если `--date` опущен, оркестратор подставит детерминированное значение из
> `local.io.default_date_prefix` (в `config/config.yaml` по умолчанию указан
> `20240101`). Для временных прогонов можно задать переменную окружения
> `CHEMBL_DA_DEFAULT_DATE_PREFIX` (доступен также алиас `CHEMBL_DA_DEFAULT_DATE`).

## Контроль качества

Каждый MR/PR обязан проходить детерминированные контрольные точки, настроенные
в CI. Быстрее всего воспроизвести их локально через обёртку тестов:

```bash
python scripts/run_tests.py
```

Перед запуском убедитесь, что установлены dev-зависимости: выполните `make init`
или `pip install -r requirements-dev.txt`, чтобы нужные плагины pytest (например,
`pytest-json-report`) были доступны. Если обязательного плагина нет, обёртка
немедленно завершит работу и выведет дружелюбный совет с той же командой
установки вместо попытки автоматической загрузки, которая всё равно провалится
в средах без доступа в интернет.

Команда запускает полный набор pytest с покрытием, формирует структурированный
отчёт `reports/test_report.json` и Markdown-выжимку `reports/test_summary.md`, а
при падении показателя `summary.success_rate` ниже **95 %** завершает выполнение
с ненулевым кодом. После правок тестов или кода пайплайна запустите скрипт
повторно и убедитесь, что отчёты идентичны за исключением временной метки —
любое расхождение указывает на потерю детерминизма, требующую исправления до
отправки MR. Не коммитьте сформированные отчёты: они считаются временными
артефактами сборки, игнорируются Git и автоматически прикладываются к CI.
Если нужно поделиться результатами локального прогона, приложите файлы к PR
или issue как артефакты, но не добавляйте их в историю репозитория. Успешный
прогон завершится с кодом `0`; нарушения порога успешности или
покрытия возвращают `1`, а ошибки генерации/валидации отчётов — `11`, позволяя
CI различать регрессии качества и проблемы окружения.

Каждый запуск также сохраняет сырое тело отчёта pytest (`reports/pytest_raw_report.json`)
рядом с агрегированными файлами. Если вы перенаправляете выходы через флаги
`--json` или `--markdown`, указывайте отдельные каталоги — перед записью новых
артефактов обёртка очищает родительскую директорию, чтобы избежать «висящих»
файлов от предыдущих запусков.

### Контракт именования логов

Все CLI-утилиты инициализируют структурированное логирование через общие
bootstrap-хелперы. По умолчанию каждая команда пишет в
`logs/<program>_<YYYYMMDD>.log`, где `<program>` — имя скрипта (например,
`get_data` → `logs/get_data_20250228.log`). Флаг `--base-path` или переменная
среды `CHEMBL_DA_BASE_PATH` переносят каталог в `<base>/logs`, сохраняя схему
`<program>_<YYYYMMDD>.log`. Тесты и операционные регламенты опираются на этот
паттерн, поэтому менять или иначе ротировать базовое имя файла нельзя.

### Порядок проверки детерминизма

Продвижение пайплайнов в прод сопровождается обязательным чек-листом:

1. Выполните `python scripts/run_tests.py`, убедившись, что суммарный успех
   pytest ≥95 %.
2. Запустите нужный CLI по чистому каталогу вывода.
3. Немедленно повторите ту же команду или вызовите
   `python scripts/check_determinism.py --no-dry-run` с нужным входным CSV,
   чтобы получить второй экспорт.
4. Сравните результаты (скрипт вычисляет SHA256 для CSV и `.meta.yaml`). Любое
   расхождение трактуется как регрессия детерминизма.
5. В средах только для чтения используйте
   `python scripts/check_determinism.py --dry-run`. Скрипт хэширует объединённые
   stdout/stderr двух последовательных запусков и выводит
   `Dry-run log hash check: matched`, если план выполнения детерминирован.

Благодаря контракту логов обе итерации пишут в один и тот же файл
`<program>_<YYYYMMDD>.log`, что упрощает анализ расхождений при диффе событий.

### Быстрая шпаргалка по CLI

| Команда | Пример запуска | Особенности |
|---------|----------------|-------------|
| Оркестратор | `get-data --base-path . --input-dir data/input --output-dir output --config config/config.yaml --date 20250228 --limit 100 --dry-run` | Запускает всю цепочку один раз, сохраняет манифесты прогонов и понимает оркестратор-специфичные флаги (`--pipeline-registry`, `--override-{input,output-stem,subcommand}`). Обёртка `python scripts/get_data.py` прокидывает только общие флаги этапов (`--limit`, `--force`, `--skip-existing`, `--dry-run`). |
| Document | `python scripts/get_document_data.py --mode all --input data/input/document.csv --final-out output/documents.csv --fallback-doi-enabled --fallback-doi-path data/input/fallback.csv --openalex-rps 2` | Поддерживает режимы `chembl`, `pubmed`, `all`, настройку размера батчей и CSV с резервными DOI. |
| Target | `python scripts/get_target_data.py all --input data/input/target.csv --final-out output/targets.csv --chembl-chunk-size 10 --uniprot-data-dir cache/uniprot --raw-out output/targets_raw.parquet --raw-format parquet` | Подкоманды (`uniprot`, `chembl`, `iuphar`, `all`) принимают префиксные оверрайды и позволяют сохранять «сырые» выгрузки. |
| Assay | `python scripts/get_assay_data.py --input data/input/assay.csv --final-out output/assay.csv --chunk-size 25 --timeout 45` | Требует словари assay, taxonomy и target в `config/dictionary` для обогащения полей `assay_group`, `assay_strain`, `year` и `accession` перед нормализацией; общие флаги плюс настройка размера пачки и таймаута запросов. |
| Test item | `python scripts/get_testitem_data.py --input data/input/testitem.csv --final-out output/testitems.csv --request-limit 500 --hierarchy-path config/dictionary/_testitem/molecule_hierarchy.csv` | Управляет обогащением родительских молекул и лимитами запросов (`--request-limit`, `--batch-size`, `--dry-run`). |
| Tissue | `python scripts/get_tissue_data.py --input data/input/tissue.csv --final-out output/tissues.csv --chunk-size 50 --xref-sources uberon,efo,bto` | Загружает метаданные тканей, объединяет онтологические кросс-ссылки и нормализует синонимы. Запускается отдельно перед `get_activity_data`, когда нужны справочники тканей. |
| Cell line | `python scripts/get_cellline_data.py --input data/input/cellline.csv --final-out output/cellline.csv --batch-size 20 --limit 100` | Выгружает данные по клеточным линиям из ChEMBL, нормализует идентификаторы и формирует стабильный CSV. |
| Activity | `python scripts/get_activity_data.py --input data/input/activity.csv --final-out output/activities.csv --column activity_id --batch-size 10 --workers 4 --dry-run` | Флаги: переопределение колонки идентификаторов (`--column activity_id`), настройка батчей и таймаутов (`--batch-size`, `--timeout`), ограничение диапазона (`--limit`, `--offset`), dry-run и количество потоков. |
| Синтетические активности | `python scripts/get_activities.py --limit 25 --dry-run` | Генерирует детерминированные тестовые строки для смоук-тестов и поддерживает те же флаги логирования, что и остальные CLI. |

Каждый пайплайн теперь оставляет детерминированный CSV вместе с
`<имя>.meta.yaml`, `<имя>_quality_report_table.csv` и
`<имя>_data_correlation_report_table.csv` в том же каталоге. Метаданные пишутся
через `io.save_metadata` и содержат параметры запуска, схему и контрольные
суммы. Наследуемые диагностические файлы (`<имя>.quality.json`,
`<имя>_failure_cases.csv` и др.) подключаются по требованию через
`--emit-legacy-artifacts`, `--debug` или `--keep-intermediate`. Поле
`generated_at` остаётся детерминированным: при наличии `--date` используется
указанная дата, иначе хэшируется нормализованный вызов CLI и вычисленный
`run_id`. Таргет-пайплайн также создаёт вспомогательные таблицы
`organism.output.target_<stamp>.csv`, `isoform.output.target_<stamp>.csv`,
`names.output.target_<stamp>.csv` и `IUPHAR.output.target_<stamp>.csv`, которые
подробно описаны в [`docs/ru/OUTPUT_TARGETS.md`](./docs/ru/OUTPUT_TARGETS.md) и
[`docs/en/OUTPUT_TARGETS.md`](./docs/en/OUTPUT_TARGETS.md). Полную спецификацию
см. в [`docs/ru/OUTPUT.md`](./docs/ru/OUTPUT.md).

#### Воспроизведение архива экспорта таргетов

Исторические примеры набора артефактов теперь находятся в каталоге
[`reports/archive/target_pipeline/`](./reports/archive/target_pipeline). Чтобы
воссоздать ту же структуру локально:

1. Активируйте виртуальное окружение и установите зависимости по инструкции из
   раздела [«Быстрый старт»](#быстрый-старт).
2. Запустите пайплайн таргетов на поставляемых примерах идентификаторов:

   ```bash
   python scripts/get_target_data.py all \
     --input data/input/target.csv \
     --final-out output/targets.csv \
     --chembl-chunk-size 10 \
     --uniprot-data-dir cache/uniprot
   ```

3. Проверьте содержимое `output/targets.csv` и QA-отчётов
   `_quality_report_table.csv`/`_data_correlation_report_table.csv`. При
   необходимости добавьте `--emit-legacy-artifacts` (или `--debug`/`--keep-intermediate`),
   чтобы восстановить исторические YAML-файлы метаданных, CSV с ошибками и
   дополнительные вспомогательные выгрузки для диагностики.

Все артефакты сохраняют описанную выше детерминированность: повторные запуски с
одинаковыми входами дают байтово идентичные файлы.

## Документация

Все руководства доступны на двух языках. Структура зеркальна:

- Протокол (DOCX, источник в Markdown): формируется по требованию скриптом
  `scripts/convert_md_to_docx.py` в файл
  `docs/ChEMBL_Data_Acquisition_Protocol_v2.1.docx` (артефакт поставки,
  не хранится в git). Для сборки перед публикацией используйте `make protocol-docx`.
  Исходники: [`docs/ru/PROTOCOL_RU.md`](./docs/ru/PROTOCOL_RU.md)
  и синхронизированная английская версия
  [`docs/en/PROTOCOL_EN.md`](./docs/en/PROTOCOL_EN.md).
- Обзор и оглавление: [`docs/ru/README.md`](./docs/ru/README.md),
  [`docs/en/README.md`](./docs/en/README.md)
- Использование и CLI: [`docs/ru/USAGE.md`](./docs/ru/USAGE.md),
  [`docs/en/USAGE.md`](./docs/en/USAGE.md)
- Руководства (расширенные сценарии, отладка, FAQ):
  [`docs/ru/guides/ADVANCED_SCENARIOS.md`](./docs/ru/guides/ADVANCED_SCENARIOS.md),
  [`docs/ru/guides/DEBUGGING.md`](./docs/ru/guides/DEBUGGING.md),
  [`docs/ru/guides/FAQ.md`](./docs/ru/guides/FAQ.md) и зеркальные английские версии в
  `docs/en/guides/`
- Руководство по постобработке: [`docs/ru/guides/POSTPROCESSING_RUNBOOK.md`](./docs/ru/guides/POSTPROCESSING_RUNBOOK.md),
  [`docs/en/guides/POSTPROCESSING_RUNBOOK.md`](./docs/en/guides/POSTPROCESSING_RUNBOOK.md)
- Конфигурация: [`docs/ru/CONFIG.md`](./docs/ru/CONFIG.md),
  [`docs/en/CONFIG.md`](./docs/en/CONFIG.md)
- Спецификация выходных данных и правила валидации:
  [`docs/ru/OUTPUT.md`](./docs/ru/OUTPUT.md), [`docs/en/OUTPUT.md`](./docs/en/OUTPUT.md)
- Архитектура и модель данных:
  [`docs/ru/architecture/ARCHITECTURE.md`](./docs/ru/architecture/ARCHITECTURE.md),
  [`docs/en/architecture/ARCHITECTURE.md`](./docs/en/architecture/ARCHITECTURE.md)
- Руководство по разработке и CI/CD:
  [`docs/ru/development/README.md`](./docs/ru/development/README.md),
  [`docs/en/development/README.md`](./docs/en/development/README.md)

## Политика тестирования

Тесты располагаются в `tests/` и запускаются через `pytest`. Локальные и CI-запуски
должны создавать (файлы добавлены в `.gitignore`, чтобы не засорять историю):

- `reports/test_report.json` — машинно читаемый протокол
- `reports/test_summary.md` — краткое Markdown-резюме

GitHub Actions публикует оба файла (а также каталог покрытия) как артефакт
`test-reports-<python-version>` для каждой записи матрицы. Откройте последний
прогон, скачайте архив и изучите JSON/Markdown, чтобы увидеть актуальное
состояние пайплайна без локальной генерации отчётов.

`test_report.json` всегда содержит три корневых секции:

```json
{
  "meta": {
    "repo": "SatoryKono/ChEMBL_data_acquisition",
    "commit": "<SHA>",
    "branch": "<branch>",
    "ts_utc": "<ISO8601>",
    "duration_sec": 0.0,
    "python": "3.11|3.12",
    "pytest": "<version>",
    "exit_code": 0
  },
  "summary": {
    "total": 0,
    "passed": 0,
    "failed": 0,
    "skipped": 0,
    "xfailed": 0,
    "xpassed": 0,
    "error": 0,
    "success_rate": 0.0
  },
  "tests": [
    {
      "nodeid": "tests/unit/test_module.py::test_case",
      "status": "passed",
      "duration_ms": 12.3,
      "stdout": "",
      "stderr": "",
      "log": [],
      "error": null
    }
  ]
}
```

`test_summary.md` дублирует агрегированные показатели и для каждого падения или
ошибки выводит значение поля `error` из JSON в виде блока кода. Этого достаточно,
чтобы диагностировать проблему, имея только Markdown-отчёт.

Для смоук-прогона подойдёт `pytest -q -k "not slow and not e2e"`, полный набор —
`pytest -q`. Детали фикстур, требований к детерминизму и целям по покрытию см. в
[`docs/ru/development/TESTING.md`](./docs/ru/development/TESTING.md).

## Лицензия

Проект распространяется по [лицензии MIT](./LICENSE).

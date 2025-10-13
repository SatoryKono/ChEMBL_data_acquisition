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

Каждый пайплайн идемпотентен и может запускаться независимо. Оркестратор
[`get-data`](./scripts/get_data.py) использует единую конфигурацию и настройки
логирования, чтобы выполнить цепочку («документы → таргеты → ассайи → тестовые
объекты → активности») автоматически и воспроизводимо. Когда нужны связи по
тканям, `get_tissue_data` запускают отдельно, чтобы обновить справочные таблицы
перед запуском пайплайна активностей.

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
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-lock.txt
pip install .[dev]
pre-commit install
```

Изучите флаги оркестратора и отдельных пайплайнов:

```bash
python scripts/get_data.py --help
python scripts/get_document_data.py --help
```

Полный прогон на образцах идентификаторов с выводом в `./output`:

```bash
python scripts/get_data.py \
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

Команда запускает полный набор pytest с покрытием, формирует структурированный
отчёт `reports/test_report.json` и Markdown-выжимку `reports/test_summary.md`, а
при падении показателя `summary.success_rate` ниже **95 %** завершает выполнение
с ненулевым кодом. После правок тестов или кода пайплайна запустите скрипт
повторно и убедитесь, что отчёты идентичны за исключением временной метки —
любое расхождение указывает на потерю детерминизма, требующую исправления до
отправки MR. Подготавливая PR, приложите полученные артефакты, чтобы ревьюеры
смогли проверить фактический процент прохождений без локального прогона.

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

Благодаря контракту логов обе итерации пишут в один и тот же файл
`<program>_<YYYYMMDD>.log`, что упрощает анализ расхождений при диффе событий.

### Быстрая шпаргалка по CLI

| Команда | Пример запуска | Особенности |
|---------|----------------|-------------|
| Оркестратор | `python scripts/get_data.py --base-path . --input-dir data/input --output-dir output --config config/config.yaml --date 20250228 --limit 100 --dry-run` | Запускает всю цепочку один раз, прокидывая `--limit`, `--force`, `--skip-existing` и `--dry-run` на отдельные этапы. Дополнительно поддерживает `--pipeline-registry` для загрузки альтернативных определений шагов и `--override-{input,output-stem,subcommand}` для точечных переопределений. |
| Document | `python scripts/get_document_data.py --mode all --input data/input/document.csv --final-out output/documents.csv --fallback-doi-enabled --fallback-doi-path data/input/fallback.csv --openalex-rps 2` | Поддерживает режимы `chembl`, `pubmed`, `all`, настройку размера батчей и CSV с резервными DOI. |
| Target | `python scripts/get_target_data.py all --input data/input/target.csv --final-out output/targets.csv --chembl-chunk-size 10 --uniprot-data-dir cache/uniprot --raw-out output/targets_raw.parquet --raw-format parquet` | Подкоманды (`uniprot`, `chembl`, `iuphar`, `all`) принимают префиксные оверрайды и позволяют сохранять «сырые» выгрузки. |
| Assay | `python scripts/get_assay_data.py --input data/input/assay.csv --final-out output/assay.csv --chunk-size 25 --timeout 45` | Требует словари assay, taxonomy и target в `config/dictionary` для обогащения полей `assay_group`, `assay_strain`, `year` и `accession` перед нормализацией; общие флаги плюс настройка размера пачки и таймаута запросов. |
| Test item | `python scripts/get_testitem_data.py --input data/input/testitem.csv --final-out output/testitems.csv --batch-size 250 --pubchem-enable --postprocess` | Включает обогащение родительских молекул с переключателями PubChem (`--pubchem-enable`/`--no-pubchem-enable`) и детерминированной постобработкой (`--postprocess`/`--no-postprocess`), объединяющей иерархии и QA-артефакты. |
| Tissue | `python scripts/get_tissue_data.py --input data/input/tissue.csv --final-out output/tissues.csv --chunk-size 50 --xref-sources uberon,efo,bto` | Загружает метаданные тканей, объединяет онтологические кросс-ссылки и нормализует синонимы. Запускается отдельно перед `get_activity_data`, когда нужны справочники тканей. |
| Cell line | `python scripts/get_cellline_data.py --input data/input/cellline.csv --final-out output/cellline.csv --batch-size 20 --limit 100` | Выгружает данные по клеточным линиям из ChEMBL, нормализует идентификаторы и формирует стабильный CSV. |
| Activity | `python scripts/get_activity_data.py --input data/input/activity.csv --final-out output/activities.csv --column activity_id --batch-size 10 --workers 4 --dry-run` | Флаги: переопределение колонки идентификаторов (`--column activity_id`), настройка батчей и таймаутов (`--batch-size`, `--timeout`), ограничение диапазона (`--limit`, `--offset`), dry-run и количество потоков. |
| Синтетические активности | `python scripts/get_activities.py --limit 25 --dry-run` | Генерирует детерминированные тестовые строки для смоук-тестов и поддерживает те же флаги логирования, что и остальные CLI. |

Каждый пайплайн теперь оставляет детерминированный CSV вместе с
`<имя>_quality_report_table.csv` и `<имя>_data_correlation_report_table.csv` в
том же каталоге. Метаданные (`<имя>.meta.yaml`), JSON-отчёты качества и CSV с
ошибками доступны по требованию через `--emit-legacy-artifacts`, `--debug` или
`--keep-intermediate`. Поле `generated_at` остаётся детерминированным: при
наличии `--date` используется указанная дата, иначе хэшируется нормализованный
вызов CLI и вычисленный `run_id`. Таргет-пайплайн также создаёт вспомогательные
таблицы `organism.output.target_<stamp>.csv`, `isoform.output.target_<stamp>.csv`,
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

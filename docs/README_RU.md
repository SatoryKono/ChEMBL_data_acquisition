# Набор инструментов ChEMBL Data Acquisition

Этот документ даёт целостное представление о проекте: какие пайплайны доступны,
как они связаны между собой и какие артефакты формируются на выходе. Используйте
его как стартовую точку перед углублением в специализированные руководства.

## Ключевые особенности

- **Пайплайны сущностей** для документов, таргетов, ассайев, активностей и
  тест-объектов. Каждый пайплайн читает идентификаторы из CSV, обращается к
  внешним API, проводит детерминированную валидацию и записывает проверенный
  экспорт.
- **Оркестратор** (`get-data`), который последовательно запускает все пайплайны
  с общими настройками и единым логированием.
- **Унифицированный CLI**: общие флаги (`--input`, `--final-out`, `--log-level`,
  `--config` и др.), а также опции стадийных снимков (`--raw-out`,
  `--raw-format`, `--id-cols`, `--no-reindex-raw`) для таргет-пайплайна. Точки
  входа объявлены в `pyproject.toml`.
- **Три источника конфигурации** — `config/config.yaml`, переменные окружения с
  префиксом `CHEMBL_DA__` и параметры командной строки. Для локальных запусков
  можно загрузить `.env` через `python-dotenv`.
- **Детерминированный вывод** — CSV-файлы записываются в фиксированном порядке,
  рядом создаются YAML-метаданные с контрольной суммой SHA-256 и отчёты качества.
- **Инструменты разработки** — строгая типизация (`mypy`), линтинг (`ruff`),
  форматирование (`black`), тесты (`pytest`), подсчёт покрытия и проверка
  детерминизма в CI.

## Структура репозитория

| Путь | Назначение |
|------|-----------|
| `scripts/` | CLI-обёртки для отдельных пайплайнов и оркестратора `get-data`. |
| `library/` | Основная логика: HTTP-клиенты, оркестрация, нормализация, валидация, пост-обработка, ввод/вывод и метаданные. |
| `library/cli/commands/` | Точки входа, регистрируемые при установке пакета. |
| `library/utils/cli_tools/` | Лёгкие утилиты (профилирование таблиц, кешированная обвязка таргетов, CSV-хелперы, маппинг). |
| `config/` | Базовый YAML-конфиг, схема и встроенные словари. |
| `dictionary/` | Справочные наборы данных (кеши UniProt, классификаторы, QA-файксы). |
| `data/` | Примеры входов и тестовые выгрузки. |
| `docs/` | Документация проекта (варианты `_EN.md` и `_RU.md`). |
| `tests/` | Юнит- и интеграционные тесты пайплайнов и утилит. |
| `Makefile` | Сценарии для форматирования, тестов, сборки и публикации. |

## Поддерживаемые команды

Установка (`pip install .` или установка wheel) регистрирует следующие консольные
скрипты. Они соответствуют модулям из `scripts/`, `library/cli/commands/` или
`library/utils/cli_tools/` и принимают те же аргументы, что и вызовы
`python -m …`.

| Консольная команда | Модуль | Описание |
|--------------------|--------|----------|
| `get-data` | `scripts.get_data:main` | Запускает все пайплайны последовательно с общими настройками. |
| `get-document-data` | `library.cli.commands.get_document_data:main` | Собирает и обогащает документы ChEMBL. |
| `get-target-data` | `library.cli.commands.get_target_data:main` | Аггрегирует таргеты из ChEMBL, UniProt и IUPHAR. |
| `get-assay-data` | `library.cli.commands.get_assay_data:main` | Выгружает метаданные ассайев. |
| `get-activity-data` | `library.cli.commands.get_activity_data:main` | Выгружает нормализованные активности. |
| `get-testitem-data` | `library.cli.commands.get_testitem_data:main` | Обогащает молекулы данными PubChem. |
| `get-document-type` | `library.utils.cli_tools.get_document_type:main` | Классифицирует публикации для QA. |
| `csv-utils` | `library.utils.cli_tools.csv_utils_main:main` | Детерминированная пересериализация CSV. |
| `mapper` | `library.utils.cli_tools.mapper_main:main` | Интерактивный маппер UniProt/ChEMBL. |
| `table-quality` | `library.utils.cli_tools.table_quality_main:main` | Формирует отчёты качества по колонкам. |
| `chunk-io` | `library.utils.cli_tools.chunk_io_main:main` | Потоковое чтение/запись CSV с сохранением порядка. |
| `get-input-initialisation` | `library.utils.cli_tools.get_input_initialisation:main` | Объединяет Excel-книги инициализации. |
| `get-activities` | `library.utils.cli_tools.get_activities:main` | Генерирует синтетические активности для смоук-тестов. |
| `check-determinism` | `library.utils.cli_tools.check_determinism:main` | Сравнивает хэши CSV между запусками. |

Подробные сценарии приведены в [`docs/USAGE_RU.md`](./USAGE_RU.md) и английской
версии [`docs/USAGE_EN.md`](./USAGE_EN.md).

## Требования и установка

| Компонент | Минимум | Последняя проверка |
|-----------|---------|--------------------|
| Python | 3.11 | 3.12 |
| numpy | 1.26 | 2.3.3 |
| pandas | 2.0 | 2.3.3 |
| requests | 2.31 | 2.32.5 |
| PyYAML | 6.0 | 6.0.3 |

Фиксированные зависимости перечислены в `requirements-lock.txt`. Пересобирайте
его только после изменения диапазонов в `pyproject.toml`, используя чистое
виртуальное окружение и `pip freeze`.

```bash
python -m pip install --upgrade pip setuptools wheel
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
git clone https://github.com/SatoryKono/ChEMBL_data_acquisition.git
cd ChEMBL_data_acquisition
pip install -r requirements-lock.txt
pre-commit install
```

При установке из wheel конфигурация копируется в пользовательские каталоги,
перечисленные в [`docs/CONFIG_RU.md`](./CONFIG_RU.md); консольные команды
регистрируются автоматически.

## Быстрый старт

1. **Подготовьте списки идентификаторов.** Используйте шаблоны из `data/input`
   или выгрузите свежие списки из своей системы. Пайплайны ожидают по одному ID
   в строке; имена колонок перечислены в [`docs/DATA_SCHEMA_RU.md`](./DATA_SCHEMA_RU.md).
2. **Проверьте конфигурацию.** Скопируйте `config/config.yaml`, если нужно
   переопределить лимиты API, каталоги вывода или стадийные флаги. Паттерн
   переменных окружения — `CHEMBL_DA__РАЗДЕЛ__КЛЮЧ`. Полное описание — в
   [`docs/CONFIG_RU.md`](./CONFIG_RU.md).
3. **Запустите пайплайн.**

   ```bash
   get-document-data all \
       --input data/input/document_ids.csv \
       --final-out output.documents_$(date +%Y%m%d).csv \
       --config config/config.yaml \
       --log-level INFO
   ```

   Для смоук-теста добавьте `--limit 10`. Таргет-пайплайн поддерживает
   стадийные опции `--raw-out` и `--raw-format parquet` для «сырых» снимков.
4. **Изучите артефакты.** Рядом с CSV появится `<имя>.meta.yaml`, отчёты качества
   и (для документов) JSON-резюме. Форматы описаны в [`docs/OUTPUT_RU.md`](./OUTPUT_RU.md).

## Тестирование и контроль качества

Перед коммитом или публикацией выполните:

```bash
pre-commit run --all-files
pip check
pytest
pytest --cov=library --cov=scripts --cov-report=term-missing
python -m library.utils.cli_tools.check_determinism --log-level DEBUG \
    --input out/latest.csv --previous out/previous.csv
```

QA-плейбук (`docs/QA_PROCESS_RU.md` / `docs/QA_PROCESS_EN.md`) описывает ворота
качества, смоук-проверки и критерии приёмки. Проверки детерминизма полагаются на
YAML-метаданные, поэтому держите их в Git для удобного сравнения запусков.

## Дополнительные материалы

- [`docs/USAGE_RU.md`](./USAGE_RU.md) / [`docs/USAGE_EN.md`](./USAGE_EN.md) —
  флаги CLI, подкоманды и типовые сценарии.
- [`docs/CONFIG_RU.md`](./CONFIG_RU.md) / [`docs/CONFIG_EN.md`](./CONFIG_EN.md) —
  источники конфигурации и пути вывода.
- [`docs/OUTPUT_RU.md`](./OUTPUT_RU.md) / [`docs/OUTPUT_EN.md`](./OUTPUT_EN.md) —
  структура артефактов и работа с «сырыми» снимками.
- [`docs/DATA_SCHEMA_RU.md`](./DATA_SCHEMA_RU.md) /
  [`docs/DATA_SCHEMA_EN.md`](./DATA_SCHEMA_EN.md) — описание колонок и схем
  валидации.
- [`docs/ETL_PROCESS_RU.md`](./ETL_PROCESS_RU.md) /
  [`docs/ETL_PROCESS_EN.md`](./ETL_PROCESS_EN.md) — поток данных end-to-end.
- [`docs/CLI_TOOLS_RU.md`](./CLI_TOOLS_RU.md) /
  [`docs/CLI_TOOLS_EN.md`](./CLI_TOOLS_EN.md) — справочник по вспомогательным
  утилитам.

Диаграммы архитектуры смотрите в [`docs/ARCHITECTURE_RU.md`](./ARCHITECTURE_RU.md)
и английской версии [`docs/ARCHITECTURE_EN.md`](./ARCHITECTURE_EN.md).

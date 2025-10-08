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
| `config/dictionary/` | Справочные наборы данных (кеши UniProt, классификаторы, QA-файксы). |
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
| `get-activities` | `library.utils.cli_tools.get_activities:main` | Генерирует синтетические активности и детерминированные CSV + `.meta.yaml` артефакты для смоук-тестов. |
| `check-determinism` | `library.utils.cli_tools.check_determinism:main` | Сравнивает хэши CSV между запусками. |

Подробные сценарии приведены в [`USAGE.md`](./USAGE.md) и английской
версии [`../en/USAGE.md`](../en/USAGE.md).

## Требования и установка

| Компонент | Поддерживаемый диапазон | Последняя проверка |
|-----------|-------------------------|--------------------|
| Python | 3.11.x-3.12.x | 3.12.3 |
| numpy | >=2.3.3,<3.0 | 2.3.3 |
| pandas | >=2.3.3,<3.0 | 2.3.3 |
| requests | >=2.32.5,<3.0 | 2.32.5 |
| PyYAML | >=6.0.3,<7.0 | 6.0.3 |
| openpyxl | >=3.1.5,<4.0 | 3.1.5 |
| pyarrow | >=17.0.0,<18.0 | 17.0.0 |
| jsonschema | >=4.25.1,<5.0 | 4.25.1 |
| pandera | >=0.26.1,<0.27 | 0.26.1 |
| pydantic | >=2.11.9,<3.0 | 2.11.9 |

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
перечисленные в [`CONFIG.md`](./CONFIG.md); консольные команды
регистрируются автоматически.

## Быстрый старт

1. **Подготовьте списки идентификаторов.** Используйте шаблоны из `data/input`
   (например, `document.csv`, `target.csv`, `assay.csv`, `activity.csv`,
   `testitem.csv`) или выгрузите свежие списки из своей системы. Расширенный
   набор примеров лежит в `data/input/full`. Пайплайны ожидают по одному ID в
   строке; имена колонок перечислены в [`DATA_SCHEMA.md`](./DATA_SCHEMA.md).
2. **Проверьте конфигурацию.** Скопируйте `config/config.yaml`, если нужно
   переопределить лимиты API, каталоги вывода или стадийные флаги. Паттерн
   переменных окружения — `CHEMBL_DA__РАЗДЕЛ__КЛЮЧ`. Полное описание — в
   [`CONFIG.md`](./CONFIG.md).
3. **Запустите пайплайн.**

   ```bash
   get-document-data all \
       --input data/input/document.csv \
       --final-out output.documents_$(date +%Y%m%d).csv \
       --config config/config.yaml \
       --log-level INFO
   ```

   Для смоук-теста добавьте `--limit 10`. Таргет-пайплайн поддерживает
   стадийные опции `--raw-out` и `--raw-format parquet` для «сырых» снимков.
4. **Изучите артефакты.** Рядом с CSV появится `<имя>.meta.yaml`, отчёты качества
   и (для документов) JSON-резюме. Форматы описаны в [`OUTPUT.md`](./OUTPUT.md).

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

QA-плейбук (`QA_PROCESS.md` / `../en/QA_PROCESS.md`) описывает ворота
качества, смоук-проверки и критерии приёмки. Проверки детерминизма полагаются на
YAML-метаданные, поэтому держите их в Git для удобного сравнения запусков.

## Дополнительные материалы

- [`USAGE.md`](./USAGE.md) / [`../en/USAGE.md`](../en/USAGE.md) —
  флаги CLI, подкоманды и типовые сценарии.
- [`CONFIG.md`](./CONFIG.md) / [`../en/CONFIG.md`](../en/CONFIG.md) —
  источники конфигурации и пути вывода.
- [`OUTPUT.md`](./OUTPUT.md) / [`../en/OUTPUT.md`](../en/OUTPUT.md) —
  структура артефактов и работа с «сырыми» снимками.
- [`DATA_SCHEMA.md`](./DATA_SCHEMA.md) /
  [`../en/DATA_SCHEMA.md`](../en/DATA_SCHEMA.md) — описание колонок и схем
  валидации.
- [`ETL_PROCESS.md`](./ETL_PROCESS.md) /
  [`../en/ETL_PROCESS.md`](../en/ETL_PROCESS.md) — поток данных end-to-end.
- [`CLI_TOOLS.md`](./CLI_TOOLS.md) /
  [`../en/CLI_TOOLS.md`](../en/CLI_TOOLS.md) — справочник по вспомогательным
  утилитам.

Диаграммы архитектуры смотрите в [`architecture/ARCHITECTURE.md`](./architecture/ARCHITECTURE.md)
и английской версии [`../en/architecture/ARCHITECTURE.md`](../en/architecture/ARCHITECTURE.md).

# ChEMBL Data Acquisition Utilities

Набор утилит и библиотек Python 3.12 для скачивания, валидации,
агрегации и экспорта биологических данных из открытых API  
(ChEMBL, PubChem, UniProt, PubMed и др.). Проект демонстрирует
типичный пайплайн обработки табличных данных: от получения
идентификаторов до сериализации нормализованных CSV/Parquet
с сопровождающими метаданными.

## Особенности

* Командные скрипты с унифицированными флагами `--input`, `--output`,
  `--log-level`, `--sep`, `--encoding`, `--column`, `--dictionary`.
* Потоковая обработка больших CSV через чанки, детерминированный вывод.
* Валидаторы схем (`schemas/`) и словари (`dictionary/`) для проверки
  типов, диапазонов и справочников.
* Конфигурация через `config.yaml`, переменные окружения и ключи CLI.
* Логирование через стандартный модуль `logging` с настраиваемым уровнем.
* Полная статическая типизация (PEP 484), линтинг `ruff`, форматирование
  `black`, проверка типов `mypy`, юнит‑тесты `pytest`.

## Требования

| Компонент     | Минимальная версия |
|---------------|-------------------|
| Python        | 3.12              |
| pandas        | 2.1               |
| requests      | 2.31              |
| PyYAML        | 6.0               |

Полный список приведён в `requirements-dev.txt` или `pyproject.toml`.

## Установка

```bash
git clone https://github.com/<org>/ChEMBL_data_acquisition.git
cd ChEMBL_data_acquisition
pip install .[dev]       # с инструментами разработки
pre-commit install       # git‑хуки: black/ruff/mypy/pytest
pre-commit run --all-files
```

## Быстрый старт

```bash
# загрузка активности по идентификаторам из тестового CSV
python -m scripts.get_activity_data \
    --input tests/data/activity_ids_small.csv \
    --output out/activities.csv \
    --limit 10 --log-level INFO

# профилирование качества таблицы
python table_quality_main.py \
    --input tests/data/activity.csv \
    --table-name activity
```

`--output` по умолчанию формируется как `output_<имя_входа>_YYYYMMDD.csv`
в каталоге, заданном `io.output_dir`.  
Для дополнительных примеров см. [`docs/USAGE.md`](docs/USAGE.md).

## Структура проекта

```
ChEMBL_data_acquisition/
├── config.yaml
├── dictionary/
├── library/
│   ├── __init__.py
│   ├── chembl_client.py
│   ├── csv_utils.py
│   ├── config.py
│   └── ...
├── scripts/
│   ├── get_activity_data.py
│   ├── get_assay_data.py
│   ├── ...
├── tests/
│   └── data/
└── docs/
    ├── CONFIG.md
    ├── OUTPUT.md
    └── USAGE.md
```

## Конфигурация

Параметры читаются из `config.yaml`, переменных окружения
(`CHEMBL_DA__SECTION__KEY`) и ключей CLI.  
Подробности в [`docs/CONFIG.md`](docs/CONFIG.md).

## Вывод и метаданные

Все сгенерированные CSV/Parquet и отчёты сохраняются в `data/output`
(см. [`docs/OUTPUT.md`](docs/OUTPUT.md)).  
Рядом создаются файлы `*.meta.yaml` с коммитом Git, параметрами запуска,
контрольной суммой SHA‑256 и статистикой строк/колонок.

## Разработка и тестирование

```bash
ruff check scripts library mapper_main.py table_quality_main.py
black scripts library mapper_main.py table_quality_main.py
mypy scripts library mapper_main.py table_quality_main.py
pytest
```

Тестовые наборы расположены в `tests/data`.  
Скрипт `scripts/check_determinism.py` проверяет повторяемость CSV‑вывода.

## Лицензия

MIT License. См. файл `LICENSE` (если присутствует).


# Структура выходных данных

## Базовый каталог

Все наборы сохраняются в `local.io.output_dir` (по умолчанию `data/output`). Если аргумент `--output` не задан, функции
используют `library.io.default_output_path` и формируют имя `output_<stem входного файла>_<YYYYMMDD>.csv` в указанном каталоге.
При `local.io.exist_ok=true` необходимые папки создаются автоматически.

Пример:

```
data/output/
└── ChEMBL/
    └── processed/
        ├── activity.csv
        ├── activity.csv.meta.yaml
        ├── activity_failure_cases.csv
        ├── activity_quality_report_table.csv
        ├── activity_data_correlation_report_table.csv
        └── ...
```

## Метаданные sidecar

Каждый CSV сопровождается `<name>.csv.meta.yaml`, который записывает `library.metadata.write_meta_yaml`. Файл включает:

* `generated_at` — отметку времени в формате ISO 8601 (UTC).
* `git_sha` — хэш коммита репозитория на момент запуска.
* `python_version`, `platform` — сведения о среде выполнения.
* `command` — точную команду запуска.
* `config` — релевантные настройки (секреты маскируются).
* `inputs` — описание исходных файлов и аргументов.
* `stats` — показатели `rows_total`, `rows_kept`, `rows_dropped`, а также контрольную сумму `output_sha256`.
* `schema` — имя схемы валидации, применённой к данным.

Если sidecar уже существует, новый набор полей аккуратно объединяется с текущим содержимым.

## Артефакты валидации

* При ошибках проверки Pandera проблемные строки сохраняются в `<stem>_failure_cases.csv` рядом с основным файлом.
* `library.table_quality.analyze_table_quality` создаёт отчёты `<stem>_quality_report_table.csv` и
  `<stem>_data_correlation_report_table.csv`. CLI-утилиты размещают их рядом с экспортом, а
  `scripts/get_input_initialisation.py` — в подкаталоге `<output>/data_validity_report/`.

Все отчёты записываются в кодировке UTF-8 и наследуют детерминированный порядок строк и колонок.

## Детерминированный экспорт

`library.io.write_csv` вызывает `library.csv_utils.write_csv_deterministic`, сортируя колонки и строки по ключевым столбцам.
Поведение следует настройкам `cfg.io.csv_sep`, `cfg.io.csv_encoding` и учитывает аргументы `key_cols`/`col_order`.

## Рекомендации по ведению архива

* Выделяйте под каждый запуск отдельную папку с датой (`YYYYMMDD/`) для удобного сравнения.
* Архивируйте устаревшие результаты — метаданные сохраняют достаточно информации для воспроизводимости.
* Следите за свободным местом при длительных или параллельных выгрузках.

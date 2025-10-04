# Руководство по использованию

Это руководство объясняет, как запускать конвейеры получения данных ChEMBL и их вспомогательные утилиты. Каждый раздел имеет английский аналог в [`../en/USAGE.md`](../en/USAGE.md).

## Общий шаблон CLI

Все инструменты командной строки используют общий набор аргументов для единообразия. Обычно они делятся на три категории:

1.  **Общие параметры ввода-вывода и конфигурации:**
    *   `--input`: Путь к входному CSV-файлу, содержащему идентификаторы.
    *   `--final-out`: Путь к конечному, очищенному выходному файлу. Устаревшие псевдонимы `--output` и `--out` все еще доступны, но будут удалены в будущей версии.
    *   `--config`: Путь к файлу конфигурации YAML. По умолчанию используется встроенный `config/config.yaml`.
    *   `--log-level`: Устанавливает уровень детализации логов (например, `INFO`, `DEBUG`).
    *   `--print-config`: Печатает итоговую конфигурацию после всех переопределений и завершает работу.

2.  **Общие параметры управления выполнением:**
    *   `--limit`: Ограничивает количество обрабатываемых записей. Установка `--limit 0` — удобный способ проверить конфигурацию без запуска полного конвейера.
    *   `--batch-size` / `--chunk-size`: Количество записей, обрабатываемых в одном вызове API или пакете.
    *   `--workers`: Количество параллельных рабочих процессов для одновременных операций.

3.  **Специфичные для конвейера параметры:**
    *   Это аргументы, уникальные для конкретного конвейера, такие как `--raw-out` для конвейера таргетов или `--fallback-doi-csv` для конвейера документов.

---

## Основные конвейеры

Это основные скрипты для получения и обработки данных по ключевым сущностям ChEMBL.

### Оркестратор (`get-data`)

Это основная точка входа для последовательного запуска всех конвейеров получения данных. Он гарантирует, что на каждый шаг передается согласованная конфигурация.

```bash
get-data --base-path /data/chembl \
    --input-dir seeds --output-dir exports \
    --config /data/chembl/config.yaml \
    --date 20250101 --limit 100 --log-level INFO
```

### Конвейер документов (`get-document-data`)

Получает и обогащает данные о публикациях из ChEMBL, PubMed, CrossRef и других источников.

```bash
get-document-data all \
    --input seeds/document_ids.csv \
    --final-out output/documents.csv \
    --limit 500
```

### Конвейер таргетов (`get-target-data`)

Получает и обогащает данные о таргетах из ChEMBL, UniProt и IUPHAR. Этот конвейер имеет уникальные флаги для поэтапной обработки для более детального контроля над выходными данными.

**Флаги поэтапной обработки:**
*   `--raw-out`: Сохраняет промежуточные, ненормализованные данные в отдельный файл.
*   `--raw-format`: Формат для "сырого" вывода (`csv` или `parquet`).
*   `--id-cols`: Колонки с идентификаторами для детерминированной сортировки.

```bash
get-target-data all \
    --input seeds/target_ids.csv \
    --final-out output/targets_final.csv \
    --raw-out output/targets_raw.parquet \
    --raw-format parquet
```

### Конвейер ассеев (`get-assay-data`)

Загружает и обрабатывает метаданные ассеев из ChEMBL.

```bash
get-assay-data --input seeds/assay_ids.csv \
    --final-out output/assays.csv \
    --limit 200
```

### Конвейер активностей (`get-activity-data`)

Извлекает данные об активностях и вычисляет нормализованные границы значений.

```bash
get-activity-data --input seeds/activity_ids.csv \
    --final-out output/activities.csv \
    --limit 500
```

### Конвейер тестовых образцов (`get-testitem-data`)

Получает данные о молекулах и обогащает их свойствами из PubChem.

```bash
get-testitem-data --input seeds/molecule_ids.csv \
    --final-out output/testitems.csv \
    --limit 400
```

---

## Вспомогательные утилиты

Это вспомогательные инструменты для диагностики, манипулирования данными и контроля качества.

| Команда | Назначение | Пример |
|---|---|---|
| `check-determinism` | Сравнивает два CSV-файла, чтобы убедиться в их идентичности. | `check-determinism --input a.csv --previous b.csv` |
| `chunk-io` | Читает и перезаписывает CSV-файл детерминированными блоками. | `chunk-io --input data.csv --final-out copy.csv` |
| `csv-utils` | Нормализует форматирование (разделители, кавычки) CSV-файла. | `csv-utils --input data.csv --final-out clean.csv` |
| `dtype-inspector` | Проверяет и сообщает о типах данных pandas, создаваемых конвейерами. | `python -m library.utils.cli_tools.dtype_inspector_main` |
| `get-activities` | Генерирует синтетические данные об активностях для тестирования. | `get-activities --limit 10 --dry-run` |
| `get-document-type` | Применяет эвристики для классификации типов публикаций. | `get-document-type --input docs.csv` |
| `get-input-initialisation` | Объединяет рабочие книги Excel в канонические входные файлы. | `get-input-initialisation --same-doc init.xlsx` |
| `mapper` | Сопоставляет идентификаторы между ChEMBL и UniProt. | `mapper --input ids.csv --final-out mapped.csv` |
| `table-quality` | Генерирует отчет о качестве и корреляции для CSV-файла. | `table-quality --input data.csv --table-name my_data` |
| `pipeline-targets-main` | Повторяет выполнение конвейера таргетов, используя кэшированные данные. | `python -m library.utils.cli_tools.pipeline_targets_main --input targets.csv` |
# Расширенные сценарии

Набор типовых шаблонов для запуска конвейеров сверх базового smoke-теста.

## Пакетная обработка релизов

Используйте `--base-path` и датированные подкаталоги, чтобы разделять выгрузки:

```bash
BASE=/data/releases/chembl-35
poetry run get-data \
  --base-path "$BASE" \
  --input-dir inbound/2025-02-01 \
  --output-dir outbound/2025-02-01 \
  --config config/config.yaml \
  --date 20250201 \
  --log-level INFO
```

Храните файлы `.meta.yaml` в системе контроля версий для отслеживания происхождения.

## Частичный перезапуск с кешами

Перезапуск только тестовых объектов:

```bash
python scripts/get_testitem_data.py \
  --input /data/inbound/testitem.csv \
  --final-out /data/outbound/testitems.csv \
  --config config/config.yaml \
  --force
```

Кэши родителей молекул находятся в `$CHEMBL_DA_BASE_PATH/cache`. Для обновления
словари очистите их:

```bash
rm -rf "$CHEMBL_DA_BASE_PATH"/cache/molecule_*
```

## Снижение лимитов для стейджинга

Временно ограничьте RPS через переменные окружения:

```bash
export CHEMBL_DA__SOURCES__CHEMBL__API__RPS=5
export CHEMBL_DA__SOURCES__CHEMBL__API__BURST=5
python scripts/get_activity_data.py --limit 200
```

Аналогично можно настраивать партнёрские API (`CHEMBL_DA_OPENALEX_RPS`,
`CHEMBL_DA_PUBMED_TIMEOUT_READ` и т.д.).

## Использование промежуточных выгрузок

`get_target_data` умеет сохранять сырые наборы данных:

```bash
python scripts/get_target_data.py all \
  --input data/input/target.csv \
  --final-out output/targets.csv \
  --raw-out output/targets_raw.parquet \
  --raw-format parquet \
  --id-cols target_chembl_id
```

Просмотрите parquet, чтобы убедиться в корректности объединений до финальной нормализации.

## Кастомный реестр пайплайнов

Для альтернативных сценариев (например, пропустить обогащение документов или
добавить собственные QA-этапы) подготовьте YAML-реестр и передайте его через
`--pipeline-registry`:

```yaml
pipelines:
  - name: document
    callable: scripts.get_document_data:main
    input: document_subset.csv
    output: documents_subset
  - name: target
    callable: scripts.get_target_data:main
    subcommand: chembl
    output: targets_chembl_only
  - name: audit
    callable: tools.audit_pipeline:main
    input: targets_chembl_only.csv
    output: audit_report
```

Запуск с адресными переопределениями:

```bash
poetry run get-data \
  --base-path /data/chembl \
  --pipeline-registry config/custom_registry.yaml \
  --override-input document=document_snapshot.csv \
  --override-subcommand target=all
```

Опция `--override-output-stem` пригодится, когда нужен стандартный реестр, но
временные имена файлов следует перенаправить без правки YAML.

## Интеграция с Makefile

Полезные цели (`Makefile`):

- `make lint` — `ruff` + `black --check`.
- `make typecheck` — `mypy` по `library/` и `scripts/`.
- `make test` — полный `pytest` с JSON-отчётом.
- `make smoke` — последовательный запуск на тестовых данных.

Комбинируйте цели в CI для воспроизводимости.

## Экспорт метрик QA

Агрегируйте JSON-отчёты при помощи `jq`:

```bash
jq -s 'map({file: input_filename, stats: .summary})' output/*.quality.json > qa_summary.json
```

Получившийся файл можно загрузить в дашборды или инструменты мониторинга для
отслеживания трендов по пропускам и нарушениям схем.

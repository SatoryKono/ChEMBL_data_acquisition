# Аномалии вывода get_target_data.py

## Ожидаемое поведение

Согласно политике стандартизации вывода каждая команда семейства `get_*_data.py` должна сохранять ровно три CSV-артефакта с помощью `save_standard_outputs`: агрегированный датасет, отчёт по корреляциям и отчёт по контролю качества. После стандартного запуска не должно оставаться никаких других файлов.

## Почему появляются лишние артефакты

### 1. Обёртка CLI автоматически активирует устаревший режим

`scripts/get_data.py` при подготовке конфигурации выставляет `args.keep_intermediate = True`, если запуск происходит с флагом `--debug`. Затем `run_cli_command` трактует любые включённые `debug`/`keep_intermediate` как сигнал выставить `args.emit_legacy_artifacts = True`. Тем самым даже стандартный прогон, выполненный через оркестратор, внезапно получает весь набор исторических файлов (raw CSV, failure cases, `.meta.yaml`). 【F:scripts/get_data.py†L188-L200】【F:library/cli_utils.py†L303-L310】

**Подзадача:** убрать автоматическую установку `keep_intermediate` либо привязать её к явному флагу, чтобы без запроса не включался устаревший писатель.

### 2. Промежуточные выгрузки создают дополнительные CSV

Комбинированная команда (`target all`) всегда вызывает `fetch_chembl`, `fetch_uniprot` и `fetch_iuphar`. Имена файлов для этих этапов вычисляются на базе финального пути: `chembl_out = final_output.stem + "_chembl.csv"`, `uniprot_out = ... + "_uniprot.csv"`, `iuphar_out = ... + "_iuphar.csv"`. Таким образом рядом с целевым `output.targets_<date>.csv` появляются промежуточные `output.targets_<date>_chembl.csv`, `output.targets_<date>_uniprot.csv` и `output.targets_<date>_iuphar.csv`. 【F:library/cli/commands/get_target_data.py†L4091-L4130】

Особенно критично, что `fetch_chembl` запускает `_run_pipeline_with_meta(..., emit_standard_outputs=True)`, а внутри `run_pipeline` это приводит к вызову `save_standard_outputs`. Последняя формирует canonical-пакет `output.<table>_<date>.csv`, где `<table>` наследует `_chembl` из имени промежуточного файла. Поэтому в каталоге вывода остаются три «лишних» артефакта: `output.targets_<date>_chembl_<date>.csv`, отчёт по корреляциям и QC-отчёт с тем же суффиксом `_chembl`. 【F:library/cli/commands/get_target_data.py†L2253-L2275】【F:library/cli_utils.py†L1012-L1067】【F:library/io/output_writer.py†L109-L135】

**Подзадача (историческая):** выключить стандартный пакет (`emit_standard_outputs=False`) для промежуточного вызова ChEMBL и удалять/перемещать черновики после слияния, чтобы в каталоге остались только три финальные таблицы.

**Обновление:** начиная с версии пайплайна `0.1.3` `run_pipeline` требует включённый стандартный пакет или легаси-артефакты. Для промежуточной выгрузки ChEMBL удерживаем `emit_standard_outputs=True`, а удаление побочных CSV/QC отчётов выполняется отдельной задачей клининга.

### 3. Шаг валидации оставляет «сырые» и failure-артефакты при активном legacy

В `validate_and_write` есть ветки, которые при `emit_legacy_artifacts=True` пишут `raw_out` (CSV/Parquet), формируют `*_failure_cases.csv` и сохраняют sidecar через `SidecarErrors`. Именно эти участки генерируют `target_raw.csv`, `target_raw.meta.yaml` и связанные файлы. 【F:library/cli/commands/get_target_data.py†L3636-L3846】

**Подзадача:** либо удалить эти ветки, либо зафиксировать их за явным флагом.

### 4. Опциональный постпроцессинг создаёт дублирующий экспорт

Если CLI вызывается с `--postprocess`, функция `run_target_postprocess_if_requested` формирует `output_postprocessed.targets.csv` рядом с основным датасетом. Позже `validate_and_write` удаляет файл только в случае, если canonical CSV оказался в другом месте; иначе артефакт остаётся. 【F:library/cli/commands/get_target_data.py†L433-L518】【F:library/cli/commands/get_target_data.py†L3914-L3922】

**Подзадача:** писать результаты постобработки напрямую через `save_standard_outputs` или очищать временный файл сразу после чтения.

### 5. Sidecar от UniProt создавался даже без legacy-режима

Этап `run_uniprot` экспортирует данные через `io.write_csv`, который внутри вызывает `write_csv_deterministic`. Последний всегда пишет `.meta.yaml` рядом с CSV, поэтому появлялись файлы вроде `output.targets_<date>_uniprot.csv.meta.yaml` даже при обычном запуске без `--emit-legacy-artifacts`. Мы удаляем sidecar сразу после чтения промежуточного CSV, чтобы в каталоге оставались только стандартизованные артефакты. 【F:library/cli/commands/get_target_data.py†L1994-L2018】【F:library/io/writers.py†L13-L44】【F:library/common/csv_utils.py†L562-L569】

## Рекомендации по исправлению

* Удалить принудительную установку `keep_intermediate` в `scripts/get_data.py` либо привязать её к явному отладочному флагу, чтобы обычные запуски не включали устаревшие артефакты.
* Передавать `emit_standard_outputs=True, emit_legacy_artifacts=False` при вызове `_run_pipeline_with_meta` для «сырых» режимов; при необходимости «сырые» CSV можно по-прежнему получать потоковой записью без сохранения лишних файлов.
* Упростить `validate_and_write`, удалив устаревшие ветки сохранения (или спрятав их за явным флагом), чтобы в обычных сценариях выполнялся только `save_standard_outputs`.
* Перенаправить вывод постпроцессинга в общий механизм записи, чтобы на диске оставались только три стандартных файла.

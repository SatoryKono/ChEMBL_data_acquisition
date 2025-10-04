# Контрольный список выгрузки по таргетам

В документе описаны шаги постобработки, форматы входных/выходных файлов и
поведение при выполнении пайплайна таргетов. Это дополнение к общей
спецификации выходных данных и ориентировано на артефакты,
которые формирует `scripts/get_target_data.py`.

## Общее описание

Пайплайн сначала нормализует объединённый набор данных ChEMBL/UniProt/IUPHAR
через `postprocess_targets`, после чего выполняет финальное выравнивание: удаляет
дубликаты, вычисляет метки клеточности и проецирует таблицу на каноническую
схему перед записью CSV.【F:library/pipelines/target/postprocessing.py†L338-L485】【F:library/pipelines/target/postprocessing.py†L535-L730】
Запись выполняется детерминированным писателем `write_csv_deterministic`, что
обеспечивает стабильный порядок строк и столбцов.【F:library/pipelines/target/postprocessing.py†L523-L532】【F:library/pipelines/target/postprocessing.py†L723-L730】

## Формат входных данных

`postprocess_targets` принимает консолидированный результат режима `get_target_data all`.
Минимально необходимы столбцы с идентификаторами ChEMBL, UniProt и
линейджем (генус, суперцарство, тип, класс), а также дополнительные признаки из
словарей — например, `multifunctional_enzyme`. Ресурс `targets_type.csv`
демонстрирует состав колонок и содержимое таксономии.【F:library/integration/input_initialisation_library.py†L520-L576】【F:config/dictionary/_target/targets_type.csv†L1-L12】 Ниже приведены
первые шесть строк (усечённые до ключевых полей):

```
target_chembl_id,uniprotkb_Id,organism_type,multifunctional_enzyme
CHEMBL1827,PDE5A,Multicellular organism,FALSE
CHEMBL1859,CACNA1H,Multicellular organism,FALSE
CHEMBL202,P00374,Multicellular organism,TRUE
CHEMBL1809,P0ABQ4,Unicellular organism,TRUE
CHEMBL1862,P00519,Multicellular organism,FALSE
CHEMBL203,P00533,Multicellular organism,FALSE
```

## Формат выходного CSV

Финальная таблица включает все поля из `TARGETS_COLUMN_ORDER`: идентификаторы,
описания белков, перекрёстные ссылки и метаданные пайплайна.【F:library/schemas/targets.py†L15-L98】
Отсутствующие значения заполняются дефисом (`-`), а в идентификаторах вместо
него остаётся пустая строка.【F:library/pipelines/target/postprocessing.py†L310-L335】 Пример (подмножество столбцов для тех же
строк):

```
target_chembl_id,uniprot_id_primary,organism,target_type
CHEMBL1827,PDE5A,Homo,Multicellular organism
CHEMBL1859,CACNA1H,Homo,Multicellular organism
CHEMBL202,P00374,Homo,Multicellular organism
CHEMBL1809,P0ABQ4,Escherichia,Unicellular organism
CHEMBL1862,P00519,Homo,Multicellular organism
CHEMBL203,P00533,Homo,Multicellular organism
```

Колонка `target_type` вычисляется по линейджу; возможные значения — `Multicellular organism`,
`Unicellular organism`, `Viral organism` и `Unknown` в соответствии с правилами
классификатора.【F:library/pipelines/target/postprocessing.py†L653-L668】【F:library/pipelines/target/organism_classification.py†L180-L270】

## Шаги постобработки

1. **Гармонизация идентификаторов.** UniProt приводятся к единому виду, при
   необходимости используются резервные значения, а генные символы и синонимы
   собираются из нескольких источников. EC-номера агрегируются в строку с
   разделителем `|`.【F:library/pipelines/target/postprocessing.py†L369-L452】
2. **Заполнение опциональных колонок.** Часто отсутствующие поля (изоформы,
   реакции и т.п.) заполняются дефисом, чтобы удовлетворить схему.【F:library/pipelines/target/postprocessing.py†L455-L484】
3. **Обогащение словарями.** Файл `targets_type.csv` добавляет иерархию IUPHAR,
   индексы генов/таксонов и флаг `multifunctional_enzyme`. Флаг приводится к
   булевому типу, а линейдж используется для первоначального расчёта
   клеточности, задействуемого дальше в конвейере.【F:library/integration/input_initialisation_library.py†L520-L593】
4. **Финальное выравнивание.** `finalise_targets` отбрасывает строки без
   устойчивого UniProt (если нет резервного значения), убирает дубликаты по
   `target_chembl_id`, приводит текстовые столбцы к строчному регистру и на базе
   линейджа вычисляет окончательную метку `target_type`.【F:library/pipelines/target/postprocessing.py†L596-L668】
5. **Проекция на схему.** `align_target_columns` упорядочивает столбцы,
   заменяет пустые значения дефисами и нормализует флаги посттрансляционных
   модификаций, обеспечивая сопоставимость с историческими выгрузками.【F:library/pipelines/target/postprocessing.py†L186-L335】

## Режим offline

Команда `python -m library.utils.cli_tools.pipeline_targets_main` воспроизводит
пайплайн без сетевых запросов: подставляется детерминированный источник ChEMBL,
данные валидируются схемой `TargetsSchema`, а запись финального CSV выполняется
тем же потоком, что и в боевом CLI.【F:library/utils/cli_tools/pipeline_targets_main.py†L467-L565】 Поскольку кэш содержит только каркас
идентификаторов, многие поля будут заполнены `-`; режим предназначен для
регрессионных проверок и контроля схемы.

## Обработка HTTP-ошибок

При реальных вызовах каждая партия данных оборачивается в перехватчик
`requests.RequestException`. Ошибки логируются с указанием размера чанка и
таймаута, трансформируются в `PipelineError` и приводят к остановке пайплайна и
записи `_failure_cases.csv`. Ошибки записи raw-снимков обрабатываются аналогично
с подробными сообщениями.【F:scripts/get_target_data.py†L1620-L1818】

## Инварианты результата

- `target_chembl_id` уникален: дубликаты удаляются перед экспортом.【F:library/pipelines/target/postprocessing.py†L624-L632】
- Плейсхолдеры детерминированы: идентификаторы очищаются до пустых строк, остальные
  пропуски заменяются на `-`, а набор текстовых колонок приводится к нижнему
  регистру.【F:library/pipelines/target/postprocessing.py†L186-L335】
- CSV формируется в исходном порядке идентификаторов и пишется детерминированным
  писателем, поэтому повторный запуск на тех же входных данных даёт битово
  идентичный файл.【F:library/pipelines/target/postprocessing.py†L523-L532】【F:library/pipelines/target/postprocessing.py†L723-L730】
- Валидатор `TargetsSchema` гарантирует присутствие обязательных колонок и
  корректные типы даже при потоковой обработке в offline-режиме.【F:library/utils/cli_tools/pipeline_targets_main.py†L369-L443】

## Дополнительные материалы

- Общая спецификация выходных данных: [`docs/ru/OUTPUT.md`](./ru/OUTPUT.md)
- Английская версия: [`docs/OUTPUT_TARGETS_EN.md`](./OUTPUT_TARGETS_EN.md)
- Руководство по CLI: [`docs/ru/guides/USAGE.md`](./ru/guides/USAGE.md)

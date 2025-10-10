# Выходные данные пайплайна целей: постобработка организмов, изоформ и имён

Это дополнение к [описанию выходных данных](./OUTPUT.md). В нём
рассматриваются вспомогательные артефакты, формируемые после завершения
`scripts/get_target_data.py`. Англоязычная версия расположена в
[`OUTPUT_TARGETS.md`](../en/OUTPUT_TARGETS.md).

## Сводка артефактов

Пайплайн целей с настройками по умолчанию публикует четыре детерминированных CSV в
директории вывода:

- `targets_<YYYYMMDD>.csv` — основной экспорт, описанный в
  [`OUTPUT.md`](./OUTPUT.md#экспорт-targets).
- `organism.output.target_<YYYYMMDD>.csv` — справочная таблица для пайплайна
  активностей и проверок QA.
- `isoform.output.target_<YYYYMMDD>.csv` — развёртка изоформ и токенов для
  контроля синонимики.
- `names.output.target_<YYYYMMDD>.csv` — компактный каталог компонент и имён
  для согласования с историческими выгрузками ChEMBL.
- `IUPHAR.output.target_<YYYYMMDD>.csv` — нормализованный справочник IUPHAR и
  синонимы, восстановленные из итогового экспорта.

Суффикс `<YYYYMMDD>` задаётся опцией `--date` либо берётся из
детерминированного значения `local.io.default_date_prefix` (его можно
переопределить переменной `CHEMBL_DA_DEFAULT_DATE_PREFIX`).

## Постобработка «организмов»

### Этапы

1. **Загрузка финальной таблицы** — `finalise_file` читает CSV из предыдущего
   шага, учитывая разделитель и кодировку из конфигурации.
2. **Нормализация таксономии** — при необходимости переименовываются исторические
   столбцы (`superkingdom` → `lineage_superkingdom`), все текстовые поля приводятся
   к строковому типу.
3. **Фильтрация идентификаторов** — строки без UniProt удаляются, дубликаты
   `target_chembl_id` отбрасываются, существующий столбец `type` переименовывается
   чтобы не конфликтовать с `target_type`.
4. **Определение клеточности** — функция
   `organism_classification.add_cellularity_smart` по роду и родословной
   присваивает значения `Multicellular organism`, `Unicellular organism` или
   `Viral`.
5. **Формирование справочника** — набор столбцов сводится к схеме
  (`target_chembl_id`, `target_type`, `unicellular_organism`,
  `multifunctional_enzyme`, `IUPHAR_class`, `IUPHAR_subclass`,
  `sortorder.target`, `gene_index`, `taxon_index`), булевы поля приводятся к
   типу `boolean`.
6. **Выгрузка `organism.output.*`** — строки сортируются по `target_chembl_id` и
   записываются с переводами строк `\n` через `write_csv_deterministic` рядом с
   основным файлом (или в путь из `--final-out`).

Полученный справочник заменяет статический
`config/dictionary/_target/targets_type.csv` и используется
`library.integration.input_initialisation_library` при обогащении активностей
метаданными об организмах.

### Формат входа/выхода

Минимальный набор входных столбцов:

```
target_chembl_id,genus,lineage_superkingdom,lineage_phylum,lineage_class,IUPHAR_class,IUPHAR_subclass,multifunctional_enzyme
CHEMBL1824,Homo,Metazoa,Chordata,Mammalia,Voltage-gated ion channels,Sodium channels,false
CHEMBL240,Homo,Metazoa,Chordata,Mammalia,Transforming growth factor beta receptors,Type I receptor,false
CHEMBL6130,Candida,Fungi,Ascomycota,Saccharomycetes,,,
CHEMBL259,Homo,Metazoa,Chordata,Mammalia,Fibroblast growth factor receptors,Type 4 receptor,true
CHEMBL1111,Homo,Metazoa,Chordata,Mammalia,Hydrolases,Serine peptidases,false
CHEMBL1922,Influenza A virus,Viruses,Negarnaviricota,Insthoviricetes,,,false
CHEMBL6138,Escherichia,Bacteria,Proteobacteria,Gammaproteobacteria,Membrane proteins,Pore-forming,false
```

Соответствующий результат (первые строки):

```
target_chembl_id,target_type,unicellular_organism,multifunctional_enzyme,IUPHAR_class,IUPHAR_subclass,sortorder.target,gene_index,taxon_index
CHEMBL1824,Multicellular organism,false,false,Voltage-gated ion channels,Sodium channels,0000012345,GENE000123,TAX000987
CHEMBL240,Multicellular organism,false,false,Transforming growth factor beta receptors,Type I receptor,0000012388,GENE000981,TAX000654
CHEMBL259,Multicellular organism,false,true,Fibroblast growth factor receptors,Type 4 receptor,0000012401,GENE000456,TAX000456
CHEMBL6130,Unicellular organism,true,false,-,-,0000012750,-,-
CHEMBL6138,Unicellular organism,true,false,Membrane proteins,Pore-forming,0000012761,-,-
CHEMBL1111,Multicellular organism,false,false,Hydrolases,Serine peptidases,0000012870,GENE000321,TAX000321
CHEMBL1922,Viral,true,false,-,-,0000012999,-,-
```

Отсутствующие значения заменяются `-`, как и в основном экспорте.

### Режим offline

Постобработка полностью локальная. При запуске оркестратора с `--offline`
(либо при пропуске загрузок) ожидается, что итоговый CSV и словари уже доступны
на диске. Если построить lookup не удаётся, CLI завершится с `FileNotFoundError`
и перечислит ожидаемые пути.

### HTTP-ошибки

Запросы к Chembl выполняются до начала постобработки. Любая
`requests.RequestException` логируется как `chembl_fetch_failed` с указанием
размера чанка и тайм-аута, после чего пробрасывается как `PipelineError`. Это
гарантирует отказ без частичных артефактов.

### Инварианты

- `target_type` принимает только значения `Multicellular organism`,
  `Unicellular organism` или `Viral`.
- `unicellular_organism` равен `true` при `target_type` из {`Unicellular organism`,
  `Viral`} и `false` в остальных случаях.
- `multifunctional_enzyme` — булево поле; пустые значения сериализуются как
  `false`.
- Строки отсортированы по `target_chembl_id`. Повторный запуск на тех же входных
  данных даёт байтово идентичный CSV.
- Кодировка UTF-8 и переводы строк `\n` используются всегда.

## Развёртка изоформ

`isoform.output.target_<stamp>.csv` повторяет Power Query-скрипт, который ранее
использовался для валидации изоформ. Алгоритм включает восемь этапов:

1. **Проекция колонок** — читаются только `isoform_synonyms`, `isoform_names`,
   `isoform_ids`, `uniprot_id_primary`, `target_chembl_id`.
2. **Приведение регистра** — `isoform_synonyms` и `isoform_names` переводятся в
   нижний регистр, идентификаторы изоформ остаются без изменений.
3. **Разделение списков** — значения разбиваются по символу `|`, пробелы
   отбрасываются, пустые токены удаляются.
4. **Выравнивание по индексам** — вспомогательная функция `MakeTriples`
   дополняет короткие списки `null`, формируя записи `{name, id, synonym}`.
5. **Токенизация** — синонимы разбиваются по `:`, затем каждый токен расширяется
   до набора `[токен, токен без "pde", токен без "pld"]` с сохранением порядка и
   удалением пустых значений.
6. **Построение таблиц** — trimmed `isoform_names` образуют первую таблицу
   (строки `""`, `"n/a"`, `"none"` исключаются), токены синонимов — вторую.
7. **Объединение и очистка** — таблицы объединяются и проходят три последовательных
   `Distinct`: (а) по (`id`, `name`, `target_chembl_id`, `uniprot_id_primary`), (б)
   после стабильной сортировки `mergesort` по (`uniprot_id_primary`, `id`) по
   набору (`id`, `target_chembl_id`, `name`), (в) финальный `Distinct` по (`id`,
   `name`).
8. **Запись** — результат сортируется в колонках `["id", "uniprot_id_primary",
   "target_chembl_id", "name"]` и сохраняется в UTF-8 без BOM.

### Запуск из Python

Воспроизвести справочник можно напрямую из Python:

```bash
python - <<'PY'
from pathlib import Path
from library.postprocessing.target import process_targets

latest = Path("data/output/output.target_20250101.csv")
process_targets(latest, verbose=True)
PY
```

Если путь не указан явно, `process_targets` найдёт последний `output.target_*.csv`
в каталоге `data/output`, повторяя поведение Power Query.

### Детерминизм и совместимость

- Реализация один-в-один повторяет Power Query (M): каждая стадия (проекция,
  нормализация, токенизация, объединение и три `Distinct`) покрыта регрессионными
  тестами.
- Перед вторым `Distinct` применяется стабильная сортировка `mergesort`, поэтому
  при коллизиях идентификаторов всегда выживает один и тот же ряд.
- CSV сохраняется в UTF-8 с переводами строк `\n` и фиксированным порядком
  колонок `id`, `uniprot_id_primary`, `target_chembl_id`, `name`.
- Поддерживаются Python 3.11+ и pandas ≥ 2.0, что совпадает с основным
  стеком пайплайна.

Пример (источник `output.target_20250101.csv`):

```
target_chembl_id,uniprot_id_primary,isoform_ids,isoform_names,isoform_synonyms
CHEMBL1,Q11111,ENSP0001|ENSP0002,Alpha|Beta,Alpha Alt|PDE3A:Alpha
CHEMBL2,Q22222,,Gamma|N/A|none,Gamma:Variant|n/a|none
CHEMBL3,Q33333,ID_UP|id_low,Theta|Lambda,PLDA:Variant
```

Результат `isoform.output.target_20250101.csv` (первые строки):

```
id,uniprot_id_primary,target_chembl_id,name
ENSP0001,Q11111,CHEMBL1,alpha
ENSP0001,Q11111,CHEMBL1,alpha alt
ENSP0002,Q11111,CHEMBL1,beta
ENSP0002,Q11111,CHEMBL1,pde3a
ENSP0002,Q11111,CHEMBL1,3a
,Q22222,CHEMBL2,gamma
,Q22222,CHEMBL2,variant
ID_UP,Q33333,CHEMBL3,theta
ID_UP,Q33333,CHEMBL3,plda
ID_UP,Q33333,CHEMBL3,a
```

Трансформация идемпотентна, не зависит от локали и не требует сетевого доступа.

## Постобработка `names.*`

`names.output.target_<stamp>.csv` формирует упрощённое представление компонент и
наименований, которые в ChEMBL разнесены по вложенным таблицам. Артефакт
поддерживает стабильность словаря имён для склейки с выгрузками активностей и
ассев, одновременно устраняя неиспользуемые поля.

### Вход и выход

- **Вход:** финальный `targets_<stamp>.csv`, собранный скриптом
  `scripts/get_target_data.py` (в каталоге `--output-dir` или `--final-out`).
  При наличии словарных оверрайдов, указанных в конфигурации, они применяются
  для фиксации порядка строк.
- **Выход:** `names.output.target_<stamp>.csv`, записанный рядом с основным
  экспортом (или в путь `--final-out`). Файл кодируется в UTF-8 с `\n` и
  сортируется по `target_chembl_id`, `active_component_type`,
  `active_component` (по убыванию), `component_id` и нормализованному
  `component_name`.

### Изменения схемы

В итоговом файле остаются только столбцы, необходимые потребителям:

```
target_chembl_id,component_id,component_name,component_synonyms,
component_accession,active_component,active_component_type
```

- Поля для отладки — `component_description`, `component_synonym_ids`,
  `component_type_raw`, `component_sequence`, промежуточные `_source_*` —
  удаляются.
- Пустые строки приводятся к `-`, булев `active_component` переводится в
  nullable-тип `boolean` из pandas.
- Отсутствующие `component_synonyms` заменяются на пустую строку с разделителем
  `|`, чтобы сохранить совместимость с историческими Power Query-процессами.

### Производные поля

`active_component_type` объединяет тип компоненты и признак активности, чтобы
потребителям не приходилось повторно открывать основной экспорт:

```
if active_component is True:
    active_component_type = (component_type or "unknown").lower()
elif active_component is False:
    active_component_type = "inactive"
else:
    active_component_type = "unassigned"
```

Значение хранится в нижнем регистре ASCII и никогда не пустует.

### Логирование

В лог добавляется событие `target_names_postprocess` с агрегированными
метриками:

- `input_rows` — количество строк во входном экспорте.
- `dropped_components` — сколько компонент отфильтровано из-за отсутствия ID или
  имени.
- `null_synonyms` — строки, где список синонимов пришлось заполнить значением по
  умолчанию.
- `written_rows` — финальное количество строк после дедупликации.

Показатели отображаются в структурированных логах и YAML с метаданными
пайплайна.

### Детерминизм

- Сортировка использует стабильный `mergesort` и фиксированные ключи, указанные
  выше.
- Исходный файл не изменяется: helper копирует данные, преобразует и записывает
  новый CSV, поэтому повторный запуск на тех же входах выдаёт байтово
  идентичный `names.output.target_<stamp>.csv`.
- Регрессионные тесты проверяют схему и значения по умолчанию, чтобы исключить
  случайное дрейфование колонок.

## Постобработка `IUPHAR.*`

`IUPHAR.output.target_<stamp>.csv` воспроизводит исторический helper, который
консолидировал цепи, семейства и синонимы Guide to PHARMACOLOGY. Файл
пересобирается при каждом запуске пайплайна таргетов, чтобы синхронизироваться с
каноническим экспортом.

### Вход и выход

- **Вход:** канонический `output.target_<stamp>.csv` с колонками обогащения
  IUPHAR, сформированный `scripts/get_target_data.py`.
- **Выход:** `IUPHAR.output.target_<stamp>.csv` в той же директории (или в путь,
  переданный через Python API). При отсутствии входного пути helper находит
  последний `output.target_*.csv` в стандартном каталоге `data/output`.
- **Схема:**

  ```
  target_chembl_id,guidetopharmacology_id,iuphar_target_id,iuphar_family_id,
  iuphar_type,iuphar_class,iuphar_subclass,iuphar_chain,iuphar_name,iuphar_synonyms
  ```

  Порядок колонок фиксируется при записи CSV, чтобы совпадать с ожиданиями
  потребителей.

### Удаляемые столбцы и переименования

Оставляются только поля, необходимые для справочника синонимов. Helper удаляет
диагностические столбцы из этапа обработки компонент:

- `component_synonym_ids`
- `component_type_raw`
- `component_sequence`
- `component_structures`

`GuidetoPHARMACOLOGY` переименовывается в `guidetopharmacology_id`, остальные
колонки сохраняются без изменений.

### Формирование синонимов

Синонимы собираются из нескольких источников и нормализуются детерминированно:

1. `gtop_synonyms`, `synonyms`, `component_description` и каноническое
   `iuphar_name` проходят через `normalise_text`, разбиваются по `|` (или по
   ключам JSON для структурированных описаний) и очищаются от пустых токенов.
2. Содержимое скобок удаляется, последовательные пробелы схлопываются, токены
   переводятся в нижний регистр.
3. Дубликаты убираются в порядке первого появления для сохранения стабильного
   порядка.
4. Результат объединяется обратно в строку `iuphar_synonyms` с разделителем `|`.

Helper фиксирует количество токенов до и после удаления дубликатов, чтобы
упростить QA-проверки.

### Логирование и детерминизм

При запуске с `verbose=True` helper выводит структурированные сообщения с путями
входа/выхода, количеством обработанных строк, числом удалённых столбцов и
статистикой по синонимам. Запись выполняется через `write_csv_deterministic` со
стабильным ключом сортировки (`target_chembl_id`), что гарантирует байтово
идентичные результаты на одинаковых входных данных и перевод строки `\n`.

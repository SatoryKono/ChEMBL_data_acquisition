# Выходные данные пайплайна целей: постобработка организмов и изоформ

Это дополнение к [описанию выходных данных](./ru/OUTPUT.md). В нём
рассматриваются вспомогательные артефакты, формируемые после завершения
`scripts/get_target_data.py`. Англоязычная версия расположена в
[`OUTPUT_TARGETS_EN.md`](./OUTPUT_TARGETS_EN.md).

## Сводка артефактов

Пайплайн целей с настройками по умолчанию публикует три детерминированных CSV в
директории вывода:

- `targets_<YYYYMMDD>.csv` — основной экспорт, описанный в
  [`OUTPUT.md`](./ru/OUTPUT.md#экспорт-targets).
- `organism.output.target_<YYYYMMDD>.csv` — справочная таблица для пайплайна
  активностей и проверок QA.
- `isoform.output.targets_<YYYYMMDD>.csv` — развёртка изоформ и токенов для
  контроля синонимики.

Суффикс `<YYYYMMDD>` задаётся явно опцией `--date` либо вычисляется по текущей
дате в UTC.

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
   `target_sort_order`, `gene_index`, `taxon_index`), булевы поля приводятся к
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
target_chembl_id,target_type,unicellular_organism,multifunctional_enzyme,IUPHAR_class,IUPHAR_subclass,target_sort_order,gene_index,taxon_index
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

Механика `isoform.output.*` сохранилась. Таблица целей читается повторно,
столбцы `isoform_*` очищаются, синонимы раскладываются в токены, дубликаты по
(`target_chembl_id`, `isoform_id`, `term`, `token`) отбрасываются, после чего
записывается детерминированный CSV для QA.

Для полноты ниже приведён исходный пример и результат:

```
target_chembl_id,isoform_ids,isoform_names,isoform_synonyms
CHEMBL1824,ENSP00000350283,Nav1.7,"Nav1.7:SCN9A isoform 3"
CHEMBL1824,ENSP00000350284,"Nav1.7 isoform 2","Nav1.7 splice variant:SCN9A-2"
CHEMBL6130,,"",""
CHEMBL240,ENSP00000456012,ALK2,"Activin receptor-like kinase 2:ACVR1"
CHEMBL259,ENSP00000263253,FGFR4,"FGFR-4 isoform:Fibroblast growth factor receptor-4"
```

```
target_chembl_id,isoform_id,isoform_name,term,token
CHEMBL1824,ENSP00000350283,Nav1.7,Nav1.7,7
CHEMBL1824,ENSP00000350283,Nav1.7,Nav1.7,nav1
CHEMBL1824,ENSP00000350283,Nav1.7,nav1.7,7
CHEMBL1824,ENSP00000350283,Nav1.7,nav1.7,nav1
CHEMBL1824,ENSP00000350283,Nav1.7,scn9a isoform 3,3
CHEMBL1824,ENSP00000350283,Nav1.7,scn9a isoform 3,isoform
CHEMBL1824,ENSP00000350283,Nav1.7,scn9a isoform 3,scn9a
CHEMBL1824,ENSP00000350284,Nav1.7 isoform 2,Nav1.7 isoform 2,2
CHEMBL1824,ENSP00000350284,Nav1.7 isoform 2,Nav1.7 isoform 2,7
CHEMBL1824,ENSP00000350284,Nav1.7 isoform 2,Nav1.7 isoform 2,isoform
```

Процесс детерминированный, не зависит от локали и не требует сетевого доступа.

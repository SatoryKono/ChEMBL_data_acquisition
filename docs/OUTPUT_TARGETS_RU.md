# Постобработка таргетов — `organism.output.target_*.csv`

Пайплайн таргетов формирует дополнительный CSV с фокусом на организмы на базе
экспорта ChEMBL/UniProt (`output.target_*.csv`). Логика Power Query (M),
использовавшаяся в отчётах по таксономии, перенесена на Python без изменений
результата по строкам и столбцам.

## Этапы преобразования

1. **Определение входа** — используется только что созданный файл
   `output.target_*.csv` (или самый свежий файл по шаблону при отдельном
   запуске).
2. **Классификация клеточности** — поля superkingdom/phylum приводятся к нижнему
   регистру. Для каждой строки вычисляется метка:
   - `ClassifyByLineage` закрывает детерминированные случаи для вирусов,
     бактерий, архей и подобранных филумов эукариот.
   - При неопределённом или пустом lineage вызывается `ClassifyByFetch`,
     выполняющий запрос к NCBI Taxonomy (`Lineage`, `ScientificName`) и
     применяющий те же правила.
3. **Флаг мультифункциональности** — колонка `reaction_ec_numbers` режется по
   символу `|`, очищается от дублей и усечённых сегментов EC; наличие более одной
   уникальной первой части кода даёт `multifunctional_enzyme = True`.
4. **Объединение и проекция** — рассчитанные признаки присоединяются обратно к
   подмножеству идентификаторов. Итоговый порядок колонок:
   `target_chembl_id`, `uniprot_id_primary`, `organism`, `taxon_id`,
   `lineage_superkingdom`, `lineage_phylum`, `lineage_class`, `cellularity`,
   `multifunctional_enzyme`.
5. **Вывод** — файл сохраняется рядом с исходным под именем с префиксом
   `organism.` в кодировке UTF-8 с фиксированным порядком строк.

## Пример

Вход (`output.target_20250301.csv`):

```csv
target_chembl_id,uniprot_id_primary,organism,taxon_id,lineage_superkingdom,lineage_phylum,lineage_class,reaction_ec_numbers
CHEMBL1,P12345,Species A,10239,Viruses,,,"1.1.1.1|2.7.11.1"
CHEMBL2,Q54321,Species B,9606,Eukaryota,Chordata,Chordata,"1.1.1.1|3.2.1.4"
CHEMBL3,R99999,Species C,111,Archaea,,,""
CHEMBL4,S88888,Species D,222,,,"","1.1.1.1|2.7.11.1|1.2.3.4"
CHEMBL5,T77777,Species E,333,Bacteria,,,"2.7.11.1"
```

Выход (`organism.output.target_20250301.csv`):

```csv
target_chembl_id,uniprot_id_primary,organism,taxon_id,lineage_superkingdom,lineage_phylum,lineage_class,cellularity,multifunctional_enzyme
CHEMBL1,P12345,Species A,10239,viruses,,,"acellular (virus)",True
CHEMBL2,Q54321,Species B,9606,eukaryota,chordata,chordata,multicellular,True
CHEMBL3,R99999,Species C,111,archaea,,,"unicellular",False
CHEMBL4,S88888,Species D,222,,,ambiguous,True
CHEMBL5,T77777,Species E,333,bacteria,,,"unicellular",False
```

Постобработанный CSV сохраняет порядок строк входного набора. Списки, полученные
из `reaction_ec_numbers`, используются только во внутренних расчётах и не
попадают в итоговый файл.

## Эксплуатационные замечания

- **Офлайн-режим** — параметр `offline=True` у
  `library.postprocessing.target.process_targets` отключает HTTP-запросы; такие
  строки получают `cellularity = "ambiguous"`.
- **Ошибки HTTP** — сетевые сбои или невалидный XML повторяют поведение M-кода:
  классификация становится неопределённой, исключения не выбрасываются.
- **Логирование** — при `verbose=True` выводятся пути входа/выхода, число строк,
  количество HTTP-запросов к NCBI и количество строк с результатом `ambiguous`.

Аналогичный материал на английском — в `docs/OUTPUT_TARGETS_EN.md`.

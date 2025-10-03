# `_testitem` reference

> Changelog
> - 2024-05-08 — добавлено описание каталога.

Каталог содержит подготовленные справочники для пайплайна `scripts/get_testitem_data.py`:

- `testitem.csv` — основной словарь тестовых веществ, включая аггрегированные названия, стереохимию и скелет InChI;
- `molecule_catalog.csv` — признаки натурального происхождения, проконвертированности и полимерности молекул;
- `../../config/dictionary/_testitem/molecule_hierarchy.csv` — связки родительских молекул, используемые при нормализации идентификаторов.

Описание колонок доступно в разделе `testitem.csv (processed export)` файла [docs/en/devel/DATA_SCHEMA.md](../../docs/en/devel/DATA_SCHEMA.md#testitemcsv-processed-export).

## Обновление
1. Подготовьте CSV с `molecule_chembl_id` (см. [docs/en/devel/DATA_SCHEMA.md](../../docs/en/devel/DATA_SCHEMA.md#input-tables)).
2. Запустите `get-testitem-data` согласно [docs/en/user/USAGE.md](../../docs/en/user/USAGE.md#test-item-pipeline-get-testitem-data), указав входной список и путь для результирующего CSV.
3. Если в набор добавлены новые молекулы, обновите словари `molecule_catalog.csv` и `molecule_hierarchy.csv`: загрузите актуальные выгрузки из внутренних источников или повторно выполните конвейер с ключом `--force`, чтобы перегенерировать кеши, заданные в `sources.chembl.molecule_catalog`.
4. Проверьте таблицы при помощи `table_quality_main` (`chembl_testitem`) и unit-тестов.
5. Зафиксируйте дату обновления и ссылки на выгрузки в разделе «История» этого README.

## Проверка качества
- Сравнивайте количество уникальных `salt_chembl_id` и флагов `natural_product` между версиями.
- При обнаружении пустых `standard_inchi_skeleton` перепроверьте upstream-выгрузку и обновите список исключений.

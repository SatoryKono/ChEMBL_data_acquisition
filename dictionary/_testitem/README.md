# `_testitem` reference

> Changelog
> - 2024-05-08 — добавлено описание каталога.

Каталог содержит подготовленные справочники для пайплайна `scripts/get_testitem_data.py`:

- `testitem.csv` — основной словарь тестовых веществ, включая аггрегированные названия, стереохимию и скелет InChI;
- `molecule_catalog.csv` — признаки натурального происхождения, проконвертированности и полимерности молекул;
- `molecule_hierarchy.csv` — связки родительских молекул, используемые при нормализации идентификаторов.

Описание колонок доступно в разделе `testitem.csv (processed export)` файла [docs/en/reference/DATA_SCHEMA.md](../../docs/en/reference/DATA_SCHEMA.md#testitemcsv-processed-export).

## Обновление
1. Подготовьте CSV с `molecule_chembl_id` (см. [DATA_SCHEMA.md](../../docs/en/reference/DATA_SCHEMA.md#input-tables)).
2. Запустите `python -m scripts.get_testitem_data chembl --input <path>/testitem.csv --output <path>/testitem.csv`.
3. Обновите вспомогательные таблицы через `python -m scripts.get_testitem_data catalog --output <path>/molecule_catalog.csv` и `python -m scripts.get_testitem_data hierarchy --output <path>/molecule_hierarchy.csv`, если изменился состав молекул.
4. Проверьте таблицы при помощи `table_quality_main` (`chembl_testitem`) и unit-тестов.
5. Зафиксируйте дату обновления и ссылки на выгрузки в разделе «История» этого README.

## Проверка качества
- Сравнивайте количество уникальных `salt_chembl_id` и флагов `natural_product` между версиями.
- При обнаружении пустых `standard_inchi_skeleton` перепроверьте upstream-выгрузку и обновите список исключений.

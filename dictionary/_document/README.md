# `_document` reference

> Changelog
> - 2024-05-08 — добавлено описание каталога.

Каталог хранит консолидированный экспорт `document.csv`, который собирает метаданные из ChEMBL, PubMed, Crossref, OpenAlex и Semantic Scholar. Файл выступает эталоном для тестирования пайплайна `scripts/get_document_data.py` и служит источником для QA-проверок. Полная схема описана в [docs/reference/en/DATA_SCHEMA.md](../../docs/reference/en/DATA_SCHEMA.md#documentscsv-processed-export).

## Файлы
| Имя | Назначение |
| --- | --- |
| `document.csv` | Нормализованная публикационная таблица с объединёнными идентификаторами, аннотациями Mesh и вычисленной классификацией публикаций. |

## Источник и обновление
1. Подготовьте CSV с `document_chembl_id` (см. раздел «Input Tables» в [DATA_SCHEMA.md](../../docs/reference/en/DATA_SCHEMA.md#input-tables)).
2. Выполните `python -m scripts.get_document_data chembl --input <path>/documents.csv --output <path>/document.csv` согласно инструкции [docs/guides/en/README.md](../../docs/guides/en/README.md#scriptsget_document_datapy).
3. Сверьте структуру с документацией и перезапустите `table_quality_main` для таблицы `chembl_document`.
4. Зафиксируйте дату обновления и внесённые изменения в этом README.

## Проверка качества
- Сверяйте счётчики публикационных типов (`publication_types_normalised`) между версиями.
- Отслеживайте ошибки аггрегаторов (`scholar.Error`, `crossref.Error`, `OpenAlex.Error`) — при появлении новых значений обновите процесс и документируйте исключения.

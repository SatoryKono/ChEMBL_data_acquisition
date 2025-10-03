# `_document` reference

> Changelog
> - 2024-05-08 — добавлено описание каталога.

Каталог хранит консолидированный экспорт `document.csv`, который собирает метаданные из ChEMBL, PubMed, Crossref, OpenAlex и Semantic Scholar. Файл выступает эталоном для тестирования пайплайна `scripts/get_document_data.py` и служит источником для QA-проверок. Полная схема описана в [docs/DATA_SCHEMA_EN.md](../../docs/DATA_SCHEMA_EN.md#documentscsv-processed-export).

## Файлы
| Имя | Назначение |
| --- | --- |
| `document.csv` | Нормализованная публикационная таблица с объединёнными идентификаторами, аннотациями Mesh и вычисленной классификацией публикаций. |

## Источник и обновление
1. Подготовьте CSV с `document_chembl_id` (см. раздел «Input Tables» в [docs/DATA_SCHEMA_EN.md](../../docs/DATA_SCHEMA_EN.md#input-tables)).
2. Выполните `get-document-data all` согласно пошаговому описанию в [docs/USAGE_EN.md](../../docs/USAGE_EN.md#document-pipeline-get-document-data), передав путь к входному списку и желаемому выходному файлу.
3. Сверьте структуру с документацией и перезапустите `table_quality_main` для таблицы `chembl_document`.
4. Зафиксируйте дату обновления и внесённые изменения в этом README.

## Проверка качества
- Сверяйте счётчики публикационных типов (`publication_types_normalised`) между версиями.
- Отслеживайте ошибки аггрегаторов (`scholar.Error`, `crossref.Error`, `OpenAlex.Error`) — при появлении новых значений обновите процесс и документируйте исключения.

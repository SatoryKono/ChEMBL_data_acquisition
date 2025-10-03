# `_target` reference

> Changelog
> - 2024-05-08 — добавлено описание каталога.

Каталог агрегирует все статические ресурсы для таргет-пайплайна:

- классификатор типов (`targets_type.csv`), который используется `finalise_targets` для определения `type` и связанных флагов;
- выгрузки IUPHAR (`_IUPHAR/*.csv`) для обогащения таргетов семействами и идентификаторами Guide to PHARMACOLOGY;
- кешированные ответы UniProt (`_uniprot/*.json`), задействованные при работе офлайн-режимов `get_target_data` и `pipeline_targets_main`.

Пути к файлам зашиты в конфигурацию (см. [docs/reference/en/CONFIG.md](../../docs/reference/en/CONFIG.md#target-data)).

## Структура
| Путь | Назначение |
| --- | --- |
| `targets_type.csv` | Таблица классификации организмов и ферментных групп, которую читает `finalise_targets` при расчёте `type`. |
| `_IUPHAR/_IUPHAR_target.csv` | Справочник таргетов IUPHAR: соответствия синонимов, SwissProt и семейств. |
| `_IUPHAR/_IUPHAR_family.csv` | Иерархия семейств IUPHAR, используемая для построения `full_name_path`. |
| `_uniprot/*.json` | Кеш ответов UniProt REST API; каждая запись соответствует UniProt accession и содержит аннотации по последовательности, таксономии и синонимам. |

## Обновление
1. Сформируйте входной CSV с `target_chembl_id` (см. [DATA_SCHEMA.md](../../docs/reference/en/DATA_SCHEMA.md#input-tables)).
2. Для онлайнового обновления запустите `python -m scripts.get_target_data all --input <path>/targets.csv --output <path>/targets.csv --raw-out <path>/raw.csv`.
3. Для проверки офлайн-режимов используйте `python -m library.utils.cli_tools.pipeline_targets_main --input <path>/targets.csv --output-dir <path>/cache`.
4. При изменении классификатора типов обновите `targets_type.csv` и протестируйте `finalise_targets` (см. [docs/processes/en/ETL_PROCESS.md](../../docs/processes/en/ETL_PROCESS.md#finalisation)).
5. Обязательно задокументируйте источник данных (ссылка на выгрузку UniProt/IUPHAR) и дату в этом README.

## Проверка качества
- После обновления запускайте `table_quality_main` с `--table-name chembl_targets`.
- Сверяйте заполненность ключевых колонок (`uniprot_id`, `IUPHAR_type`, `type`) и дифференцируйте изменения по числу уникальных организмов.
- Если менялась `_uniprot` директория, убедитесь, что кэш содержит только JSON с валидной схемой UniProt и их количество соответствует числу используемых акцессий.

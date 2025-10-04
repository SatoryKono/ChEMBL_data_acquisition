# Справочные наборы данных

Конвейеры используют статические ресурсы из `config/dictionary`. Таблица ниже
описывает назначение файлов и затрагиваемые пайплайны.

| Путь | Используется в | Назначение |
|------|----------------|-----------|
| `_target/_IUPHAR/_IUPHAR_target.csv` | Target (`iuphar`, `all`) | Сопоставление UniProt↔Guide to PHARMACOLOGY. |
| `_target/_IUPHAR/_IUPHAR_family.csv` | Target (`iuphar`, `all`) | Иерархия семейств IUPHAR для построения `iuphar_full_*`. |
| `_target/_uniprot/*.json` | Target (`uniprot`, `all`) | Кэш ответов UniProt KB по акцессиям (минимальный набор для тестов). |
| `_target/targets_type.csv` | QA таргетов | Карта типов таргетов ChEMBL в агрегированные категории. |
| `_testitem/molecule_catalog.csv` | Test item | Справочник `molecule_chembl_id`→`parent_molecule_chembl_id` и связанных признаков. |
| `_testitem/molecule_hierarchy.csv` | Test item | Фолбэк-иерархия родителей при отсутствии данных в API. |
| `_document/fallback_doi_template.csv` | Document | Шаблон структуры CSV для DOI-фолбэков. |

## Правила сопровождения

- Коммитьте изменения атомарно; smoke-тесты проверяют хэши файлов.
- Добавляя новые поля, фиксируйте их назначение здесь и проверяйте обработку
  `null` в коде.
- JSON-ответы UniProt должны содержать только необходимый минимум; полноценные
  данные будут подтягиваться в продакшене.

## Обогащение родительских молекул

Справочники `_testitem` — единственный источник истины для `parent_molecule_chembl_id`.
Обязательные колонки:

| Файл | Колонки |
|------|---------|
| `_testitem/molecule_catalog.csv` | `molecule_chembl_id`, `parent_molecule_chembl_id`, дополнительные поля PubChem |
| `_testitem/molecule_hierarchy.csv` | `molecule_chembl_id`, `parent_molecule_chembl_id`, `pref_name`, `level` |

Отсутствующие значения логируются; при включённом `parent_fallback` они
подставляются из `molecule_hierarchy.csv`.

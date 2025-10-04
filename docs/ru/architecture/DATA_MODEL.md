# Модель данных

Экспортируемые CSV формируют звёздную схему с таблицей фактов активностей.

```mermaid
erDiagram
    ACTIVITY ||--o{ ASSAY : assay_chembl_id
    ACTIVITY ||--o{ DOCUMENT : document_chembl_id
    ACTIVITY ||--o{ TARGET : target_chembl_id
    ACTIVITY ||--o{ TESTITEM : molecule_chembl_id

    DOCUMENT {
        string document_chembl_id PK
        string doi
        string publication_class
        string pipeline_version
        string timestamp_utc
    }
    TARGET {
        string target_chembl_id PK
        string uniprot_id_primary
        string pref_name
        string pipeline_version
        string timestamp_utc
    }
    ASSAY {
        string assay_chembl_id PK
        string target_chembl_id FK
        string bao_format
        string pipeline_version
        string timestamp_utc
    }
    TESTITEM {
        string molecule_chembl_id PK
        string parent_molecule_chembl_id
        string pref_name
        string pipeline_version
        string timestamp_utc
    }
    ACTIVITY {
        string activity_id PK
        string assay_chembl_id FK
        string document_chembl_id FK
        string molecule_chembl_id FK
        float standard_value
        string action_type
        string pipeline_version
        string timestamp_utc
    }
```

## Использование измерений

- **Документы** — библиографические сведения для трассировки источников.
- **Таргеты** — белковые и генетические аннотации для дальнейшего обогащения.
- **Ассайи** — экспериментальный контекст, BAO и QA-флаги.
- **Тестовые объекты** — идентичность молекул, способы введения, данные PubChem,
  связи с родительскими молекулами.
- **Факт активностей** — нормализованные значения измерений, класс `action_type`,
  вычисленные границы и хэш дополнительного свойства.

## Ключи и временные метки

Каждая таблица содержит `pipeline_version` и `timestamp_utc`, что позволяет
отслеживать происхождение и строить инкрементальные загрузки. Активности ссылаются
на измерения по натуральным ключам; при необходимости downstream-хранилище может
создавать суррогатные ключи.

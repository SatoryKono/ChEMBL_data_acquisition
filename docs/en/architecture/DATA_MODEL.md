# Analytical data model

The exported CSV files form a star schema centred on the activity fact table.

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

## Dimensional usage

- **Document dimension** — bibliographic metadata used for evidence tracking and
  literature analysis.
- **Target dimension** — protein and gene annotations for downstream enrichment.
- **Assay dimension** — experimental context, BAO metadata and QA flags.
- **Test item dimension** — molecular identity, administration routes, PubChem
  enrichments and parent relationships.
- **Activity fact** — measurement values, action type classification, derived
  bounds and hash of auxiliary properties.

## Surrogate keys and timestamps

Each table records `pipeline_version` and `timestamp_utc`, enabling incremental
loading and provenance checks. Activity rows reference the other tables via
natural keys; surrogate keys can be introduced downstream if necessary.

# Output specification

Each pipeline writes a deterministic CSV together with metadata (`.meta.yaml`),
table-quality profiles (`_quality_report_table.csv` and `.quality.json`) and log
records. This document describes the CSV schemas, column meanings, data types
and validation rules enforced by Pandera (`library/schemas/*.py`).

## Metadata sidecars

For an output named `output.targets_20250228.csv` the sidecars are:

- `output.targets_20250228.meta.yaml` — pipeline name, CLI parameters, effective
  configuration, schema version and SHA-256 hash of the CSV.
- `output.targets_20250228_quality_report_table.csv` — tabular profile covering
  counts, missing values, numeric statistics and pattern coverage for selected
  columns.
- `output.targets_20250228.quality.json` — JSON summary with the same metrics
  plus warnings raised during profiling.

These artefacts support reproducibility and downstream QA pipelines.

## Document export (`documents`)

Schema: [`library/schemas/documents.py`](../../library/schemas/documents.py).

### ChEMBL base columns

| Column | Type | Required | Notes |
|--------|------|----------|-------|
| `document_chembl_id` | string | Yes | Primary ChEMBL identifier. |
| `title` | string | No | Document title. |
| `abstract` | string | No | Abstract text. |
| `doi` | string | No | DOI provided by ChEMBL. |
| `year` | string/int | No | Publication year (coerced to numeric when possible). |
| `journal` | string | No | Journal name. |
| `journal_abbrev` | string | No | Journal abbreviation. |
| `volume` | string/int | No | Volume information. |
| `issue` | string/int | No | Issue information. |
| `first_page` | string/int | No | First page. |
| `last_page` | string/int | No | Last page. |
| `pubmed_id` | string/int | No | PubMed identifier from ChEMBL. |
| `authors` | string | No | Author list provided by ChEMBL. |
| `source` | string | No | Origin flag from ChEMBL dataset. |

### Derived enrichment fields

| Column | Type | Notes |
|--------|------|-------|
| `doi_normalised` | string | DOI normalised to lowercase without whitespace. |
| `publication_types_normalised` | string | Canonicalised publication types (semicolon separated). |
| `publication_type_score_review` | integer | Heuristic score for review classification. |
| `publication_type_score_experimental` | integer | Heuristic score for experimental articles. |
| `publication_type_score_unknown` | integer | Heuristic score for unknown types. |
| `publication_class` | string | Final class label (`review`, `experimental`, `unknown`). |
| `fetch_status` | string | Fetch status per DOI/PMID. |
| `error_source` | string | Subsystem raising an error (PubMed, CrossRef, etc.). |
| `date_code` | string | Year-month identifier emitted by post-processing. |
| `Index` | string/int | Original positional index before merges. |
| `pipeline_version` | string | Library version recorded during export. |
| `timestamp_utc` | string | ISO 8601 timestamp of export. |

### PubMed enrichment (`PubMed.*`)

All columns are optional strings unless noted. Values originate from the PubMed
E-utilities payload and may contain structured lists encoded as strings.

| Column | Description |
|--------|-------------|
| `PubMed.PMID` | PubMed identifier returned by E-utilities. |
| `PubMed.DOI` | DOI returned by PubMed. |
| `PubMed.ArticleTitle` | Article title. |
| `PubMed.Abstract` | Abstract text. |
| `PubMed.JournalTitle` | Journal title. |
| `PubMed.Volume` / `PubMed.Issue` | Volume/issue values. |
| `PubMed.StartPage` / `PubMed.EndPage` | Page range. |
| `PubMed.PublicationType` | Publication type list (semicolon separated). |
| `PubMed.MeSH_Descriptors` / `PubMed.MeSH_Qualifiers` | MeSH descriptors and qualifiers. |
| `PubMed.ChemicalList` | Chemical substance descriptors. |
| `PubMed.DayRevised` / `PubMed.MonthRevised` / `PubMed.YearRevised` | Date of last revision. |
| `PubMed.YearCompleted` / `MonthCompleted` / `DayCompleted` | Indexing completion date. |
| `PubMed.Error` | Error code if the record failed to download. |
| `PubMed.ISSN` | ISSN value. |
| `PubMed.is_review` | Boolean-like flag set during post-processing. |

### Semantic Scholar enrichment (`scholar.*`)

| Column | Description |
|--------|-------------|
| `scholar.PMID` | PubMed identifier resolved by Semantic Scholar. |
| `scholar.Venue` | Venue metadata. |
| `scholar.PublicationTypes` | Publication type list. |
| `scholar.SemanticScholarId` | Internal Semantic Scholar identifier. |
| `scholar.ExternalIds` | JSON-style mapping of external IDs. |
| `scholar.DOI` | DOI reported by Semantic Scholar. |
| `scholar.Error` | Error message for failed lookups. |
| `scholar.is_review` | Derived review flag. |

### OpenAlex enrichment (`OpenAlex.*`)

| Column | Description |
|--------|-------------|
| `OpenAlex.PublicationTypes` | Publication types from OpenAlex. |
| `OpenAlex.TypeCrossref` | Crossref document type. |
| `OpenAlex.Genre` | Genre/format classification. |
| `OpenAlex.Id` | OpenAlex identifier. |
| `OpenAlex.Venue` | Venue metadata. |
| `OpenAlex.MeshDescriptors` / `OpenAlex.MeshQualifiers` | MeSH terms via OpenAlex. |
| `OpenAlex.Error` | Error text for failed lookups. |
| `OpenAlex.is_review` | Review flag derived from OpenAlex metadata. |

### CrossRef enrichment (`crossref.*`)

| Column | Description |
|--------|-------------|
| `crossref.Type` | Primary type from CrossRef. |
| `crossref.Subtype` | Subtype classification. |
| `crossref.Title` | Title as reported by CrossRef. |
| `crossref.Subtitle` | Subtitle component. |
| `crossref.Subject` | Subject categories. |
| `crossref.Error` | Error message when the request fails. |

## Target export (`targets`)

Schema: [`library/schemas/targets.py`](../../library/schemas/targets.py). Column
order follows `TARGETS_COLUMN_ORDER`.

### Core identifiers and naming

| Column | Type | Notes |
|--------|------|-------|
| `target_chembl_id` | string | Primary identifier (nullable for incomplete merges but expected in final exports). |
| `pref_name` | string | Preferred ChEMBL name. |
| `target_type` | string | ChEMBL target type. |
| `gene_symbol` | string | Primary gene symbol. |
| `gene_symbol_list` | string | Additional gene symbols (comma separated). |
| `protein_name_canonical` | string | Canonical protein name. |
| `protein_name_alt` | string | Alternative names. |
| `protein_synonym_list` | string | Synonyms aggregated from sources. |
| `species_group_flag` | string | Species grouping flag from ChEMBL. |
| `organism` | string | Organism name. |
| `tax_id` | string | NCBI taxonomy ID from ChEMBL. |
| `taxon_id` | string/int | Taxonomy ID from UniProt payload. |

### Optional classification summary

The modular post-processing pipeline for targets
(`library.postprocessing.targets.run_target_pipeline`) validates an aggregated view
that powers QA metrics and downstream consumers needing concise taxonomy
labels. The helper is invoked automatically after the canonical CSV is written
and enforces the schema defined in
`library/postprocessing/targets/schema.py`.

| Column | Description |
|--------|-------------|
| `target_class` | Primary classification extracted from the ChEMBL protein classes. |
| `protein_family` | First-level protein family description surfaced during aggregation. |
| `synonyms` | Deterministically ordered synonym list combining preferred names, component descriptions and alternative names. |
| `pipeline_version` | Version recorded by the post-processing runner (may differ from the export when executed independently). |

### UniProt accession data

| Column | Type | Description |
|--------|------|-------------|
| `uniprot_id_primary` | string | Primary UniProt accession. |
| `uniprot_ids_all` | string | All UniProt accessions combined. |
| `isoform_ids` | string | List of isoform accessions. |
| `isoform_names` | string | Isoform names. |
| `isoform_synonyms` | string | Synonyms per isoform. |
| `uniprot_last_update` | string | Timestamp of the UniProt record. |
| `uniprot_version` | string | UniProt sequence version. |
| `secondaryAccessions` | string | Secondary accessions. |
| `recommendedName` | string | Recommended UniProt name. |
| `geneName` | string | Gene name from UniProt entry. |
| `secondaryAccessionNames` | string | Names associated with secondary accessions. |
| `uniProtkbId` | string | UniProt KB identifier. |

### Functional annotations

| Column | Description |
|--------|-------------|
| `molecular_function` | GO molecular function annotations. |
| `cellular_component` | GO cellular component annotations. |
| `subcellular_location` | Free-text location description. |
| `topology`, `transmembrane`, `intramembrane` | Topology and membrane features. |
| `signal_peptide`, `propeptide` | Signal/propeptide presence flags. |
| `features_signal_peptide`, `features_transmembrane`, `features_topology` | Parsed feature flags from UniProt. |
| `ptm_glycosylation`, `ptm_lipidation`, `ptm_disulfide_bond`, `ptm_modified_residue` | Post-translational modification lists. |
| `glycosylation`, `lipidation`, `disulfide_bond`, `modified_residue`, `phosphorylation`, `acetylation`, `ubiquitination` | Expanded PTM data. |

### Cross-references

| Column | Description |
|--------|-------------|
| `xref_chembl`, `xref_uniprot`, `xref_ensembl`, `xref_iuphar` | Source cross references. |
| `xref_pdb`, `xref_alphafold` | Structural database links. |
| `GuidetoPHARMACOLOGY`, `gtop_target_id`, `gtop_synonyms`, `gtop_natural_ligands_n`, `gtop_interactions_n`, `gtop_function_text_short` | Guide to PHARMACOLOGY metadata. |
| `cross_references` | Aggregated cross references from UniProt. |

### Protein family classification

| Column | Description |
|--------|-------------|
| `pfam`, `Pfam` | Pfam domain identifiers (both ChEMBL and UniProt perspectives). |
| `interpro`, `InterPro` | InterPro annotations. |
| `SUPFAM`, `PROSITE`, `PRINTS`, `TCDB` | Domain/family classification datasets. |
| `protein_classifications` | ChEMBL protein class hierarchy. |
| `protein_class_pred_L1` … `protein_class_pred_confidence` | Predicted classes, rule identifiers, evidence and confidence scores. |

### Reaction annotations

| Column | Description |
|--------|-------------|
| `target_components` | Raw component list. |
| `reactions` | Reaction identifiers aggregated from UniProt/ChEMBL. |
| `reaction_ec_numbers` | EC numbers from reaction annotations. |

### IUPHAR mapping

| Column | Description |
|--------|-------------|
| `iuphar_target_id` | Guide to PHARMACOLOGY target ID. |
| `iuphar_family_id` | Family identifier. |
| `iuphar_type` / `iuphar_class` / `iuphar_subclass` | Classification levels. |
| `iuphar_chain` | Chain identifier if available. |
| `iuphar_name` | Human-readable name. |
| `iuphar_full_id_path` / `iuphar_full_name_path` | Hierarchical path concatenated with `>` separators. |

### Audit columns

| Column | Description |
|--------|-------------|
| `pipeline_version` | Package version used for the export. |
| `timestamp_utc` | UTC timestamp of export. |

## Assay export (`assays`)

Schema: [`library/schemas/assays.py`](../../library/schemas/assays.py).

| Column | Type | Notes |
|--------|------|-------|
| `assay_chembl_id` | string | Primary assay identifier. |
| `accession` | string | UniProt accession for the annotated target. |
| `assay_cell_type` | string | Reported cell type. |
| `assay_subcellular_fraction` | string | Sub-cellular fraction tested in the assay. |
| `assay_group` | string | Group classification delivered with the export. |
| `assay_tissue` | string | Tissue context. |
| `assay_strain` | string | Strain information provided by ChEMBL. |
| `bao_format` | string | BAO format code. |
| `description` | string | Assay description. |
| `document_chembl_id` | string | Source document identifier. |
| `isoform` | string/nullable | Isoform identifier (`None` disables dtype coercion). |
| `mutation` | string | Target mutation details. |
| `target_chembl_id` | string | Linked target identifier. |
| `year` | numeric/string | Publication year, accepted as string or numeric. |
| `pipeline_version` | string | Package version written during export. |
| `timestamp_utc` | string | ISO 8601 export timestamp. |

> **Post-processing note:** Columns listed in
> [`_ASSAY_OUTPUT_DROP_COLUMNS`](../../scripts/get_assay_data.py) are removed
> after download to prevent legacy fields (for example `ASSAY_ID`, `Target TYPE`
> or `acts_per_assay_step5`) from resurfacing in published CSVs.

## Activity export (`activities`)

Schema: [`library/schemas/activities.py`](../../library/schemas/activities.py).

| Column | Type | Validation |
|--------|------|-----------|
| `activity_id` | string/int | Identifier from ChEMBL. |
| `molecule_chembl_id` | string | Molecule identifier. |
| `assay_chembl_id` | string | Assay identifier. |
| `activity_comment` | string | Free-text comment. |
| `assay_description` | string | Description of linked assay. |
| `assay_variant_accession` | string | Variant accession (if provided). |
| `assay_variant_mutation` | string | Mutation description. |
| `bao_format` / `bao_label` | string | BAO annotations. |
| `data_validity_comment` / `data_validity_description` | string | Validity notes. |
| `document_chembl_id` | string | Source document. |
| `pchembl_value` | numeric/string | pChEMBL value. |
| `potential_duplicate` | boolean/string | Flag for duplicates. |
| `qudt_units` | string | Units in QUDT vocabulary. |
| `record_id` | string/int | Record identifier. |
| `relation` | string | Raw relational operator. |
| `src_assay_id` / `src_id` | string/int | Source identifiers. |
| `standard_relation` | string | Normalised relation. |
| `standard_units` | string | Normalised units. |
| `type` / `units` / `value` | string / string / numeric | Original measurement triplet. |
| `standard_type` | string | Must be in allowlist derived from configuration (`metrics` keys). |
| `standard_value` | float | **Validated:** non-negative, coerced to float. |
| `lower_value` / `upper_value` | float | Optional bounds. |
| `activity_properties` | string | Raw JSON-like property list. |
| `action_type` | string | **Validated:** must be one of `PAM`, `NAM`, `activation`, `inhibition`, `binding`, `triaged`, `unknown`. |
| `properties_hash` | string | Hash of parsed properties. |
| `pipeline_version` | string | Package version. |
| `timestamp_utc` | string | Export timestamp. |

## Test item export (`testitems`)

Schema: [`library/schemas/testitems.py`](../../library/schemas/testitems.py).

| Column | Type | Notes |
|--------|------|-------|
| `molecule_chembl_id` | string | Primary molecule identifier. |
| `parent_molecule_chembl_id` | string | Parent molecule identifier. |
| `salt_chembl_id` | string | Salt identifier when present. |
| `natural_product` / `prodrug` / `polymer_flag` | boolean | Boolean flags from ChEMBL (nullable). |
| `black_box_warning` | string | Black box warning status. |
| `first_approval` | string | Year of first approval. |
| `max_phase` | string | Clinical development phase. |
| `canonical_smiles` | string | ChEMBL canonical SMILES. |
| `standard_inchi` | string | Standard InChI. |
| `standard_inchi_key` | string | Standard InChI key. |
| `molecule_type` | string | Molecule type classification. |
| `oral` / `parenteral` / `topical` | string/bool | Administration flags. |
| `pref_name` | string | Preferred name. |
| `pubchem_canonical_smiles` | string | PubChem canonical SMILES. |
| `pubchem_cid` | string/int | PubChem CID. |
| `pubchem_inchi` / `pubchem_inchikey` / `pubchem_isomeric_smiles` / `pubchem_iupac_name` / `pubchem_molecular_formula` | string | Fields retrieved from PubChem. |
| `structure_type` | string | Structural type (e.g., `SMALL MOLECULE`). |
| `pipeline_version` | string | Package version. |
| `timestamp_utc` | string | Export timestamp. |

## Tissue export (`tissues`)

Schema: [`library/schemas/tissues.py`](../../library/schemas/tissues.py). Column
ordering follows `TISSUES_COLUMN_ORDER`.

### Columns

| Column | Type | Notes |
|--------|------|-------|
| `tissue_chembl_id` | string | Primary identifier returned by the `/tissue` endpoint. Placeholders for missing records preserve the requested identifier even if ChEMBL omits it. |
| `pref_name` | string | Preferred tissue name where available. |
| `uberon_id` | string | UBERON cross-reference supplied by ChEMBL. |
| `efo_id` | string | EFO cross-reference supplied by ChEMBL. |
| `bto_id` | string | BRENDA Tissue Ontology cross-reference supplied by ChEMBL. |
| `caloha_id` | string | Caloha anatomical ontology identifier supplied by ChEMBL. |
| `pipeline_version` | string | Package version inserted by :func:`library.pipelines.common.add_pipeline_metadata`. |
| `timestamp_utc` | string | UTC timestamp of the pipeline run in ISO 8601 format. |

### Sorting, nulls and companion artefacts

- Rows are sorted by `tissue_chembl_id` (ascending) to guarantee deterministic
  ordering. Ties are broken deterministically by the CSV writer.
- Cross-reference columns are emitted using the nullable pandas string type and
  serialised as empty fields (`<NA>` when reloaded with `dtype="string"`).
- Metadata columns are populated for every record, including placeholders for
  missing identifiers.
- When `--final-out` is omitted the CLI writes
  `output.tissue_<YYYYMMDD>.csv` alongside
  `output.tissue_<YYYYMMDD>.meta.yaml`,
  `output.tissue_<YYYYMMDD>_quality_report_table.csv` and
  `output.tissue_<YYYYMMDD>.quality.json`.

## Cell line export (`cellline`)

Schema: [`library/schemas/celllines.py`](../../library/schemas/celllines.py). The
columns follow `CELL_LINE_COLUMN_ORDER` and include the pipeline metadata
(`pipeline_version`, `timestamp_utc`).

| Column | Type | Notes |
|--------|------|-------|
| `cell_chembl_id` | string | Primary ChEMBL identifier (required, unique). |
| `cell_name` | string | Preferred cell line name. |
| `cell_description` | string | Free-text description from ChEMBL. |
| `cell_id` | integer | Internal numeric identifier provided by ChEMBL (nullable). |
| `cell_source_organism` | string | Organism name for the originating tissue. |
| `cell_source_tax_id` | integer | NCBI taxonomy identifier for the source organism (nullable). |
| `cell_source_tissue` | string | Tissue or organ where the cell line originates. |
| `cellosaurus_id` | string | Cellosaurus accession when available. |
| `cl_lincs_id` | string | LINCS identifier (nullable). |
| `clo_id` | string | CLO (Cell Line Ontology) identifier. |
| `efo_id` | string | EFO cross-reference (nullable). |
| `pipeline_version` | string | Package version recorded during export. |
| `timestamp_utc` | string | UTC timestamp of the export. |

Values are normalised using `normalize_cell_lines`: identifier columns are
trimmed and coerced to `string`, numeric fields use pandas nullable integer
types (`Int64`). Missing values are exported as empty strings to preserve CSV
compatibility. Rows are sorted by `cell_chembl_id` to guarantee deterministic
ordering.

## Quality reports

The CSV quality report contains the following columns:

| Column | Description |
|--------|-------------|
| `column` | Column name analysed. |
| `row_count` | Total rows inspected. |
| `non_null` | Count of non-null entries. |
| `non_empty_ratio` | Fraction of values passing `_non_empty_mask`. |
| `distinct_count` | Number of distinct values. |
| `numeric_min` / `numeric_mean` / `numeric_max` | Numeric summary (for coercible columns). |
| `pattern_email_ratio`, `pattern_doi_ratio`, `pattern_url_ratio`, `pattern_issn_ratio` | Pattern coverage metrics. |
| `bool_like_ratio` | Share of boolean-like values. |

The JSON summary reproduces the same metrics and adds severity tags (`info`,
`warn`, `error`) reflecting threshold breaches configured in `system.doc_quality`.

# Data Schema

## Input Tables

### activity.csv
- **Purpose:** seed list of ChEMBL activity identifiers consumed by the `get_activity_data` CLI before querying the API.【F:config.yaml†L25-L31】【356235†L1-L5】

| Column | Data type | Source | Description |
| --- | --- | --- | --- |
| activity_chembl_id | string | External CSV of identifiers | Activity identifier in ChEMBL; used to request `/activity` records.|

### assay.csv
- **Purpose:** list of ChEMBL assay identifiers processed by the `get_assay_data` pipeline.【F:config.yaml†L32-L36】【29bb29†L1-L5】

| Column | Data type | Source | Description |
| --- | --- | --- | --- |
| assay_chembl_id | string | External CSV of identifiers | ChEMBL assay ID used to call `/assay`.|

### documents.csv
- **Purpose:** collection of ChEMBL document identifiers that feeds the `get_document_data` CLI (`chembl`/`all` modes).【F:config.yaml†L49-L56】【fca225†L1-L5】

| Column | Data type | Source | Description |
| --- | --- | --- | --- |
| document_chembl_id | string | External CSV of identifiers | Primary publication identifier in ChEMBL.|

### targets.csv
- **Purpose:** set of ChEMBL target identifiers used by the `get_target_data` pipeline.【F:config.yaml†L62-L80】【086d4d†L1-L5】

| Column | Data type | Source | Description |
| --- | --- | --- | --- |
| target_chembl_id | string | External CSV of identifiers | ChEMBL target ID passed to the `/target` endpoint.|

### testitem.csv
- **Purpose:** list of ChEMBL molecule identifiers that `get_testitem_data` enriches with structural and PubChem attributes.【F:config.yaml†L37-L41】【85074e†L1-L5】

| Column | Data type | Source | Description |
| --- | --- | --- | --- |
| molecule_chembl_id | string | External CSV of identifiers | ChEMBL molecule ID used when calling `/molecule`.|

## Output Tables

### activity.csv (processed export)
- **Purpose:** normalized set of experimental activity measurements retrieved from the ChEMBL API, extended with derived bounds, structured annotations, and pipeline metadata.【F:scripts/get_activity_data.py†L63-L357】【F:library/chembl_assay.py†L62-L111】【F:schemas/activities.py†L16-L78】【F:library/pipeline_metadata.py†L60-L84】

| Column | Data type | Source | Description |
| --- | --- | --- | --- |
| activity_id | string/number | ChEMBL `/activity` | Internal activity record identifier mirroring `ACTIVITY_ID` from the API.|
| molecule_chembl_id | string | ChEMBL `/activity` | Identifier of the molecule whose activity was measured.|
| assay_chembl_id | string | ChEMBL `/activity` | Identifier of the assay that produced the activity.|
| activity_comment | string | ChEMBL `/activity` | Free-text comment describing the experiment or result.|
| assay_description | string | ChEMBL `/activity` | Description of the linked assay.|
| assay_variant_accession | string | ChEMBL `/activity` | UniProt accession of the assayed protein variant when provided.|
| assay_variant_mutation | string | ChEMBL `/activity` | Mutation details for the assayed protein.|
| bao_format | string | ChEMBL `/activity` | BioAssay Ontology format identifier.|
| bao_label | string | ChEMBL `/activity` | BAO label describing the result type.|
| data_validity_comment | string | ChEMBL `/activity` | Validation comment supplied by the source.|
| data_validity_description | string | ChEMBL `/activity` | Detailed description of the validation status.|
| document_chembl_id | string | ChEMBL `/activity` | Publication identifier linked to the experiment.|
| pchembl_value | string/number | ChEMBL `/activity` | pChEMBL logarithmic value if computed by ChEMBL.|
| potential_duplicate | boolean/string | ChEMBL `/activity` | Flag indicating that the activity may duplicate another record.|
| qudt_units | string | ChEMBL `/activity` | Units expressed using the QUDT vocabulary when available.|
| record_id | string/number | ChEMBL `/activity` | Identifier of the originating record in the source database.|
| relation | string | ChEMBL `/activity` | Raw comparison operator (`<`, `<=`, `=`, `>=`, etc.) attached to the value.|
| src_assay_id | string/number | ChEMBL `/activity` | Source-specific assay identifier.|
| src_id | string/number | ChEMBL `/activity` | Identifier of the contributing data source.|
| standard_relation | string | ChEMBL `/activity` | Normalized comparison operator applied during processing.|
| standard_units | string | ChEMBL `/activity` | Normalized measurement units.|
| type | string | ChEMBL `/activity` | Original measurement type returned by the API.|
| units | string | ChEMBL `/activity` | Original units associated with `value`.|
| value | string/number | ChEMBL `/activity` | Raw measurement reported by the API.|
| standard_type | string | ChEMBL `/activity` | Normalized measurement type constrained by configuration (e.g., `IC50`, `EC50`, `Ki`, `KD`).|
| standard_value | number (float) | ChEMBL `/activity` | Normalized numeric value in molar units; guaranteed non-negative.|
| standard_lower_value | string/number | ChEMBL `/activity` | Lower bound supplied by ChEMBL when the measurement is a range.|
| standard_upper_value | string/number | ChEMBL `/activity` | Upper bound supplied by ChEMBL when available.|
| lower_value | number (float) | Derived field | Consolidated lower limit inferred from ChEMBL bounds, relations or uncertainty metadata.|
| upper_value | number (float) | Derived field | Consolidated upper limit inferred from ChEMBL bounds, relations or uncertainty metadata.|
| activity_properties | JSON string | Derived field | Canonical JSON payload describing the measurement, assay context and enrichment flags.|
| action_type | string | Derived field | Normalized action label (e.g., `PAM`, `NAM`) detected from annotations and configured lookups.|
| properties_hash | string | Derived field | SHA-256 checksum of `activity_properties` for change tracking.|
| pipeline_version | string | Pipeline metadata | Version of the `chembl-data-acquisition` package stamped onto the export.|
| timestamp_utc | ISO 8601 string | Pipeline metadata | UTC timestamp indicating when the export was created.|

### assay.csv (processed export)
- **Purpose:** aggregated assay descriptions with a per-target count of co-occurring assays and pipeline metadata.【F:scripts/get_assay_data.py†L47-L167】【F:library/chembl_assay.py†L22-L59】【F:library/assay_postprocessing.py†L1-L41】【F:schemas/assays.py†L1-L85】【F:library/pipeline_metadata.py†L60-L84】

| Column | Data type | Source | Description |
| --- | --- | --- | --- |
| assay_chembl_id | string | ChEMBL `/assay` | Primary assay identifier.|
| ASSAY_ID | string | ChEMBL `/assay` | Secondary identifier taken from the contributing source file.|
| Target TYPE | string | ChEMBL `/assay` | Target category (cell-based, protein, etc.).|
| accession | string | ChEMBL `/assay` | UniProt accession of the target when provided.|
| acts_per_assay_step5 | string/number | ChEMBL `/assay` | Number of associated activities for step 5; preserved as returned.|
| assay_cell_type | string | ChEMBL `/assay` | Cell type used during the experiment.|
| assay_subcellular_fraction | string | ChEMBL `/assay` | Subcellular fraction tested in the assay.|
| assay_tissue | string | ChEMBL `/assay` | Tissue or organ from which the material was derived.|
| bao_format | string | ChEMBL `/assay` | BioAssay Ontology format identifier.|
| cited_assay_corr | string/boolean | ChEMBL `/assay` | Whether the assay is cited as correlated.|
| description | string | ChEMBL `/assay` | Human-readable description of the assay.|
| document_chembl_id | string | ChEMBL `/assay` | Identifier of the source publication.|
| error_assay_corr | string/boolean | ChEMBL `/assay` | Error indicator for the correlation citation flag.|
| higly_correlated_cit | string/boolean | ChEMBL `/assay` | Flag for highly correlated citations.|
| isoform | string | ChEMBL `/assay` | Protein isoform reference.|
| month | string/number | ChEMBL `/assay` | Month of publication as delivered by the source.|
| mutation | string | ChEMBL `/assay` | Mutation affecting the target protein.|
| shuffled_cit | string/boolean | ChEMBL `/assay` | Indicator of a shuffled citation.|
| shuffled_target_assay | string/boolean | ChEMBL `/assay` | Indicator of a shuffled target/assay pair.|
| substrate_name | string | ChEMBL `/assay` | Name of the substrate used in the assay.|
| target_chembl_id | string | ChEMBL `/assay` | Linked ChEMBL target identifier.|
| target_name | string | ChEMBL `/assay` | Human-readable target name.|
| version | string/number | ChEMBL `/assay` | Version of the assay record.|
| year | string/number | ChEMBL `/assay` | Publication year.|
| assay_with_same_target | number (int) | `postprocess_assays` | Count of assays sharing the same `document_chembl_id` and `target_chembl_id`.|
| pipeline_version | string | Pipeline metadata | `chembl-data-acquisition` package version.|
| timestamp_utc | ISO 8601 string | Pipeline metadata | Export timestamp in UTC.|

### documents.csv (processed export)
- **Purpose:** consolidated bibliographic metadata combining ChEMBL, PubMed, Semantic Scholar, OpenAlex and Crossref outputs with computed publication classes.【F:scripts/get_document_data.py†L18-L884】【F:library/document_pipeline.py†L20-L119】【F:schemas/documents.py†L14-L119】【F:library/document_postprocessing.py†L18-L154】【F:library/pipeline_metadata.py†L60-L84】

| Column | Data type | Source | Description |
| --- | --- | --- | --- |
| document_chembl_id | string | ChEMBL `/document` | Primary identifier of the publication.|
| title | string | ChEMBL `/document` | Publication title.|
| abstract | string | ChEMBL `/document` | Abstract text.|
| doi | string | ChEMBL `/document` and aggregators | DOI in the original format before normalization.|
| year | number/string | ChEMBL `/document` | Publication year.|
| journal | string | ChEMBL `/document` | Full journal title.|
| journal_abbrev | string | ChEMBL `/document` | Journal abbreviation.|
| volume | number/string | ChEMBL `/document` | Volume number.|
| issue | number/string | ChEMBL `/document` | Issue number.|
| first_page | number/string | ChEMBL `/document` | First page of the article.|
| last_page | number/string | ChEMBL `/document` | Last page of the article.|
| pubmed_id | number/string | ChEMBL `/document` | PubMed ID captured by ChEMBL when available.|
| authors | string | ChEMBL `/document` | List of authors.|
| source | string | ChEMBL `/document` | Origin of the record (typically `ChEMBL`).|
| doi_normalised | string | Derived field | Canonical DOI computed across sources.|
| publication_types_normalised | string | Derived field | Unified publication-type terms collected from all sources.|
| publication_type_score_review | number (int) | Derived field | Score indicating evidence for a review article.|
| publication_type_score_experimental | number (int) | Derived field | Score indicating experimental content.|
| publication_type_score_unknown | number (int) | Derived field | Score assigned to unknown classes.|
| publication_class | string | Derived field | Final publication label (`review`, `experimental`, etc.).|
| PubMed.PMID | number/string | PubMed API | PubMed identifier returned by the API.|
| PubMed.DOI | string | PubMed API | DOI as normalized from PubMed.|
| PubMed.ArticleTitle | string | PubMed API | Article title from PubMed.|
| PubMed.Abstract | string | PubMed API | Abstract text from PubMed.|
| PubMed.JournalTitle | string | PubMed API | Journal title from PubMed.|
| PubMed.Volume | number/string | PubMed API | Volume reported by PubMed.|
| PubMed.Issue | number/string | PubMed API | Issue reported by PubMed.|
| PubMed.StartPage | number/string | PubMed API | First page according to PubMed.|
| PubMed.EndPage | number/string | PubMed API | Last page according to PubMed.|
| PubMed.PublicationType | string | PubMed API | Raw publication-type terms.|
| PubMed.MeSH_Descriptors | string | PubMed API | MeSH descriptors linked to the record.|
| PubMed.MeSH_Qualifiers | string | PubMed API | MeSH qualifiers associated with the record.|
| PubMed.ChemicalList | string | PubMed API | List of chemical entities mentioned.|
| PubMed.DayRevised | number/string | PubMed API | Day of the last revision.|
| PubMed.MonthRevised | number/string | PubMed API | Month of the last revision.|
| PubMed.YearRevised | number/string | PubMed API | Year of the last revision.|
| PubMed.YearCompleted | number/string | PubMed API | Year when indexing was completed.|
| PubMed.MonthCompleted | number/string | PubMed API | Month when indexing was completed.|
| PubMed.DayCompleted | number/string | PubMed API | Day when indexing was completed.|
| PubMed.Error | string | PubMed API | Error message returned by the PubMed API.|
| PubMed.ISSN | string | PubMed API | Journal ISSN as reported by PubMed.|
| scholar.PMID | number/string | Semantic Scholar | PubMed ID linked by Semantic Scholar.|
| scholar.Venue | string | Semantic Scholar | Publication venue recorded by Semantic Scholar.|
| scholar.PublicationTypes | string | Semantic Scholar | Publication-type terms from Semantic Scholar.|
| scholar.SemanticScholarId | string | Semantic Scholar | Semantic Scholar internal identifier.|
| scholar.ExternalIds | string | Semantic Scholar | External identifiers serialized as JSON.|
| scholar.DOI | string | Semantic Scholar | DOI value returned by Semantic Scholar.|
| scholar.Error | string | Semantic Scholar | Error diagnostics from the Semantic Scholar API.|
| OpenAlex.PublicationTypes | string | OpenAlex | Publication-type list returned by OpenAlex.|
| OpenAlex.TypeCrossref | string | OpenAlex/Crossref | Crossref type relayed by OpenAlex.|
| OpenAlex.Genre | string | OpenAlex | OpenAlex genre classification.|
| OpenAlex.Id | string | OpenAlex | OpenAlex record identifier.|
| OpenAlex.Venue | string | OpenAlex | Publication venue according to OpenAlex.|
| OpenAlex.MeshDescriptors | string | OpenAlex | MeSH descriptors supplied by OpenAlex.|
| OpenAlex.MeshQualifiers | string | OpenAlex | MeSH qualifiers supplied by OpenAlex.|
| OpenAlex.Error | string | OpenAlex | Error message returned by OpenAlex.|
| crossref.Type | string | Crossref | Crossref article type.|
| crossref.Subtype | string | Crossref | Crossref subtype.|
| crossref.Title | string | Crossref | Title from Crossref metadata.|
| crossref.Subtitle | string | Crossref | Subtitle from Crossref metadata.|
| crossref.Subject | string | Crossref | Subject categories from Crossref.|
| crossref.Error | string | Crossref | Error message returned by Crossref.|
| date_code | string | Post-processing | Auxiliary date code generated during normalization.|
| Index | number/string | Post-processing | Row index preserved from the source table.|
| PubMed.is_review | boolean/string | Post-processing | Derived review flag based on PubMed terminology.|
| scholar.is_review | boolean/string | Post-processing | Review flag derived from Semantic Scholar.|
| OpenAlex.is_review | boolean/string | Post-processing | Review flag derived from OpenAlex metadata.|
| pipeline_version | string | Pipeline metadata | `chembl-data-acquisition` package version.|
| timestamp_utc | ISO 8601 string | Pipeline metadata | Export timestamp in UTC.|

### target.csv (final export)
- **Purpose:** unified target table combining ChEMBL, UniProt and IUPHAR attributes while preserving the deterministic column order expected by downstream BI processes.【F:scripts/get_target_data.py†L31-L1148】【F:library/chembl_target.py†L185-L315】【F:library/target_postprocessing.py†L181-L442】【F:schemas/targets.py†L17-L214】【F:library/pipeline_metadata.py†L60-L84】

| Column | Data type | Source | Description |
| --- | --- | --- | --- |
| target_chembl_id | string | ChEMBL `/target` | Primary identifier of the target.|
| uniprot_id_primary | string | UniProt post-processing | Canonical UniProt accession (normalized `uniProtkbId`).|
| uniprot_ids_all | string | UniProt post-processing | Pipe-delimited list of all known UniProt accessions (primary, secondary, mapping).|
| isoform_ids | string | UniProt merge | List of isoform identifiers.|
| isoform_names | string | UniProt merge | Names of the isoforms.|
| isoform_synonyms | string | UniProt merge | Isoform synonyms converted to lowercase.|
| hgnc_id | string | HGNC via cross references | HGNC identifier stripped of the `HGNC:` prefix.|
| gene_symbol | string | Post-processing | Final gene symbol rendered in upper case.|
| protein_name_canonical | string | UniProt/ChEMBL | Canonical protein name.|
| protein_name_alt | string | UniProt/ChEMBL | Pipe-delimited list of alternative protein names.|
| organism | string | ChEMBL/UniProt | Genus name representing the organism.|
| taxon_id | string | ChEMBL | NCBI taxonomy identifier.|
| lineage_superkingdom | string | UniProt taxonomy merge | Taxonomic superkingdom.|
| lineage_phylum | string | UniProt taxonomy merge | Taxonomic phylum.|
| lineage_class | string | UniProt taxonomy merge | Taxonomic class.|
| sequence_length | string | UniProt | Amino-acid sequence length.|
| features_signal_peptide | string | UniProt features | Presence of signal peptide annotations (lowercase).|
| features_transmembrane | string | UniProt features | Presence of transmembrane segments (lowercase).|
| features_topology | string | UniProt features | Reported membrane topology.|
| ptm_glycosylation | string | UniProt PTM | Glycosylation annotations.|
| ptm_lipidation | string | UniProt PTM | Lipidation annotations.|
| ptm_disulfide_bond | string | UniProt PTM | Disulfide bond annotations.|
| ptm_modified_residue | string | UniProt PTM | Other modified residues.|
| xref_chembl | string | Post-processing | Duplicate of `target_chembl_id` for compatibility.|
| xref_uniprot | string | Post-processing | Duplicate of `uniprot_id_primary`.|
| xref_ensembl | string | ChEMBL/UniProt | Ensembl cross references serialized as a string.|
| xref_iuphar | string | IUPHAR | IUPHAR target identifier.|
| gtop_target_id | string | Guide to PHARMACOLOGY | Identifier used in the Guide to PHARMACOLOGY database.|
| gtop_synonyms | string | Guide to PHARMACOLOGY | Synonyms reported by GToP (lowercase).|
| gtop_natural_ligands_n | string | Guide to PHARMACOLOGY | Number of natural ligands as reported by GToP.|
| gtop_interactions_n | string | Guide to PHARMACOLOGY | Number of recorded interactions.|
| gtop_function_text_short | string | Guide to PHARMACOLOGY | Short functional description.|
| type | string | Taxonomy classifier | Final organism class assigned by the taxonomy module (e.g., `Multicellular organism`).|
| uniprot_last_update | string | UniProt | Last update timestamp of the UniProt record.|
| uniprot_version | string | UniProt | UniProt record version.|
| pipeline_version | string | Pipeline metadata | `chembl-data-acquisition` package version.|
| timestamp_utc | ISO 8601 string | Pipeline metadata | Export timestamp in UTC.|
| pfam | string | UniProt | Pfam domains reported by UniProt.|
| interpro | string | UniProt | InterPro domains reported by UniProt.|
| xref_pdb | string | UniProt | PDB structure references.|
| xref_alphafold | string | UniProt/AlphaFold | AlphaFold model references.|
| hgnc_name | string | HGNC | HGNC gene name.|
| uniProtkbId | string | UniProt | Original UniProt identifier before normalization.|
| secondaryAccessions | string | UniProt | List of secondary UniProt accessions.|
| recommendedName | string | UniProt | Recommended protein name from UniProt.|
| geneName | string | UniProt | Primary gene symbol from UniProt.|
| secondaryAccessionNames | string | UniProt | Names associated with secondary accessions.|
| molecular_function | string | UniProt GO | Gene Ontology molecular function annotations.|
| cellular_component | string | UniProt GO | Gene Ontology cellular component annotations.|
| subcellular_location | string | UniProt | Reported subcellular localization.|
| topology | string | UniProt | Transmembrane topology classification.|
| transmembrane | string | UniProt | Transmembrane segment indicator (lowercase).|
| intramembrane | string | UniProt | Intramembrane segment indicator.|
| glycosylation | string | UniProt | Glycosylation annotation (lowercase).|
| lipidation | string | UniProt | Lipidation annotation.|
| disulfide_bond | string | UniProt | Disulfide bond annotation.|
| modified_residue | string | UniProt | Miscellaneous residue modifications.|
| phosphorylation | string | UniProt | Phosphorylation annotation.|
| acetylation | string | UniProt | Acetylation annotation.|
| ubiquitination | string | UniProt | Ubiquitination annotation.|
| signal_peptide | string | UniProt | Signal peptide annotation (lowercase).|
| propeptide | string | UniProt | Propeptide annotation.|
| GuidetoPHARMACOLOGY | string | IUPHAR/GToP | Primary Guide to PHARMACOLOGY identifier.|
| family | string | IUPHAR | IUPHAR family name.|
| SUPFAM | string | UniProt/IUPHAR | Superfamily descriptor preserved from the source.|
| PROSITE | string | UniProt | PROSITE motif annotations.|
| InterPro | string | UniProt | InterPro domain annotations (alternate field).|
| Pfam | string | UniProt | Pfam domain annotations (alternate field).|
| PRINTS | string | UniProt | PRINTS motif annotations.|
| TCDB | string | UniProt | Transporter Classification Database annotation.|
| pref_name | string | ChEMBL | Preferred target name from ChEMBL.|
| target_type | string | ChEMBL | Target type reported by ChEMBL (`SINGLE PROTEIN`, `PROTEIN FAMILY`, etc.).|
| tax_id | string | ChEMBL | Taxonomy identifier as stored in ChEMBL.|
| species_group_flag | string | ChEMBL | Species group flag from ChEMBL.|
| target_sort_order | string | Taxonomy classifier | Deterministic sort key derived from the classifier for analytics pivots.|
| gene_index | string | Taxonomy classifier | Classifier-provided gene ordering indicator for curation worksheets.|
| taxon_index | string | Taxonomy classifier | Derived taxonomy ordering indicator used by downstream reports.|
| multifunctional_enzyme | boolean/string | Taxonomy classifier | Flag raised when taxonomy rules classify the target as multifunctional.|
| unicellular_organism | boolean/string | Taxonomy classifier | Classifier flag indicating unicellular or viral organisms.|
| target_components | JSON string | ChEMBL | Serialized representation of target components.|
| protein_classifications | JSON string | ChEMBL | Serialized protein classification hierarchy.|
| cross_references | JSON string | ChEMBL | Serialized list of external cross references.|
| gene_symbol_list | string | Post-processing | Pipe-delimited list of gene symbols (upper case).|
| protein_synonym_list | string | Post-processing | Pipe-delimited list of protein synonyms (lower case).|
| reactions | string | ChEMBL | Reactions associated with the target, when provided.|
| reaction_ec_numbers | string | Post-processing | Pipe-delimited list of EC numbers aggregated from synonyms and cross references.|
| protein_class_pred_L1 | string | IUPHAR supplemental data | Predicted protein class level 1.|
| protein_class_pred_L2 | string | IUPHAR supplemental data | Predicted protein class level 2.|
| protein_class_pred_L3 | string | IUPHAR supplemental data | Predicted protein class level 3.|
| protein_class_pred_rule_id | string | IUPHAR supplemental data | Identifier of the classification rule.|
| protein_class_pred_evidence | string | IUPHAR supplemental data | Evidence type backing the prediction.|
| protein_class_pred_confidence | string | IUPHAR supplemental data | Confidence score for the prediction.|
| iuphar_target_id | string | IUPHAR | IUPHAR target identifier.|
| iuphar_family_id | string | IUPHAR | IUPHAR family identifier.|
| iuphar_type | string | IUPHAR | IUPHAR object type.|
| iuphar_class | string | IUPHAR | IUPHAR class.|
| iuphar_subclass | string | IUPHAR | IUPHAR subclass.|
| iuphar_chain | string | IUPHAR | Chain or subunit identifier.|
| iuphar_name | string | IUPHAR | Official IUPHAR name.|
| iuphar_full_id_path | string | IUPHAR | Full identifier path within the IUPHAR hierarchy.|
| iuphar_full_name_path | string | IUPHAR | Full name path within the IUPHAR hierarchy.|

> The taxonomy classifier consumes `genus`, `lineage_superkingdom`, `lineage_phylum`,
> `lineage_class`, `taxon_id`, and `species_group_flag` to populate `type` and the
> derived flags listed above.

### testitem.csv (processed export)
- **Purpose:** enriched description of ChEMBL compounds combining parent hierarchy, structural attributes, PubChem augmentation, catalog flags, and pipeline metadata.【F:scripts/get_testitem_data.py†L36-L356】【F:library/testitem_enrichment.py†L151-L239】【F:schemas/testitems.py†L14-L47】【F:library/pipeline_metadata.py†L60-L84】

| Column | Data type | Source | Description |
| --- | --- | --- | --- |
| molecule_chembl_id | string | ChEMBL `/molecule` | Primary molecule identifier.|
| parent_molecule_chembl_id | string | ChEMBL hierarchy/catalog | Identifier of the parent molecule used for salt roll-ups.|
| salt_chembl_id | string | Derived field | Copies the child identifier when the molecule is a salt; otherwise blank or `-` per configuration.|
| natural_product | boolean (nullable) | Molecule catalog | Normalized flag indicating natural-product origin.|
| prodrug | boolean (nullable) | Molecule catalog | Normalized flag identifying prodrug records.|
| polymer_flag | boolean (nullable) | Molecule catalog | Normalized polymer indicator from the catalog.|
| black_box_warning | string/boolean | ChEMBL `/molecule` | Black-box warning flag.|
| first_approval | string/date | ChEMBL `/molecule` | Year or date of the first regulatory approval (as provided).|
| max_phase | string | ChEMBL `/molecule` | Maximum clinical phase reached.|
| canonical_smiles | string | ChEMBL `/molecule` | Canonical SMILES reported by ChEMBL.|
| standard_inchi | string | ChEMBL `/molecule` | Standard InChI.|
| standard_inchi_key | string | ChEMBL `/molecule` | InChIKey derived from the standard InChI.|
| molecule_type | string | ChEMBL `/molecule` | Molecule category (e.g., small molecule, antibody).|
| oral | string/boolean | ChEMBL `/molecule` | Oral dosage form indicator.|
| parenteral | string/boolean | ChEMBL `/molecule` | Parenteral administration indicator.|
| pref_name | string | ChEMBL `/molecule` | Preferred name of the molecule.|
| pubchem_canonical_smiles | string | PubChem | Canonical SMILES retrieved from PubChem.|
| pubchem_cid | string/number | PubChem | PubChem Compound ID resolved from SMILES.|
| pubchem_inchi | string | PubChem | Standard InChI from PubChem.|
| pubchem_inchikey | string | PubChem | InChIKey from PubChem.|
| pubchem_isomeric_smiles | string | PubChem | Isomeric SMILES from PubChem.|
| pubchem_iupac_name | string | PubChem | IUPAC name returned by PubChem.|
| pubchem_molecular_formula | string | PubChem | Molecular formula reported by PubChem.|
| structure_type | string | ChEMBL `/molecule` | Structure representation type (e.g., `MOLFILE`).|
| topical | string/boolean | ChEMBL `/molecule` | Topical administration indicator.|
| pipeline_version | string | Pipeline metadata | `chembl-data-acquisition` package version.|
| timestamp_utc | ISO 8601 string | Pipeline metadata | Export timestamp in UTC.|


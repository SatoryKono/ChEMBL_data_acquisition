# Missing canonical SMILES and InChI identifiers in `get_testitem_data`

## Summary

The ChEMBL extraction step **does** download the structural identifiers (SMILES
and InChI) for every requested molecule, but these attributes never make it past
the initial normalisation step in `get_testitem`. The DataFrame returned by the
ChEMBL API contains nested column names such as
`molecule_structures.canonical_smiles`, yet the pipeline only renames the
`pubchem.*` subset before reindexing to the canonical schema. As a result the
`molecule_structures.*` columns are silently dropped during `reindex`, leaving
empty `canonical_smiles`, `standard_inchi`, and `standard_inchi_key` fields in
`raw.testitem.csv` and, consequently, in `output.testitem.csv`.

## Detailed findings

1. **Fields requested from ChEMBL**  \
   Both the configuration (`Config.testitem.fields`) and the default fallback
   (`TESTITEM_FIELD_DEFAULTS`, see
   `library/config/models.py`) explicitly request
   `molecule_structures.canonical_smiles`,
   `molecule_structures.standard_inchi`, and
   `molecule_structures.standard_inchi_key` alongside the `pubchem.*`
   enrichment attributes, so the API response does contain these values.

2. **Data loss during `get_testitem` normalisation**  \
   The helper (`library/pipelines/assay/chembl_assay.py`) renames raw columns
   using `TESTITEM_STRUCTURE_COLUMNS` and then calls
   `df.reindex(columns=TESTITEM_COLUMNS)`. `TESTITEM_STRUCTURE_COLUMNS` only
   defines mappings for `pubchem.*`, so the `molecule_structures.*` columns keep
   their fully qualified names. The final `reindex` discards those unmatched
   columns and creates new, all-null `canonical_smiles`, `standard_inchi`, and
   `standard_inchi_key` columns to fit the schema.

3. **Downstream stages keep the empty identifiers**  \
   Subsequent pipeline stages (parent enrichment, schema validation, CSV
   writers, and post-processing) operate exclusively on the `TESTITEM_COLUMNS`
   layout. Because the identifiers were lost earlier, these steps simply pass
   through the empty columns into `raw.testitem.csv` and `output.testitem.csv`.

4. **Raw vs. final output**  \
   Comparing `raw.testitem.csv` and `output.testitem.csv` shows that both files
   already contain empty structural identifier columns, confirming that the loss
   happened before the export logic.

## Suggested fix

Extend `TESTITEM_STRUCTURE_COLUMNS` (or add an explicit rename/flattening step)
to include the three `molecule_structures.*` attributes **before** reindexing.
Once those columns are renamed to `canonical_smiles`, `standard_inchi`, and
`standard_inchi_key`, they will persist through aggregation, schema validation,
and export without further changes.

# Missing canonical SMILES and InChI identifiers in `get_testitem_data`

## Summary

An investigation of the `get_testitem_data` pipeline shows that the ChEMBL API
response contains the requested structural identifiers, but they are dropped
during normalisation of the `/molecule` payload. The pipeline asks ChEMBL for
`molecule_structures.canonical_smiles`, `molecule_structures.standard_inchi`,
and `molecule_structures.standard_inchi_key`, yet the renaming step only maps
`pubchem.*` fields. When the frame is reindexed to the canonical output schema,
the untouched `molecule_structures.*` columns are discarded and replaced with
empty columns, which explains the missing data in `output.testitem.csv`.

## Detailed findings

1. **Fields requested from ChEMBL**  \
   Both the configuration (`Config.testitem.fields`) and the default fallback
   (`TESTITEM_FIELD_DEFAULTS`) explicitly request the structural attributes from
   the ChEMBL API alongside PubChem enrichment fields.  
   Relevant references: `library/config/models.py` (`TESTITEM_FIELD_DEFAULTS`).

2. **Data loss during normalisation**  \
   `get_testitem` normalises each JSON chunk and renames a handful of columns
   via `TESTITEM_STRUCTURE_COLUMNS` before reindexing to `TESTITEM_COLUMNS`.  
   The rename mapping only covers `pubchem.*` keys, so the columns named
   `molecule_structures.canonical_smiles`, `molecule_structures.standard_inchi`,
   and `molecule_structures.standard_inchi_key` are not converted to their
   flattened counterparts. When the DataFrame is reindexed, pandas creates empty
   `canonical_smiles`, `standard_inchi`, and `standard_inchi_key` columns while
   dropping the original names.  
   Relevant references: `library/pipelines/assay/chembl_assay.py` (definitions
   of `TESTITEM_STRUCTURE_COLUMNS`, `TESTITEM_COLUMNS`, and the `get_testitem`
   implementation).

3. **Downstream processing preserves column structure**  \
   Later pipeline stages (`finalize_output`, schema normalisation, post-
   processing) operate on the already-flattened column list and do not reintroduce
   or remove these identifiers, so the loss occurs before aggregation/output.

## Suggested fix

Extend `TESTITEM_STRUCTURE_COLUMNS` (or add a dedicated rename step) so that
`molecule_structures.canonical_smiles`, `molecule_structures.standard_inchi`,
and `molecule_structures.standard_inchi_key` are mapped to the canonical output
column names before the frame is reindexed. This ensures the structural data
survives the transformation into `raw.testitem.csv` and `output.testitem.csv`.

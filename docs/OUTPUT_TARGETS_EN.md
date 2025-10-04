# Target isoform post-processing

The target acquisition pipeline now emits an additional isoform summary derived
from the canonical `output.target_<date>.csv` export. The logic mirrors the
Power Query (M) workbook used in analytics dashboards which makes the Python
implementation a drop-in replacement with byte-identical results.

## Workflow overview

1. Select the required input columns:
   `isoform_synonyms`, `isoform_names`, `isoform_ids`, `uniprot_id_primary`,
   `target_chembl_id`.
2. Normalise case for `isoform_synonyms` and `isoform_names` only (IDs preserve
   their original casing).
3. Split pipe-delimited strings (`|`) into trimmed lists and drop empty tokens.
4. Align name, identifier and synonym lists by index using null padding.
5. Expand each synonym token by splitting on `:` and generating additional
   variants with the substrings `pde` and `pld` removed.
6. Build two tables:
   - Names: trimmed isoform names excluding empty, `n/a` and `none` values.
   - Synonym tokens: flattened variants filtered with the same rules.
7. Combine the two tables, deduplicate in three passes and apply a stable sort:
   - Distinct on four columns (`id`, `name`, `target_chembl_id`,
     `uniprot_id_primary`).
   - Stable `mergesort` by `uniprot_id_primary`, then `id` (nulls first).
   - Distinct on three columns (`id`, `target_chembl_id`, `name`).
   - Final distinct on two columns (`id`, `name`).
8. Export the table with columns ordered exactly as
   `['id', 'uniprot_id_primary', 'target_chembl_id', 'name']` using UTF-8
   encoding and comma separators.

Every step is deterministic and mirrors the behaviour of the reference M script.

## Example

### Input (`output.target_20250301.csv`)

| target_chembl_id | uniprot_id_primary | isoform_names  | isoform_ids         | isoform_synonyms        |
|------------------|--------------------|----------------|---------------------|-------------------------|
| CHEMBL1000       | Q11111             | Alpha\|A       | IsoA\|IsoA-Variant  | PDE3A:alpha\|Alpha variant |
| CHEMBL1001       | Q22222             |  n/a \| none   | IsoB                |                         |
| CHEMBL1003       |                    |                | MixedID             | PLDALPHA\|pdeBeta        |

### Output (`isoform.output.target_20250301.csv`)

| id            | uniprot_id_primary | target_chembl_id | name          |
|---------------|--------------------|------------------|---------------|
|               |                    | CHEMBL1003       | pdebeta       |
|               |                    | CHEMBL1003       | beta          |
| MixedID       |                    | CHEMBL1003       | pldalpha      |
| MixedID       |                    | CHEMBL1003       | alpha         |
| IsoA          | Q11111             | CHEMBL1000       | alpha         |
| IsoA          | Q11111             | CHEMBL1000       | pde3a         |
| IsoA          | Q11111             | CHEMBL1000       | 3a            |
| IsoA-Variant  | Q11111             | CHEMBL1000       | a             |
| IsoA-Variant  | Q11111             | CHEMBL1000       | alpha variant |

Rows are shown in export order. Empty identifiers are preserved when the source
contained names or tokens without a matching isoform ID.

## Usage

The post-processing step runs automatically at the end of
`python scripts/get_target_data.py all ...`. A manual invocation is also
available:

```python
from library.postprocessing import target

# produces isoform.output.target_20250301.csv next to the source file
isoform_path = target.process_targets("output.target_20250301.csv")
```

The loader retries common encodings (`utf-8`, `utf-8-sig`, `cp1252`) and emits a
concise progress summary. All deduplication and sorting steps are stable which
ensures the output is identical across repeated runs given the same input.

## Compatibility and guarantees

- Python ≥ 3.11 and pandas ≥ 2.3.
- No additional dependencies beyond the standard library and pandas.
- Deterministic CSV export in UTF-8 without BOM and comma separators.
- Exact equivalence with the legacy Power Query transformation down to the
  column order and row ordering.

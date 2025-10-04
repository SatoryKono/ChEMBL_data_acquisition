# Target pipeline output: isoform post-processing

This addendum complements the [Output specification](./en/OUTPUT.md) by detailing
the optional `isoform.*` post-processing stage executed after the target pipeline.
The same content is available in Russian in
[`OUTPUT_TARGETS_RU.md`](./OUTPUT_TARGETS_RU.md).

## Post-processing `isoform.*`

### Stages
1. **Read aggregated targets** — load the pipeline CSV using the encoding fallbacks
   described in `library.postprocessing.target.DEFAULT_ENCODINGS`.
2. **Normalise name fields** — trim whitespace, collapse sentinels (`"None"`,
   `"n/a"`), and lowercase synonym columns to mirror the legacy Power Query rules.
3. **Expand isoforms** — split `isoform_ids`, `isoform_names`, and
   `isoform_synonyms` on the pipe delimiter and align entries with `zip_longest`.
4. **Tokenise synonyms** — split colon-separated synonym groups, drop empty
   variants, preserve order, and generate per-token records via `_syn_expand`.
5. **Deduplicate & sort** — drop duplicate
   (`target_chembl_id`, `isoform_id`, `term`, `token`) combinations and sort with
   a stable mergesort before emission.
6. **Emit artefact** — write `isoform.output.<original>.csv` encoded as UTF-8 with
   Unix newlines alongside the source file (or into the chosen output folder).

### Invariants
- Output columns are exactly: `target_chembl_id`, `isoform_id`, `isoform_name`,
  `term`, `token`.
- `token` values are lowercase ASCII without punctuation; `term` keeps the
  normalised original text.
- Rows are globally sorted by (`target_chembl_id`, `isoform_id`, `term`, `token`).
- Duplicate tokens per isoform are removed while retaining the first occurrence.
- Empty isoform identifiers are permitted only when at least one synonym/token is
  present.
- The writer always uses `\n` line endings to guarantee cross-platform diffs.

### Example (input → output)

Input snapshot:
```
target_chembl_id,isoform_ids,isoform_names,isoform_synonyms
CHEMBL1824,ENSP00000350283,Nav1.7,"Nav1.7:SCN9A isoform 3"
CHEMBL1824,ENSP00000350284,"Nav1.7 isoform 2","Nav1.7 splice variant:SCN9A-2"
CHEMBL6130,,"",""
CHEMBL240,ENSP00000456012,ALK2,"Activin receptor-like kinase 2:ACVR1"
CHEMBL259,ENSP00000263253,FGFR4,"FGFR-4 isoform:Fibroblast growth factor receptor-4"
```

Generated artefact (first 10 lines):
```
target_chembl_id,isoform_id,isoform_name,term,token
CHEMBL1824,ENSP00000350283,Nav1.7,Nav1.7,7
CHEMBL1824,ENSP00000350283,Nav1.7,Nav1.7,nav1
CHEMBL1824,ENSP00000350283,Nav1.7,nav1.7,7
CHEMBL1824,ENSP00000350283,Nav1.7,nav1.7,nav1
CHEMBL1824,ENSP00000350283,Nav1.7,scn9a isoform 3,3
CHEMBL1824,ENSP00000350283,Nav1.7,scn9a isoform 3,isoform
CHEMBL1824,ENSP00000350283,Nav1.7,scn9a isoform 3,scn9a
CHEMBL1824,ENSP00000350284,Nav1.7 isoform 2,Nav1.7 isoform 2,2
CHEMBL1824,ENSP00000350284,Nav1.7 isoform 2,Nav1.7 isoform 2,7
CHEMBL1824,ENSP00000350284,Nav1.7 isoform 2,Nav1.7 isoform 2,isoform
```

### Determinism
- Stable mergesort ordering, explicit newline handling, and deterministic token
  extraction guarantee repeatable CSVs for identical inputs.
- The function avoids locale-dependent transformations; all string operations are
  ASCII and Unicode normalisation free.
- Randomness, timestamps, and external services are not involved.

### Version requirements
- `chembl-data-acquisition` **0.1.3** or newer (provides
  `library.postprocessing.target.process_targets`).
- Python **3.11** runtime (matching the project constraint).
- `pandas` **2.3.3** to ensure availability of the dataframe operations used.
- When invoked from the CLI helper, ensure the main output schema follows
  [`docs/en/OUTPUT.md`](./en/OUTPUT.md) so that column expectations match.


    # Target pipeline outputs: isoform post-processing

    This document extends the [output reference](./en/OUTPUT.md) and describes the
    optional `isoform.*` stage executed after the target pipeline has produced
    `output.target_YYYYMMDD.csv`. The Russian translation is available in
    [`OUTPUT_TARGETS_RU.md`](./OUTPUT_TARGETS_RU.md).

    ## `isoform.*` post-processing

    ### Steps
    1. **Load the aggregated CSV.** The file is read by trying `utf-8`,
       `utf-8-sig`, then `cp1252`. Only the columns `isoform_synonyms`,
       `isoform_names`, `isoform_ids`, `uniprot_id_primary` and
       `target_chembl_id` are retained.
    2. **Normalise text.** `isoform_synonyms` and `isoform_names` are converted to
       lowercase, split by `|`, trimmed and stripped of empty entries. The
        `isoform_ids` column is split by `|` without changing the letter case.
    3. **Align records.** Name, identifier and synonym lists are aligned by index
       with `null` padding when one of the lists is shorter.
    4. **Tokenise synonyms.** Every synonym is split on `":"`; empty segments are
       discarded. For each token the variants `[token, token without "pde",
       token without "pld"]` are generated and deduplicated while keeping the
       original order.
    5. **Build tables.** Two tables are produced: `names` from
       `isoform_names` and `synonyms` from the expanded tokens (renaming the
       `tokens` column to `name`). The tables are concatenated.
    6. **Deduplicate.** Three passes are executed in sequence:
       `distinct(id, name, target_chembl_id, uniprot_id_primary)`, stable sorting
       by `(uniprot_id_primary, id)` with `mergesort`, followed by
       `distinct(id, target_chembl_id, name)` and a final `distinct(id, name)`.
    7. **Export.** The resulting frame is written as
       `isoform.output.<basename>.csv` next to the input file (or the provided
       path) using UTF-8 without BOM.

    ### Invariants
    - Output columns are exactly `id`, `uniprot_id_primary`,
      `target_chembl_id`, `name` in that order.
    - `name` values are lowercase; `id` preserves the original `isoform_ids`
      casing and may be empty when no identifier was available.
    - Placeholders `"", "n/a", "none"` are filtered before combining the
      tables.
    - The deduplication sequence guarantees deterministic results for identical
      inputs.

    ### Example (input → output)

    Input excerpt `output.target_20250228.csv`:

    ```csv
    target_chembl_id,uniprot_id_primary,isoform_ids,isoform_names,isoform_synonyms
    CHEMBL3587,P12345,"ENSP0001|ENSP0002","Alpha|Beta","PDE3A:Alpha|Alpha variant|Beta"
    CHEMBL1250,P67890,"ENSP0003","Gamma","Gamma isoform 1|PLD1A"
    CHEMBL3135,Q11111,"","",""
    CHEMBL2205,Q22222,"ENSP0004","None","None"
    ```

    Resulting `isoform.output.target_20250228.csv`:

    ```csv
    id,uniprot_id_primary,target_chembl_id,name
    ,P12345,CHEMBL3587,beta
    ENSP0001,P12345,CHEMBL3587,alpha
    ENSP0001,P12345,CHEMBL3587,pde3a
    ENSP0001,P12345,CHEMBL3587,3a
    ENSP0002,P12345,CHEMBL3587,beta
    ENSP0002,P12345,CHEMBL3587,alpha variant
    ,P67890,CHEMBL1250,pld1a
    ,P67890,CHEMBL1250,1a
    ENSP0003,P67890,CHEMBL1250,gamma
    ENSP0003,P67890,CHEMBL1250,gamma isoform 1
    ```

    ### Determinism
    - Sorting uses the stable `mergesort` algorithm, ensuring identical ordering
      on repeated runs.
    - The transformation is pure: no randomness, timestamps or locale-dependent
      operations are involved.
    - Output encoding is UTF-8 without BOM and uses Unix (`
`) line endings.

    ### Compatibility
    - `library.postprocessing.target.process_targets` is available starting from
      `chembl-data-acquisition` 0.1.3.
    - Python 3.11+ and `pandas` 2.3.x are required.
    - The input CSV must follow the targets schema described in
      [`docs/en/OUTPUT.md`](./en/OUTPUT.md).

    ### Usage example

    ```bash
    python - <<'PY'
    from library.postprocessing import target
    target.process_targets("data/output/output.target_20250228.csv")
    PY
    ```

    The command produces `data/output/isoform.output.target_20250228.csv` with
    isoform names aligned to the Power Query logic.

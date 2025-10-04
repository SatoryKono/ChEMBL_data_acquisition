# Target pipeline output checklist

This note summarises the post-processing logic, file formats and runtime
behaviour of the target acquisition pipeline. It complements the generic output
specification by focusing on the organism-specific artefacts produced by
`scripts/get_target_data.py`.

## Overview

The pipeline first normalises the merged ChEMBL/UniProt/IUPHAR payload via
`postprocess_targets`, then applies a final clean-up stage that deduplicates
rows, derives lineage-based labels and projects the table onto the canonical
schema before writing the CSV.【F:library/pipelines/target/postprocessing.py†L338-L485】【F:library/pipelines/target/postprocessing.py†L535-L730】 The resulting
artefact is saved alongside metadata and quality sidecars using the
`write_csv_deterministic` helper to guarantee stable ordering.【F:library/pipelines/target/postprocessing.py†L523-L532】【F:library/pipelines/target/postprocessing.py†L723-L730】

## Expected input

`postprocess_targets` expects the consolidated target export emitted by the
`get_target_data all` mode. At minimum the frame must include ChEMBL identifiers,
UniProt accessions and the lineage columns used for cellularity inference, plus
optional dictionary enrichments such as `multifunctional_enzyme`. The shared
resource `targets_type.csv` illustrates the required fields and their taxonomy
content.【F:library/integration/input_initialisation_library.py†L520-L576】【F:config/dictionary/_target/targets_type.csv†L1-L12】 A trimmed example (first six rows) is shown below:

```
target_chembl_id,uniprotkb_Id,organism_type,multifunctional_enzyme
CHEMBL1827,PDE5A,Multicellular organism,FALSE
CHEMBL1859,CACNA1H,Multicellular organism,FALSE
CHEMBL202,P00374,Multicellular organism,TRUE
CHEMBL1809,P0ABQ4,Unicellular organism,TRUE
CHEMBL1862,P00519,Multicellular organism,FALSE
CHEMBL203,P00533,Multicellular organism,FALSE
```

## Output format

The final CSV contains all fields from `TARGETS_COLUMN_ORDER`, including
identifier aliases, protein descriptors, cross references and pipeline
metadata.【F:library/schemas/targets.py†L15-L98】 Columns missing in the upstream
payload are populated with `-`, while identifier slots prefer empty strings over
the placeholder.【F:library/pipelines/target/postprocessing.py†L310-L335】 A
condensed sample of the export (subset of columns for the same rows) looks as
follows:

```
target_chembl_id,uniprot_id_primary,organism,target_type
CHEMBL1827,PDE5A,Homo,Multicellular organism
CHEMBL1859,CACNA1H,Homo,Multicellular organism
CHEMBL202,P00374,Homo,Multicellular organism
CHEMBL1809,P0ABQ4,Escherichia,Unicellular organism
CHEMBL1862,P00519,Homo,Multicellular organism
CHEMBL203,P00533,Homo,Multicellular organism
```

The `target_type` column is inferred from the lineage columns; values fall into
`Multicellular organism`, `Unicellular organism`, `Viral organism` or
`Unknown` depending on the rules baked into the classifier.【F:library/pipelines/target/postprocessing.py†L653-L668】【F:library/pipelines/target/organism_classification.py†L180-L270】

## Post-processing steps

1. **Identifier harmonisation** – UniProt accessions are normalised, fallback
   IDs are applied and gene symbols are reconstructed from multiple sources.
   Synonyms and EC numbers are aggregated into pipe-delimited strings to avoid
   duplicates.【F:library/pipelines/target/postprocessing.py†L369-L452】
2. **Optional field padding** – frequently absent columns (isoforms, reaction
   paths, etc.) are filled with `-` to maintain schema compatibility.【F:library/pipelines/target/postprocessing.py†L455-L484】
3. **Dictionary enrichments** – `targets_type.csv` is joined to propagate IUPHAR
   hierarchy, genomic indexes and the `multifunctional_enzyme` flag. The flag is
   coerced to boolean and the lineage fields feed the first pass of cellularity
   inference used downstream.【F:library/integration/input_initialisation_library.py†L520-L593】
4. **Final normalisation** – `finalise_targets` drops rows without stable
   UniProt IDs (unless a fallback is available), removes duplicate ChEMBL IDs,
   lowercases text-heavy columns and reuses the lineage columns to compute the
   definitive `target_type` label.【F:library/pipelines/target/postprocessing.py†L596-L668】
5. **Schema alignment** – `align_target_columns` reorders every column,
   replaces residual blanks with `-` and lowercases the post-translational
   modification flags to match historical extracts.【F:library/pipelines/target/postprocessing.py†L186-L335】

## Offline replay

`python -m library.utils.cli_tools.pipeline_targets_main` replays the pipeline
without contacting external APIs. The command wires a cached fetcher that emits
Deterministic ChEMBL frames, re-validates them against `TargetsSchema` and
writes both the raw dump (optional) and the final CSV through the same alignment
routine as the production CLI.【F:library/utils/cli_tools/pipeline_targets_main.py†L467-L565】 Because the cached payload only carries
identifier scaffolding, enrichment columns may be filled with the `-` placeholder;
the mode is intended for regression testing and schema verification, not for
producing publication-grade exports.

## HTTP error handling

During live runs the `chembl` fetcher guards every batch against
`requests.RequestException` and schema coercion failures. Errors are logged with
the offending chunk size/timeout, converted into `PipelineError` instances and
propagate to the orchestrator so the run halts and a `_failure_cases.csv`
sidecar is generated. Raw snapshot writing also guards against filesystem
errors and reports them consistently.【F:scripts/get_target_data.py†L1620-L1818】

## Result invariants

- Each `target_chembl_id` appears at most once; duplicates are removed before
  export.【F:library/pipelines/target/postprocessing.py†L624-L632】
- Placeholder handling is deterministic: identifier fields prefer empty strings,
  every other missing value becomes `-`, and selected columns are lowercased to
  ensure diff-friendly comparisons.【F:library/pipelines/target/postprocessing.py†L186-L335】
- The CSV is sorted by the original order of identifiers and written through the
  deterministic writer, so reruns with the same inputs produce byte-identical
  files and matching hashes.【F:library/pipelines/target/postprocessing.py†L523-L532】【F:library/pipelines/target/postprocessing.py†L723-L730】
- Validation against `TargetsSchema` guarantees that required columns are
  present and typed correctly even when streaming chunks in offline mode.【F:library/utils/cli_tools/pipeline_targets_main.py†L369-L443】

## Related references

- Generic output specification: [`docs/en/OUTPUT.md`](./en/OUTPUT.md)
- Russian counterpart: [`docs/OUTPUT_TARGETS_RU.md`](./OUTPUT_TARGETS_RU.md)
- CLI usage guide: [`docs/en/guides/USAGE.md`](./en/guides/USAGE.md)

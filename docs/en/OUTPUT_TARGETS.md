# Target pipeline output: organism, isoform, and names post-processing

This addendum complements the [Output specification](./OUTPUT.md) and
describes the auxiliary artefacts created after `scripts/get_target_data.py`
finalises the canonical `targets_*.csv` export. A Russian translation is
available in [`OUTPUT_TARGETS.md`](../ru/OUTPUT_TARGETS.md).

## Artefact overview

Running the target pipeline with the default writer produces four deterministic
CSV artefacts under the chosen output directory:

- `targets_<YYYYMMDD>.csv` — the main table described in
  [`OUTPUT.md`](./OUTPUT.md#target-export-targets).
- `organism.output.target_<YYYYMMDD>.csv` — organism level helpers used by the
  activity pipeline and QA checks (detailed below).
- `isoform.output.target_<YYYYMMDD>.csv` — per-isoform token expansion used for
  synonym coverage analysis.
- `names.output.target_<YYYYMMDD>.csv` — flattened component/name catalogue used
  to reconcile legacy exports and ChEMBL component groupings.
- `IUPHAR.output.target_<YYYYMMDD>.csv` — normalised IUPHAR classifications and
  synonym lists rebuilt from the canonical target export.

The `<YYYYMMDD>` suffix comes from the CLI `--date` argument or the deterministic
default configured via `local.io.default_date_prefix` (overridable with
`CHEMBL_DA_DEFAULT_DATE_PREFIX`).

## Organism enrichment post-processing

### Processing steps

1. **Read the final targets table** — `finalise_file` loads the CSV emitted by
   `postprocess_file`, applying the configured separator and encoding fallbacks.
2. **Normalise taxonomy fields** — missing lineage columns are reconstructed
   from legacy names (e.g. `superkingdom` → `lineage_superkingdom`) and text
   columns are coerced to string dtype.
3. **Filter invalid identifiers** — rows lacking UniProt accessions are dropped,
   duplicates on `target_chembl_id` are removed, and existing `type` columns are
   renamed to avoid clashes with the derived output.
4. **Infer cellularity** — `organism_classification.add_cellularity_smart`
   consumes genus and lineage information to populate `target_type` with the
   canonical labels `Multicellular organism`, `Unicellular organism`, or
   `Viral`.
5. **Project helper columns** — the frame is reduced to the lookup schema
   (`target_chembl_id`, `target_type`, `unicellular_organism`,
   `multifunctional_enzyme`, `IUPHAR_class`, `IUPHAR_subclass`, `sortorder.target`,
   `gene_index`, `taxon_index`). Boolean columns are normalised to the Pandas
   nullable `boolean` dtype.
6. **Emit `organism.output.*`** — the result is sorted by `target_chembl_id`,
   written with Unix line endings via `write_csv_deterministic`, and stored next
   to the main export (or under `--final-out` when provided).

The generated lookup replaces the static `config/dictionary/_target/targets_type.csv`
bundle and is consumed by `library.integration.input_initialisation_library`
when activities are enriched with species metadata.

### Input/output format

Minimal input columns used by the enrichment stage:

```
target_chembl_id,genus,lineage_superkingdom,lineage_phylum,lineage_class,IUPHAR_class,IUPHAR_subclass,multifunctional_enzyme
CHEMBL1824,Homo,Metazoa,Chordata,Mammalia,Voltage-gated ion channels,Sodium channels,false
CHEMBL240,Homo,Metazoa,Chordata,Mammalia,Transforming growth factor beta receptors,Type I receptor,false
CHEMBL6130,Candida,Fungi,Ascomycota,Saccharomycetes,,,
CHEMBL259,Homo,Metazoa,Chordata,Mammalia,Fibroblast growth factor receptors,Type 4 receptor,true
CHEMBL1111,Homo,Metazoa,Chordata,Mammalia,Hydrolases,Serine peptidases,false
CHEMBL1922,Influenza A virus,Viruses,Negarnaviricota,Insthoviricetes,,,false
CHEMBL6138,Escherichia,Bacteria,Proteobacteria,Gammaproteobacteria,Membrane proteins,Pore-forming,false
```

Corresponding lookup (first rows):

```
target_chembl_id,target_type,unicellular_organism,multifunctional_enzyme,IUPHAR_class,IUPHAR_subclass,sortorder.target,gene_index,taxon_index
CHEMBL1824,Multicellular organism,false,false,Voltage-gated ion channels,Sodium channels,0000012345,GENE000123,TAX000987
CHEMBL240,Multicellular organism,false,false,Transforming growth factor beta receptors,Type I receptor,0000012388,GENE000981,TAX000654
CHEMBL259,Multicellular organism,false,true,Fibroblast growth factor receptors,Type 4 receptor,0000012401,GENE000456,TAX000456
CHEMBL6130,Unicellular organism,true,false,-,-,0000012750,-,-
CHEMBL6138,Unicellular organism,true,false,Membrane proteins,Pore-forming,0000012761,-,-
CHEMBL1111,Multicellular organism,false,false,Hydrolases,Serine peptidases,0000012870,GENE000321,TAX000321
CHEMBL1922,Viral,true,false,-,-,0000012999,-,-
```

Missing dictionary attributes are rendered as `-` to match the main export.

### Offline mode

The enrichment is a pure local transformation. When the orchestrator runs with
`--offline` (or when upstream downloads are skipped), the stage expects the
canonical targets CSV and the optional dictionary overrides to exist locally. If
the lookup cannot be produced, the CLI aborts with a `FileNotFoundError` that
lists the expected locations for the dictionary bundle.

### HTTP failures

HTTP interactions occur before post-processing while `ChemblClient` streams raw
payloads. Network errors bubble up as `PipelineError` instances after being
logged under the `chembl_fetch_failed` key with the chunk size and timeout
metadata, ensuring the orchestrator fails fast instead of yielding partial
artefacts.

### Result invariants

- `target_type` is always one of `Multicellular organism`, `Unicellular organism`,
  or `Viral`.
- `unicellular_organism` equals `true` when `target_type` is `Unicellular
  organism` or `Viral`, and `false` otherwise.
- `multifunctional_enzyme` is a nullable boolean sourced from the lookup. Empty
  cells are serialised as `false`.
- Rows are sorted by `target_chembl_id` using a stable strategy. Re-running the
  pipeline on identical inputs yields byte-identical CSV files.
- All strings use UTF-8 with `\n` line endings.

## Isoform token expansion

`isoform.output.target_<stamp>.csv` reproduces the Power Query workbook that was
previously used to sanity-check isoform coverage. The pipeline performs the
following deterministic stages:

1. **Project required columns** — only `isoform_synonyms`, `isoform_names`,
   `isoform_ids`, `uniprot_id_primary`, and `target_chembl_id` are loaded.
2. **Normalise case** — `isoform_synonyms` and `isoform_names` are converted to
   lowercase text. Isoform identifiers keep their original case.
3. **Split pipe-separated lists** — the three isoform columns are split on `|`,
   whitespace is trimmed, and empty tokens are discarded.
4. **Align values by index** — `MakeTriples` pads the name/id/synonym lists with
   `null` so that each position forms a `{name, id, synonym}` record.
5. **Tokenise synonyms** — every synonym is split on `:` and each token is
   expanded into the variants `[token, token without "pde", token without
   "pld"]`, with duplicates removed while preserving order.
6. **Build name and synonym tables** — trimmed isoform names feed the first
   table (filtering out `""`, `"n/a"`, and `"none"`), while the expanded
   synonym tokens form the second table.
7. **Union and deduplicate** — the two tables are concatenated and then deduped
   in three passes: (a) on
   (`id`, `name`, `target_chembl_id`, `uniprot_id_primary`), (b) sorted with
   `mergesort` by (`uniprot_id_primary`, `id`) and deduped on (`id`,
   `target_chembl_id`, `name`), and (c) a final dedupe on (`id`, `name`). The
   stable sort guarantees deterministic survivors.
8. **Emit artefact** — the resulting columns are ordered as
   `["id", "uniprot_id_primary", "target_chembl_id", "name"]` and written as a
   UTF-8 CSV without BOM.

### Running the helper

Invoke the helper directly to regenerate the isoform table from any canonical
target export:

```bash
python - <<'PY'
from pathlib import Path
from library.postprocessing.target import process_targets

latest = Path("data/output/output.target_20250101.csv")
process_targets(latest, verbose=True)
PY
```

The script discovers the output path automatically when `process_targets` is
called without arguments, mirroring the behaviour of the original Power Query
workbook.

### Determinism and compatibility

- The Python implementation mirrors the Power Query (M) workbook used by the
  data team. Each transformation stage (projection, normalisation, tokenisation,
  union, and deduplication) has regression tests that assert byte-identical
  results.
- Stable `mergesort` ordering is enforced before the second deduplication pass,
  guaranteeing deterministic survivors when identifiers collide.
- Outputs are encoded as UTF-8 with Unix newlines and preserve the canonical
  column order `id`, `uniprot_id_primary`, `target_chembl_id`, `name`.
- The helper is compatible with Python 3.11+ and pandas ≥ 2.0, matching the
  baseline versions used across the pipeline.

Example input snapshot (from `output.target_20250101.csv`):

```
target_chembl_id,uniprot_id_primary,isoform_ids,isoform_names,isoform_synonyms
CHEMBL1,Q11111,ENSP0001|ENSP0002,Alpha|Beta,Alpha Alt|PDE3A:Alpha
CHEMBL2,Q22222,,Gamma|N/A|none,Gamma:Variant|n/a|none
CHEMBL3,Q33333,ID_UP|id_low,Theta|Lambda,PLDA:Variant
```

Corresponding `isoform.output.target_20250101.csv` (first rows):

```
id,uniprot_id_primary,target_chembl_id,name
ENSP0001,Q11111,CHEMBL1,alpha
ENSP0001,Q11111,CHEMBL1,alpha alt
ENSP0002,Q11111,CHEMBL1,beta
ENSP0002,Q11111,CHEMBL1,pde3a
ENSP0002,Q11111,CHEMBL1,3a
,Q22222,CHEMBL2,gamma
,Q22222,CHEMBL2,variant
ID_UP,Q33333,CHEMBL3,theta
ID_UP,Q33333,CHEMBL3,plda
ID_UP,Q33333,CHEMBL3,a
```

The helper is idempotent, locale independent, and produces byte-identical
outputs when rerun on the same input CSV.

## Post-processing `names.*`

`names.output.target_<stamp>.csv` materialises a compact view of the component
and nomenclature fields that ChEMBL historically exposed via multiple nested
tables. The helper keeps the nomenclature feed stable for downstream matching
against assay and activity exports while eliminating redundant columns.

### Input and output

- **Input:** the canonical `targets_<stamp>.csv` produced by
  `scripts/get_target_data.py` (either from `--output-dir` or from
  `--final-out`). Optional dictionary overrides referenced in the configuration
  are honoured when present to keep sort orders deterministic.
- **Output:** `names.output.target_<stamp>.csv` written next to the main export
  (or under `--final-out`). The file is encoded as UTF-8 with Unix line endings
  and sorted by `target_chembl_id`, `active_component_type`, `active_component`
  (descending), `component_id`, and the normalised `component_name`.

### Schema adjustments

The transformation keeps just the columns required by the consumer lookups:

```
target_chembl_id,component_id,component_name,component_synonyms,
component_accession,active_component,active_component_type
```

- Columns that only served for debugging in the canonical export—such as
  `component_description`, `component_synonym_ids`, `component_type_raw`,
  `component_sequence`, and the intermediate `_source_*` markers—are dropped.
- Empty strings are normalised to `-` for text columns, while boolean
  `active_component` values are cast to the pandas nullable `boolean` dtype.
- Null `component_synonyms` become empty pipe-delimited strings to keep the file
  compatible with legacy Power Query automations.

### Derived fields

`active_component_type` blends the raw component type with the activation flag
so that downstream joins do not need to re-open the main targets export:

```
if active_component is True:
    active_component_type = (component_type or "unknown").lower()
elif active_component is False:
    active_component_type = "inactive"
else:
    active_component_type = "unassigned"
```

The resulting value is stored in lowercase ASCII and is never left empty.

### Logging summary

The CLI logger records a single `target_names_postprocess` event with
aggregated counters:

- `input_rows` — rows read from the canonical export.
- `dropped_components` — components removed due to missing IDs or names.
- `null_synonyms` — rows where the synonym list had to be defaulted.
- `written_rows` — final row count after deduplication.

These metrics surface in the structured logs and the pipeline metadata YAML.

### Determinism

- Sorting uses `mergesort` for stability and a fixed tie-breaker order as noted
  above.
- The helper never touches the source file in place: it copies, reshapes, and
  writes to a brand new file, ensuring that repeated executions over the same
  input yield byte-identical `names.output.target_<stamp>.csv` artefacts.
- The output schema and default values are covered by regression tests to guard
  against accidental column drift.

## Post-processing `IUPHAR.*`

`IUPHAR.output.target_<stamp>.csv` reproduces the legacy helper that consolidated
Guide to PHARMACOLOGY chains, families, and synonym tokens. It is regenerated on
every target run to keep the helper aligned with the canonical export.

### Input/output

- **Input:** canonical `output.target_<stamp>.csv` with the IUPHAR enrichment
  columns produced by `scripts/get_target_data.py`.
- **Output:** helper written as `IUPHAR.output.target_<stamp>.csv` in the same
  directory (or to a custom path supplied via the Python API). When no input is
  provided, the helper discovers the most recent `output.target_*.csv` under the
  default `data/output` location.
- **Schema:**

  ```
  target_chembl_id,guidetopharmacology_id,iuphar_target_id,iuphar_family_id,
  iuphar_type,iuphar_class,iuphar_subclass,iuphar_chain,iuphar_name,iuphar_synonyms
  ```

  The column order is enforced when writing the CSV to match downstream
  consumers.

### Dropped columns and renames

Only the fields required by the synonym catalogue are kept. The helper removes
diagnostic columns inherited from the component processing stage:

- `component_synonym_ids`
- `component_type_raw`
- `component_sequence`
- `component_structures`

`GuidetoPHARMACOLOGY` is renamed to `guidetopharmacology_id`, while all other
columns are preserved as-is.

### Synonym construction

Synonyms aggregate multiple sources and are rewritten deterministically:

1. `gtop_synonyms`, `synonyms`, `component_description`, and the canonical
   `iuphar_name` are normalised via `normalise_text`, split on `|` (or JSON keys
   for structured descriptions), and stripped of empty tokens.
2. Parenthetical or bracketed fragments are removed, consecutive whitespace is
   collapsed, and tokens are lower-cased.
3. Tokens are deduplicated in first-seen order to preserve deterministic
   ordering.
4. The resulting list joins into a pipe-delimited `iuphar_synonyms` column.

The helper records the number of tokens observed before and after deduplication
to simplify QA checks.

### Logging and determinism

When called with `verbose=True`, the helper emits structured lines with the
input/output paths, counts of processed rows, the number of dropped columns, and
synonym statistics. Outputs are written through `write_csv_deterministic` with a
stable sort key (`target_chembl_id`), ensuring byte-identical results for
identical inputs and Unix line endings.

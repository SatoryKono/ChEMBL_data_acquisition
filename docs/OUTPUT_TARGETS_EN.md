# Target pipeline output: organism and isoform post-processing

This addendum complements the [Output specification](./en/OUTPUT.md) and
describes the auxiliary artefacts created after `scripts/get_target_data.py`
finalises the canonical `targets_*.csv` export. A Russian translation is
available in [`OUTPUT_TARGETS_RU.md`](./OUTPUT_TARGETS_RU.md).

## Artefact overview

Running the target pipeline with the default writer produces three deterministic
CSV artefacts under the chosen output directory:

- `targets_<YYYYMMDD>.csv` — the main table described in
  [`OUTPUT.md`](./en/OUTPUT.md#target-export-targets).
- `organism.output.target_<YYYYMMDD>.csv` — organism level helpers used by the
  activity pipeline and QA checks (detailed below).
- `isoform.output.target_<YYYYMMDD>.csv` — per-isoform token expansion used for
  synonym coverage analysis.

The `<YYYYMMDD>` suffix is inherited from the CLI `--date` argument or the
automatic stamp derived from the execution date.

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
   `multifunctional_enzyme`, `IUPHAR_class`, `IUPHAR_subclass`, `target_sort_order`,
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
target_chembl_id,target_type,unicellular_organism,multifunctional_enzyme,IUPHAR_class,IUPHAR_subclass,target_sort_order,gene_index,taxon_index
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

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
- `isoform.output.targets_<YYYYMMDD>.csv` — per-isoform token expansion used for
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

The isoform helper retains the semantics previously documented for
`isoform.output.*`. The stage reads the canonical targets CSV, cleans the
`isoform_*` columns, expands synonyms into tokens, removes duplicates on
(`target_chembl_id`, `isoform_id`, `term`, `token`), and emits a deterministic
CSV for downstream QA tooling.

The input snapshot and invariants remain unchanged from earlier revisions. An
end-to-end example is included below for completeness.

```
target_chembl_id,isoform_ids,isoform_names,isoform_synonyms
CHEMBL1824,ENSP00000350283,Nav1.7,"Nav1.7:SCN9A isoform 3"
CHEMBL1824,ENSP00000350284,"Nav1.7 isoform 2","Nav1.7 splice variant:SCN9A-2"
CHEMBL6130,,"",""
CHEMBL240,ENSP00000456012,ALK2,"Activin receptor-like kinase 2:ACVR1"
CHEMBL259,ENSP00000263253,FGFR4,"FGFR-4 isoform:Fibroblast growth factor receptor-4"
```

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

The transformation is deterministic, locale independent, and requires no
network connectivity.

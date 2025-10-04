# Target post-processing — `organism.output.target_*.csv`

The target pipeline now emits an additional organism-focused CSV derived from
the ChEMBL/UniProt export (`output.target_*.csv`). The Power Query (M) logic
previously used for taxonomy dashboards has been ported to Python without
changing the resulting rows or columns.

## Transformation steps

1. **Input detection** — the post-processor consumes the freshly written
   `output.target_*.csv` file (or the most recent file matching the pattern
   when invoked standalone).
2. **Cellularity classification** — superkingdom/phylum fields are normalised
   to lower case. Each row is labelled using the taxonomy rules:
   - `ClassifyByLineage` covers deterministic cases for viruses, bacteria,
     archaea and curated eukaryotic phyla.
   - Ambiguous or empty lineage entries trigger `ClassifyByFetch`, which reads
     NCBI Taxonomy XML (`Lineage`, `ScientificName`) and applies the same
     rulebook.
3. **Multifunctional enzyme flag** — `reaction_ec_numbers` is split on `|`,
   deduplicated, and truncated to the first EC segment. Targets with more than
   one distinct EC prefix are marked as multifunctional.
4. **Join and projection** — cellularity and multifunctional flags are merged
   back into the identifier subset, producing the final column order:
   `target_chembl_id`, `uniprot_id_primary`, `organism`, `taxon_id`,
   `lineage_superkingdom`, `lineage_phylum`, `lineage_class`, `cellularity`,
   `multifunctional_enzyme`.
5. **Output** — results are saved next to the input file with the prefix
   `organism.` while preserving UTF-8 encoding and deterministic ordering.

## Example

Input (`output.target_20250301.csv`):

```csv
target_chembl_id,uniprot_id_primary,organism,taxon_id,lineage_superkingdom,lineage_phylum,lineage_class,reaction_ec_numbers
CHEMBL1,P12345,Species A,10239,Viruses,,,"1.1.1.1|2.7.11.1"
CHEMBL2,Q54321,Species B,9606,Eukaryota,Chordata,Chordata,"1.1.1.1|3.2.1.4"
CHEMBL3,R99999,Species C,111,Archaea,,,""
CHEMBL4,S88888,Species D,222,,,"","1.1.1.1|2.7.11.1|1.2.3.4"
CHEMBL5,T77777,Species E,333,Bacteria,,,"2.7.11.1"
```

Output (`organism.output.target_20250301.csv`):

```csv
target_chembl_id,uniprot_id_primary,organism,taxon_id,lineage_superkingdom,lineage_phylum,lineage_class,cellularity,multifunctional_enzyme
CHEMBL1,P12345,Species A,10239,viruses,,,"acellular (virus)",True
CHEMBL2,Q54321,Species B,9606,eukaryota,chordata,chordata,multicellular,True
CHEMBL3,R99999,Species C,111,archaea,,,"unicellular",False
CHEMBL4,S88888,Species D,222,,,ambiguous,True
CHEMBL5,T77777,Species E,333,bacteria,,,"unicellular",False
```

The organism CSV preserves the row order of the input file. Lists produced
from `reaction_ec_numbers` remain internal to the computation and are not
exposed in the final export.

## Operational notes

- **Offline mode** — pass `offline=True` to
  `library.postprocessing.target.process_targets` to disable HTTP calls. Rows
  that would otherwise rely on NCBI fall back to `cellularity = "ambiguous"`.
- **HTTP errors** — network failures or invalid XML payloads mirror the M
  implementation and result in an ambiguous classification without raising.
- **Logging** — the processor reports the input/output paths, number of rows,
  HTTP requests made to NCBI and the count of ambiguous classifications when
  `verbose=True`.

Refer to the Russian counterpart (`docs/OUTPUT_TARGETS_RU.md`) for the same
information in Russian.

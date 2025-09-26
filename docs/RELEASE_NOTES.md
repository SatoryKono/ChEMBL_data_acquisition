# Release Notes

## Unreleased

- Target taxonomy no longer requires the legacy `dictionary/_Target/organism.csv` lookup. Remove any `organism_csv` overrides from
  custom configurations; the pipeline now derives organism group and lineage data from UniProt metadata automatically.

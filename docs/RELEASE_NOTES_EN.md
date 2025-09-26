# Release Notes

## Unreleased

- **Migration:** Target taxonomy classification is now derived directly from UniProt lineage data. The legacy `dictionary/_Target/organism.csv` lookup and the `organism_csv` configuration key were removed; drop any overrides pointing to that file when upgrading.

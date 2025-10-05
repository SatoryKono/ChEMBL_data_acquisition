# Test resources

## `tissue_pipeline/`

Recorded fixtures for the tissue pipeline integration tests.

- `tissue_ids.csv` — minimal input identifiers used to seed the pipeline run with a mixture of present and missing tissues.
- `chembl_tissue_page1.json` — first page of a captured ChEMBL API response containing tissues with null and empty string fields.
- `chembl_tissue_page2.json` — second page continuing the response with additional tissues and terminating the pagination chain.
- `integration_payload.json` — representative payload used for integration tests covering deterministic CSV ordering and metadata enrichment.
- `integration_expected_output.csv` — expected CSV produced by the integration scenario, including pipeline metadata columns.

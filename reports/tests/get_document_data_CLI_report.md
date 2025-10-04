# Document CLI Test Report

## Summary
- ✅ `pytest tests/unit/test_cli_get_document_data.py`
- ✅ `pytest tests/integration/test_fallback_doi.py`

## CLI Flag Matrix

| Scenario | Flags | Expected Outcome |
|----------|-------|------------------|
| PubMed defaults | `get-document-data pubmed --input <file>` | Parser applies default batch size, offsets, and disables DOI fallback. |
| PubMed CLI override | `get-document-data pubmed --batch-size 77 --final-out <out>` | CLI flag overrides configuration and default values for batch sizing. |
| DOI fallback CSV | `get-document-data pubmed --fallback-doi-csv fallback.csv` | PMID→DOI overrides are loaded and handed to the fetch pipeline, logging fallback metrics. |
| Skip existing output | `get-document-data pubmed --skip-existing [--force] --final-out <out>` | Without `--force` the run short-circuits with `pipeline_skip_existing`; with `--force` it reprocesses the input. |

## Notes
- Unit tests cover parser defaults, CLI/config/default precedence, and boundary validation for `--limit`/`--offset`.
- Integration tests stub network fetches, replay deterministic CSV fixtures, and assert structured log events for fallback DOI handling.

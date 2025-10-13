# Quick start

Follow these steps to run the full pipeline locally using the bundled sample
inputs.

## 1. Prepare the environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-lock.txt
pre-commit install
```

Optional extras:

- `pip install -r requirements-dev.txt` to run linting and type checks.
- `export PYTHONPATH=$(pwd)` when invoking modules directly.

## 2. Verify configuration

The default configuration lives in `config/config.yaml`. Create
`config/config.local.yaml` for overrides (e.g. custom API keys):

```yaml
sources:
  chembl:
    api:
      user_agent: "ChEMBL-ETL/2.1 (mailto:your-team@example.org)"
```

Confirm the final configuration:

```bash
python scripts/get_document_data.py --mode chembl --print-config | less
```

## 3. Run smoke pipelines

The `data/input` directory contains minimal CSVs covering each entity. Execute
pipelines individually or via the orchestrator:

```bash
# Run a single pipeline
python scripts/get_document_data.py --mode all \
  --input data/input/document.csv \
  --final-out output/documents.csv

# Run the full chain
poetry run get-data \
  --base-path . \
  --input-dir data/input \
  --output-dir output \
  --config config/config.yaml \
  --date $(date -u +%Y%m%d)
```

Outputs are written to `output/` (create the directory beforehand or rely on
`exist_ok=true`). Inspect `<name>.meta.yaml` for provenance.

## 4. Execute tests

```bash
pytest -q --disable-warnings
pytest -q --json-report --json-report-file=reports/test_report.json
python tools/make_md_summary.py --input reports/test_report.json --output reports/test_summary.md
# arguments are optional when using the default locations:
# python tools/make_md_summary.py
# make-md-summary
```

Expect `reports/test_summary.md` to report ≥95 % success rate. Upload the JSON and
Markdown artefacts to CI when filing issues.

## 5. Optional: determinism check

To confirm reproducibility, run the full pipeline twice and compare artefacts:

```bash
poetry run get-data --output-dir output/run1
poetry run get-data --output-dir output/run2
check-determinism --baseline output/run1 --candidate output/run2
```

Hash mismatches indicate non-deterministic behaviour or configuration drift.
Consult [`QA_PROCESS.md`](../QA_PROCESS.md) for troubleshooting hints.

# ChEMBL data acquisition

Utilities for collecting target, assay, and document information from the
[ChEMBL](https://www.ebi.ac.uk/chembl/) API and related services.  The
repository provides reusable library functions along with small command line
wrappers for bulk data retrieval.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

### Assay data

```bash
python csv/python/get_assay_data.py assays.csv output.csv
```

The input CSV must contain a column named `assay_chembl_id` with ChEMBL assay
identifiers.  Additional options such as delimiter or encoding can be adjusted
via command line flags (`--sep`, `--encoding`).

### Document data

```bash
python get_document_data.py documents.csv output.csv
```

By default the script reads the `document_chembl_id` column from the input
file.  Use `--column` to override the field name.

## Development

Run the unit tests with:

```bash
pytest
```

Recommended quality tools:

```bash
black .
ruff .
mypy .
```

## License

This project is provided under the MIT license.  See `LICENSE` for details.

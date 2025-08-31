# ChEMBL_data_acquisition

## Review classification

This repository now includes a small, reproducible pipeline for classifying
publications as **review**, **non_review**, or **uncertain**.

### Installation

```bash
pip install -r requirements.txt
```

### Usage

```bash
python main.py --input tests/data/sample_records.json --output output
# explicit encoding can be supplied if needed
python main.py --input non_utf8.json --output output --encoding cp1252
```

The same CLI is also available as ``get_document_classification.py`` for
compatibility with existing scripts.

Artifacts are written to the specified output directory using the structure
outlined in the project documentation.

### Quality tools

The following tools are recommended:

```bash
black .
ruff .
mypy .
pytest
```

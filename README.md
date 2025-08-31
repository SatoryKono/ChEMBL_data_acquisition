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
```

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

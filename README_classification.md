# Document Classification

This tool classifies documents as `review`, `non_review`, or `uncertain`
using publication type (PT) and MeSH annotations from multiple sources.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Quality tools

```bash
ruff .
black .
mypy .
pytest
```

## Usage

```bash
python main.py --input INPUT.csv --output-dir output
```

Encoding is auto-detected, but you may specify it explicitly:

```bash
python main.py --input INPUT.csv --output-dir output --encoding cp1251
```

Input must contain the required columns. Outputs include normalized PT
and MeSH features, scores, decisions, logs, and `run_meta.json` with the
rules version.

Example:

```bash
python main.py --input tests/data/sample.csv --output-dir demo
```

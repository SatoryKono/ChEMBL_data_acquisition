# ChEMBL data acquisition

Utilities for working with ChEMBL related documents. This repository now
includes a document classifier that assigns records to `review`,
`non-review` or `unknown` categories based on PublicationType and MeSH
information from multiple sources.

## Installation

```bash
pip install -r requirements.txt
```

Optional development tools:

```bash
pip install black ruff mypy
```

## Usage

```bash
python classify.py --input /mnt/data/all_documents.csv --output out.csv --log logs.jsonl
```

### Command line options

- `--input` – path to the source CSV/TSV file.
- `--output` – destination CSV with classification results.
- `--log` – optional JSONL log file with detailed decision traces.
- `--chunk-size` – number of rows to process at once (default 1000).
- `--threshold` – decision margin threshold (default 1.0).
- `--log-level` – logging verbosity (default INFO).

## Testing

```bash
pytest
```

## Code quality

The project follows PEP 8 and uses type hints. Recommended tooling:

```bash
black .
ruff .
mypy .
```

## Algorithm Notes

The classifier normalises text to lower case, replaces known synonyms and
extracts signals indicating review or non-review nature. Each source and
signal type has configurable weights. Scores are aggregated and compared
with a configurable threshold. Details for each processed record can be
logged in JSONL format for full traceability.

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
python classify_documents.py \
  --input-all all_documents.csv \
  --input-mesh experimental_mesh.csv \
  --output result.csv
```

### Command line options

- `--input-all` – path to the `all_documents.csv` file.
- `--input-mesh` – path to `experimental_mesh.csv` with term probabilities.
- `--output` – destination path for results (CSV or Parquet).
- `--delta` – MeSH score difference required for a confident decision (default 0.5).
- `--k-min` – minimum number of MeSH terms for refinement (default 3).
- `--sep-list` – separators for list-like fields (default `|;,/`).
- `--unknown-mode` – output `unknown` when MeSH signal is weak.
- `--prefer-pubmed-epsilon` – optional extra weight when only PubMed votes review.
- `--chunksize` – read input in chunks of this size.
- `--output-format` – `csv` (default) or `parquet`.
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

1. Publication types from PubMed, OpenAlex and Google Scholar are
   normalised and mapped to canonical labels. If at least two sources
   mark the record as a review, the final label is ``review``; if none do,
   it is ``non-review``.
2. Records with exactly one review vote or missing data are refined using
   MeSH terms. Experimental probabilities for each term are summed to
   produce scores for experimental vs review articles. A difference above
   ``delta`` determines the label; otherwise the record defaults to
   ``non-review`` or ``unknown`` depending on the CLI flags.
3. The output contains provenance columns with individual source votes,
   MeSH terms used and top contributing terms.

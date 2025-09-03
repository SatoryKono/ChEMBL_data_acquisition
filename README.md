# ChEMBL data acquisition

Utilities for working with ChEMBL related documents. This repository now
includes a publication classifier that assigns records to `review`,
`experimental` or `unknown` categories based on publication type metadata from
multiple sources.

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
python classify_publications.py --input all_documents.csv --output out.csv
```

### Command line options

- `--input` – path to the source CSV file.
- `--output` – destination CSV with classification results.
- `--min-review-score` – minimum score for the `review` class (default 1).
- `--min-unknown-score` – minimum score for the `unknown` class (default 2).
- `--min-experimental-score` – minimum score for the `experimental` class (default 1).
- `--encoding-fallbacks` – comma separated list of encodings to try (default `utf-8-sig,cp1251,latin1`).
- `--delimiters` – delimiters to try when reading the input (default `,;|\t`).
- `--log-level` – logging verbosity (default `INFO`).

The script prints class distribution to stdout and logs diagnostics to stderr.

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

Publication type fields from PubMed, Google Scholar and OpenAlex are
normalised: text is lower‑cased, split on common delimiters and mapped to
canonical forms (e.g. “mini review” → “review”). Weighted signals from the three
sources are accumulated for `review`, `experimental` and `unknown` term sets. If
`review` has the highest score it is selected. Otherwise `unknown` is chosen
when its score is highest and ≥2, followed by `experimental` when its score is
highest and ≥1. Records without signals fall back to `unknown`.

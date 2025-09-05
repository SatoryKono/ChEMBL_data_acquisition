# ChEMBL Data Acquisition Utilities

Utilities for downloading and processing biological data from public APIs.
The project demonstrates a typical Python 3.12 data pipeline including
parsing, validation, aggregation and export of tabular data.

## Installation

```bash
pip install -r requirements.txt
```

## Command line interface

Individual scripts provide specialised data retrieval utilities:

* ``get_activity_data.py`` – fetch ChEMBL activity information.
* ``get_assay_data.py`` – retrieve assay descriptions from ChEMBL.
* ``get_document_data.py`` – gather publication metadata.
* ``get_target_data.py`` – combine ChEMBL, UniProt and IUPHAR target data.
* ``get_testitem_data.py`` – download compound data and enrich with PubChem.

All scripts share a common set of flags:

* ``--input`` – input CSV file (default ``input.csv``)
* ``--output`` – destination CSV file (default: auto-generated next to the input)
* ``--log-level`` – logging verbosity (default ``INFO``)
* ``--sep`` – CSV delimiter (default ``,``)
* ``--encoding`` – file encoding (default ``utf8``)
* ``--column`` – column containing identifiers (script specific)

Example fetching assay data::

    python get_assay_data.py --input assays.csv --output assays_out.csv \
        --column assay_chembl_id

Each command validates required columns before querying external APIs and
writes the resulting table to the specified output file.

## Development

Formatting, linting and type checking are handled by *black*, *ruff* and
*mypy* respectively:

```bash
black get_*.py library/io.py library/validation.py tests
ruff get_*.py library/io.py library/validation.py tests
mypy get_*.py library/io.py library/validation.py
pytest
```

Test data live in ``tests/data`` and provide coverage for utility
functions in the library modules.

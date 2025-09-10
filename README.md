# ChEMBL Data Acquisition Utilities

Utilities for downloading and processing biological data from public APIs.
The project demonstrates a typical Python 3.12 data pipeline including
parsing, validation, aggregation and export of tabular data.

## Installation

Clone the repository and install the runtime and development dependencies:

```bash
git clone https://example.com/ChEMBL_data_acquisition.git
cd ChEMBL_data_acquisition
pip install -r requirements.txt
```

## Project structure

```
data/             Example input and output files
dictionary/       Lookup tables used during processing
library/          Reusable data-processing modules
tests/            Pytest suite and sample datasets
get_*.py          Command-line utilities for specific tasks
mapper_main.py    Mapping CLI
table_quality_main.py  CSV profiling CLI
```

## Command line interface

Individual scripts provide specialised data retrieval utilities:

* ``get_activity_data.py`` – fetch ChEMBL activity information.
* ``get_assay_data.py`` – retrieve assay descriptions from ChEMBL.
* ``get_document_data.py`` – gather publication metadata.
* ``get_target_data.py`` – combine ChEMBL, UniProt and IUPHAR target data.
* ``get_testitem_data.py`` – download compound data and enrich with PubChem.
* ``get_input_initialisation.py`` – merge ChEMBL initialisation workbooks.

Detailed command line examples using the bundled smoke datasets can be found in
``docs/USAGE.md``.

### Table quality analysis

``table_quality_main.py`` profiles arbitrary CSV files and reports column
statistics along with correlations between numeric fields. Example usage:

```python
import pandas as pd
from library.table_quality import analyze_table_quality

df = pd.read_csv("data.csv", encoding="utf-8-sig")
quality, corr = analyze_table_quality(df, table_name="data")
```

Running the CLI saves ``data_quality_report_table.csv`` and
``data_data_correlation_report_table.csv`` in the current working directory::

    python table_quality_main.py --input data.csv --table-name data

All scripts share a common set of flags:

## Configuration

Default settings live in `config.yaml`. Values can be overridden via environment variables using the `CHEMBL_DA__SECTION__KEY` pattern or short aliases like `CHEMBL_DA_RPS`.

Example usage:

```bash
python get_activity_data.py --config config.yaml
CHEMBL_DA__API__RPS=2 python get_assay_data.py
```

| Script | Config sections |
| --- | --- |
| get_activity_data.py | api, io, jobs, quality, log |
| get_assay_data.py | api, io, jobs, quality, log |
| get_document_data.py | api, io, jobs, quality, log |
| get_document_type.py | io, log |
| get_target_data.py | api, uniprot, iuphar, io, jobs, quality, log |
| get_testitem_data.py | api, pubchem, io, jobs, quality, log |
| get_input_initialisation.py | init, io, quality, log |
| mapper_main.py | mapper, io, jobs, log |
| table_quality_main.py | io, quality, log |


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

Example merging initialisation tables::

    python get_input_initialisation.py --config config.yaml

The ``same_doc`` and ``all_doc`` workbook paths default to values from
``config.yaml`` but can be overridden on the command line::

    python get_input_initialisation.py \
      --same-doc path/to/ChEMBL_same_document_20_05.xlsx \
      --all-doc  path/to/ChEMBL_all_10_05_step5.xlsx \
      --out-dir  ./out

The script also profiles each exported table and writes
``<name>_quality_report_table.csv`` and
``<name>_data_correlation_report_table.csv`` alongside the original CSVs
in ``--out-dir``.

## Development

Formatting, linting and type checking are handled by *black*, *ruff* and
*mypy* respectively. Run the following before committing changes:

```bash
black get_*.py library mapper_main.py table_quality_main.py
ruff check get_*.py library mapper_main.py table_quality_main.py
mypy get_*.py library mapper_main.py table_quality_main.py
pytest
```

Test data live in ``tests/data`` and provide coverage for utility
functions in the library modules.

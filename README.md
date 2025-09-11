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
config.yaml       Global configuration defaults
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


Default settings live in ``config.yaml`` and are split into sections for each
API (``api``, ``openalex``, ``crossref``, ``uniprot``, ``iuphar``, ``pubchem``),
I/O and processing (``io``, ``jobs``, ``batch``, ``quality``, ``mapper``) and
general infrastructure (``init``, ``rate``, ``retry``, ``log``). A minimal
configuration looks like::


    api:
      rps: 5
    io:
      output_dir: data/output
    jobs:
      concurrency: 8

Environment variables override values from the YAML file. Variables use the
``CHEMBL_DA__SECTION__KEY`` pattern and also support short aliases:

* ``CHEMBL_DA__API__CHEMBL_BASE`` / ``CHEMBL_DA_BASE``
* ``CHEMBL_DA__API__TIMEOUT_CONNECT`` / ``CHEMBL_DA_TIMEOUT_CONNECT``
* ``CHEMBL_DA__API__TIMEOUT_READ`` / ``CHEMBL_DA_TIMEOUT_READ``
* ``CHEMBL_DA__API__RPS`` / ``CHEMBL_DA_RPS``

* ``CHEMBL_DA__OPENALEX__TIMEOUT_CONNECT`` / ``CHEMBL_DA_OPENALEX_TIMEOUT_CONNECT``
* ``CHEMBL_DA__OPENALEX__TIMEOUT_READ`` / ``CHEMBL_DA_OPENALEX_TIMEOUT_READ``
* ``CHEMBL_DA__OPENALEX__RPS`` / ``CHEMBL_DA_OPENALEX_RPS``
* ``CHEMBL_DA__CROSSREF__TIMEOUT_CONNECT`` / ``CHEMBL_DA_CROSSREF_TIMEOUT_CONNECT``
* ``CHEMBL_DA__CROSSREF__TIMEOUT_READ`` / ``CHEMBL_DA_CROSSREF_TIMEOUT_READ``
* ``CHEMBL_DA__CROSSREF__RPS`` / ``CHEMBL_DA_CROSSREF_RPS``
* ``CHEMBL_DA__UNIPROT__TIMEOUT_CONNECT`` / ``CHEMBL_DA_UNIPROT_TIMEOUT_CONNECT``
* ``CHEMBL_DA__UNIPROT__TIMEOUT_READ`` / ``CHEMBL_DA_UNIPROT_TIMEOUT_READ``
* ``CHEMBL_DA__UNIPROT__RPS`` / ``CHEMBL_DA_UNIPROT_RPS``
* ``CHEMBL_DA__IUPHAR__TIMEOUT_CONNECT`` / ``CHEMBL_DA_IUPHAR_TIMEOUT_CONNECT``
* ``CHEMBL_DA__IUPHAR__TIMEOUT_READ`` / ``CHEMBL_DA_IUPHAR_TIMEOUT_READ``
* ``CHEMBL_DA__IUPHAR__RPS`` / ``CHEMBL_DA_IUPHAR_RPS``
* ``CHEMBL_DA__PUBCHEM__TIMEOUT_CONNECT`` / ``CHEMBL_DA_PUBCHEM_TIMEOUT_CONNECT``
* ``CHEMBL_DA__PUBCHEM__TIMEOUT_READ`` / ``CHEMBL_DA_PUBCHEM_TIMEOUT_READ``
* ``CHEMBL_DA__PUBCHEM__RPS`` / ``CHEMBL_DA_PUBCHEM_RPS``

* ``CHEMBL_DA__IO__OUTPUT_DIR`` / ``CHEMBL_DA_OUTDIR``
* ``CHEMBL_DA__JOBS__CONCURRENCY`` / ``CHEMBL_DA_CONCURRENCY``
* ``CHEMBL_DA__JOBS__CHUNK_SIZE`` / ``CHEMBL_DA_CHUNK_SIZE``
* ``CHEMBL_DA__RETRY__MAX_ATTEMPTS`` / ``CHEMBL_DA_RETRY_MAX_ATTEMPTS``
* ``CHEMBL_DA__RETRY__BACKOFF_FACTOR`` / ``CHEMBL_DA_RETRY_BACKOFF_FACTOR``
* ``CHEMBL_DA__LOG__LEVEL`` / ``CHEMBL_DA_LOG_LEVEL``
* ``CHEMBL_DA__LOG__FORMAT`` / ``CHEMBL_DA_LOG_FORMAT``

### Schema validation

Configuration values are validated against a JSON Schema via the
``jsonschema`` package. The schema mirrors the dataclass structure and checks
types and value ranges, producing helpful error messages for nested fields.

Command line flags have the highest priority. All utilities accept ``--config``
to point at a configuration file and ``--print-config`` to show the effective
values after all overrides have been applied. The final precedence is::

    YAML < environment variables < CLI options

Only the top-level command line scripts read the configuration file. Modules
under ``library/`` expect a :class:`Config` (or one of its subsections) to be
passed explicitly, making dependencies clear and avoiding hidden global state.
The directories referenced by ``io.output_dir`` and ``io.cache_dir`` are checked
but not created when loading the configuration. Scripts that need these paths
can call :func:`library.config.ensure_dirs` after :func:`load_config` to create
them if they are missing and ``io.exist_ok`` permits it.

Path values such as ``io.output_dir``, ``io.cache_dir`` and the ``init``
workbook paths are exposed as :class:`pathlib.Path` objects. String values in
``config.yaml`` or overrides from the environment and command line are
automatically converted.


Common flags shared by scripts include:

* ``--input`` – input CSV file (default ``input.csv``)
* ``--output`` – destination CSV file (default: auto-generated next to the input)
* ``--log-level`` – logging verbosity (default ``INFO``)
* ``--sep`` – CSV delimiter (default taken from configuration)
* ``--encoding`` – file encoding (default taken from configuration)
* ``--column`` – column containing identifiers (script specific)

Example fetching assay data::

    python get_assay_data.py --input assays.csv --output assays_out.csv \
        --column assay_chembl_id

Each command validates required columns before querying external APIs and
writes the resulting table to the specified output file.

## Data contracts

Each output table is validated with ``pandera`` to guarantee a consistent
layout. Columns must satisfy the following contracts.

### Activities

Required columns

* ``activity_id`` (*int*, ``>= 0``)
* ``testitem_id`` (*str*)
* ``standard_value`` (*float*, ``>= 0``)

Optional columns

* ``target_id`` (*str*)
* ``standard_type`` (*str*, one of ``IC50``, ``EC50``, ``Ki``, ``Kd``)
* ``pA_value`` (*float*, ``0–14``)

Valid row

```csv
activity_id,testitem_id,target_id,standard_type,standard_value,pA_value
1,TST1,TGT1,IC50,50,9
```

Invalid row (``standard_type`` outside enum, ``pA_value`` > 14)

```csv
activity_id,testitem_id,target_id,standard_type,standard_value,pA_value
2,TST2,TGT2,IC90,100,20
```

### Assays

Required columns

* ``assay_chembl_id`` (*str*)
* ``document_chembl_id`` (*str*)
* ``year`` (*int*, ``1900–2100``)
* ``month`` (*int*, ``1–12``)

Optional columns

* ``target_chembl_id`` (*str*)

Valid row

```csv
assay_chembl_id,document_chembl_id,target_chembl_id,year,month
A1,D1,T1,2023,5
```

Invalid row (``month`` > 12)

```csv
assay_chembl_id,document_chembl_id,target_chembl_id,year,month
A2,D2,T2,2023,13
```

### Documents

Required columns

* ``document_chembl_id`` (*str*)
* ``title`` (*str*)
* ``year`` (*int*, ``1900–2100``)
* ``month`` (*int*, ``1–12``)

Optional columns

* ``doi`` (*str*)
* ``day`` (*int*, ``1–31``)
* ``citation`` (*int*, ``>= 0``)

Valid row

```csv
document_chembl_id,doi,title,year,month,day,citation
D1,10.1000/test,A study,2022,7,15,3
```

Invalid row (``day`` > 31, ``citation`` < 0)

```csv
document_chembl_id,doi,title,year,month,day,citation
D2,10.1000/test2,Another study,2022,7,45,-1
```

### Targets

Required columns

* ``target_chembl_id`` (*str*)
* ``organism`` (*str*)

Optional columns

* ``target_uniprot_id`` (*str*)
* ``pH_dependence`` (*float*, ``0–14``)
* ``isoforms`` (*float*, ``>= 0``)

Valid row

```csv
target_chembl_id,organism,target_uniprot_id,pH_dependence,isoforms
T1,Homo sapiens,P12345,7.4,2
```

Invalid row (``pH_dependence`` > 14)

```csv
target_chembl_id,organism,target_uniprot_id,pH_dependence,isoforms
T2,Mus musculus,P67890,15,1
```

### Testitems

Required columns

* ``salt_chembl_id`` (*str*)
* ``molecule_chembl_id`` (*str*)
* ``molecule_type`` (*str*, ``Small molecule``, ``Biopolymer``,
  ``Oligosaccharide``, ``Unknown``)
* ``mw_freebase`` (*float*, ``0–2000``)

Optional columns

* ``chirality`` (*int*, ``-1``, ``0``, ``1``, ``2``)
* ``num_ro5_violations`` (*float*, ``0–5``)
* ``is_radical`` (*bool*)

Valid row

```csv
salt_chembl_id,molecule_chembl_id,molecule_type,chirality,mw_freebase,num_ro5_violations,is_radical
S1,M1,Small molecule,1,350.5,0,false
```

Invalid row (``molecule_type`` outside enum, ``mw_freebase`` > 2000)

```csv
salt_chembl_id,molecule_chembl_id,molecule_type,chirality,mw_freebase,num_ro5_violations,is_radical
S2,M2,Peptide,0,2500,1,true
```

## Configuration

Default settings such as API endpoints, network timeouts, rate limits and
output directories live in `config.yaml` at the repository root. Each field is
documented inside the file and has a sensible fallback that the utilities use
if the entry is missing. See `docs/CONFIG.md` for a detailed description of all
available options.

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

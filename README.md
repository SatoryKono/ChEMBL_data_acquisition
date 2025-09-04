# ChEMBL_data_acquisition
Utilities for downloading and integrating target information from
[ChEMBL](https://www.ebi.ac.uk/chembl/),
[UniProt](https://www.uniprot.org/) and the
[IUPHAR Guide to Pharmacology](https://www.guidetopharmacology.org/).
It also provides helpers for collecting publication metadata from
PubMed, Semantic Scholar, OpenAlex and CrossRef.

Three command line tools are available:
`get_activity_data.py`
    Fetch bioactivity of small molecyle compounds  from the ChEMBL API for a list of activity IDs.
    
`get_assay_data.py`
    Fetch assay information from the ChEMBL API for a list of assay IDs.

`get_testitem_data.py`
    Fetch small molecule compounds information from the ChEMBL API for a list of molecule IDs.
    
`get_target_data.py`
    Query biological data sources ([ChEMBL], [UniProt] and [IUPHAR Guide to Pharmacology]) in the combined
    pipeline to provide target data acquisition in a unified CSV table. 

`get_document_data.py`
    Retrieve document information from PubMed, Semantic Scholar, OpenAlex and CrossRef. Also classify 
    documents as Review, non-Review and Unknown.
    
## Installation

Install the runtime dependencies:

```bash
pip install -r requirements.txt
```

Optional type stubs for development:

```bash
pip install pandas-stubs types-requests
```

## Usage
Each sub-command reads an input CSV and writes a new CSV with the
requested annotations. Delimiters and encodings can be customised with
`--sep` and `--encoding`.

Fetch ChEMBL targets for the identifiers in `targets.csv`:

Parse target related data (id, protein and gene names, uniprot id, etc) from chembl db 
python get_target_data.py chembl targets.csv chembl_results.csv
```

Parse UniProt JSON files in `uniprot/` and enrich the accessions listed. 
in `ids.csv`:

```bash
python get_target_data.py uniprot ids.csv uniprot_results.csv --data-dir uniprot
```

Map UniProt IDs to IUPHAR classifications:

```bash
python get_target_data.py iuphar uniprot_results.csv iuphar_results.csv \
    --target-csv data/_IUPHAR_target.csv \
    --family-csv data/_IUPHAR_family.csv
```

### Assay metadata

Retrieve assay information from the ChEMBL API for identifiers listed in
`assays.csv`:

```bash
python get_assay_data.py assays.csv assay_results.csv
```

### Document metadata

Fetch PubMed, Semantic Scholar, OpenAlex and CrossRef records for PMIDs
listed in `pmids.csv`:

```bash
python get_document_data.py pubmed pmids.csv document_data.csv
```

Retrieve document information from the ChEMBL API for the identifiers in
`docs.csv`:

```bash
python get_document_data.py chembl docs.csv chembl_docs.csv
```

Run both pipelines and merge the outputs:

```bash
python get_document_data.py all docs.csv merged_docs.csv
```

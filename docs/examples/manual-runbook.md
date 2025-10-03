# Manual execution runbook

This example collects the one-off command-line invocations used to generate CSV exports from local ChEMBL data files. Adjust the output filenames to match the desired snapshot date before running.

## Activity export

```bash
python scripts/get_activity_data.py \
  --output data/output/activity/output.activity_20250928.csv \
  --input data/input/full/activity3.csv
```

## Assay export

```bash
python scripts/get_assay_data.py \
  --output data/output/assay/output.assay_20250928.csv \
  --input data/input/full/assay.csv
```

## Document export

```bash
python scripts/get_document_data.py \
  --output data/output/document/output.document_20250928.csv \
  --input data/input/full/document.csv
```

## Target export

```bash
python scripts/get_target_data.py \
  --output data/output/target/output.target_20250928.csv \
  --input data/input/full/target.csv
```

## Test item export

```bash
python scripts/get_testitem_data.py \
  --output data/output/testitem/output.testitem_20250928-3.csv \
  --input data/input/full/testitem.csv \
  --limit 1000 \
  --offset 1000
```

> **Tip:** Replace the snapshot date (for example, `20250928`) and chunk suffixes in the output filenames as needed for your run.

# Manual execution runbook

This example collects the one-off command-line invocations used to generate CSV exports from local ChEMBL data files. Adjust the output filenames to match the desired snapshot date before running.

## Activity export

```bash
get-activity-data \
  --input data/input/full/activity3.csv \
  --final-out data/output/activity/output.activity_20250928.csv
```

## Assay export

```bash
get-assay-data \
  --input data/input/full/assay.csv \
  --final-out data/output/assay/output.assay_20250928.csv
```

## Document export

```bash
get-document-data all \
  --input data/input/full/document.csv \
  --final-out data/output/document/output.document_20250928.csv
```

## Target export

```bash
get-target-data all \
  --input data/input/full/target.csv \
  --final-out data/output/target/output.target_20250928.csv
```

## Test item export

```bash
get-testitem-data \
  --input data/input/full/testitem.csv \
  --final-out data/output/testitem/output.testitem_20250928-3.csv \
  --limit 1000 \
  --offset 1000
```

> **Tip:** Replace the snapshot date (for example, `20250928`) and chunk suffixes in the output filenames as needed for your run.

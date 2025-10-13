# Advanced scenarios

This guide collects repeatable patterns for running the pipelines beyond the
basic smoke flow.

## Batch processing multiple drops

Use `--base-path` with date-stamped subdirectories to isolate runs:

```bash
BASE=/data/releases/chembl-35
poetry run get-data \
  --base-path "$BASE" \
  --input-dir inbound/2025-02-01 \
  --output-dir outbound/2025-02-01 \
  --config config/config.yaml \
  --date 20250201 \
  --log-level INFO
```

Store the metadata sidecars in version control to track provenance.

## Partial re-runs with caching

To reprocess only test items:

```bash
python scripts/get_testitem_data.py \
  --input /data/inbound/testitem.csv \
  --final-out /data/outbound/testitems.csv \
  --config config/config.yaml \
  --force  # overwrite the previous export
```

Parent molecule lookups are cached under `$CHEMBL_DA_BASE_PATH/cache`. Delete the
cache when refreshing dictionary data:

```bash
rm -rf "$CHEMBL_DA_BASE_PATH"/cache/molecule_*
```

## Custom rate limiting for staging environments

Override limits temporarily via environment variables:

```bash
export CHEMBL_DA__SOURCES__CHEMBL__API__RPS=5
export CHEMBL_DA__SOURCES__CHEMBL__API__BURST=5
python scripts/get_activity_data.py --limit 200
```

Similarly adjust partner APIs (`CHEMBL_DA_OPENALEX_RPS`, `CHEMBL_DA_PUBMED_TIMEOUT_READ`, …).

## Using raw intermediate outputs

`get_target_data` supports writing raw intermediate datasets:

```bash
python scripts/get_target_data.py all \
  --input data/input/target.csv \
  --final-out output/targets.csv \
  --raw-out output/targets_raw.parquet \
  --raw-format parquet \
  --id-cols target_chembl_id
```

You can inspect the parquet file to debug merges before the final normalisation.

## Loading a custom pipeline registry

When orchestrating alternative workflows (e.g. skipping document enrichment or
adding bespoke QA steps) create a registry YAML and pass it via
`--pipeline-registry`:

```yaml
pipelines:
  - name: document
    callable: scripts.get_document_data:main
    input: document_subset.csv
    output: documents_subset
  - name: target
    callable: scripts.get_target_data:main
    subcommand: chembl
    output: targets_chembl_only
  - name: audit
    callable: tools.audit_pipeline:main
    input: targets_chembl_only.csv
    output: audit_report
```

Invoke the orchestrator with targeted overrides:

```bash
poetry run get-data \
  --base-path /data/chembl \
  --pipeline-registry config/custom_registry.yaml \
  --override-input document=document_snapshot.csv \
  --override-subcommand target=all
```

`--override-output-stem` is useful when reusing the default registry but routing
results to temporary filenames without editing YAML.

## Integrating with Makefile targets

The repository ships with handy targets (see `Makefile`):

- `make lint` — run `ruff` and `black --check`.
- `make typecheck` — run `mypy` against `library/` and `scripts/`.
- `make test` — execute the full pytest suite with JSON reporting.
- `make smoke` — orchestrated smoke run using bundled inputs.

Combine targets in CI pipelines to keep jobs reproducible.

## Exporting QA metrics for dashboards

Quality JSON files can be aggregated with `jq`:

```bash
jq -s 'map({file: input_filename, stats: .summary})' output/*.quality.json > qa_summary.json
```

Feed the result into dashboards or regression monitoring tools to detect trends
in missing data or schema violations.

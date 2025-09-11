# Output Directory Structure

Generated datasets are stored beneath `data/output`. The directory mirrors the
configured `io.output_dir` and may contain project-specific subfolders such as
`ChEMBL/processed`.

Each CSV written by the pipeline is accompanied by a sidecar metadata file with
the same name and the additional suffix `.meta.yaml`. These files capture the
Git commit, exact command line invocation and the configuration values that
produced the dataset so results can be audited and reproduced.

Example layout:

```
data/output/
└── ChEMBL/
    └── processed/
        ├── activity.csv
        ├── activity.csv.meta.yaml
        ├── assay.csv
        ├── assay.csv.meta.yaml
        └── ...
```

The `.meta.yaml` sidecars also record basic table statistics such as the number
of rows written and the SHA-256 checksum of the CSV file. Downstream consumers
can read this metadata to verify integrity or to trace the provenance of a
generated dataset.

# Pipeline configuration schema

All pipeline CLIs share declarative metadata stored under `config/pipeline/`.
Every `<table>.yaml` file describes how the corresponding CLI behaves: which
steps are enabled, the configuration parameters that may be overridden from the
command line, the expected inputs/outputs, and the metadata columns stamped on
exports.  The files are consumed by
`library.pipelines.configuration.load_pipeline_config`, which returns
`PipelineConfig` dataclasses for downstream tooling.

```pycon
>>> from library.pipelines.configuration import load_pipeline_config
>>> cfg = load_pipeline_config("activity")
>>> cfg.name
'activity'
>>> cfg.parameters.mapping_for()
{'column': 'activity.column', 'batch_size': 'activity.batch_size', ...}
```

The loader automatically resolves the installed package version when a file uses
`pipeline_version: value: auto` and exposes convenience helpers to merge
parameter mappings per mode/command.

## Top-level keys

| Key                | Type                          | Description |
|--------------------|-------------------------------|-------------|
| `name`             | string                        | Canonical pipeline identifier. Defaults to the filename stem when omitted. |
| `summary`          | string (optional)             | Human-readable description used in docs and diagnostics. |
| `schemas`          | mapping                       | Mapping of variant → Pandera schema name applied to the export. Use `default` when a pipeline has a single schema. |
| `pipeline_version` | string or mapping             | Describes the metadata column stamped on exports. Accepts either a column name (`pipeline_version: pipeline_version`) or a mapping with `column`, optional `description`, and optional `value`. The loader resolves `value: auto` to the installed package version. |
| `inputs`           | mapping (optional)            | Declares the input artefacts. Recognised keys: `primary` (string) and `dependencies` (list of strings). Additional keys are preserved verbatim in `PipelineIO.variants`. |
| `outputs`          | mapping (optional)            | Declares the output artefacts using the same structure as `inputs`. Common variants include `modes`/`commands` for multi-stage pipelines and `artefacts` for sidecars. |
| `parameters`       | mapping (optional)            | Maps CLI options to configuration paths. Accepts inline key/value pairs for defaults and nested mappings under `default`, `shared`, `modes`, and `commands`. |
| `steps`            | sequence (optional)           | Ordered list describing logical pipeline stages. Each entry is a mapping documented below. |

### Parameter mappings

The `parameters` section supports four buckets:

- Inline key/value pairs or a `default` mapping for options shared by all
  variants.
- `shared`: applied in addition to the defaults for every variant.
- `modes`: mapping of mode name → parameters. Used by the document pipeline.
- `commands`: mapping of sub-command name → parameters. Used by the target
  pipeline.

The loader exposes `PipelineParameters.mapping_for(mode="all")` and
`PipelineParameters.mapping_for(command="chembl")` helpers to merge the
appropriate dictionaries.

### Step entries

Each `steps` entry can define the following keys:

| Key           | Type                | Description |
|---------------|---------------------|-------------|
| `name`        | string              | Identifier for the stage. |
| `description` | string (optional)   | Human-readable summary shown in docs and manifests. |
| `enabled`     | bool/string (optional) | When present, overrides the default `True`. Values other than explicit `false`/`0`/`no` are treated as enabled. |
| `depends_on`  | list of strings (optional) | Declares upstream stages that must finish before this stage runs. |
| `produces`    | list of strings (optional) | Declares artefacts created by the stage. |
| `applies_to`  | string or list of strings (optional) | Restricts the step to specific modes/commands (e.g. `chembl`, `pubmed`, `all`). |

Steps allow downstream tooling (for example the orchestrator) to render clearer
execution plans without hard-coding descriptions.

## Adding a new pipeline definition

1. Create `config/pipeline/<name>.yaml` following the structure above.
2. Populate `schemas` with the Pandera schema(s) used by the CLI.
3. Declare the `pipeline_version` column. Use `value: auto` to resolve the
   current package version at runtime.
4. List default and variant-specific CLI parameters under `parameters` so that
   automation can build accurate config override mappings.
5. Describe the logical steps executed by the pipeline under `steps`.
6. Update or add unit tests covering `library.pipelines.configuration` if you
   introduce new keys so that schema regressions are caught early.

Refer to existing files such as `config/pipeline/activity.yaml` and
`config/pipeline/document.yaml` for concrete examples.

# Variant 1 Post-processing Pipeline — Current Architecture and Roadmap

The Variant 1 modular post-processing pipeline is fully integrated into the ETL
scripts and replaces the historical ad-hoc helpers. This note summarises the
current state, highlights the implementation anchors, and captures the remaining
ideas that have not yet been prioritised.

## Implemented scope

- **Core execution layer.** Transformation steps are represented by
  `StepDefinition` dataclasses and executed sequentially through
  `run_steps`, which provides schema validation hooks, defensive argument
  filtering, and per-step timing/shape metrics. The implementation lives in
  [`library/postprocessing/pipeline/common/types.py`](../library/postprocessing/pipeline/common/types.py)
  and [`library/postprocessing/pipeline/common/runner.py`](../library/postprocessing/pipeline/common/runner.py).
- **Declarative configuration.** Pipelines are described via YAML files under
  [`config/pipeline/*.yaml`](../config/pipeline/), resolved through
  [`library/postprocessing/pipeline/common/config.py`](../library/postprocessing/pipeline/common/config.py).
  The loader expands environment placeholders, validates callable signatures, and
  exposes ordered `StepDefinition` collections for the runner.
- **Schema governance.** Canonical output schemas, column ordering rules, and
  deterministic sorting helpers are defined per domain in the `schema.py` files
  (for example, [`library/postprocessing/pipeline/activities/schema.py`](../library/postprocessing/pipeline/activities/schema.py)).
  The shared helpers in [`library/postprocessing/pipeline/common/schema.py`](../library/postprocessing/pipeline/common/schema.py)
  enforce the contract before/after each pipeline run.
- **Domain step libraries.** Each entity exposes a `steps.py` module that houses
  pure DataFrame transforms and a lazily loaded `PIPELINE_STEPS` tuple sourced
  from the declarative configuration. Examples include
  [`library/postprocessing/pipeline/documents/steps.py`](../library/postprocessing/pipeline/documents/steps.py),
  [`library/postprocessing/pipeline/assays/steps.py`](../library/postprocessing/pipeline/assays/steps.py),
  [`library/postprocessing/pipeline/activities/steps.py`](../library/postprocessing/pipeline/activities/steps.py),
  and [`library/postprocessing/pipeline/targets/steps.py`](../library/postprocessing/pipeline/targets/steps.py).
- **Metrics and reporting.** The execution layer captures structured per-step
  telemetry and final run summaries through
  [`library/postprocessing/pipeline/common/logging.py`](../library/postprocessing/pipeline/common/logging.py).
  CLI entry points collect these metrics and persist JSON reports by delegating to
  [`library/postprocessing/pipeline/common/utils.py`](../library/postprocessing/pipeline/common/utils.py).
- **Integration with ETL scripts.** Entity scripts such as
  [`scripts/get_document_data.py`](../scripts/get_document_data.py),
  [`scripts/get_assay_data.py`](../scripts/get_assay_data.py),
  [`scripts/get_activity_data.py`](../scripts/get_activity_data.py), and
  [`scripts/get_target_data.py`](../scripts/get_target_data.py)
  load the exported CSVs, execute the Variant 1 pipeline, and attach the metrics
  artefacts to the existing logging/QA flow.

Together these components cover Milestones 0–4 of the original plan: discovery,
core infrastructure, declarative configuration, migration of legacy logic, and
QA automation have been delivered by the referenced modules.

## Operational model

1. **Config resolution.** The ETL script determines which pipeline to run and
   resolves the matching YAML file from `config/pipeline`. Environment variables
   (for example `CHEMBL_DOCUMENT_PIPELINE_VERSION`) can override defaults without
   editing the repository configuration.
2. **Runner bootstrap.** The loader materialises an ordered list of steps and
   passes it to `run_steps` together with optional pre/post `DataFrameSchema`
   contracts supplied by each domain module.
3. **Step execution.** Each step receives a defensive copy of the DataFrame and
   only the parameters declared in the callable signature. Unsupported keyword
   arguments are logged and dropped to keep runs deterministic.
4. **Schema enforcement.** Pre- and post-run validations coerce column order,
   ensure required fields, and detect dtype mismatches. Failures raise
   `SchemaValidationError`, which is surfaced to the CLI and reported via the
   metrics payload.
5. **Reporting.** On success the runner stores the final frame with an optional
   `pipeline_version` attribute, records timing/shape deltas per step, and builds
   a JSON summary that is persisted under `reports/` alongside QA artefacts.

## Further improvements

The following ideas remain open for future iterations:

- Extend `config.schema.json` with a section dedicated to post-processing
  pipelines so that declarative configs can be validated outside of runtime
  loading.
- Provide a lightweight CLI/Make wrapper (e.g. `make postprocess TABLE=targets`)
  to execute Variant 1 pipelines independently from the full ETL scripts.
- Publish contributor templates/snippets for new post-processing steps and tests
  to simplify onboarding.
- Document ownership and on-call expectations for the post-processing layer in
  the development handbook.
- Add richer end-to-end fixtures that cover degraded inputs (empty files,
  malformed headers) to harden the integration tests beyond the current happy
  paths.

These items replace the unfinished portions of Milestones 5–6 from the original
backlog and can be pulled in as the pipeline evolves.

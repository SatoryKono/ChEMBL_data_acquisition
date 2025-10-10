# Proposed architecture improvements

1. **Introduce a declarative pipeline registry**
   - **Description.** Centralise pipeline orchestration metadata in `library/pipelines/registry.py` so `scripts/get_data.py` can enumerate stages dynamically instead of hard-coding them.
   - **Problem.** `_PIPELINE_STEPS` is a hand-maintained tuple with step-specific options embedded directly in the orchestrator, which couples the CLI to every pipeline signature and forces manual edits for new stages or flags.
   - **Expected effect.** Adding or reordering steps becomes a data/config change, enabling reuse from tests and alternative entry points while improving maintainability of the ETL chain.

2. **Refactor CLI-heavy pipelines into reusable services**
   - **Description.** Move the fetching, concurrency and enrichment routines from `scripts/get_document_data.py` into a `DocumentPipeline` class under `library/pipelines/document/service.py` that the CLI merely instantiates.
   - **Problem.** The script currently mixes argument parsing, fallback DOI handling, network concurrency and Pandera validation in one file, making it difficult to test and reuse from other orchestrators.
   - **Expected effect.** A service layer would let automated tests and orchestrators call the document pipeline without re-executing CLI boilerplate, improving modularity and readability.

3. **Standardise CLI bootstrap logic**
   - **Description.** Provide a shared `library/cli/bootstrap.py` that encapsulates `ensure_project_root` handling and logging setup so each `scripts/get_*` module can import a single helper.
   - **Problem.** Multiple scripts embed ad-hoc sys.path manipulation and duplicate `_option` helpers, leading to inconsistent startup behaviour and harder maintenance.
   - **Expected effect.** Reduced duplication and a consistent CLI surface that is easier to reason about and extend.

4. **Split configuration loading into a package**
   - **Description.** (Completed) Configuration code now lives in `library/config/models.py`, `library/config/loader.py` and `library/config/runtime.py`, separating schema definitions, loading logic and runtime helpers.
   - **Problem.** The current module interleaves path normalisation, schema expansion, HTTP session factories and environment overrides in one 1k+ LOC file, reducing readability and testability.
   - **Expected effect.** Clear separation of concerns with focused units that can be unit-tested and reused independently across pipelines.

5. **Modularise the JSON schema**
   - **Description.** Replace the monolithic `config.schema.json` with per-pipeline schema fragments (e.g. `config/schema/activity.json`) composed at build time.
   - **Problem.** A single giant schema couples unrelated sections and makes schema reviews or overrides harder because every change touches the same file.
   - **Expected effect.** Pipeline owners can evolve their sections independently while keeping validation reusable across tools.

6. **Inject deterministic retry strategies**
   - **Description.** Allow HTTP clients to accept a deterministic random generator or jitter policy so retry delays can be reproduced under test.
   - **Problem.** Clients compute backoff jitter via `random.uniform`, making retries non-deterministic and hindering reproducibility guarantees.
   - **Expected effect.** Tests and CI runs can assert exact retry timing and logs, supporting deterministic ETL pipelines.

7. **Version dictionary resources via manifests**
   - **Description.** Introduce a manifest file (e.g. `config/dictionary/manifest.yaml`) listing checksums and semantic versions for enrichment CSVs, with helper loaders in `library/resources`.
   - **Problem.** Pipelines consume raw paths from the config without provenance, so updating a dictionary silently changes enrichment outputs without traceability.
   - **Expected effect.** Reproducible resource management with explicit versioning and integrity checks integrated into metadata sidecars.

8. **Derive schema metadata from declarative specs**
   - **Description.** Move column group declarations (document, PubMed, OpenAlex, etc.) into structured metadata (YAML/JSON) that generates both the Pandera schema and CSV ordering logic.
   - **Problem.** Column lists are maintained as Python literals, so any addition requires touching multiple modules and risks divergence between schema and ordering.
   - **Expected effect.** Single-source-of-truth definitions reduce duplication and keep schema validation, ordering and documentation in sync.

9. **Provide pipeline-level dependency injection**
   - **Description.** Pass service factories (ChemblClient, rate limiters, etc.) via a `PipelineContext` object so scripts can be replaced with programmatic runners or mocks.
   - **Problem.** CLI scripts instantiate clients directly, making it hard to substitute mock clients or adjust retry policies without editing the CLI layer.
   - **Expected effect.** Easier testing, clearer separation between infrastructure and business logic, and better support for offline execution.

10. **Align test layout with pipeline layers**
    - **Description.** Collapse helper suites under `tests/unit`, `tests/integration` and `tests/e2e` (e.g. move `tests/integration/postprocessing/` into `tests/integration/`) and document fixtures per stage.
    - **Problem.** The current structure spreads assertions across extra top-level folders, obscuring the unit/integration/e2e boundaries described in the testing policy.
    - **Expected effect.** A predictable test taxonomy that mirrors the ETL layering and simplifies coverage tracking in CI.

## Keeping this document bilingual

- Update this English version and [`../ru/architecture/IMPROVEMENT_PROPOSALS.md`](../../ru/architecture/IMPROVEMENT_PROPOSALS.md) in the same pull request.
- Preserve the numbering, headings and emphasis so the two files stay diff-friendly.
- When adding new proposals, mirror them in both languages and note translation status in the PR description if any section is pending review.


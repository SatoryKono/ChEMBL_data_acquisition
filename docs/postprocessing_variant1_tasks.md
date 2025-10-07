# Variant 1 Post-processing Pipeline — Implementation Tasks

## Overview
This document enumerates the implementation backlog required to introduce the modular post-processing pipeline (Variant 1) to the ChEMBL ETL project. The work is split into milestones that can be delivered incrementally while keeping the current scripts operational.

## Milestone 0 — Discovery & Alignment
- [ ] Confirm current post-processing entry points in `scripts/get_*.py` and existing helpers under `library/postprocessing/`.
- [ ] Align with data consumers on the list of mandatory outputs per entity (`target`, `activity`, etc.).
- [ ] Collect configuration requirements (environment-specific overrides, secrets handling, default locations).

## Milestone 1 — Core Infrastructure
- [ ] Implement `library/postprocessing/base.py` with `PostprocessingStep`, `PipelineContext`, and result typing (DataFrame in/out contract).
- [ ] Add registration decorator (`register_step`) and central registry with uniqueness validation per table/step name.
- [ ] Provide `PipelineRunner` capable of reading YAML definitions, loading registered steps, and executing them sequentially.
- [ ] Define structured error handling and logging hooks (warn vs error) aligned with deterministic logging policy.

## Milestone 2 — Configuration Layer
- [ ] Design YAML schema for step sequences under `config/postprocessing/<table>.yaml` (incl. optional kwargs, gating flags).
- [ ] Extend `config.schema.json` with a section for the post-processing pipeline; generate JSON schema for validation.
- [ ] Implement loader/validator that merges defaults and table-specific overrides while checking schema compliance.
- [ ] Provide tooling to diff/validate configs (CLI command or `make validate-config-postprocessing`).

## Milestone 3 — Migration of Existing Logic
- [ ] Move existing helper functions into dedicated `library/postprocessing/<table>/steps/` modules.
- [ ] Refactor each helper into a `PostprocessingStep` subclass with deterministic I/O and context usage.
- [ ] Update scripts (`scripts/get_target_data.py`, etc.) to invoke `PipelineRunner` after main CSV export.
- [ ] Maintain backward compatibility by ensuring generated files match current naming and sorting conventions.

## Milestone 4 — Testing & Quality Gates
- [ ] Create unit tests for `PostprocessingStep` base behaviors and registry validation under `tests/unit/postprocessing/`.
- [ ] Add integration tests covering YAML-driven pipelines with fixture CSVs (`tests/integration/postprocessing/`).
- [ ] Introduce deterministic e2e tests for at least one entity (e.g., `target`) that execute CLI path and assert outputs and logging.
- [ ] Ensure reports (`reports/test_report.json`, `reports/test_summary.md`) capture new suites and maintain ≥95% success rate.

## Milestone 5 — Tooling & Documentation
- [ ] Update developer docs (`docs/` and README) with instructions for authoring new steps and configuring YAML.
- [ ] Provide code templates/snippets for new step classes and testing patterns.
- [ ] Add `make` or CLI commands for running only post-processing pipelines (e.g., `make postprocess TARGET=target`).
- [ ] Prepare onboarding checklist for future step contributions (coding standards, deterministic behaviors, review checklist).

## Milestone 6 — Rollout & Maintenance
- [ ] Plan phased rollout by table (toggle via YAML) with rollback strategy.
- [ ] Instrument logging/metrics to track pipeline duration and failures per step after rollout.
- [ ] Schedule periodic config reviews to prune unused steps and ensure registry hygiene.
- [ ] Document maintenance responsibilities and escalation path for post-processing failures.

---

### Traceability to Requirements
- Modular step classes with YAML-driven orchestration → Milestones 1–3.
- Deterministic, testable pipeline with ≥95% success rate → Milestone 4.
- Configuration and developer enablement via docs/tooling → Milestone 5.
- Sustainable operations and monitoring → Milestone 6.

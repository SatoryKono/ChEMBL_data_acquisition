# Performance regression analysis — 2025-10-06

## Methodology

- Profiled `scripts/get_activity_data.py` before (`cb23b2d`) and after the CLI refactor on the current `work` branch using a deterministic harness (`tools/profile_activity_pipeline.py`).
- The harness reuses the existing integration fixtures, patches the ChEMBL client with in-memory stubs, and executes the pipeline on 300 identifiers to stress the chunked fetch/write loop while isolating network variance. 【F:tools/profile_activity_pipeline.py†L1-L227】
- Stage timings (load → fetch → process → postprocess → save) were captured by wrapping the corresponding functions and are reported as wall-clock seconds.

## Stage-level timings

| Stage | Pre-refactor (cb23b2d) | Post-refactor (work) | Δ (post − pre) | Relative change |
|---|---|---|---|---|
| Load input identifiers | 0.0035 s | 0.0044 s | +0.0009 s | +25.7% |
| Fetch chunks (stubbed API) | 0.1783 s | 0.1682 s | −0.0101 s | −5.7% |
| Processing (normalise + metadata + validation) | 3.6161 s | 3.3766 s | −0.2395 s | −6.6% |
| Post-processing hook | ≈0.0000 s | ≈0.0000 s | −0.0000 s | −4.0% |
| Deterministic CSV write | 18.5863 s | 17.2740 s | −1.3123 s | −7.1% |
| **Total** | **22.3841 s** | **20.8220 s** | **−1.5621 s** | **−7.0%** |

_Source data:_ post-refactor timings【5c5e8b†L1-L34】, pre-refactor timings【65c7c3†L1-L33】.

## Findings

1. **The deterministic writer dominates runtime.** Both revisions spend ~80–90 % of their time inside `write_csv_chunks_deterministic`, which serialises hundreds of small chunks to temporary files before merging. With the default batch size (3 records per chunk in the harness), this results in 100+ temp files, repeated sorting and repeated YAML sidecar writes. Even modest overhead in this stage amplifies the overall runtime.【5c5e8b†L1-L34】【65c7c3†L1-L33】
2. **Processing and fetch phases improved slightly post-refactor.** The shared `ETLContext` avoids repeatedly constructing `ChemblClient` instances, so the cumulative processing time dropped by ~6 %.【5c5e8b†L1-L34】【65c7c3†L1-L33】
3. **Metadata lookups were retried for every pipeline run.** When dictionary resources are unavailable the pipeline now emits warnings for each metadata lookup attempt (`metadata_dictionary_lookup_failed`), which compounds logging cost and conceals the actual manifest error. Ensuring the allow-list override is configured prevents needless retries.【5c5e8b†L17-L25】
4. **No extra network calls or validation loops were introduced.** The refactor keeps the same chunking algorithm (`ChunkedFetchConfig`) and validator list; cProfile shows no new iterative hot spots (no additional `map/apply` usage inside loops) beyond the deterministic writer.

## Hypotheses on reported 10× slowdown

- **I/O amplification in real runs.** On full datasets the deterministic writer’s merge pass scales with the number of chunk files. If the refactor changed default batch sizes or disabled streaming (`cfg.activity.batch_size` sourced from config), the number of temp files could jump by an order of magnitude, producing the observed 10× slowdown.
- **Dictionary metadata verification.** In production the new pipeline enriches metadata with dictionary checksums. If the manifest is large (thousands of entries) each run now hashes and deserialises the manifest once per execution. Lacking caching across runs could multiply runtime when the manifest resides on network storage.【F:tools/profile_activity_pipeline.py†L53-L101】【5c5e8b†L17-L25】
- **Logging overhead.** The CLI base class enables INFO logging by default. When ETL runs stream millions of records, the added structured log events (per chunk and per validation pass) can increase CPU consumption and I/O, particularly with synchronous appenders.

## Recommendations

1. **Tame deterministic writer overhead.**
   - Increase `activity.batch_size` (e.g., 100–250) so fewer temporary chunk files are emitted.【F:tools/profile_activity_pipeline.py†L107-L119】
   - Allow `write_csv_chunks_deterministic` to stream directly to the final file when chunk order is already deterministic; skip the temp file merge where feasible.
   - Cache per-run metadata (dtypes, column order) so the expensive reconciliation step isn’t repeated for each chunk.

2. **Cache dictionary metadata aggressively.** Persist the parsed manifest and checksum map in a module-level LRU cache (or reuse `_dictionary_manifest_metadata`) to avoid repeated YAML parsing across invocations. Provide a fast-path for allow-listed resources to bypass checksum recomputation.【F:library/common/metadata.py†L134-L200】

3. **Reduce logging verbosity for steady-state loops.** Gating chunk-level INFO logs behind `--verbose` or throttling them to every *n* chunks will cut disk writes during long exports.

4. **Add regression guardrails.** Extend the harness into an automated benchmark (e.g., `pytest --benchmark-only`) that runs in CI after structural refactors. This will flag future increases in the deterministic writer or metadata stages early.

## Next steps

- Confirm production configuration (batch size, dictionary manifest size) to reproduce the reported 10× slowdown and validate the hypotheses above.
- Prototype a streaming writer variant and measure its impact on wall-clock time using the harness.
- Ensure dictionary checksum allowlists are documented to prevent repeated warning floods.

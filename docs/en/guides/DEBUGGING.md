# Debugging playbook

Common issues and remediation steps when running the pipelines.

## 1. CLI argument mistakes

- `get_document_data` requires an explicit mode—pass `--mode <chembl|pubmed|all>`
  or the legacy positional command. `--mode all` (or the `all` positional
  command) runs the combined workflow. `get_target_data` still requires an
  explicit subcommand (`chembl`, `uniprot`, `iuphar`, `all`).
- When combining orchestrator flags, remember that `--limit 0` disables a
  pipeline entirely.
- Use absolute paths or `--base-path` with `--input-dir`/`--output-dir` to avoid
  dependence on the current working directory.

## 2. Encoding / delimiter errors

Symptoms: Pandera complains about missing columns or `csv.Error: iterator should
return strings, not bytes`.

- Check the encoding with `file -bi data/input/document.csv`. Override with
  `--encoding` if necessary.
- Specify `--sep '\t'` for TSV inputs. The reader falls back to tab/semicolon but
  explicit flags are faster and more reliable.

## 3. API failures

- Inspect log entries with `event=api_request` and `status_code`. Persistent 404
  errors usually mean the identifier is invalid or the dictionary is outdated.
- Retryable errors (429/500) should show `retry_attempt` counters; tune
  `sources.*.rps`, `backoff_factor` or `batch_size`.
- For PubChem lookups ensure `sources.pubchem.user_agent` contains a valid email
  address.

## 4. Missing dictionary data

- Targets rely on `config/dictionary/_target` files. Run
  `python scripts/get_target_data.py all --print-config` to confirm the paths.
- For test items check the molecule hierarchy CSV. Use
  `python -m library.integration.molecule_catalog --help` to rebuild the cache if
  necessary.

## 5. Activity enrichment surprises

- `action_type=unknown`: check whether the metric name is absent from
  `activity_enrichment.action_type.metrics` or whether the raw value is missing.
- Negative `standard_value`: inspect the raw `value`/`relation` combination; the
  bounds module clamps values to zero by default.

## 6. Determinism failures

- Ensure the same Python version and dependency lock file are used across runs.
- Clean caches when debugging: remove `$CHEMBL_DA_BASE_PATH/cache` and rerun.
- Use `check-determinism --explain` to list differing files and hashes.

## 7. CI-specific tips

- When running inside containers set `CHEMBL_DA_BASE_PATH=/tmp/chembl-da` to avoid
  permission issues.
- Publish `reports/test_report.json` and `reports/test_summary.md` as artefacts
  for auditability.

If issues persist, gather logs (`data/logs/*.log` by default or `<base>/logs` when
`CHEMBL_DA_BASE_PATH` is set), the metadata sidecar and failing rows, then raise
an issue referencing the troubleshooting steps already taken.

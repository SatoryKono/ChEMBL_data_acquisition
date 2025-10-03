# Quickstart

Follow this checklist to go from a clean checkout to a validated pipeline run. Each step mirrors the development defaults so analysts and engineers land on the same tooling.

## 1. Prepare the environment

1. **Match the Python toolchain.** Install Python `3.11.12` (see `.python-version`) and ensure `pip`, `setuptools`, and `wheel` are current.

   ```bash
   python -m pip install --upgrade pip setuptools wheel
   ```

2. **Create and activate a virtual environment.** Keeping dependencies isolated avoids polluting the system interpreter.

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
   ```

3. **Install the locked dependencies and development extras.** The lock file aligns local and CI environments; the editable install exposes the console scripts used throughout the docs.

   ```bash
   pip install -r requirements-lock.txt
   pip install -e .[dev]
   ```

4. **Enable repository automation.** Pre-commit hooks guarantee formatting, linting, type checks, and tests fire before commits.

   ```bash
   pre-commit install
   ```

## 2. Run a minimal pipeline

1. **Create a scratch output directory.**

   ```bash
   mkdir -p tmp/quickstart
   ```

2. **Execute the target pipeline against the bundled smoke input.** The run fetches a handful of ChEMBL targets, writes both raw and cleaned exports, and exercises the staged logging.

   ```bash
   get-target-data chembl \
     --input tests/data/chembl_targets_min.csv \
     --column target_chembl_id \
     --limit 5 \
     --raw-out tmp/quickstart/targets.raw.csv \
     --final-out tmp/quickstart/targets.final.csv
   ```

   *Use `--log-level DEBUG` if you need verbose diagnostics or `--print-config` to inspect the resolved configuration before the run.*

## 3. Validate the setup

Run the quality gates below to mirror the continuous-integration workflow. All commands assume the virtual environment is still active.

```bash
make lint
make test
make smoke
```

`make lint` expands to Ruff, Black (check-only), and MyPy. `make test` executes the full PyTest suite, while `make smoke` replays the networked smoke tests with in-repo fixtures (skipping the heavyweight test-item pipeline by default).

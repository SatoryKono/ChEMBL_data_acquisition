# Activity log missing investigation

| Potential cause | Confirmed | Comment | Fix |
|-----------------|-----------|---------|-----|
| Log directory resolved relative to the working directory when `CHEMBL_DA_BASE_PATH` is a relative path | ✅ | Reproduced by running `scripts/get_activity_data.py` from `scripts/` with `CHEMBL_DA_BASE_PATH=data`, which wrote logs to `scripts/logs` instead of the repository-level `logs`. | Normalise the log directory against the project root so relative paths are anchored consistently. |
| Logger misconfiguration or missing handlers | ❌ | `setup_cli_logging` attaches a `FileHandler` and `logger.handlers` contains the expected entry during runs. | No change required. |
| File permissions or missing directories | ❌ | `setup_cli_logging` creates the directory with `mkdir(parents=True, exist_ok=True)`; failures would raise an exception. | No change required. |

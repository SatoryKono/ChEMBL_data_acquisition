# ChEMBL Data Acquisition Utilities

This repository contains a suite of Python-based utilities for fetching, processing, and enriching data from the ChEMBL database and other related life sciences resources like UniProt, IUPHAR, PubMed, and PubChem.

| Language | Documentation Root |
|----------|--------------------|
| English  | [`docs/en/README.md`](./docs/en/README.md) |
| Русский  | [`docs/ru/README.md`](./docs/ru/README.md) |

All documentation is maintained in synchronized English and Russian variants.

## Feature highlights

- **Data Pipelines:** End-to-end pipelines for documents, targets, assays, activities, and test items.
- **Orchestration:** A top-level orchestrator (`get-data`) to run all pipelines sequentially with shared configuration.
- **Unified CLI:** A consistent command-line interface with shared flags for I/O, configuration, and logging.
- **Flexible Configuration:** Configure via YAML files, environment variables, and CLI arguments.
- **Quality Gates:** Built-in schema validation, deterministic CSV writing, static analysis, and unit tests to ensure data quality and reproducibility.


## Quick Start

1.  **Create a virtual environment and install dependencies:**

    ```bash
    python -m venv .venv
    source .venv/bin/activate  # On Windows, use: .venv\Scripts\activate
    pip install -r requirements-lock.txt
    ```

2.  **Initialize pre-commit hooks:**

    ```bash
    pre-commit install
    ```

3.  **Explore the main orchestrator script:**

    Run the main `get-data` script with `--help` to see all available options.

    ```bash
    get-data --help
    ```

## Full Documentation

For detailed information on usage, configuration, and output formats, please refer to the full documentation:

- **[Usage Guide](./docs/en/USAGE.md)**: A comprehensive guide to all command-line tools, their arguments, and examples.
  - *[Руководство по использованию](./docs/ru/USAGE.md)*
- **[Configuration Reference](./docs/en/CONFIG.md)**: A complete reference for all configuration options available in `config.yaml` and via environment variables.
  - *[Справка по конфигурации](./docs/ru/CONFIG.md)*
- **[Output Reference](./docs/en/OUTPUT.md)**: A detailed description of the generated files, including CSV formats and metadata sidecars.
  - *[Справка по выходным данным](./docs/ru/OUTPUT.md)*
- **[Developer Guide](./docs/en/DEVELOPER.md)**: Information for developers contributing to the project.
  - *[Руководство для разработчика](./docs/ru/DEVELOPER.md)*
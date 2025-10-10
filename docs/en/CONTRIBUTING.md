# Contributing

> **Languages:** English · [Русский](../../ru/CONTRIBUTING.md)

## Dependency pinning policy

- Runtime dependencies (`[project.dependencies]`) should use compatible release specifiers (e.g. `~=`, `==`) that refer only to versions published on PyPI.
- Development extras (`[project.optional-dependencies].dev`) must not reference unreleased or future-dated placeholder versions.
- Typing stub packages should reference currently released builds (e.g. `package==YYYY.MM.DD` or a lower bound such as `package>=YYYY.MM.DD`) and must never point to placeholder future versions.
- When updating dependencies, run `pip install .[dev]` locally or in CI to verify that all declared extras can be resolved.
- After dependency refreshes land on the main branch, run `pip install -U .[dev]` to pick up the latest compatible toolchain ranges.
- Keep `.pre-commit-config.yaml` `additional_dependencies` in sync with the versions defined in `pyproject.toml` and `requirements-lock.txt`. If you need a newer toolchain, update the version bounds in `pyproject.toml`, regenerate `requirements-lock.txt`, and only then adjust the pre-commit hook pins.
- Mirror updates between `[project.optional-dependencies].dev`, `requirements-dev.txt`, and `requirements-lock.txt` so the development environment remains reproducible. Regenerate the lock file after editing the extras to propagate the new pins everywhere.


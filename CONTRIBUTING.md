# Contributing

## Dependency pinning policy

- Runtime dependencies (`[project.dependencies]`) should use compatible release specifiers (e.g. `~=`, `==`) that refer only to versions published on PyPI.
- Development extras (`[project.optional-dependencies].dev`) must not reference unreleased or future-dated placeholder versions.
- Typing stub packages should reference currently released builds (e.g. `package==YYYY.MM.DD` or a lower bound such as `package>=YYYY.MM.DD`) and must never point to placeholder future versions.
- When updating dependencies, run `pip install .[dev]` locally or in CI to verify that all declared extras can be resolved.
- After dependency refreshes land on the main branch, run `pip install -U .[dev]` to pick up the latest compatible toolchain ranges.


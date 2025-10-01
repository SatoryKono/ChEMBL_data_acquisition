# Release Notes

## Unreleased

- Clarified the minimum supported Python version as 3.11 across documentation and runtime checks.
- Documented the requirement to configure `api.user_agent` with a real contact, including validation behaviour and the
  `CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT` environment override.
- Downgraded ChEMBL cache hit/miss logging to DEBUG to reduce noise in default INFO logs.


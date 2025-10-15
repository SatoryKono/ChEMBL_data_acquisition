# Historical `project/` namespace

This repository previously exposed a top-level `project/` package that
hosted analytical utilities such as the duplicate table profiling tools.
During the migration to the modular `library/` layout these modules were
consolidated under the `library.qa` package.  Some downstream notebooks and
legacy scripts, however, still import from the original `project.*`
locations.

To preserve backwards compatibility we provide a thin compatibility layer in
this directory.  The modules lazily redirect imports to their new home under
`library/`.  Consumers are encouraged to update their imports to point at the
`library` package directly, but the shim keeps existing tooling operational
while coordination with stakeholders continues.

For duplicate analysis specifically, see `project/duplicate_analysis`, which
now proxies to `library.qa.table_quality`, `library.qa.reporting` and
`library.qa.validation`.

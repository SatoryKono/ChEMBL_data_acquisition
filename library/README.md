# `library` package

Utility helpers for ChEMBL data acquisition and processing pipelines. The
package bundles frequently used submodules (CLI helpers, IO helpers, schemas)
and keeps them importable from the top level while remaining usable as a
standalone package inside larger projects.

Historically, some scripts relied on absolute imports which caused import
failures when the package was executed from different working directories. The
package now uses relative imports within submodules, and the top-level
``__init__`` keeps only lightweight re-exports to avoid pulling heavy
dependencies during import.

The following names are intentionally exposed at the top level for downstream
callers:

- Submodules: `cli`, `cli_utils`, `io`, `offline`, `postprocess`, `qa`,
  `schemas`, `testitem_pipeline`, `validation`.
- Sidecar helpers: `SidecarErrors`, `resolve_failure_chunk_size`.

If additional helpers are required upstream, prefer importing them directly
from their defining modules to keep the package boundary narrow.

# Test item pipeline module migration

The test item pipeline modules now live under `library.pipelines.testitem`. The
former `library.testitem_pipeline` package remains available as a compatibility
shim that proxies attribute access to the new location so existing automation
continues to run. New code should import from `library.pipelines.testitem` to
benefit from the consolidated namespace that aligns the test item pipeline with
the other ETL packages.

## Updating imports

Replace legacy imports with the consolidated path:

```python
# Before
from library.testitem_pipeline import cli
from library.testitem_pipeline import pubchem as testitem_pubchem

# After
from library.pipelines.testitem import cli
from library.pipelines.testitem import pubchem as testitem_pubchem
```

The same pattern applies to the catalog and enrichment helpers:

```python
from library.pipelines.testitem import catalog
from library.pipelines.testitem import enrichment
```

Scripts that previously imported the entire package can switch to the
namespaced variant without losing functionality:

```python
from library.pipelines import testitem as pipeline
```

## Legacy shim behaviour

The `library.testitem_pipeline` compatibility package lazily forwards attribute
access to `library.pipelines.testitem` (and its submodules) to minimise upgrade
risk. The shim will continue to exist for a transitional period, but projects
should update to the new import path to avoid future deprecation work.

## Related documentation

- [`docs/en/ARCHITECTURE.md`](../ARCHITECTURE.md)
- [`docs/en/architecture/ARCHITECTURE.md`](../architecture/ARCHITECTURE.md)

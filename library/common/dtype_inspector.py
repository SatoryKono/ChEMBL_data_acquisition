"""Inspect data retrieval dtypes.

This module executes the various ``get_*_data`` retrieval functions on small
sample inputs and reports the resulting :class:`pandas.DataFrame` dtypes.  The
utility aids schema development by showing how external APIs currently
serialise their payloads.  It intentionally keeps the number of records small
so the calls remain inexpensive.

Example
-------
Run interactively from a Python shell::

    from library.common.dtype_inspector import inspect_dtypes
    inspect_dtypes()

The function logs the dtypes for activities, assays, documents, targets and
special test items.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from library.orchestration import ETLContext

from ..config import Config
from ..integration import chembl_library as cl
from .log import logger

# Default sample identifiers for each dataset.  These are deliberately minimal
# and can be overridden via the ``samples`` argument to :func:`inspect_dtypes`.
DEFAULT_SAMPLES: Mapping[str, Sequence[str]] = {
    "activities": ["1"],
    "assays": ["1"],
    "documents": ["1"],
    "targets": ["1"],
    "testitems": ["1"],
}


def _dtype_map(df: pd.DataFrame) -> dict[str, str]:
    """Return column dtypes for ``df`` as plain strings."""

    return {str(col): str(dtype) for col, dtype in df.dtypes.items()}


def inspect_dtypes(
    samples: Mapping[str, Sequence[str]] | None = None,
) -> dict[str, dict[str, str]]:
    """Log and return dtypes produced by retrieval helpers.

    Parameters
    ----------
    samples:
        Optional mapping of dataset name to sample identifiers.  Missing keys
        fall back to :data:`DEFAULT_SAMPLES`.

    Returns
    -------
    dict[str, dict[str, str]]
        Nested mapping of dataset names to column dtype mappings.
    """

    cfg = Config()
    combined = {**DEFAULT_SAMPLES, **(samples or {})}
    results: dict[str, dict[str, str]] = {}

    # ``ChemblClient`` manages HTTP connections and retries.  The context
    # manager ensures the underlying ``requests.Session`` is properly closed
    # even if one of the retrieval functions fails.
    with ETLContext(cfg) as context:
        client = context.chembl_client
        df = cl.get_activities(combined["activities"], cfg=cfg.api, client=client)
        results["activities"] = _dtype_map(df)
        logger.info("dtypes", dataset="activities", dtypes=results["activities"])

        df = cl.get_assays(combined["assays"], cfg=cfg.api, client=client)
        results["assays"] = _dtype_map(df)
        logger.info("dtypes", dataset="assays", dtypes=results["assays"])

        df = cl.get_documents(combined["documents"], cfg=cfg.api, client=client)
        results["documents"] = _dtype_map(df)
        logger.info("dtypes", dataset="documents", dtypes=results["documents"])

        df = cl.get_targets(
            combined["targets"],
            cfg=cfg.api,
            client=client,
            mapping_cfg=cfg.uniprot_mapping,
        )
        results["targets"] = _dtype_map(df)
        logger.info("dtypes", dataset="targets", dtypes=results["targets"])

        df = cl.get_testitem(combined["testitems"], cfg=cfg.api, client=client)
        results["testitems"] = _dtype_map(df)
        logger.info("dtypes", dataset="testitems", dtypes=results["testitems"])

    return results


__all__ = ["inspect_dtypes"]

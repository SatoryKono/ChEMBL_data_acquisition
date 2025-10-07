"""Functions to map ChEMBL target IDs to UniProt accessions using the UniProt ID Mapping API."""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from collections.abc import Callable
from functools import lru_cache
from typing import Any
from urllib.error import HTTPError

from ..common.log import logger
from ..common.rate_limiter import get_limiter
from ..config import UniprotMappingCfg

_UNIPROT_MAPPING_CACHE_MAXSIZE = 128


def _map_chembl_to_uniprot(
    chembl_target_id: str,
    cfg: UniprotMappingCfg,
    opener: Callable[..., Any] | None = None,
) -> str | None:
    """Return the UniProt accession for ``chembl_target_id`` without caching."""
    if opener is None:
        opener = urllib.request.urlopen

    def _open_json(
        url: str, data: bytes | None = None, timeout: float | None = None
    ) -> Any:
        """Open ``url`` and parse the JSON response."""
        try:
            # Use the configured timeout for network requests to avoid hangs.
            if opener is urllib.request.urlopen:
                with opener(url, data=data, timeout=timeout) as response:
                    return json.load(response)
            else:  # For mock openers used in tests
                with opener(url, data=data) as response:
                    return json.load(response)
        except HTTPError as exc:  # pragma: no cover - network failure simulation
            body = ""
            if exc.fp is not None:
                try:
                    body = exc.fp.read().decode()
                except Exception:  # pragma: no cover - fallback if decode fails
                    body = ""
            raise ValueError(
                f"UniProt API request to {url} failed with status {exc.code}: {body or exc.reason}"
            ) from exc

    # Submit the mapping job
    data = urllib.parse.urlencode(
        {"from": "ChEMBL", "to": "UniProtKB", "ids": chembl_target_id}
    ).encode()
    logger.debug("Submitting ID mapping job for %s", chembl_target_id)
    base = cfg.base.rstrip("/")
    run_data = _open_json(f"{base}/run", data=data, timeout=cfg.timeout)
    job_id = run_data.get("jobId")
    if not job_id:
        raise ValueError("UniProt ID Mapping API did not return a job ID")

    # Poll the job status until it finishes
    status_url = f"{base}/status/{job_id}"
    start = time.time()
    result_data = {}  # Initialize result_data
    limiter = get_limiter(
        "uniprot_mapping", 1 / cfg.poll_interval if cfg.poll_interval else 0
    )
    while True:
        limiter.acquire()
        status_data = _open_json(status_url, timeout=cfg.timeout)

        # Check for results in the response, which indicates a redirect has occurred
        if "results" in status_data:
            result_data = status_data
            status = "FINISHED"
        else:
            status = status_data.get("jobStatus") or status_data.get("status")

        logger.debug("Job %s status: %s", job_id, status)
        if status == "FINISHED":
            break
        if status == "FAILED":
            raise ValueError("UniProt ID mapping job failed")
        if time.time() - start > cfg.timeout:
            raise TimeoutError("UniProt ID mapping job timed out")

    # Retrieve the results
    # If the last status check was a redirect, status_data already contains the results
    if not result_data:
        result_url = f"{base}/uniprotkb/results/{job_id}?format=json"
        result_data = _open_json(result_url, timeout=cfg.timeout)

    results = result_data.get("results", [])
    if not results:
        return None

    first = results[0]
    to = first.get("to", {})
    accession = to.get("primaryAccession")
    if not accession:
        raise ValueError("Unexpected response format from UniProt ID mapping API")

    return str(accession)


@lru_cache(maxsize=_UNIPROT_MAPPING_CACHE_MAXSIZE)
def _map_chembl_to_uniprot_cached(
    chembl_target_id: str,
    cfg: UniprotMappingCfg,
    opener: Callable[..., Any] | None,
    _ttl_hash: int | None,
) -> str | None:
    """Cached wrapper for :func:`_map_chembl_to_uniprot`."""

    return _map_chembl_to_uniprot(chembl_target_id, cfg, opener)


def map_chembl_to_uniprot(
    chembl_target_id: str,
    cfg: UniprotMappingCfg,
    opener: Callable[..., Any] | None = None,
) -> str | None:
    """Map a ChEMBL target identifier to a UniProt accession.

    Parameters
    ----------
    chembl_target_id:
        ChEMBL target identifier (e.g., ``"CHEMBL204"``).
    cfg:
        Configuration for the UniProt ID Mapping API including base URL,
        polling interval, timeout and cache settings.
    opener:
        Optional callable matching :func:`urllib.request.urlopen` used to
        perform HTTP requests. Primarily intended for testing.

    Returns
    -------
    str or None
        UniProt accession corresponding to ``chembl_target_id`` or ``None`` if
        no mapping is found.

    Notes
    -----
    Results are cached using :func:`functools.lru_cache` with a maximum size
    of :data:`_UNIPROT_MAPPING_CACHE_MAXSIZE` entries. If ``cfg.cache_ttl`` is
    provided, cache entries expire after the given number of seconds.
    """

    ttl_hash = int(time.time() // cfg.cache_ttl) if cfg.cache_ttl else None
    return _map_chembl_to_uniprot_cached(chembl_target_id, cfg, opener, ttl_hash)

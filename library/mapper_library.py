"""Functions to map ChEMBL target IDs to UniProt accessions using the UniProt ID Mapping API."""

from __future__ import annotations

import time
from collections.abc import Callable
from typing import Any

import requests
from requests import Session
from requests.adapters import HTTPAdapter
from requests.exceptions import HTTPError
from urllib3.util.retry import Retry

from .config import UniprotMappingCfg
from .log import logger
from .rate_limiter import get_limiter

# Default session with basic retry logic and placeholder user agent.
_session: Session = Session()
_adapter = HTTPAdapter(
    max_retries=Retry(
        total=3,
        backoff_factor=0.5,
        status_forcelist=[500, 502, 503, 504],
        allowed_methods=None,
        raise_on_status=False,
    )
)
_session.mount("http://", _adapter)
_session.mount("https://", _adapter)


def map_chembl_to_uniprot(
    chembl_target_id: str,
    cfg: UniprotMappingCfg,
    session: Session | None = None,
) -> str | None:
    """Map a ChEMBL target identifier to a UniProt accession.

    Parameters
    ----------
    chembl_target_id:
        ChEMBL target identifier (e.g., ``"CHEMBL204"``).
    cfg:
        Configuration for the UniProt ID Mapping API including base URL,
        poll interval and timeout settings.
    session:
        Optional :class:`requests.Session` instance used for HTTP requests.
        Defaults to a module-level session with retry logic applied. Primarily
        intended for testing.

    Returns
    -------
    str or None
        UniProt accession corresponding to ``chembl_target_id``, or ``None``
        if no mapping is found.

    Raises
    ------
    ValueError
        If the API reports failure or returns an unexpected response format.
    TimeoutError
        If the mapping job does not complete within ``cfg.timeout`` seconds.
    requests.RequestException
        If a network-related error occurs.
    """

    sess = session or _session

    def _open_json(
        method: Callable[..., requests.Response],
        url: str,
        *,
        data: dict[str, Any] | None = None,
        timeout: float | None = None,
    ) -> Any:
        """Execute an HTTP request and return the decoded JSON body."""

        try:
            kwargs: dict[str, Any] = {"timeout": timeout}
            if data is not None:
                kwargs["data"] = data
            with method(url, **kwargs) as response:
                response.raise_for_status()
                return response.json()
        except HTTPError as exc:  # pragma: no cover - network failure simulation
            body = ""
            if exc.response is not None:
                try:
                    body = exc.response.text
                except Exception:  # pragma: no cover - fallback if decode fails
                    body = ""
            raise ValueError(
                f"UniProt API request to {url} failed with status {exc.response.status_code}: {body or exc.response.reason}"
            ) from exc

    # Submit the mapping job
    data = {"from": "ChEMBL", "to": "UniProtKB", "ids": chembl_target_id}
    logger.debug("Submitting ID mapping job for %s", chembl_target_id)
    base = cfg.base.rstrip("/")
    run_data = _open_json(sess.post, f"{base}/run", data=data, timeout=cfg.timeout)
    job_id = run_data.get("jobId")
    if not job_id:
        raise ValueError("UniProt ID Mapping API did not return a job ID")

    # Poll the job status until it finishes
    status_url = f"{base}/status/{job_id}"
    start = time.time()
    result_data: dict[str, Any] = {}
    limiter = get_limiter(
        "uniprot_mapping", 1 / cfg.poll_interval if cfg.poll_interval else 0
    )
    while True:
        limiter.acquire()
        status_data = _open_json(sess.get, status_url, timeout=cfg.timeout)

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
        result_data = _open_json(sess.get, result_url, timeout=cfg.timeout)

    results = result_data.get("results", [])
    if not results:
        return None

    first = results[0]
    to = first.get("to", {})
    accession = to.get("primaryAccession")
    if not accession:
        raise ValueError("Unexpected response format from UniProt ID mapping API")

    return str(accession)

"""Functions to map ChEMBL target IDs to UniProt accessions using the UniProt ID Mapping API."""

from __future__ import annotations

import json
import logging
import time
import urllib.parse
import urllib.request
from urllib.error import HTTPError
from typing import Optional, Callable, Any
import argparse
import pandas as pd

API_BASE = "https://rest.uniprot.org/idmapping"

# Configure basic logging
logging.basicConfig(level=logging.INFO)

def map_chembl_to_uniprot(
    chembl_target_id: str,
    poll_interval: float = 0.5,
    timeout: float = 300.0,  # Reduced timeout for faster failure if something is wrong
    opener: Optional[Callable[..., Any]] = None,
) -> str | None:
    """Map a ChEMBL target identifier to a UniProt accession.

    Parameters
    ----------
    chembl_target_id:
        ChEMBL target identifier (e.g., ``"CHEMBL204"``).
    poll_interval:
        Seconds to wait between polling the UniProt API for job completion.
    timeout:
        Maximum number of seconds to wait for the mapping job to finish.
    opener:
        Optional callable with the same signature as :func:`urllib.request.urlopen`
        used to perform HTTP requests. Primarily intended for testing.

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
        If the mapping job does not complete within `timeout` seconds.
    URLError
        If a network-related error occurs.
    """

    if opener is None:
        opener = urllib.request.urlopen

    def _open_json(url: str, data: bytes | None = None) -> Any:
        """Open `url` and parse the JSON response."""
        try:
            # Add a 30-second timeout to the underlying urlopen call
            # to prevent hangs on network requests.
            if opener is urllib.request.urlopen:
                with opener(url, data=data, timeout=30) as response:
                    return json.load(response)
            else: # For mock openers used in tests
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
    logging.debug("Submitting ID mapping job for %s", chembl_target_id)
    run_data = _open_json(f"{API_BASE}/run", data=data)
    job_id = run_data.get("jobId")
    if not job_id:
        raise ValueError("UniProt ID Mapping API did not return a job ID")

    # Poll the job status until it finishes
    status_url = f"{API_BASE}/status/{job_id}"
    start = time.time()
    result_data = {} # Initialize result_data
    while True:
        status_data = _open_json(status_url)
        
        # Check for results in the response, which indicates a redirect has occurred
        if "results" in status_data:
            result_data = status_data
            status = "FINISHED"
        else:
            status = status_data.get("jobStatus") or status_data.get("status")

        logging.debug("Job %s status: %s", job_id, status)
        if status == "FINISHED":
            break
        if status == "FAILED":
            raise ValueError("UniProt ID mapping job failed")
        if time.time() - start > timeout:
            raise TimeoutError("UniProt ID mapping job timed out")
        time.sleep(poll_interval)

    # Retrieve the results
    # If the last status check was a redirect, status_data already contains the results
    if not result_data:
        result_url = f"{API_BASE}/uniprotkb/results/{job_id}?format=json"
        result_data = _open_json(result_url)
    
    results = result_data.get("results", [])
    if not results:
        return None

    first = results[0]
    to = first.get("to", {})
    accession = to.get("primaryAccession")
    if not accession:
        raise ValueError("Unexpected response format from UniProt ID mapping API")

    return accession


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map ChEMBL IDs to UniProt accessions in a CSV file.")
    parser.add_argument("csv_file", help="Path to the CSV file to process.")
    parser.add_argument("chembl_id_column", help="Name of the column containing ChEMBL IDs.")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.csv_file)
    except FileNotFoundError:
        logging.error(f"Error: The file '{args.csv_file}' was not found.")
        exit(1)
    except Exception as e:
        logging.error(f"Error reading CSV file: {e}")
        exit(1)

    if args.chembl_id_column not in df.columns:
        logging.error(f"Error: Column '{args.chembl_id_column}' not found in the CSV file.")
        exit(1)

    uniprot_ids = []
    for chembl_id in df[args.chembl_id_column]:
        if pd.isna(chembl_id) or not str(chembl_id).strip():
            uniprot_ids.append(None)
            continue
        try:
            uniprot_id = map_chembl_to_uniprot(str(chembl_id))
            uniprot_ids.append(uniprot_id)
            if uniprot_id:
                logging.info(f"Successfully mapped {chembl_id} to {uniprot_id}")
            else:
                logging.warning(f"No UniProt ID found for {chembl_id}")
        except (ValueError, TimeoutError) as e:
            logging.warning(f"Failed to map {chembl_id}: {e}")
            uniprot_ids.append(None)

    df["mapping_uniprot_id"] = uniprot_ids
    
    try:
        df.to_csv(args.csv_file, index=False)
        logging.info(f"Successfully updated '{args.csv_file}' with UniProt IDs.")
    except Exception as e:
        logging.error(f"Error writing to CSV file: {e}")
        exit(1)
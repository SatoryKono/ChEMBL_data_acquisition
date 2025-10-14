"""Utility helpers orchestrating assay postprocessing steps."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd
import requests

from library.io.config_loader import Config, load_config
from library.schemas.assay_schema import AssayDataSchema
from library.utils.logging import Logger, get_logger
from library.utils.retry import DEFAULT_MAX_TRIES, DEFAULT_TIMEOUT, with_retry

DEFAULT_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
DEFAULT_PAGE_SIZE = 200


@dataclass(frozen=True, slots=True)
class FetchConfig:
    """Runtime parameters controlling assay retrieval."""

    base_url: str = DEFAULT_BASE_URL
    timeout: float = DEFAULT_TIMEOUT
    max_tries: int = DEFAULT_MAX_TRIES
    page_size: int = DEFAULT_PAGE_SIZE


def _build_endpoint(base_url: str, endpoint: str) -> str:
    """Return the absolute URL for ``endpoint`` relative to ``base_url``."""

    return f"{base_url.rstrip('/')}/{endpoint.lstrip('/')}"


def _normalise_limit(limit: int | None) -> int | None:
    """Ensure that ``limit`` is either ``None`` or a non-negative integer."""

    if limit is None:
        return None
    if limit < 0:
        raise ValueError("limit must be non-negative")
    return limit


class AssayFetcher:
    """Iterate over raw assay payloads produced by the ChEMBL API."""

    def __init__(
        self,
        config: FetchConfig,
        *,
        session: requests.Session | None = None,
        logger: Logger | None = None,
    ) -> None:
        self._config = config
        self._session = session or requests.Session()
        self._logger = get_logger(__name__) if logger is None else logger
        self._request = with_retry(
            max_tries=config.max_tries,
            timeout=config.timeout,
            logger=self._logger,
            log_event="chembl_assay_request",
        )(self._session.get)

    def fetch(self, limit: int | None = None) -> list[dict[str, object]]:
        """Return assay payload dictionaries respecting ``limit``."""

        effective_limit = _normalise_limit(limit)
        collected: list[dict[str, object]] = []
        offset = 0
        remaining = effective_limit
        while True:
            page_limit = self._config.page_size
            if remaining is not None:
                if remaining <= 0:
                    break
                page_limit = min(page_limit, remaining)
            payload = self._request_page(offset=offset, limit=page_limit)
            items = payload.get("assays")
            if not isinstance(items, list):
                break
            for item in items:
                if not isinstance(item, dict):
                    continue
                collected.append(item)
                if remaining is not None:
                    remaining -= 1
                    if remaining <= 0:
                        break
            if remaining is not None and remaining <= 0:
                break
            page_meta = payload.get("page_meta")
            next_link = None
            if isinstance(page_meta, Mapping):
                next_link = page_meta.get("next")
            if not next_link:
                break
            offset += page_limit
        return collected

    def _request_page(self, *, offset: int, limit: int) -> Mapping[str, object]:
        endpoint = _build_endpoint(self._config.base_url, "assay.json")
        params = {"offset": offset, "limit": limit}
        response = self._request(endpoint, params=params)
        response.raise_for_status()
        try:
            payload = response.json()
        except ValueError:
            return {}
        if not isinstance(payload, Mapping):
            return {}
        return payload


def _prepare_dataframe(records: Sequence[Mapping[str, object]]) -> pd.DataFrame:
    """Return a normalised DataFrame built from ``records``."""

    if not records:
        columns = [
            "assay_chembl_id",
            "assay_type",
            "assay_test_type",
            "target_chembl_id",
            "created_on",
            "updated_on",
            "assay_strain",
            "assay_group",
            "accession",
        ]
        return pd.DataFrame(columns=columns)

    frame = pd.DataFrame(records)
    expected_columns = {
        "assay_chembl_id",
        "assay_type",
        "assay_test_type",
        "target_chembl_id",
        "created_on",
        "updated_on",
        "assay_strain",
        "assay_group",
        "accession",
    }
    missing = expected_columns - set(frame.columns)
    for column in missing:
        frame[column] = pd.NA
    frame = frame.loc[:, sorted(expected_columns)]
    frame["assay_chembl_id"] = frame["assay_chembl_id"].astype("string")
    for column in [
        "assay_type",
        "assay_test_type",
        "target_chembl_id",
        "assay_strain",
        "assay_group",
        "accession",
    ]:
        frame[column] = frame[column].astype("string")
    return frame


def _derive_timestamp(frame: pd.DataFrame) -> pd.Series:
    """Return ``timestamp_utc`` derived from created/updated columns."""

    created = frame.get("created_on")
    updated = frame.get("updated_on")
    combined = pd.Series(pd.NA, index=frame.index, dtype="string")
    if created is not None:
        combined = created.astype("string")
    if updated is not None:
        combined = combined.fillna(updated.astype("string"))
    return pd.to_datetime(combined, utc=True, errors="coerce")


def _apply_timestamp(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach ``timestamp_utc`` and ``year`` columns to ``frame``."""

    enriched = frame.copy()
    timestamp = _derive_timestamp(enriched)
    enriched["timestamp_utc"] = timestamp
    enriched["year"] = timestamp.dt.year.astype("Int64")
    return enriched


def _load_dictionary_frames(dictionary_dir: Path) -> list[pd.DataFrame]:
    """Load dictionary CSV files located under ``dictionary_dir``."""

    frames: list[pd.DataFrame] = []
    if not dictionary_dir.exists():
        return frames
    for csv_path in sorted(dictionary_dir.glob("*.csv")):
        try:
            frame = pd.read_csv(csv_path, dtype="string")
        except (OSError, ValueError, pd.errors.EmptyDataError):
            continue
        if "assay_chembl_id" not in frame.columns:
            continue
        frames.append(frame)
    return frames


def _merge_dictionary(frame: pd.DataFrame, dictionary: pd.DataFrame) -> pd.DataFrame:
    """Left-join ``dictionary`` rows into ``frame`` on ``assay_chembl_id``."""

    merged = frame.merge(
        dictionary.drop_duplicates(subset=["assay_chembl_id"]),
        how="left",
        on="assay_chembl_id",
        suffixes=("", "_dict"),
    )
    for column in dictionary.columns:
        if column == "assay_chembl_id":
            continue
        dict_column = f"{column}_dict"
        if dict_column in merged.columns:
            merged[column] = merged[column].combine_first(merged[dict_column])
            merged = merged.drop(columns=[dict_column])
    return merged


def enrich_with_dictionaries(frame: pd.DataFrame, dictionary_dir: Path) -> pd.DataFrame:
    """Return ``frame`` enriched with lookup tables stored in ``dictionary_dir``."""

    enriched = frame.copy()
    frames = _load_dictionary_frames(dictionary_dir)
    if not frames:
        return enriched
    enriched["assay_chembl_id"] = enriched["assay_chembl_id"].astype("string")
    for dictionary in frames:
        dictionary["assay_chembl_id"] = dictionary["assay_chembl_id"].astype("string")
        enriched = _merge_dictionary(enriched, dictionary)
    return enriched


def fetch_normalize_assay(
    limit: int | None = None,
    *,
    config: FetchConfig | None = None,
    session: requests.Session | None = None,
    dictionary_dir: Path | None = None,
) -> pd.DataFrame:
    """Fetch assays from ChEMBL and return a validated DataFrame."""

    fetch_config = config or FetchConfig()
    fetcher = AssayFetcher(fetch_config, session=session)
    records = fetcher.fetch(limit)
    frame = _prepare_dataframe(records)
    frame = _apply_timestamp(frame)
    dict_dir = dictionary_dir or Path("dictionary") / "_assay"
    frame = enrich_with_dictionaries(frame, dict_dir)
    validated = AssayDataSchema.validate(frame, lazy=True)
    return validated


def load_assay_data(
    *,
    limit: int | None = None,
    dictionary_dir: Path | None = None,
    config_path: str | Path | None = None,
) -> pd.DataFrame:
    """Return assay data using configuration derived from ``config_path``."""

    cfg = load_config(config_path)
    fetch_cfg = _build_fetch_config(cfg)
    return fetch_normalize_assay(limit, config=fetch_cfg, dictionary_dir=dictionary_dir)


def _build_fetch_config(cfg: Config) -> FetchConfig:
    """Construct :class:`FetchConfig` from ``cfg`` contents."""

    section = cfg.get("assay", {})
    if isinstance(section, Mapping):
        base_url = str(section.get("base_url", DEFAULT_BASE_URL))
        timeout = float(section.get("timeout", DEFAULT_TIMEOUT))
        max_tries = int(section.get("max_tries", DEFAULT_MAX_TRIES))
        page_size = int(section.get("page_size", DEFAULT_PAGE_SIZE))
    else:
        base_url = DEFAULT_BASE_URL
        timeout = DEFAULT_TIMEOUT
        max_tries = DEFAULT_MAX_TRIES
        page_size = DEFAULT_PAGE_SIZE
    return FetchConfig(
        base_url=base_url,
        timeout=timeout,
        max_tries=max_tries,
        page_size=page_size,
    )


__all__ = [
    "AssayFetcher",
    "FetchConfig",
    "enrich_with_dictionaries",
    "fetch_normalize_assay",
    "load_assay_data",
]

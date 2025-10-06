"""Target retrieval helpers for the ChEMBL API."""

from __future__ import annotations

import json
import math
import re
from collections.abc import Iterable, Iterator, Sequence
from typing import Any, Callable

import pandas as pd
import requests

from library.clients import ChemblClient, _chunked
from ...config import ApiCfg, TargetChemblBatchRetryCfg, UniprotMappingCfg
from ...common.log import logger

TARGET_FIELDS = [
    "pref_name",
    "target_chembl_id",
    "component_description",
    "component_id",
    "relationship",
    "gene",
    "uniprot_id",
    "mapping_uniprot_id",
    "chembl_alternative_name",
    "ec_code",
    "hgnc_name",
    "hgnc_id",
    "target_type",
    "tax_id",
    "species_group_flag",
    "target_components",
    "protein_classifications",
    "cross_references",
    "reaction_ec_numbers",
]

EMPTY_TARGET: dict[str, str] = {field: "" for field in TARGET_FIELDS}

TARGET_INCLUDE_PARAMS = "protein_classifications,cross_references"

REACTION_EC_EXCLUDED_XREFS = {
    "REACTOME",
    "RHEA",
    "METACYC",
    "EC_REACTION",
}

_EC_TOKEN_SPLIT = re.compile(r"[|;,/\\\s]+")
_EC_FULL_PATTERN = re.compile(r"^\d+(?:\.(?:\d+|-)){3}$")


def normalize_reaction_ec_numbers(values: Iterable[str | None]) -> str:
    """Return a pipe-delimited string of sanitized EC numbers from ``values``."""

    numbers = _collect_normalized_ec_tokens(values)
    return "|".join(sorted(numbers))


def _collect_normalized_ec_tokens(values: Iterable[str | None]) -> set[str]:
    """Return a set of EC numbers extracted from ``values``."""

    numbers: set[str] = set()
    for value in values:
        if not isinstance(value, str) or not value:
            continue
        for token in _split_ec_candidates(value):
            cleaned = _normalise_ec_token(token)
            if cleaned and _EC_FULL_PATTERN.fullmatch(cleaned):
                numbers.add(cleaned)
    return numbers


def _split_ec_candidates(value: str) -> list[str]:
    """Return tokens parsed from ``value`` using standard separators."""

    return [token for token in _EC_TOKEN_SPLIT.split(value) if token]


def _normalise_ec_token(token: str) -> str:
    """Return a cleaned EC token stripped of prefixes and whitespace."""

    token = token.strip()
    if not token:
        return ""
    upper = token.upper()
    if upper.startswith("EC"):
        token = token[2:]
        token = token.lstrip(":._- ")
    return token.strip()


def _parse_gene_synonyms(synonyms: list[dict[str, str]]) -> str:
    """Return a sorted, pipe separated list of gene synonyms."""
    names = {
        s["component_synonym"]
        for s in synonyms
        if s.get("syn_type") in {"GENE_SYMBOL", "GENE_SYMBOL_OTHER"}
    }
    return "|".join(sorted(names))


def _parse_ec_codes(synonyms: list[dict[str, str]]) -> str:
    """Return a sorted, pipe separated list of EC numbers."""
    codes = {
        s["component_synonym"] for s in synonyms if s.get("syn_type") == "EC_NUMBER"
    }
    return "|".join(sorted(codes))


def _parse_alt_names(synonyms: list[dict[str, str]]) -> str:
    """Return a sorted, pipe separated list of UniProt alternative names."""
    names = {s["component_synonym"] for s in synonyms if s.get("syn_type") == "UNIPROT"}
    return "|".join(sorted(names))


def _parse_uniprot_id(
    xrefs: list[dict[str, str]], chembl_id: str, mapping_cfg: UniprotMappingCfg
) -> tuple[str, str]:
    """Return UniProt IDs from cross references and mapping."""
    uniprot_id = ""
    for x in xrefs:
        src = (x.get("xref_src_db") or "").upper()
        if src in {"UNIPROT", "UNIPROT ACCESSION", "UNIPROT ACC", "UNIPROTKB"}:
            ident = x.get("xref_id", "")
            if ident:
                uniprot_id = ident
                break
    try:
        mapping_uniprot_id = _map_to_uniprot(chembl_id, mapping_cfg) or ""
    except Exception as exc:  # pragma: no cover - network failure paths
        logger.warning(
            "uniprot_mapping_error",
            chembl_id=str(chembl_id),
            error=str(exc),
        )
        mapping_uniprot_id = ""
    return uniprot_id, mapping_uniprot_id


def _parse_hgnc(xrefs: list[dict[str, str]]) -> tuple[str, str]:
    """Extract HGNC name and identifier from a list of cross references."""
    for x in xrefs:
        if x.get("xref_src_db") == "HGNC":
            name = x.get("xref_name", "")
            ident = x.get("xref_id", "")
            hgnc_id = ident.split(":")[-1] if ident else ""
            return name, hgnc_id
    return "", ""


def _stringify(value: Any) -> str:
    """Return a string representation of ``value`` avoiding ``None``."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _serialize_structure(value: Any) -> str:
    """Return a JSON string for ``value`` or ``""`` when empty."""

    if value in (None, "", [], {}):
        return ""
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    except (TypeError, ValueError):
        return _stringify(value)


def _collect_reaction_ec_numbers(components: list[dict[str, Any]]) -> str:
    """Return pipe-delimited reaction EC numbers discovered in ``components``."""

    candidates: list[str] = []
    for component in components:
        synonyms = _get_items(
            component.get("target_component_synonyms"), "target_component_synonym"
        )
        for synonym in synonyms:
            syn_type = (synonym.get("syn_type") or "").upper()
            if syn_type in {"REACTION", "REACTION_NUMBER", "EC_REACTION_NUMBER"}:
                value = synonym.get("component_synonym", "")
                if isinstance(value, str) and value:
                    candidates.append(value)
        xrefs = _get_items(component.get("target_component_xrefs"), "target")
        for xref in xrefs:
            src = (xref.get("xref_src_db") or "").upper()
            if src in REACTION_EC_EXCLUDED_XREFS:
                continue
            value = xref.get("xref_id")
            if isinstance(value, str) and value:
                candidates.append(value)
    return normalize_reaction_ec_numbers(candidates)


def _get_items(container: Any, key: str) -> list[Any]:
    """Return a list of items from a container that may be a dict or list."""
    if isinstance(container, dict):
        items = container.get(key, [])
    else:
        items = container or []
    if isinstance(items, dict):
        return [items]
    if isinstance(items, list):
        return items
    return []


def _extract_target_payload(data: Any) -> dict[str, Any]:
    """Return a target record from ``data`` regardless of envelope shape."""

    if isinstance(data, dict):
        if "target" in data:
            items = data["target"]
            if isinstance(items, list):
                return items[0] if items else {}
            if isinstance(items, dict):
                return items
        if "targets" in data:
            items = data["targets"]
            if isinstance(items, list):
                return items[0] if items else {}
            if isinstance(items, dict):
                return items
        return data
    if isinstance(data, list):
        return data[0] if data else {}
    return {}


def _parse_target_record(
    data: dict[str, Any], mapping_cfg: UniprotMappingCfg
) -> dict[str, Any]:
    """Transform a raw target record into a flat dictionary."""
    components = _get_items(data.get("target_components"), "target_component")
    if not components:
        logger.debug("No components found in target record: %s", data)
        components = []
    comp = components[0] if components else {}
    synonyms = _get_items(
        comp.get("target_component_synonyms"), "target_component_synonym"
    )
    xrefs = _get_items(comp.get("target_component_xrefs"), "target")

    gene_syn = _parse_gene_synonyms(synonyms)
    ec_code = _parse_ec_codes(synonyms)
    alt_name = _parse_alt_names(synonyms)
    uniprot_id, mapping_uniprot_id = _parse_uniprot_id(
        xrefs, data.get("target_chembl_id", ""), mapping_cfg
    )
    hgnc_name, hgnc_id = _parse_hgnc(xrefs)

    res = dict(EMPTY_TARGET)
    res.update(
        {
            "pref_name": data.get("pref_name", ""),
            "target_chembl_id": data.get("target_chembl_id", ""),
            "component_description": comp.get("component_description", ""),
            "component_id": str(comp.get("component_id", "")),
            "relationship": data.get("target_type", ""),
            "gene": gene_syn,
            "uniprot_id": uniprot_id,
            "mapping_uniprot_id": mapping_uniprot_id,
            "chembl_alternative_name": alt_name,
            "ec_code": ec_code,
            "hgnc_name": hgnc_name,
            "hgnc_id": hgnc_id,
            "target_type": _stringify(data.get("target_type")),
            "tax_id": _stringify(data.get("tax_id")),
            "species_group_flag": _stringify(data.get("species_group_flag")),
            "target_components": _serialize_structure(components),
            "protein_classifications": _serialize_structure(
                data.get("protein_classifications")
            ),
            "cross_references": _serialize_structure(data.get("cross_references")),
            "reaction_ec_numbers": _collect_reaction_ec_numbers(components),
        }
    )
    return res


def _normalise_target_payloads(payloads: Sequence[dict[str, Any]]) -> pd.DataFrame:
    """Return a DataFrame created from raw ``payloads`` using :func:`json_normalize`."""

    if not payloads:
        return pd.DataFrame()
    frame = pd.json_normalize(payloads, sep=".")
    return frame


def iter_target_batches(
    ids: Sequence[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    chunk_size: int = 5,
    timeout: float | None = None,
) -> Iterator[tuple[list[dict[str, Any]], pd.DataFrame, pd.DataFrame]]:
    """Yield payloads, raw and parsed target data frames for ``ids``."""

    if not ids:
        return

    base = (
        f"{cfg.chembl_base.rstrip('/')}/target.json?format=json"
        f"&include={TARGET_INCLUDE_PARAMS}"
    )
    effective_timeout = timeout if timeout is not None else cfg.timeout_read

    for chunk in _chunked(ids, chunk_size):
        yield from _iter_target_chunk_with_fallback(
            chunk,
            base_url=base,
            cfg=cfg,
            client=client,
            mapping_cfg=mapping_cfg,
            timeout=effective_timeout,
        )


def _iter_target_chunk_with_fallback(
    chunk: Sequence[str],
    *,
    base_url: str,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    timeout: float,
) -> Iterator[tuple[list[dict[str, Any]], pd.DataFrame, pd.DataFrame]]:
    """Yield processed records for ``chunk`` with timeout-aware retries."""

    if not chunk:
        return

    url = f"{base_url}&target_chembl_id__in={','.join(chunk)}"
    handled_exc: requests.RequestException | None
    try:
        data = client.request_json(url, cfg=cfg, timeout=timeout)
    except requests.ReadTimeout as exc:
        event = "chembl_timeout_split"
        handled_exc = exc
    except requests.RequestException as exc:
        event = "chembl_request_split"
        handled_exc = exc
    else:
        event = ""
        handled_exc = None

    if event:
        if len(chunk) <= 1 or handled_exc is None:
            raise handled_exc
        logger.warning(
            event,
            extra={
                "chunk_size": len(chunk),
                "ids": list(chunk),
                "timeout": timeout,
                "error": str(handled_exc),
            },
        )
        for identifier in chunk:
            yield from _iter_target_chunk_with_fallback(
                [identifier],
                base_url=base_url,
                cfg=cfg,
                client=client,
                mapping_cfg=mapping_cfg,
                timeout=timeout,
            )
        return

    items = data.get("targets") or data.get("target") or []
    payloads = [
        _extract_target_payload(item)
        for item in items
        if isinstance(item, (dict, list))
    ]
    payloads = [payload for payload in payloads if payload]
    if not payloads:
        return
    raw_frame = _normalise_target_payloads(payloads)
    records = [_parse_target_record(payload, mapping_cfg) for payload in payloads]
    parsed_frame = pd.DataFrame(records)
    if parsed_frame.empty:
        parsed_frame = pd.DataFrame(columns=TARGET_FIELDS)
    else:
        parsed_frame = parsed_frame.reindex(columns=TARGET_FIELDS)
    yield payloads, raw_frame, parsed_frame


def iter_target_batches_with_retry(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    chunk_size: int = 5,
    timeout: float | None = None,
    retry_cfg: TargetChemblBatchRetryCfg | None = None,
    log: Any | None = None,
    on_attempt: Callable[[], None] | None = None,
) -> Iterator[tuple[list[dict[str, Any]], pd.DataFrame, pd.DataFrame]]:
    """Yield payloads, raw and parsed frames with adaptive chunk sizing."""

    effective_logger = log or logger
    base_chunk_size = max(1, chunk_size)

    if retry_cfg is None or not getattr(retry_cfg, "enable", False):
        for batch in iter_target_batches(
            ids,
            cfg=cfg,
            client=client,
            mapping_cfg=mapping_cfg,
            chunk_size=base_chunk_size,
            timeout=timeout,
        ):
            if on_attempt is not None:
                on_attempt()
            yield batch
        return

    shrink_factor = retry_cfg.shrink_factor
    min_size = max(1, retry_cfg.min_size)

    buffer: list[str] = []

    def _drain_buffer(batch: list[str]) -> Iterator[tuple[list[dict[str, Any]], pd.DataFrame, pd.DataFrame]]:
        queue = list(batch)
        current_size = min(len(queue), base_chunk_size)

        while queue:
            attempt_size = min(current_size, len(queue))
            if attempt_size <= 0:
                break
            attempt_ids = queue[:attempt_size]
            if on_attempt is not None:
                on_attempt()
            try:
                for result in iter_target_batches(
                    attempt_ids,
                    cfg=cfg,
                    client=client,
                    mapping_cfg=mapping_cfg,
                    chunk_size=attempt_size,
                    timeout=timeout,
                ):
                    yield result
            except (requests.RequestException, ValueError) as exc:
                if attempt_size <= min_size:
                    raise
                next_size = int(math.floor(attempt_size * shrink_factor))
                if next_size < min_size:
                    next_size = min_size
                if next_size >= attempt_size:
                    next_size = max(min_size, attempt_size - 1)
                effective_logger.warning(
                    "chembl_chunk_retry",
                    chunk_size=attempt_size,
                    next_chunk_size=next_size,
                    remaining=len(queue),
                    error=str(exc),
                )
                current_size = next_size
                continue
            del queue[:attempt_size]
            current_size = min(base_chunk_size, len(queue)) if queue else 0

        batch.clear()

    for target_id in ids:
        buffer.append(target_id)
        if len(buffer) >= base_chunk_size:
            yield from _drain_buffer(buffer)
    if buffer:
        yield from _drain_buffer(buffer)

def get_target(
    chembl_target_id: str,
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    timeout: float | None = None,
) -> dict[str, str]:
    """Fetch target information for a single ChEMBL identifier.

    Parameters
    ----------
    chembl_target_id:
        ChEMBL target identifier (e.g., ``"CHEMBL203"``).
    cfg:
        API configuration providing base URL and timeouts.
    client:
        ChemblClient used for HTTP requests and caching.
    mapping_cfg:
        Configuration for the UniProt mapping service.
    timeout:
        Optional override for the HTTP request timeout.
    """
    if chembl_target_id in {"", "#N/A"}:
        return dict(EMPTY_TARGET)
    base = cfg.chembl_base.rstrip("/")
    url = f"{base}/target/{chembl_target_id}.json?include={TARGET_INCLUDE_PARAMS}"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    data = client.request_json(url, cfg=cfg, timeout=effective_timeout)
    payload = _extract_target_payload(data)
    if not payload:
        return dict(EMPTY_TARGET)
    return _parse_target_record(payload, mapping_cfg)


def get_target_payload(
    chembl_target_id: str,
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    timeout: float | None = None,
) -> dict[str, Any]:
    """Return the raw payload for ``chembl_target_id`` without post-processing."""

    if chembl_target_id in {"", "#N/A"}:
        return {}
    base = cfg.chembl_base.rstrip("/")
    url = f"{base}/target/{chembl_target_id}.json?include={TARGET_INCLUDE_PARAMS}"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read
    data = client.request_json(url, cfg=cfg, timeout=effective_timeout)
    return _extract_target_payload(data)


def get_targets(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    chunk_size: int = 5,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Fetch target records for ``ids``.

    Parameters
    ----------
    ids:
        ChEMBL target identifiers to retrieve.
    cfg:
        API configuration with base URL and timeouts.
    client:
        ChemblClient used for HTTP requests and caching.
    mapping_cfg:
        Settings for the UniProt mapping service.
    chunk_size:
        Number of identifiers to request per HTTP call.
    timeout:
        Optional override for the HTTP request timeout.
    """
    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=TARGET_FIELDS)

    parsed_frames: list[pd.DataFrame] = []
    for _, _, parsed in iter_target_batches(
        valid,
        cfg=cfg,
        client=client,
        mapping_cfg=mapping_cfg,
        chunk_size=chunk_size,
        timeout=timeout,
    ):
        if parsed.empty:
            continue
        parsed_frames.append(parsed)
    if not parsed_frames:
        return pd.DataFrame(columns=TARGET_FIELDS)
    df = pd.concat(parsed_frames, ignore_index=True)
    return df.reindex(columns=TARGET_FIELDS)


def get_targets_payloads(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    chunk_size: int = 5,
    timeout: float | None = None,
) -> list[dict[str, Any]]:
    """Return a list of raw payloads for ``ids`` without flattening."""

    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return []

    payloads: list[dict[str, Any]] = []
    for batch_payloads, _, _ in iter_target_batches(
        valid,
        cfg=cfg,
        client=client,
        mapping_cfg=mapping_cfg,
        chunk_size=chunk_size,
        timeout=timeout,
    ):
        payloads.extend(batch_payloads)
    return payloads


def get_targets_raw_frame(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    chunk_size: int = 5,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Return a DataFrame created from the raw payloads for ``ids``."""

    valid = [i for i in ids if i not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame()

    frames: list[pd.DataFrame] = []
    for _, raw_frame, _ in iter_target_batches(
        valid,
        cfg=cfg,
        client=client,
        mapping_cfg=mapping_cfg,
        chunk_size=chunk_size,
        timeout=timeout,
    ):
        if raw_frame.empty:
            continue
        frames.append(raw_frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def extend_target(
    df: pd.DataFrame,
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    mapping_cfg: UniprotMappingCfg,
    id_column: str = "target_chembl_id",
) -> pd.DataFrame:
    """Augment ``df`` with columns returned from :func:`get_target`.

    Parameters
    ----------
    df:
        DataFrame containing a column with ChEMBL target identifiers.
    cfg:
        API configuration providing base URL and timeouts.
    client:
        ChemblClient used for HTTP requests and caching.
    mapping_cfg:
        Settings for the UniProt mapping service.
    id_column:
        Name of the column holding the identifiers.

    """
    if id_column not in df.columns:
        raise ValueError(f"missing required column: {id_column}")
    targets = [
        get_target(i, cfg=cfg, client=client, mapping_cfg=mapping_cfg)
        for i in df[id_column].fillna("")
    ]
    extra = pd.DataFrame(targets)
    return pd.concat([df.reset_index(drop=True), extra], axis=1)


__all__ = [
    "get_target",
    "get_target_payload",
    "get_targets",
    "get_targets_payloads",
    "get_targets_raw_frame",
    "iter_target_batches",
    "extend_target",
    "TARGET_FIELDS",
]
def _map_to_uniprot(chembl_id: str, mapping_cfg: UniprotMappingCfg) -> str | None:
    from ...integration.mapper_library import map_chembl_to_uniprot

    return map_chembl_to_uniprot(chembl_id, mapping_cfg)


"""Low-level PubChem client utilities."""

from __future__ import annotations

import json
import random
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from email.utils import parsedate_to_datetime
from numbers import Real
from threading import Lock, local
from time import monotonic
from typing import Any, cast
from urllib.parse import quote

import requests
from cachetools import TTLCache
from requests import Session

from ..common.log import logger
from ..common.rate_limiter import get_limiter, sleep
from ..config.models import ApiCfg, PubChemCfg, RetryCfg
from ..config.runtime import session_with_retry

__all__ = [
    "Properties",
    "PubChemServiceUnavailable",
    "SERVICE_UNAVAILABLE_OUTCOMES",
    "get_all_cid",
    "get_cid",
    "get_cid_from_inchi",
    "get_cid_from_inchikey",
    "get_cid_from_smiles",
    "get_properties",
    "get_session",
    "get_standard_name",
    "init_session",
    "last_request_outcome",
    "make_request",
    "url_encode",
    "validate_cid",
]


def _config_signature(api: ApiCfg, retry: RetryCfg) -> str:
    """Return a stable signature for the session configuration."""

    api_dump = api.model_dump(mode="json")
    retry_dump = retry.model_dump(mode="json")
    return json.dumps({"api": api_dump, "retry": retry_dump}, sort_keys=True)


@dataclass(frozen=True)
class _CacheEntry:
    """Cache entry capturing both payload and outcome."""

    payload: dict[str, Any] | None
    outcome: str
    details: dict[str, Any] | None = None

    @property
    def is_hit(self) -> bool:
        return self.payload is not None


def _retry_after_seconds(
    value: str | None, *, now: datetime | None = None
) -> float | None:
    """Return the delay in seconds represented by ``Retry-After`` header ``value``."""

    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        seconds = float(text)
    except ValueError:
        try:
            parsed = parsedate_to_datetime(text)
        except (TypeError, ValueError):
            return None
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=UTC)
        reference = now or datetime.now(UTC)
        seconds = (parsed - reference).total_seconds()
    if seconds < 0:
        return 0.0
    return seconds


# Cache is initialised lazily to allow configuration of the TTL via
# :class:`PubChemCfg`. The cache is recreated when the TTL changes.
_CACHE: TTLCache[str, _CacheEntry] | None = None
_CACHE_LOCK = Lock()

_SERVICE_OUTAGE_LOCK = Lock()
_SERVICE_OUTAGE_UNTIL: float | None = None
_SERVICE_OUTAGE_REASON: str | None = None
_SERVICE_OUTAGE_DETAILS: dict[str, Any] | None = None

_SESSION_LOCK = Lock()
_DEFAULT_API_CFG = ApiCfg(user_agent="ChEMBL-ETL/2.1 (mailto:chembl-data@ebi.ac.uk)")
_DEFAULT_RETRY_CFG = RetryCfg()
_SESSION_CFG: tuple[ApiCfg, RetryCfg] = (_DEFAULT_API_CFG, _DEFAULT_RETRY_CFG)
_SESSION_SIGNATURE = _config_signature(*_SESSION_CFG)
_SESSION_INITIALISED = False
_session: Session | None = session_with_retry(*_SESSION_CFG)

_SERVICE_UNAVAILABLE_UNTIL: float | None = None
_SERVICE_UNAVAILABLE_DETAILS: dict[str, Any] | None = None
_SERVICE_UNAVAILABLE_LOCK = Lock()

_PLACEHOLDER_CONTACT = "contact@example.org"
_USER_AGENT_ERROR = (
    "PubChem client requires a custom User-Agent; "
    "set sources.pubchem.user_agent to include valid contact details before making requests."
)

_THREAD_STATE = local()

_JITTER_LOCK = Lock()
_JITTER_GENERATORS: dict[int, random.Random] = {}


def _set_last_outcome(outcome: str | None, details: Mapping[str, Any] | None) -> None:
    if details is not None:
        details = dict(details)
    _THREAD_STATE.last_outcome = (outcome, details)


def _ensure_default_headers(session: Session) -> None:
    """Ensure the shared session advertises JSON support."""

    headers = getattr(session, "headers", None)
    if isinstance(headers, Mapping):
        if "Accept" not in headers:
            try:
                session.headers["Accept"] = "application/json"  # type: ignore[index]
            except Exception:  # pragma: no cover - defensive for custom sessions
                pass
        return
    if hasattr(session, "headers"):
        try:
            session.headers.setdefault("Accept", "application/json")  # type: ignore[attr-defined]
        except AttributeError:
            pass


def _compute_retry_sleep(delay: float, cfg: PubChemCfg) -> tuple[float, float | None]:
    """Return the sleep duration and applied jitter for retry back-offs."""

    jitter_max = getattr(cfg, "retry_jitter_seconds", 0.0) or 0.0
    if delay <= 0 or jitter_max <= 0:
        return delay, None
    seed = getattr(cfg, "retry_jitter_seed", 0)
    if seed is None:
        jitter = random.uniform(0.0, jitter_max)
    else:
        with _JITTER_LOCK:
            rng = _JITTER_GENERATORS.get(seed)
            if rng is None:
                rng = random.Random(seed)
                _JITTER_GENERATORS[seed] = rng
            jitter = rng.uniform(0.0, jitter_max)
    applied = delay + jitter
    if applied < 0:
        applied = 0.0
    return applied, jitter


def _service_outage_remaining(
    *, now: float | None = None,
) -> tuple[float | None, str | None, dict[str, Any] | None]:
    """Return remaining global cooldown seconds, outcome and details if active."""

    current = monotonic() if now is None else now
    with _SERVICE_OUTAGE_LOCK:
        global _SERVICE_OUTAGE_UNTIL, _SERVICE_OUTAGE_REASON, _SERVICE_OUTAGE_DETAILS
        if _SERVICE_OUTAGE_UNTIL is None:
            return None, None, None
        remaining = _SERVICE_OUTAGE_UNTIL - current
        if remaining <= 0:
            _SERVICE_OUTAGE_UNTIL = None
            _SERVICE_OUTAGE_REASON = None
            _SERVICE_OUTAGE_DETAILS = None
            return None, None, None
        details = (
            dict(_SERVICE_OUTAGE_DETAILS)
            if _SERVICE_OUTAGE_DETAILS is not None
            else None
        )
        if details is not None:
            details.setdefault("cache", True)
        return remaining, _SERVICE_OUTAGE_REASON, details


def _start_service_outage(
    outcome: str | None,
    details: Mapping[str, Any] | None,
    cfg: PubChemCfg,
) -> None:
    """Begin a global service cooldown when *outcome* denotes unavailability."""

    if outcome not in SERVICE_UNAVAILABLE_OUTCOMES:
        return

    retry_after: float | None = None
    source: str | None = None
    if details:
        retry_after_value = details.get("retry_after")
        if isinstance(retry_after_value, Real) and float(retry_after_value) > 0:
            retry_after = float(retry_after_value)
            retry_after_source = details.get("retry_after_source")
            if isinstance(retry_after_source, str) and retry_after_source:
                source = retry_after_source
        if retry_after is None:
            timeout_retry = details.get("timeout_retry_after")
            if isinstance(timeout_retry, Real) and float(timeout_retry) > 0:
                retry_after = float(timeout_retry)
                source = "timeout_cache"
            elif isinstance(details.get("timeout_retry_in"), Real):
                pending_retry = float(details["timeout_retry_in"])
                if pending_retry > 0:
                    retry_after = pending_retry
                    source = "timeout_cache"

    if retry_after is None or retry_after <= 0:
        fallback_candidates: tuple[float | None, ...] = (
            cfg.backoff_initial_seconds,
            cfg.delay,
        )
        for candidate in fallback_candidates:
            if isinstance(candidate, Real) and float(candidate) > 0:
                retry_after = float(candidate)
                source = "fallback"
                break

    if retry_after is None or retry_after <= 0:
        return

    now_monotonic = monotonic()
    until = now_monotonic + retry_after
    available_at = datetime.now(UTC) + timedelta(seconds=retry_after)
    stored_details = dict(details) if details else {}
    stored_details.setdefault("reason", outcome)
    stored_details.setdefault("retry_after", retry_after)
    if source and not stored_details.get("retry_after_source"):
        stored_details["retry_after_source"] = source
    stored_details["cooldown_started_at"] = now_monotonic
    stored_details["cooldown_until"] = until
    stored_details["cooldown_available_at"] = available_at.isoformat()

    with _SERVICE_OUTAGE_LOCK:
        global _SERVICE_OUTAGE_UNTIL, _SERVICE_OUTAGE_REASON, _SERVICE_OUTAGE_DETAILS
        if _SERVICE_OUTAGE_UNTIL is not None and until <= _SERVICE_OUTAGE_UNTIL:
            return
        _SERVICE_OUTAGE_UNTIL = until
        _SERVICE_OUTAGE_REASON = outcome
        _SERVICE_OUTAGE_DETAILS = stored_details


def last_request_outcome() -> tuple[str | None, dict[str, Any] | None]:
    state = getattr(_THREAD_STATE, "last_outcome", (None, None))
    outcome, details = state
    if details is None:
        return outcome, None
    return outcome, dict(details)


def _service_unavailable_remaining() -> tuple[float | None, dict[str, Any] | None]:
    """Return remaining Retry-After seconds from the cached service outage."""

    global _SERVICE_UNAVAILABLE_UNTIL, _SERVICE_UNAVAILABLE_DETAILS

    with _SERVICE_UNAVAILABLE_LOCK:
        if _SERVICE_UNAVAILABLE_UNTIL is None:
            return None, None
        now = monotonic()
        if now >= _SERVICE_UNAVAILABLE_UNTIL:
            _SERVICE_UNAVAILABLE_UNTIL = None
            _SERVICE_UNAVAILABLE_DETAILS = None
            return None, None
        remaining = _SERVICE_UNAVAILABLE_UNTIL - now
        details = (
            dict(_SERVICE_UNAVAILABLE_DETAILS)
            if _SERVICE_UNAVAILABLE_DETAILS is not None
            else None
        )
    return remaining, details


def _remember_service_unavailable(
    delay: float,
    details: Mapping[str, Any] | None,
    *,
    source: str,
) -> None:
    """Remember ``delay`` seconds of service outage for future requests."""

    if delay <= 0:
        return

    global _SERVICE_UNAVAILABLE_UNTIL, _SERVICE_UNAVAILABLE_DETAILS

    expires = monotonic() + delay
    stored_details = dict(details or {})
    stored_details.setdefault("reason", "server_error")
    stored_details["retry_after"] = delay
    stored_details["retry_after_source"] = source

    with _SERVICE_UNAVAILABLE_LOCK:
        current_until = _SERVICE_UNAVAILABLE_UNTIL
        if current_until is None or expires > current_until:
            _SERVICE_UNAVAILABLE_UNTIL = expires
            _SERVICE_UNAVAILABLE_DETAILS = stored_details


def _clear_service_unavailable() -> None:
    """Clear cached service outage hints after a successful request."""

    global _SERVICE_UNAVAILABLE_UNTIL, _SERVICE_UNAVAILABLE_DETAILS
    with _SERVICE_UNAVAILABLE_LOCK:
        _SERVICE_UNAVAILABLE_UNTIL = None
        _SERVICE_UNAVAILABLE_DETAILS = None


SERVICE_UNAVAILABLE_OUTCOMES: frozenset[str] = frozenset(
    {"network_error", "rate_limited", "response_error", "server_error", "timeout", "unexpected_status"}
)


class PubChemServiceUnavailable(RuntimeError):
    """Raised when PubChem temporarily cannot satisfy a request."""

    def __init__(
        self,
        outcome: str | None,
        details: Mapping[str, Any] | None = None,
    ) -> None:
        detail_map = dict(details) if details else {}
        status = detail_map.get("status")
        if isinstance(status, Real):
            status_value: int | None = int(status)
        else:
            status_value = None
        message = "PubChem service unavailable"
        if outcome:
            message = f"{message} ({outcome})"
        if status_value is not None:
            message = f"{message}; status={status_value}"
        super().__init__(message)
        self.outcome: str | None = outcome
        self.details: dict[str, Any] = detail_map
        self.status: int | None = status_value


def _raise_for_service_unavailable() -> None:
    outcome, details = last_request_outcome()
    if outcome in SERVICE_UNAVAILABLE_OUTCOMES:
        raise PubChemServiceUnavailable(outcome, details)


def _has_contact_details(user_agent: str | None) -> bool:
    """Return ``True`` when *user_agent* contains usable contact details."""

    if not user_agent:
        return False
    lowered = user_agent.casefold()
    if _PLACEHOLDER_CONTACT in lowered:
        return False
    return "@" in lowered


def _assert_user_agent(user_agent: str | None) -> None:
    """Raise :class:`ValueError` when *user_agent* lacks contact details."""

    if not _has_contact_details(user_agent):
        raise ValueError(_USER_AGENT_ERROR)


def init_session(api: ApiCfg, retry: RetryCfg) -> None:
    """Initialise the shared HTTP session."""

    global _SESSION_CFG, _SESSION_SIGNATURE, _SESSION_INITIALISED, _session
    _assert_user_agent(api.user_agent)
    signature = _config_signature(api, retry)
    new_session = session_with_retry(api, retry)
    _ensure_default_headers(new_session)
    old_session: Session | None = None
    with _SESSION_LOCK:
        old_session = _session
        _session = new_session
        _SESSION_CFG = (api, retry)
        _SESSION_SIGNATURE = signature
        _SESSION_INITIALISED = True
    if old_session is not None:
        old_session.close()


def get_session(cfg: ApiCfg | None = None) -> Session:
    """Return the shared HTTP session, creating or refreshing as needed."""

    global _SESSION_CFG, _SESSION_SIGNATURE, _SESSION_INITIALISED, _session
    old_session: Session | None = None
    with _SESSION_LOCK:
        current_api, current_retry = _SESSION_CFG
        target_api = cfg or current_api
        signature = _config_signature(target_api, current_retry)
        needs_refresh = signature != _SESSION_SIGNATURE or _session is None
        _SESSION_CFG = (target_api, current_retry)
        if needs_refresh:
            _assert_user_agent(target_api.user_agent)
            new_session = session_with_retry(target_api, current_retry)
            old_session = _session
            _session = new_session
            _SESSION_SIGNATURE = signature
        session = _session
    if old_session is not None:
        old_session.close()
    if session is None:
        raise RuntimeError("Failed to initialise PubChem session")
    user_agent = session.headers.get("User-Agent")
    _assert_user_agent(user_agent)
    _ensure_default_headers(session)
    _SESSION_INITIALISED = True
    return session


def url_encode(text: str) -> str:
    """Return ``text`` URL-encoded for safe usage in HTTP requests."""

    return quote(text, safe="")


def _cids_from_identifier_list(data: dict[str, Any]) -> list[str]:
    """Extract CIDs from a JSON ``IdentifierList`` structure."""

    return [str(cid) for cid in data.get("IdentifierList", {}).get("CID", [])]


def get_cid_from_smiles(smiles: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for a SMILES string."""

    base = cfg.base.rstrip("/")
    contains_path_separators = "/" in smiles or "\\" in smiles

    if not contains_path_separators:
        safe_smiles = url_encode(smiles)
        url = f"{base}/compound/smiles/{safe_smiles}/cids/JSON"
        response = make_request(url, cfg)
        if not response:
            _raise_for_service_unavailable()
            return None
    else:
        url = f"{base}/compound/smiles/cids/JSON"
        response = make_request(url, cfg, method="POST", payload={"smiles": smiles})
        if not response:
            _raise_for_service_unavailable()
            return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_cid_from_inchi(inchi: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for an InChI string."""

    safe_inchi = url_encode(inchi)
    base = cfg.base.rstrip("/")
    url = f"{base}/compound/inchi/{safe_inchi}/cids/JSON"
    response = make_request(url, cfg)
    if not response:
        _raise_for_service_unavailable()
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_cid_from_inchikey(inchikey: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for an InChIKey."""

    safe_inchikey = url_encode(inchikey)
    base = cfg.base.rstrip("/")
    url = f"{base}/compound/inchikey/{safe_inchikey}/cids/JSON"
    response = make_request(url, cfg)
    if not response:
        _raise_for_service_unavailable()
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def _build_cache_key(
    method: str, url: str, data: Mapping[str, Any] | None = None
) -> str:
    """Return a stable cache key for ``method``/``url``/``data`` combinations."""

    method_upper = method.upper()
    if not data:
        return f"{method_upper} {url}"
    try:
        serialised = json.dumps(
            data, sort_keys=True, separators=(",", ":"), ensure_ascii=False
        )
    except TypeError:
        serialised = repr(sorted(data.items()))
    return f"{method_upper} {url}?{serialised}"


def make_request(
    url: str,
    cfg: PubChemCfg,
    *,
    method: str = "GET",
    payload: Mapping[str, Any] | None = None,
) -> dict[str, Any] | None:
    """Make an HTTP request and return parsed JSON."""

    global _CACHE
    _set_last_outcome(None, None)
    method_upper = method.upper()
    cache_key = _build_cache_key(method_upper, url, payload)
    timeout_retry_in: float | None = None
    with _CACHE_LOCK:
        cache = _ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
        cached = cache.get(cache_key) if cache is not None else None
        if (
            cached is not None
            and not cached.is_hit
            and cached.outcome == "timeout"
            and cached.details
        ):
            stored_at = cached.details.get("timeout_stored_at")
            retry_after = cached.details.get("timeout_retry_after")
            if isinstance(stored_at, Real) and isinstance(retry_after, Real):
                elapsed = monotonic() - float(stored_at)
                if elapsed < float(retry_after):
                    timeout_retry_in = float(retry_after) - elapsed
                else:
                    cache.pop(cache_key, None)
                    cached = None
    if cached is not None:
        if cached.is_hit:
            logger.debug(
                "cache_hit", url=url, rps=cfg.rps, status="hit", method=method_upper
            )
            _set_last_outcome("hit", {"cache": True})
            return cast(dict[str, Any], cached.payload)
        miss_details: dict[str, Any] = {}
        for key, value in (cached.details or {}).items():
            if key == "timeout_stored_at":
                continue
            if key == "status":
                miss_details["http_status"] = value
            else:
                miss_details[key] = value
        if timeout_retry_in is not None:
            miss_details["timeout_retry_in"] = timeout_retry_in
        logger.debug(
            "cache_hit",
            url=url,
            rps=cfg.rps,
            status="miss",
            method=method_upper,
            outcome=cached.outcome,
            **miss_details,
        )
        cache_details = dict(miss_details)
        cache_details["cache"] = True
        _set_last_outcome(cached.outcome, cache_details)
        return None
    logger.debug("cache_miss", url=url, rps=cfg.rps, status="miss", method=method_upper)

    outage_remaining, outage_outcome, outage_details = _service_outage_remaining()
    if outage_remaining is not None:
        skip_context: dict[str, Any] = {"url": url, "remaining": outage_remaining}
        if outage_outcome:
            skip_context["outcome"] = outage_outcome
        logger.debug("request_service_unavailable_skip", **skip_context)
        details = dict(outage_details) if outage_details else {}
        details.setdefault("reason", outage_outcome)
        details["cooldown_remaining"] = outage_remaining
        _set_last_outcome(outage_outcome, details)
        return None

    api_cfg = ApiCfg(user_agent=cfg.user_agent)

    total_attempts = cfg.retries + 1
    if total_attempts <= 0:
        total_attempts = 1
    backoff_delay = cfg.backoff_initial_seconds
    start_time = monotonic()
    deadline: float | None = None
    deadline_limit: float | None = None
    if cfg.timeout_seconds > 0:
        deadline = start_time + cfg.timeout_seconds
        grace = getattr(cfg, "retry_after_grace_seconds", 0)
        if grace and grace > 0:
            deadline_limit = deadline + grace
    last_failure_details: dict[str, Any] | None = None

    def details_excluding(*keys: str) -> dict[str, Any]:
        if not last_failure_details:
            return {}
        return {
            key: value for key, value in last_failure_details.items() if key not in keys
        }

    for attempt in range(1, total_attempts + 1):
        remaining_retry, cached_details = _service_unavailable_remaining()
        if remaining_retry is not None:
            outcome_details: dict[str, Any] = dict(cached_details or {})
            outcome_details.setdefault("reason", "server_error")
            outcome_details["retry_after"] = remaining_retry
            outcome_details.setdefault("retry_after_source", "cached")
            outcome_details["cache"] = True
            last_failure_details = outcome_details
            log_payload: dict[str, Any] = {
                "url": url,
                "method": method_upper,
                "attempt": attempt,
                "total_attempts": total_attempts,
                "rps": cfg.rps,
                "retry_after": remaining_retry,
            }
            status_hint = outcome_details.get("status")
            if status_hint is not None:
                log_payload["status"] = status_hint
            reason_hint = outcome_details.get("reason")
            if reason_hint:
                log_payload["reason"] = reason_hint
            logger.warning("request_service_unavailable_cached", **log_payload)
            _store_cache_miss(
                cache_key,
                cfg,
                outcome_details["reason"],
                outcome_details,
                url=url,
            )
            _set_last_outcome(outcome_details["reason"], outcome_details)
            return None
        if deadline is not None and monotonic() >= deadline:
            logger.warning(
                "request_timeout",
                url=url,
                method=method_upper,
                attempt=attempt,
                total_attempts=total_attempts,
                rps=cfg.rps,
                **(last_failure_details or {}),
            )
            _store_cache_miss(
                cache_key,
                cfg,
                "timeout",
                last_failure_details,
                url=url,
            )
            logger.debug(
                "request_fail",
                url=url,
                status=None,
                total_attempts=total_attempts,
                method=method_upper,
                rps=cfg.rps,
                **details_excluding("status"),
            )
            outcome_details = dict(last_failure_details or {})
            outcome_details.setdefault("reason", "timeout")
            _start_service_outage("timeout", outcome_details, cfg)
            _set_last_outcome("timeout", outcome_details)
            return None
        event = "request_start" if attempt == 1 else "request_retry"
        logger.debug(
            event,
            url=url,
            attempt=attempt,
            total_attempts=total_attempts,
            rps=cfg.rps,
        )
        get_limiter("pubchem", cfg.rps, cfg.burst).acquire()
        try:
            session = get_session(api_cfg)
            if getattr(session, "verify", True) != cfg.verify:
                session.verify = cfg.verify
            request_kwargs: dict[str, Any] = {
                "timeout": (cfg.timeout_connect, cfg.timeout_read),
            }
            if method_upper != "GET" and payload is not None:
                request_kwargs["data"] = payload
            with session.request(method_upper, url, **request_kwargs) as response:
                status = response.status_code
                if status == 404:
                    last_failure_details = {"reason": "not_found", "status": status}
                    logger.info(
                        "request_not_found",
                        url=url,
                        method=method_upper,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    _store_cache_miss(
                        cache_key,
                        cfg,
                        "not_found",
                        last_failure_details,
                        url=url,
                    )
                    logger.debug(
                        "request_fail",
                        url=url,
                        total_attempts=total_attempts,
                        method=method_upper,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    _set_last_outcome("not_found", last_failure_details)
                    return None
                if status == 429 or 500 <= status < 600:
                    retry_after_header = response.headers.get("Retry-After")
                    retry_after = _retry_after_seconds(retry_after_header)
                    reason = "rate_limited" if status == 429 else "server_error"
                    last_failure_details = {"reason": reason, "status": status}
                    if retry_after is not None:
                        last_failure_details["retry_after"] = retry_after
                        last_failure_details["retry_after_source"] = "header"
                    event_name = (
                        "request_rate_limited"
                        if status == 429
                        else "request_server_error"
                    )
                    warning_context: dict[str, Any] = {
                        "url": url,
                        "method": method_upper,
                        "status": status,
                        "attempt": attempt,
                        "total_attempts": total_attempts,
                        "rps": cfg.rps,
                    }
                    if retry_after is not None:
                        warning_context["retry_after"] = retry_after
                    logger.warning(
                        event_name,
                        **warning_context,
                    )
                    if retry_after is not None:
                        now_retry = monotonic()
                        can_extend_deadline = (
                            deadline is not None
                            and deadline_limit is not None
                            and deadline < deadline_limit
                        )
                        if (
                            deadline is not None
                            and now_retry + retry_after >= deadline
                            and not can_extend_deadline
                        ):
                            if retry_after > 0:
                                _remember_service_unavailable(
                                    retry_after,
                                    last_failure_details,
                                    source=last_failure_details.get(
                                        "retry_after_source", "header"
                                    ),
                                )
                            timeout_details = dict(last_failure_details)
                            timeout_details.setdefault("retry_after", retry_after)
                            timeout_details.setdefault("retry_after_source", "header")
                            timeout_details["timeout_reason"] = "retry_after_exceeds_deadline"
                            logger.warning(
                                "request_timeout",
                                url=url,
                                method=method_upper,
                                attempt=attempt,
                                total_attempts=total_attempts,
                                rps=cfg.rps,
                                **timeout_details,
                            )
                            _store_cache_miss(
                                cache_key,
                                cfg,
                                "timeout",
                                timeout_details,
                                url=url,
                            )
                            status_value = last_failure_details.get("status")
                            logger.debug(
                                "request_fail",
                                url=url,
                                status=status_value,
                                total_attempts=total_attempts,
                                method=method_upper,
                                rps=cfg.rps,
                                **details_excluding("status"),
                            )
                            timeout_info = dict(timeout_details)
                            timeout_info["reason"] = "timeout"
                            _set_last_outcome("timeout", timeout_info)
                            return None
                        if retry_after > 0:
                            _remember_service_unavailable(
                                retry_after,
                                last_failure_details,
                                source=last_failure_details.get(
                                    "retry_after_source", "header"
                                ),
                            )
                        _store_cache_miss(
                            cache_key,
                            cfg,
                            reason,
                            last_failure_details,
                            url=url,
                        )
                        if reason in SERVICE_UNAVAILABLE_OUTCOMES:
                            last_failure_details.setdefault("cache", True)
                        logger.debug(
                            "request_fail",
                            url=url,
                            total_attempts=total_attempts,
                            method=method_upper,
                            rps=cfg.rps,
                            status=status,
                            **details_excluding("status"),
                        )
                        _start_service_outage(reason, last_failure_details, cfg)
                        _set_last_outcome(reason, last_failure_details)
                        return None
                    if attempt >= total_attempts:
                        delay_for_cache = backoff_delay if backoff_delay > 0 else cfg.delay
                        if (
                            delay_for_cache
                            and reason in SERVICE_UNAVAILABLE_OUTCOMES
                        ):
                            last_failure_details.setdefault(
                                "retry_after", delay_for_cache
                            )
                            source_for_cache = (
                                "backoff" if backoff_delay > 0 else "config"
                            )
                            last_failure_details.setdefault(
                                "retry_after_source", source_for_cache
                            )
                            _remember_service_unavailable(
                                delay_for_cache,
                                last_failure_details,
                                source=last_failure_details["retry_after_source"],
                            )
                        _store_cache_miss(
                            cache_key,
                            cfg,
                            reason,
                            last_failure_details,
                            url=url,
                        )
                        if reason in SERVICE_UNAVAILABLE_OUTCOMES:
                            last_failure_details.setdefault("cache", True)
                        logger.debug(
                            "request_fail",
                            url=url,
                            total_attempts=total_attempts,
                            method=method_upper,
                            rps=cfg.rps,
                            status=status,
                            **details_excluding("status"),
                        )
                        _start_service_outage(reason, last_failure_details, cfg)
                        _set_last_outcome(reason, last_failure_details)
                        return None
                    if retry_after is not None:
                        delay = retry_after
                        backoff_delay = delay * 2 if delay > 0 else backoff_delay
                        delay_source = "header"
                        sleep_delay = delay
                        jitter_value: float | None = None
                    else:
                        delay_source = "config"
                        delay = backoff_delay if backoff_delay > 0 else cfg.delay
                        if backoff_delay > 0:
                            delay_source = "backoff"
                        if delay > 0:
                            backoff_delay = delay * 2
                        sleep_delay, jitter_value = _compute_retry_sleep(delay, cfg)
                        if jitter_value is not None and delay > 0:
                            logger.debug(
                                "request_delay_jitter",
                                url=url,
                                method=method_upper,
                                rps=cfg.rps,
                                base_delay=delay,
                                jitter=jitter_value,
                                applied_delay=sleep_delay,
                            )
                        else:
                            sleep_delay = delay
                    if delay > 0:
                        backoff_delay = delay * 2
                        if reason in SERVICE_UNAVAILABLE_OUTCOMES:
                            last_failure_details.setdefault("retry_after", delay)
                            last_failure_details.setdefault(
                                "retry_after_source", delay_source
                            )
                            _remember_service_unavailable(
                                delay,
                                last_failure_details,
                                source=last_failure_details["retry_after_source"],
                            )
                        now = monotonic()
                        if (
                            delay_source == "header"
                            and deadline is not None
                            and deadline_limit is not None
                            and deadline < deadline_limit
                        ):
                            extension_target = max(deadline, now + delay + cfg.timeout_read)
                            capped_extension = min(extension_target, deadline_limit)
                            if capped_extension > deadline:
                                logger.debug(
                                    "request_deadline_extended",
                                    url=url,
                                    method=method_upper,
                                    delay=delay,
                                    previous_remaining=max(deadline - now, 0.0),
                                    new_remaining=max(capped_extension - now, 0.0),
                                    limit_remaining=max(deadline_limit - now, 0.0),
                                )
                                deadline = capped_extension
                        if deadline is not None and now + delay >= deadline:
                            timeout_details = dict(last_failure_details or {})
                            timeout_details.setdefault(
                                "retry_after", delay
                            )
                            timeout_details.setdefault(
                                "retry_after_source",
                                "header" if retry_after is not None else "backoff",
                            )
                            timeout_details["timeout_reason"] = (
                                "retry_after_exceeds_deadline"
                            )
                            logger.warning(
                                "request_timeout",
                                url=url,
                                method=method_upper,
                                attempt=attempt,
                                total_attempts=total_attempts,
                                rps=cfg.rps,
                                **timeout_details,
                            )
                            _store_cache_miss(
                                cache_key,
                                cfg,
                                "timeout",
                                timeout_details,
                                url=url,
                            )
                            status_value = (
                                last_failure_details.get("status")
                                if last_failure_details
                                else None
                            )
                            logger.debug(
                                "request_fail",
                                url=url,
                                status=status_value,
                                total_attempts=total_attempts,
                                method=method_upper,
                                rps=cfg.rps,
                                **details_excluding("status"),
                            )
                            timeout_info = dict(timeout_details)
                            timeout_info["reason"] = "timeout"
                            _set_last_outcome("timeout", timeout_info)
                            return None
                    sleep(sleep_delay)
                    continue
                if status >= 400:
                    if status == 400:
                        last_failure_details = {
                            "reason": "invalid_identifier",
                            "status": status,
                        }
                        logger.info(
                            "request_invalid_identifier",
                            url=url,
                            method=method_upper,
                            rps=cfg.rps,
                            status=status,
                            **details_excluding("status"),
                        )
                        _store_cache_miss(
                            cache_key,
                            cfg,
                            "invalid_identifier",
                            last_failure_details,
                            url=url,
                        )
                    else:
                        last_failure_details = {
                            "reason": "unexpected_status",
                            "status": status,
                        }
                        logger.warning(
                            "request_unexpected_status",
                            url=url,
                            method=method_upper,
                            rps=cfg.rps,
                            status=status,
                            **details_excluding("status"),
                        )
                    logger.debug(
                        "request_fail",
                        url=url,
                        total_attempts=total_attempts,
                        method=method_upper,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    outcome_name = last_failure_details.get("reason") if last_failure_details else None
                    _start_service_outage(outcome_name, last_failure_details, cfg)
                    _set_last_outcome(outcome_name, last_failure_details)
                    return None

                try:
                    response.raise_for_status()
                    response_data = cast(dict[str, Any], response.json())
                except requests.RequestException as exc:  # pragma: no cover - network
                    last_failure_details = {
                        "reason": "response_error",
                        "status": status,
                        "error": str(exc),
                    }
                    if attempt >= total_attempts:
                        logger.error(
                            "request_error",
                            url=url,
                            error=str(exc),
                            attempt=attempt,
                            total_attempts=total_attempts,
                            method=method_upper,
                            rps=cfg.rps,
                        )
                        logger.debug(
                            "request_fail",
                            url=url,
                            total_attempts=total_attempts,
                            method=method_upper,
                            rps=cfg.rps,
                            status=status,
                            **details_excluding("status"),
                        )
                        _start_service_outage("response_error", last_failure_details, cfg)
                        _set_last_outcome("response_error", last_failure_details)
                        return None
                    if cfg.delay > 0:
                        sleep_delay, jitter_value = _compute_retry_sleep(cfg.delay, cfg)
                        if jitter_value is not None:
                            logger.debug(
                                "request_delay_jitter",
                                url=url,
                                method=method_upper,
                                rps=cfg.rps,
                                base_delay=cfg.delay,
                                jitter=jitter_value,
                                applied_delay=sleep_delay,
                            )
                        sleep(sleep_delay)
                    continue
                except ValueError:
                    last_failure_details = {
                        "reason": "response_not_json",
                        "status": status,
                    }
                    logger.warning(
                        "response_not_json",
                        url=url,
                        method=method_upper,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    logger.debug(
                        "request_fail",
                        url=url,
                        total_attempts=total_attempts,
                        method=method_upper,
                        rps=cfg.rps,
                        status=status,
                        **details_excluding("status"),
                    )
                    _set_last_outcome("response_not_json", last_failure_details)
                    return None

                logger.debug(
                    "request_ok",
                    url=url,
                    status=status,
                    method=method_upper,
                    rps=cfg.rps,
                )
                _clear_service_unavailable()
                with _CACHE_LOCK:
                    cache = _ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)
                    cache[cache_key] = _CacheEntry(payload=response_data, outcome="hit")
                logger.debug(
                    "cache_set",
                    url=url,
                    rps=cfg.rps,
                    status="hit",
                    method=method_upper,
                )
                _set_last_outcome("hit", {"status": status})
                return response_data
        except requests.RequestException as exc:  # pragma: no cover - network
            last_failure_details = {
                "reason": "network_error",
                "status": None,
                "error": str(exc),
            }
            if attempt >= total_attempts:
                logger.error(
                    "request_error",
                    url=url,
                    error=str(exc),
                    attempt=attempt,
                    total_attempts=total_attempts,
                    method=method_upper,
                    rps=cfg.rps,
                )
                logger.debug(
                    "request_fail",
                    url=url,
                    status=None,
                    total_attempts=total_attempts,
                    method=method_upper,
                    rps=cfg.rps,
                    **details_excluding("status"),
                )
                _start_service_outage("network_error", last_failure_details, cfg)
                _set_last_outcome("network_error", last_failure_details)
                return None
            if cfg.delay > 0:
                sleep_delay, jitter_value = _compute_retry_sleep(cfg.delay, cfg)
                if jitter_value is not None:
                    logger.debug(
                        "request_delay_jitter",
                        url=url,
                        method=method_upper,
                        rps=cfg.rps,
                        base_delay=cfg.delay,
                        jitter=jitter_value,
                        applied_delay=sleep_delay,
                    )
                sleep(sleep_delay)
            continue

    outcome_name = (last_failure_details or {}).get("reason")
    _start_service_outage(outcome_name, last_failure_details, cfg)
    _set_last_outcome(outcome_name, last_failure_details)
    return None


def _ensure_cache(ttl: float, maxsize: int) -> TTLCache[str, _CacheEntry]:
    """Return the shared cache instance, recreating it when TTL or maxsize change."""

    global _CACHE
    cache = _CACHE
    if cache is None or cache.ttl != ttl or cache.maxsize != maxsize:
        cache = _CACHE = TTLCache(maxsize=maxsize, ttl=ttl)
    return cache


def _store_cache_miss(
    cache_key: str,
    cfg: PubChemCfg,
    outcome: str,
    details: dict[str, Any] | None = None,
    *,
    url: str | None = None,
) -> None:
    """Persist a cached miss outcome for ``url`` including optional details."""

    details_data = dict(details) if details else {}
    cached = False
    log_details: dict[str, Any] = {}

    service_unavailable_outcome = (
        outcome in SERVICE_UNAVAILABLE_OUTCOMES and outcome != "timeout"
    )

    with _CACHE_LOCK:
        cache = _ensure_cache(cfg.cache_ttl, cfg.cache_maxsize)

        if service_unavailable_outcome:
            cache_details = details_data.copy()
            cache_details.setdefault("cache", True)
            cache[cache_key] = _CacheEntry(
                payload=None,
                outcome=outcome,
                details=cache_details,
            )
            cached = True
            log_details = cache_details.copy()
        elif outcome == "timeout":
            base_backoff = (
                cfg.backoff_initial_seconds
                if cfg.backoff_initial_seconds > 0
                else cfg.delay
            )
            retry_after_hint = details_data.get("retry_after")
            if isinstance(retry_after_hint, Real) and float(retry_after_hint) > 0:
                hint_value = float(retry_after_hint)
                if base_backoff and base_backoff > 0:
                    base_backoff = max(base_backoff, hint_value)
                else:
                    base_backoff = hint_value
            existing = cache.get(cache_key)
            if (
                existing is not None
                and not existing.is_hit
                and existing.outcome == "timeout"
                and existing.details
            ):
                previous_retry = existing.details.get("timeout_retry_after")
                if isinstance(previous_retry, Real):
                    doubled = float(previous_retry) * 2
                    if base_backoff and base_backoff > 0:
                        base_backoff = max(base_backoff, doubled)
                    else:
                        base_backoff = doubled
            max_backoff = (
                cfg.timeout_seconds
                if cfg.timeout_seconds and cfg.timeout_seconds > 0
                else None
            )
            retry_after_source = details_data.get("retry_after_source")
            effective_backoff = (
                base_backoff if base_backoff and base_backoff > 0 else None
            )
            if effective_backoff is not None:
                if max_backoff is not None and retry_after_source != "header":
                    effective_backoff = min(effective_backoff, max_backoff)
                stored_at = monotonic()
                details_data.update(
                    {
                        "timeout": True,
                        "timeout_retry_after": effective_backoff,
                        "timeout_stored_at": stored_at,
                    }
                )
                cache[cache_key] = _CacheEntry(
                    payload=None,
                    outcome=outcome,
                    details=details_data.copy(),
                )
                cached = True
                log_details = {
                    key: value
                    for key, value in details_data.items()
                    if key != "timeout_stored_at"
                }
            else:
                cache.pop(cache_key, None)
                log_details = details_data.copy()
        else:
            cache[cache_key] = _CacheEntry(
                payload=None,
                outcome=outcome,
                details=details_data.copy() if details_data else None,
            )
            cached = True
            log_details = details_data.copy()
    log_data: dict[str, Any] = {
        "url": url or cache_key,
        "rps": cfg.rps,
        "status": "miss",
        "outcome": outcome,
    }
    if log_details:
        log_data.update(log_details)
    event = "cache_set" if cached else "cache_skip"
    logger.debug(event, **log_data)


def _extract_cids(bindings: list[dict[str, Any]]) -> list[str]:
    """Extract CIDs from API bindings."""

    cids: list[str] = []
    for item in bindings:
        cid_field = item.get("cid")
        if isinstance(cid_field, dict):
            cid_value = cid_field.get("value", "")
        else:
            cid_value = str(cid_field)
        cid_value = cid_value.replace(
            "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID", ""
        )
        if cid_value:
            cids.append(cid_value)
    return cids


def _dedup_join(cids: list[str]) -> str | None:
    """Return a ``|``-joined, sorted representation of ``cids`` or ``None``."""

    unique_cids = sorted({cid for cid in cids if cid})
    return "|".join(unique_cids) if unique_cids else None


def _get_cids_from_name_via_rdf(
    compound_name: str, cfg: PubChemCfg, *, partial: bool
) -> list[str]:
    safe_name = url_encode(compound_name)
    rdf_base = cfg.base.rstrip("/").rsplit("/", 1)[0] + "/rdf"
    contain_suffix = "&contain=true" if partial else ""
    url = (
        f"{rdf_base}/query?graph=synonym&return=cid&format=json&name="
        f"{safe_name}{contain_suffix}"
    )
    response = make_request(url, cfg)
    if not response:
        return []
    bindings = response.get("results", {}).get("bindings", [])
    return _extract_cids(bindings)


def _get_cids_from_name_via_pug(
    compound_name: str, cfg: PubChemCfg, *, partial: bool
) -> list[str]:
    safe_name = url_encode(compound_name)
    base = cfg.base.rstrip("/")
    url = f"{base}/compound/name/{safe_name}/cids/JSON"
    if partial:
        separator = "?" if "?" not in url else "&"
        url = f"{url}{separator}name_type=word"
    response = make_request(url, cfg)
    if not response:
        return []
    return _cids_from_identifier_list(response)


def get_cid(compound_name: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for *compound_name* (exact match)."""

    cids = _get_cids_from_name_via_rdf(compound_name, cfg, partial=False)
    if not cids:
        cids = _get_cids_from_name_via_pug(compound_name, cfg, partial=False)
    return _dedup_join(cids)


def get_all_cid(compound_name: str, cfg: PubChemCfg) -> str | None:
    """Retrieve PubChem CID(s) for *compound_name* (partial match)."""

    cids = _get_cids_from_name_via_rdf(compound_name, cfg, partial=True)
    if not cids:
        cids = _get_cids_from_name_via_pug(compound_name, cfg, partial=True)
    return _dedup_join(cids)


def validate_cid(cid: str) -> str | None:
    """Validate PubChem CID."""

    if cid in {"", "0", "-1"}:
        return None
    return cid


def get_standard_name(cid: str, cfg: PubChemCfg) -> str | None:
    """Retrieve the standard compound name for a given CID."""

    validated = validate_cid(cid)
    if not validated:
        return None
    base = cfg.base.rstrip("/")
    url = f"{base}/compound/cid/{validated}/description/JSON"
    response = make_request(url, cfg)
    if not response:
        _raise_for_service_unavailable()
        return None
    info = response.get("InformationList", {}).get("Information", [])
    if not info:
        return None
    return cast(str | None, info[0].get("Title"))


@dataclass
class Properties:
    """Chemical properties for a PubChem compound."""

    IUPACName: str | None
    MolecularFormula: str | None
    iSMILES: str | None
    cSMILES: str | None
    InChI: str | None
    InChIKey: str | None


def get_properties(cid: str, cfg: PubChemCfg) -> Properties:
    """Retrieve chemical properties for a compound by CID."""

    validated = validate_cid(cid)
    if not validated:
        return Properties(None, None, None, None, None, None)
    base = cfg.base.rstrip("/")
    url = (
        f"{base}/compound/cid/{validated}/property/MolecularFormula,IUPACName,SMILES,"
        "ConnectivitySMILES,InChI,InChIKey/JSON"
    )
    response = make_request(url, cfg)
    if not response:
        _raise_for_service_unavailable()
        return Properties(None, None, None, None, None, None)
    props = response.get("PropertyTable", {}).get("Properties", [])
    if not props:
        return Properties(None, None, None, None, None, None)
    item = props[0]
    return Properties(
        cast(str | None, item.get("IUPACName")),
        cast(str | None, item.get("MolecularFormula")),
        cast(str | None, item.get("SMILES")),
        cast(str | None, item.get("ConnectivitySMILES")),
        cast(str | None, item.get("InChI")),
        cast(str | None, item.get("InChIKey")),
    )

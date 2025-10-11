"""Path helper functions for pipeline outputs."""

from __future__ import annotations

import re
from datetime import UTC, datetime
from pathlib import Path

from ..config import IoCfg

_DATE_TOKEN_RE = re.compile(r"^\d{8}$", flags=re.ASCII)


def _strip_output_prefix(value: str) -> str:
    """Return ``value`` without redundant ``output.`` prefixes or leading dots."""

    candidate = value.lstrip(".")
    prefix = "output."
    while candidate.lower().startswith(prefix) and candidate:
        candidate = candidate[len(prefix) :].lstrip(".")
    return candidate


def _strip_csv_suffix(value: str) -> str:
    """Remove trailing ``.csv`` tokens and surrounding dots from ``value``."""

    candidate = value
    lowered = candidate.lower()
    while lowered.endswith(".csv") and candidate:
        candidate = candidate[: -len(".csv")].rstrip(".")
        lowered = candidate.lower()
    return candidate


def _strip_date_suffix(value: str) -> tuple[str, str | None]:
    """Return ``value`` without ``.csv_<YYYYMMDD>`` tail and the extracted date."""

    match = re.search(r"\.csv_(\d{8})$", value, flags=re.IGNORECASE)
    if match and _DATE_TOKEN_RE.match(match.group(1)):
        stripped = value[: match.start()].rstrip("_" )
        return stripped, match.group(1)
    return value, None


def _normalise_table_name(raw: str, default: str) -> str:
    """Return a canonical table label derived from ``raw``."""

    candidate = _strip_output_prefix(_strip_csv_suffix(raw.strip()))
    if len(candidate) >= 9 and candidate[-9] == "_" and _DATE_TOKEN_RE.match(
        candidate[-8:]
    ):
        candidate = candidate[:-9]
    candidate = candidate.strip("._")
    return candidate or default


def _coerce_date_token(value: str | None) -> str | None:
    """Return ``value`` when it encodes an eight digit date, otherwise ``None``."""

    if value is None:
        return None
    candidate = value.strip().strip("._")
    candidate = _strip_csv_suffix(candidate)
    if _DATE_TOKEN_RE.match(candidate):
        return candidate
    return None


def default_output_path(
    input_path: str | Path,
    cfg: IoCfg,
    *,
    date: str | None = None,
    stamp_mode: str | None = None,
) -> Path:
    """Construct a deterministic output path based on the input file name.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Location of the CSV supplied to the pipeline. Only the stem is used
        when generating the output file name.
    cfg : IoCfg
        IO configuration containing the destination directory for derived
        artefacts.
    date : str, optional keyword-only
        Date prefix explicitly supplied by the caller. When omitted, the
        function falls back to ``cfg.default_date_prefix`` or generates a
        ``YYYYMMDD`` string for the current UTC date when the effective
        ``stamp_mode`` resolves to ``"date"``.
    stamp_mode : {"omit", "date", "require"}, optional keyword-only
        Overrides the configured strategy controlling how date components are
        appended to the output file name. When ``None`` (the default) the
        function consults ``cfg.output_stamp_mode`` and uses ``"date"`` as the
        ultimate fallback for backwards compatibility.

    Returns
    -------
    pathlib.Path
        Path pointing to ``cfg.output_dir / f"output_<stem>_<YYYYMMDD>.csv"``.

    Raises
    ------
    TypeError
        Raised when ``cfg.output_dir`` cannot be coerced into a
        :class:`pathlib.Path` instance.
    """

    inp = Path(input_path)
    stem = inp.stem
    suffixes = inp.suffixes

    # When the source file carries multiple extensions (for example
    # ``.csv.tmp``), :attr:`Path.stem` only removes the last suffix.  This
    # leaves artefacts such as ``<stem>.csv`` in the name, which would then
    # result in outputs like ``output.<stem>.csv.csv``.  Trim the additional
    # ``.csv`` segment proactively to keep the generated filenames stable.
    if (
        len(suffixes) >= 2
        and suffixes[-2].lower() == ".csv"
        and stem.lower().endswith(".csv")
    ):
        stem = stem[: -len(".csv")]

    # Normalise hidden temporary files that start with ``.output.`` to avoid
    # producing names with double dots such as ``output..output.<name>``.
    if stem.startswith("."):
        stripped = stem.lstrip(".")
        stem = stripped or stem
    prefix = "output."
    if stem.lower().startswith(prefix):
        stripped = stem
        # Avoid generating names like ``output.output.<stem>`` when the
        # source filename already carries the canonical prefix. The slicing is
        # case-preserving to keep any meaningful capitalisation while still
        # matching against lower-cased prefixes.
        while stripped.lower().startswith(prefix) and stripped:
            candidate = stripped[len(prefix) :]
            # ``output..foo`` style stems introduce a leading ``.`` after the
            # first iteration. Normalise it here so the next loop iteration can
            # drop the remaining redundant ``output.`` prefix.
            candidate = candidate.lstrip(".")
            if not candidate:
                # If the filename is literally ``output`` (with or without the
                # trailing dot), fall back to the previous value to keep the
                # stem non-empty.
                break
            stripped = candidate
        stem = stripped or stem

    # Determine the stamp mode
    mode = stamp_mode or getattr(cfg, "output_stamp_mode", None) or "date"
    if isinstance(mode, str):
        normalized_mode = mode.strip().lower()
    else:
        normalized_mode = "date"

    if normalized_mode not in {"omit", "date", "require"}:
        raise ValueError(
            "stamp_mode must be 'omit', 'date' or 'require'",
        )

    # Try provided date string first
    date_str: str | None = None
    if isinstance(date, str):
        candidate = date.strip()
        if candidate:
            date_str = candidate

    # Short-circuit for "omit" mode (no date in filename)
    if normalized_mode == "omit":
        return Path(cfg.output_dir) / f"output.{stem}.csv"

    # For "require" mode, a date must be present, from caller
    if normalized_mode == "require" and date_str is None:
        raise ValueError("date must be provided when stamp_mode is 'require'")

    # Try to use default_date_prefix from cfg if available and not already found
    if date_str is None:
        cfg_prefix = getattr(cfg, "default_date_prefix", None)
        if isinstance(cfg_prefix, str):
            candidate = cfg_prefix.strip()
            if candidate:
                date_str = candidate

    # Fallback to current utc date string if nothing else
    if date_str is None:
        date_str = datetime.now(UTC).strftime("%Y%m%d")

    return Path(cfg.output_dir) / f"output.{stem}_{date_str}.csv"


def derive_output_labels(
    source: str | Path,
    *,
    default_table: str,
    fallback_date: str | None = None,
) -> tuple[str, str]:
    """Return the logical table name and date tag inferred from ``source``.

    Parameters
    ----------
    source:
        Path pointing to a dataset generated by the pipeline. Only the filename
        is inspected; parent directories do not affect the derived labels.
    default_table:
        Table identifier to use when ``source`` does not encode one explicitly.
    fallback_date:
        Optional date token to return when ``source`` does not include a
        timestamp. When omitted the current UTC date is used, matching the
        behaviour of :func:`default_output_path`.

    Returns
    -------
    (str, str)
        Tuple of ``(table_name, date_tag)`` suitable for constructing canonical
        output filenames.
    """

    path = Path(source)
    base = _strip_output_prefix(path.stem)
    base, trailing_date = _strip_date_suffix(base)
    base = _strip_csv_suffix(base)

    table_candidate = base or default_table

    date_candidates: list[str | None] = []
    if "_" in base:
        maybe_table, maybe_date = base.rsplit("_", 1)
        if maybe_table:
            table_candidate = maybe_table
        date_candidates.append(maybe_date)
    date_candidates.append(trailing_date)
    date_candidates.append(fallback_date)

    table_name = _normalise_table_name(table_candidate, default_table)

    for candidate in date_candidates:
        token = _coerce_date_token(candidate)
        if token is not None:
            date_tag = token
            break
    else:
        date_tag = datetime.now(UTC).strftime("%Y%m%d")

    return table_name, date_tag

"""Path helper functions for pipeline outputs."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from ..config import IoCfg


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
    base = path.stem

    if base.endswith(".csv"):
        base = base[: -len(".csv")]

    # Normalise leading dots introduced by temporary files and repeated
    # ``output.`` prefixes to ensure the inferred table name and timestamp do
    # not retain placeholder markers. Use a loop to gracefully handle chained
    # prefixes such as ``.output.output.activities``.
    base = base.lstrip(".")
    prefix = "output."
    while base.lower().startswith(prefix):
        base = base[len(prefix) :]
        base = base.lstrip(".")

    table_candidate = base or default_table
    date_candidate = fallback_date
    if not isinstance(date_candidate, str) or not date_candidate.strip():
        date_candidate = datetime.now(UTC).strftime("%Y%m%d")

    if "_" in base:
        maybe_table, maybe_date = base.rsplit("_", 1)
        if maybe_table:
            table_candidate = maybe_table
        if maybe_date:
            date_candidate = maybe_date

    table_name = table_candidate or default_table
    date_tag = date_candidate or datetime.now(UTC).strftime("%Y%m%d")
    return table_name, date_tag

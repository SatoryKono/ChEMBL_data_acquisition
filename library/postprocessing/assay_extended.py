"""Assay metadata enrichment mirroring the legacy Power Query workbook."""

from __future__ import annotations

import functools
import json
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from library.common.log import logger
from config.paths import DICTIONARY_DIR

from . import helpers

__all__ = ["AssayExtendedError", "enrich_assay_metadata"]

_DEFAULT_DICTIONARY_DIR = DICTIONARY_DIR


class AssayExtendedError(RuntimeError):
    """Raised when the assay enrichment workflow cannot proceed."""


def _resolve_dictionary_root(dictionary_dir: Path | str | None) -> Path:
    if dictionary_dir is None:
        return _DEFAULT_DICTIONARY_DIR
    return Path(dictionary_dir)


def _normalise_column_key(name: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in str(name)).strip("_")


def _lookup_column_name(frame: pd.DataFrame, *candidates: str) -> str | None:
    normalised = {_normalise_column_key(column): column for column in frame.columns}
    for candidate in candidates:
        key = _normalise_column_key(candidate)
        if key in normalised:
            return normalised[key]
    return None


@functools.lru_cache(maxsize=None)
def _load_assay_lookup_cached(dictionary_root: str) -> pd.DataFrame:
    root = Path(dictionary_root)
    candidate = root / "_assay" / "assay.csv"
    if not candidate.exists():
        raise AssayExtendedError(
            "assay.csv not found; expected at "
            f"'{candidate}'. Provide dictionary_dir pointing to the bundled dictionaries."
    )

    frame: pd.DataFrame | None = None
    errors: list[str] = []
    for sep in helpers.CSV_SEPARATORS:
        try:
            candidate_frame = helpers.read_csv_with_fallbacks(candidate, sep=sep)
        except Exception as exc:  # pragma: no cover - defensive
            errors.append(f"sep={sep!r}: {exc!s}")
            continue
        if "assay_chembl_id" in candidate_frame.columns:
            frame = candidate_frame
            break
        errors.append(
            "sep={}: available columns -> {}".format(
                repr(sep), ", ".join(candidate_frame.columns) or "<none>"
            )
        )
    if frame is None:
        detail = "; ".join(errors)
        message = "assay.csv missing required 'assay_chembl_id' column; unable to join assay metadata"
        if detail:
            message += f". {detail}"
        raise AssayExtendedError(message)

    expected = {
        "assay_chembl_id",
        "assay_group",
        "assay_strain",
        "assay_tax_id",
        "document_chembl_id",
        "target_chembl_id",
    }
    available = expected & set(frame.columns)
    subset = frame.loc[:, ["assay_chembl_id", *sorted(available - {"assay_chembl_id"})]].copy()
    subset = helpers.sort_power_query(subset, ("assay_chembl_id",))
    rename_map = {
        column: f"{column}_lookup"
        for column in subset.columns
        if column != "assay_chembl_id"
    }
    subset = subset.rename(columns=rename_map)
    for column in subset.columns:
        subset[column] = subset[column].astype("string")
    subset = subset.drop_duplicates(subset=["assay_chembl_id"], keep="first")
    return subset


def _load_assay_lookup(dictionary_root: Path) -> pd.DataFrame:
    return _load_assay_lookup_cached(str(dictionary_root.resolve()))


@functools.lru_cache(maxsize=None)
def _load_taxonomy_lookup_cached(dictionary_root: str) -> pd.DataFrame:
    root = Path(dictionary_root)
    candidate = root / "_taxonomy" / "taxonomy.csv"
    if not candidate.exists():
        logger.warning(
            "assay_extended_missing_taxonomy_dictionary",
            path=str(candidate),
            dictionary_root=str(root),
        )
        empty = pd.DataFrame(
            {
                "assay_tax_id": pd.Series(dtype="string"),
                "assay_group_taxonomy": pd.Series(dtype="string"),
                "assay_strain_taxonomy": pd.Series(dtype="string"),
            }
        )
        return empty

    frame = helpers.read_csv_with_fallbacks(candidate)
    aliases = {
        "tax_id": "assay_tax_id",
        "taxonomy_id": "assay_tax_id",
        "assay_tax_id": "assay_tax_id",
        "assay_group": "assay_group_taxonomy",
        "group": "assay_group_taxonomy",
        "assay_strain": "assay_strain_taxonomy",
        "strain": "assay_strain_taxonomy",
    }
    resolved: dict[str, str] = {}
    for column in frame.columns:
        key = _normalise_column_key(column)
        for alias, target in aliases.items():
            if key == _normalise_column_key(alias) and target not in resolved.values():
                resolved[column] = target
                break
    if "assay_tax_id" not in resolved.values():
        raise AssayExtendedError("taxonomy.csv missing 'assay_tax_id' column")

    renamed = frame.rename(columns=resolved)
    required = {"assay_tax_id", "assay_group_taxonomy", "assay_strain_taxonomy"}
    missing = [column for column in required if column not in renamed.columns and column != "assay_tax_id"]
    for column in missing:
        renamed[column] = pd.Series(dtype="string")
    subset = renamed.loc[:, ["assay_tax_id", "assay_group_taxonomy", "assay_strain_taxonomy"]].copy()
    for column in subset.columns:
        subset[column] = subset[column].astype("string")
    subset = helpers.sort_power_query(subset, ("assay_tax_id", "assay_group_taxonomy", "assay_strain_taxonomy"))
    subset = subset.drop_duplicates(subset=["assay_tax_id"], keep="first")
    return subset


def _load_taxonomy_lookup(dictionary_root: Path) -> pd.DataFrame:
    return _load_taxonomy_lookup_cached(str(dictionary_root.resolve()))


@functools.lru_cache(maxsize=None)
def _load_document_year_lookup_cached(dictionary_root: str) -> pd.DataFrame:
    root = Path(dictionary_root)
    candidate = root / "_document" / "document.csv"
    if not candidate.exists():
        raise AssayExtendedError(
            "document.csv not found; expected at "
            f"'{candidate}'. Provide dictionary_dir pointing to the bundled dictionaries."
    )

    frame: pd.DataFrame | None = None
    errors: list[str] = []
    for sep in helpers.CSV_SEPARATORS:
        try:
            candidate_frame = helpers.read_csv_with_fallbacks(candidate, sep=sep)
        except Exception as exc:  # pragma: no cover - defensive
            errors.append(f"sep={sep!r}: {exc!s}")
            continue
        column = _lookup_column_name(candidate_frame, "document_chembl_id", "chembl.document_chembl_id")
        if column is not None:
            frame = candidate_frame.rename(columns={column: "document_chembl_id"})
            break
        errors.append(
            "sep={}: available columns -> {}".format(
                repr(sep), ", ".join(candidate_frame.columns) or "<none>"
            )
        )
    if frame is None:
        detail = "; ".join(errors)
        raise AssayExtendedError(
            "document.csv missing required 'document_chembl_id' column"
            + (f". {detail}" if detail else "")
        )

    year_column = _lookup_column_name(
        frame,
        "year",
        "chembl.year",
        "pubmed.yearcompleted",
        "pubmed.year",
    )
    if year_column is None:
        raise AssayExtendedError("document.csv missing a year column")

    subset = frame.loc[:, ["document_chembl_id", year_column]].rename(columns={year_column: "year_document"})
    subset = helpers.sort_power_query(subset, ("document_chembl_id", "year_document"))
    subset = subset.drop_duplicates(subset=["document_chembl_id"], keep="first")
    subset["document_chembl_id"] = subset["document_chembl_id"].astype("string")
    subset["year_document"] = subset["year_document"].astype("string")
    return subset


def _load_document_year_lookup(dictionary_root: Path) -> pd.DataFrame:
    return _load_document_year_lookup_cached(str(dictionary_root.resolve()))


_LEGACY_TARGET_EXPORT_FILENAMES: tuple[str, ...] = (
    "output.target.csv",
    "output.targets.csv",
)


def _latest_target_export(target_dir: Path) -> Path:
    candidates = sorted(target_dir.glob("output.target_*.csv"))
    plural_candidates = sorted(target_dir.glob("output.targets_*.csv"))
    candidates = sorted({*candidates, *plural_candidates})
    if candidates:
        return candidates[-1]

    for legacy_name in _LEGACY_TARGET_EXPORT_FILENAMES:
        legacy_path = (target_dir / legacy_name).resolve()
        if legacy_path.exists():
            return legacy_path

    expected_names = "output.target_YYYYMMDD.csv"
    if _LEGACY_TARGET_EXPORT_FILENAMES:
        legacy_variants = ", ".join(_LEGACY_TARGET_EXPORT_FILENAMES)
        expected_names = f"{expected_names} or one of {{{legacy_variants}}}"

    raise AssayExtendedError(
        "No target exports matching "
        f"'{expected_names}' found in '{target_dir}'. Provide the target export in the dictionary bundle."
    )


def _parse_accessions(payload: object) -> list[str]:
    if payload is None:
        return []
    if helpers.is_missing_scalar(payload):
        return []
    text = str(payload).strip()
    if not text:
        return []
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        logger.warning("assay_extended_invalid_components", payload=text[:200])
        return []
    items: Iterable[object]
    if isinstance(data, dict):
        items = [data]
    elif isinstance(data, list):
        items = data
    else:
        return []
    tokens: list[str] = []
    for item in items:
        if not isinstance(item, dict):
            continue
        accession = str(item.get("accession", "")).strip()
        if accession and accession not in tokens:
            tokens.append(accession)
    return tokens


def _aggregate_accessions(frame: pd.DataFrame) -> pd.DataFrame:
    records: list[tuple[str, str]] = []
    for row in frame.itertuples(index=False):
        target_id = str(getattr(row, "target_chembl_id", "") or "").strip()
        if not target_id:
            continue
        for accession in _parse_accessions(getattr(row, "target_components", None)):
            records.append((target_id, accession))
    if not records:
        return pd.DataFrame({"target_chembl_id": pd.Series(dtype="string"), "accession_lookup": pd.Series(dtype="string")})
    accessions = pd.DataFrame(records, columns=["target_chembl_id", "accession_lookup"])
    accessions = helpers.sort_power_query(accessions, ("target_chembl_id", "accession_lookup"))
    accessions = accessions.drop_duplicates(subset=["target_chembl_id", "accession_lookup"], keep="first")

    def _join(series: pd.Series) -> str:
        tokens = [token for token in series.astype("string") if token and token != "<NA>"]
        return "|".join(dict.fromkeys(tokens))

    aggregated = accessions.groupby("target_chembl_id", as_index=False)["accession_lookup"].agg(_join)
    aggregated["target_chembl_id"] = aggregated["target_chembl_id"].astype("string")
    aggregated["accession_lookup"] = aggregated["accession_lookup"].astype("string")
    return aggregated


@functools.lru_cache(maxsize=None)
def _load_target_accession_lookup_cached(dictionary_root: str) -> pd.DataFrame:
    root = Path(dictionary_root)
    target_dir = root / "_target"
    if not target_dir.exists():
        raise AssayExtendedError(
            "Target dictionary directory not found; expected at "
            f"'{target_dir}'. Provide dictionary_dir pointing to the bundled dictionaries."
    )
    export_path = _latest_target_export(target_dir)
    frame = helpers.read_csv_with_fallbacks(export_path)
    if "target_chembl_id" not in frame.columns:
        raise AssayExtendedError(
            "Target export missing 'target_chembl_id' column; unable to join accessions"
        )
    if "target_components" not in frame.columns:
        logger.warning(
            "assay_extended_missing_target_components",
            path=str(export_path),
        )
        result = frame.loc[:, ["target_chembl_id"]].copy()
        result["accession_lookup"] = pd.Series(dtype="string")
        return helpers.sort_power_query(result, ("target_chembl_id",))
    subset = frame.loc[:, ["target_chembl_id", "target_components"]].copy()
    subset = helpers.sort_power_query(subset, ("target_chembl_id",))
    subset["target_chembl_id"] = subset["target_chembl_id"].astype("string")
    return _aggregate_accessions(subset)


def _load_target_accession_lookup(dictionary_root: Path) -> pd.DataFrame:
    return _load_target_accession_lookup_cached(str(dictionary_root.resolve()))


def _coalesce_strings(series_list: Sequence[pd.Series]) -> pd.Series:
    if not series_list:
        return pd.Series(dtype="string")
    result = pd.Series(pd.NA, index=series_list[0].index, dtype="string")
    for series in series_list:
        candidate = series.astype("string") if series is not None else None
        if candidate is None:
            continue
        cleaned = candidate.where(candidate.str.strip().ne(""), pd.NA)
        mask = result.isna()
        result = result.mask(~mask, result).mask(mask, cleaned)
    return result.astype("string")


def _coalesce_years(series_list: Sequence[pd.Series]) -> pd.Series:
    if not series_list:
        return pd.Series(dtype="Int64")
    result = pd.Series(pd.NA, index=series_list[0].index, dtype="Int64")
    for series in series_list:
        if series is None:
            continue
        numeric = pd.to_numeric(series, errors="coerce").astype("Int64")
        mask = result.isna()
        result = result.mask(~mask, result).mask(mask, numeric)
    return result


def enrich_assay_metadata(
    df: pd.DataFrame,
    *,
    dictionary_dir: Path | str | None = None,
) -> pd.DataFrame:
    """Return ``df`` enriched with dictionary-driven assay metadata."""

    dictionary_root = _resolve_dictionary_root(dictionary_dir)
    if not isinstance(dictionary_root, Path):  # pragma: no cover - defensive
        dictionary_root = Path(dictionary_root)

    result = df.copy()
    if result.empty:
        for column, dtype in (
            ("assay_group", "string"),
            ("assay_strain", "string"),
            ("accession", "string"),
        ):
            if column not in result.columns:
                result[column] = pd.Series(dtype=dtype)
        if "year" not in result.columns:
            result["year"] = pd.Series(dtype="Int64")
        return result

    index = result.index
    ensure_columns = {
        "assay_chembl_id": "string",
        "assay_tax_id": "string",
        "assay_group": "string",
        "assay_strain": "string",
        "target_chembl_id": "string",
        "document_chembl_id": "string",
    }
    for column, dtype in ensure_columns.items():
        if column not in result.columns:
            result[column] = pd.Series(pd.NA, index=index, dtype="string")
        else:
            result[column] = result[column].astype("string")

    assay_lookup = _load_assay_lookup(dictionary_root)
    taxonomy_lookup = _load_taxonomy_lookup(dictionary_root)
    document_lookup = _load_document_year_lookup(dictionary_root)
    target_lookup = _load_target_accession_lookup(dictionary_root)

    merged = result.merge(assay_lookup, on="assay_chembl_id", how="left")
    if "assay_tax_id_lookup" in merged.columns:
        merged["assay_tax_id"] = _coalesce_strings([merged["assay_tax_id"], merged["assay_tax_id_lookup"]])
    if "document_chembl_id_lookup" in merged.columns:
        merged["document_chembl_id"] = _coalesce_strings(
            [merged["document_chembl_id"], merged["document_chembl_id_lookup"]]
        )
    if "target_chembl_id_lookup" in merged.columns:
        merged["target_chembl_id"] = _coalesce_strings([
            merged["target_chembl_id"], merged["target_chembl_id_lookup"]
        ])

    merged = merged.merge(taxonomy_lookup, on="assay_tax_id", how="left")
    merged = merged.merge(document_lookup, on="document_chembl_id", how="left")
    merged = merged.merge(target_lookup, on="target_chembl_id", how="left")

    strain_candidates = [
        merged["assay_strain"],
        merged.get("assay_strain_taxonomy", pd.Series(pd.NA, index=index, dtype="string")),
        merged.get("assay_strain_lookup", pd.Series(pd.NA, index=index, dtype="string")),
    ]
    group_candidates = [
        merged["assay_group"],
        merged.get("assay_group_taxonomy", pd.Series(pd.NA, index=index, dtype="string")),
        merged.get("assay_group_lookup", pd.Series(pd.NA, index=index, dtype="string")),
    ]
    merged["assay_strain"] = _coalesce_strings(strain_candidates)
    merged["assay_group"] = _coalesce_strings(group_candidates)

    year_candidates = []
    if "year" in merged.columns:
        year_candidates.append(merged["year"])
    year_candidates.append(merged.get("year_document", pd.Series(pd.NA, index=index, dtype="string")))
    merged["year"] = _coalesce_years(year_candidates)

    accession_candidates = [
        merged.get("accession", pd.Series(pd.NA, index=index, dtype="string")),
        merged.get("accession_lookup", pd.Series(pd.NA, index=index, dtype="string")),
    ]
    merged["accession"] = _coalesce_strings(accession_candidates)

    drop_columns = [
        column
        for column in merged.columns
        if column.endswith("_lookup") or column.endswith("_taxonomy") or column == "year_document"
    ]
    merged = merged.drop(columns=[column for column in drop_columns if column in merged.columns])

    merged["assay_group"] = merged["assay_group"].astype("string")
    merged["assay_strain"] = merged["assay_strain"].astype("string")
    merged["accession"] = merged["accession"].astype("string")
    merged["year"] = merged["year"].astype("Int64")

    return merged

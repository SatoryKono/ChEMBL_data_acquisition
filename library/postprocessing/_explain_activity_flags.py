"""Explain provenance for activity flag columns.

The functions in this module mirror the production logic implemented in
``activity_extended`` and annotate each flag with a textual trace describing the
rule path that produced its value.  They are intentionally verbose to aid QA
and parity validation against the historical Power Query workbook.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import pandas as pd

from .activity_extended import (
    _DEFAULT_DICTIONARY_DIR,
    _GROUP_KEY_COLUMNS,
    _annotate_high_citation,
    _compute_citation_flags,
    _load_citation_fraction,
    _safe_to_bool,
    _safe_to_int,
)


def _coerce_boolean(value: object, column: str) -> object:
    series = pd.Series([value])
    coerced = _safe_to_bool(series, column).iloc[0]
    if coerced is pd.NA:
        return pd.NA
    return bool(coerced)


def _format_bool(value: object) -> str:
    if value is pd.NA or pd.isna(value):
        return "<NA>"
    return "True" if bool(value) else "False"


def _boolean_from_sources(
    df: pd.DataFrame,
    target: str,
    sources: Sequence[str],
    *,
    default: bool = False,
) -> tuple[pd.Series, pd.Series]:
    index = df.index
    resolved: list[bool | pd.NA] = []
    explanations: list[str] = []

    for row_idx in index:
        source_column: str | None = None
        raw_value: object | None = None
        coerced: object | None = None
        for candidate in sources:
            if candidate not in df.columns:
                continue
            value = df.at[row_idx, candidate]
            if pd.isna(value):
                continue
            source_column = candidate
            raw_value = value
            coerced = _coerce_boolean(value, candidate)
            break

        if source_column is None:
            resolved.append(default)
            explanations.append(
                f"{target} defaulted to {default} because none of {', '.join(sources)} provided data"
            )
            continue

        if coerced is pd.NA:
            resolved.append(default)
            explanations.append(
                f"{source_column}={raw_value!r} -> not coercible; defaulted to {default}"
            )
            continue

        resolved.append(bool(coerced))
        explanations.append(
            f"{source_column}={raw_value!r} -> {_format_bool(coerced)}"
        )

    result_series = pd.Series(resolved, index=index, dtype="boolean")
    explanation_series = pd.Series(
        [f"{explanations[i]} => {target}={_format_bool(resolved[i])}" for i in range(len(explanations))],
        index=index,
        dtype="string",
    )
    return result_series, explanation_series


def explain_unknown_chirality(df: pd.DataFrame) -> pd.DataFrame:
    index = df.index
    if "nstereo" in df.columns:
        nstereo = _safe_to_int(df["nstereo"], "nstereo")
        values: list[bool] = []
        reasons: list[str] = []
        for raw in nstereo:
            if pd.isna(raw):
                values.append(True)
                reasons.append("nstereo missing -> default True to mirror Power Query")
            elif int(raw) == 1:
                values.append(False)
                reasons.append("nstereo == 1 -> chirality known -> False")
            else:
                values.append(True)
                reasons.append(f"nstereo == {int(raw)} != 1 -> True")
        series = pd.Series(values, index=index, dtype="boolean")
        explanations = pd.Series(
            [f"{reasons[i]} => unknown_chirality={_format_bool(values[i])}" for i in range(len(reasons))],
            index=index,
            dtype="string",
        )
        return pd.DataFrame(
            {
                "unknown_chirality": series,
                "__explain_unknown_chirality": explanations,
            }
        )

    if "unknown_chirality" in df.columns:
        series = _safe_to_bool(df["unknown_chirality"], "unknown_chirality")
        explanations = pd.Series(
            [
                f"unknown_chirality column provided -> {_format_bool(value)}"
                for value in series
            ],
            index=index,
            dtype="string",
        )
        return pd.DataFrame(
            {
                "unknown_chirality": series,
                "__explain_unknown_chirality": explanations,
            }
        )

    series = pd.Series(True, index=index, dtype="boolean")
    explanations = pd.Series(
        "nstereo column absent -> forced True",
        index=index,
        dtype="string",
    )
    return pd.DataFrame(
        {
            "unknown_chirality": series,
            "__explain_unknown_chirality": explanations,
        }
    )


def explain_multmol_assay(df: pd.DataFrame) -> pd.DataFrame:
    index = df.index
    if "multmol_assay" in df.columns:
        raw_series = df["multmol_assay"]
    else:
        raw_series = pd.Series(pd.NA, index=index)
    source_bool = _safe_to_bool(raw_series, "multmol_assay")
    source_bool = source_bool.astype("boolean")

    counts = (
        df.groupby(list(_GROUP_KEY_COLUMNS), dropna=False)
        .size()
        .rename("Count")
        .reset_index()
    )
    work = df.copy()
    work["__row_id"] = range(len(work))
    merged = work.merge(counts, on=list(_GROUP_KEY_COLUMNS), how="left")
    merged = merged.sort_values("__row_id")
    merged.index = index

    unknown = (
        _safe_to_bool(df.get("unknown_chirality", pd.Series(pd.NA, index=index)), "unknown_chirality")
        if "unknown_chirality" in df.columns
        else pd.Series(pd.NA, index=index, dtype="boolean")
    )

    duplicate_mask = (
        unknown.fillna(True).eq(False)
        & raw_series.isna()
        & merged["Count"].fillna(0).gt(1)
    )
    duplicated_assays = set(
        merged.loc[duplicate_mask, "assay_chembl_id"].astype("string").dropna()
    )
    derived = merged["assay_chembl_id"].astype("string").isin(duplicated_assays)

    final_series = _safe_to_bool(source_bool.fillna(False) | derived, "multmol_assay")
    explanations: list[str] = []
    for idx in index:
        raw_value = raw_series.at[idx] if idx in raw_series.index else pd.NA
        base = source_bool.at[idx] if idx in source_bool.index else pd.NA
        duplicates = bool(derived.loc[idx]) if idx in derived.index else False
        count_value = merged.loc[idx, "Count"] if idx in merged.index else pd.NA
        pieces: list[str] = []
        if pd.isna(raw_value):
            pieces.append("source missing -> defaults to False")
        else:
            pieces.append(f"raw {raw_value!r} -> {_format_bool(base)}")
        if duplicates:
            pieces.append(
                f"duplicate group Count={int(count_value)} with known chirality triggers True"
            )
        else:
            pieces.append(
                f"duplicate rule not triggered (Count={count_value if not pd.isna(count_value) else 'NA'}, unknown={_format_bool(unknown.at[idx] if idx in unknown.index else pd.NA)})"
            )
        final_value = final_series.at[idx]
        explanations.append("; ".join(pieces) + f" => multmol_assay={_format_bool(final_value)}")

    return pd.DataFrame(
        {
            "multmol_assay": final_series,
            "__explain_multmol_assay": pd.Series(explanations, index=index, dtype="string"),
        }
    )


def explain_exact_data_citation(df: pd.DataFrame) -> pd.DataFrame:
    series, explanations = _boolean_from_sources(
        df,
        "exact_data_citation",
        ("exact_cited_activity", "exact_data_citation"),
    )
    return pd.DataFrame(
        {
            "exact_data_citation": series,
            "__explain_exact_data_citation": explanations,
        }
    )


def explain_higly_correlated_assay(df: pd.DataFrame) -> pd.DataFrame:
    series, explanations = _boolean_from_sources(
        df,
        "higly_correlated_assay",
        ("higly_correlated_cit", "higly_correlated_assay", "highly_correlated_assay"),
    )
    alias = series.copy()
    alias_explain = pd.Series(
        [exp.replace("higly_correlated_assay", "highly_correlated_assay") for exp in explanations],
        index=series.index,
        dtype="string",
    )
    return pd.DataFrame(
        {
            "higly_correlated_assay": series,
            "__explain_higly_correlated_assay": explanations,
            "highly_correlated_assay": alias,
            "__explain_highly_correlated_assay": alias_explain,
        }
    )


def explain_shuffled_assay(df: pd.DataFrame) -> pd.DataFrame:
    series, explanations = _boolean_from_sources(
        df,
        "shuffled_assay",
        ("shuffled_cit", "shuffled_assay"),
    )
    return pd.DataFrame(
        {
            "shuffled_assay": series,
            "__explain_shuffled_assay": explanations,
        }
    )


def explain_review(df: pd.DataFrame) -> pd.DataFrame:
    series, explanations = _boolean_from_sources(
        df,
        "review",
        ("review_doc", "review"),
    )
    return pd.DataFrame(
        {
            "review": series,
            "__explain_review": explanations,
        }
    )


def explain_rounded_data_citation(df: pd.DataFrame) -> pd.DataFrame:
    series, explanations = _boolean_from_sources(
        df,
        "rounded_data_citation",
        ("approx_cited_activity", "rounded_data_citation"),
    )
    return pd.DataFrame(
        {
            "rounded_data_citation": series,
            "__explain_rounded_data_citation": explanations,
        }
    )


def explain_high_citation_rate(
    df: pd.DataFrame,
    *,
    dictionary_root: Path | None = None,
) -> pd.DataFrame:
    root = Path(dictionary_root) if dictionary_root is not None else _DEFAULT_DICTIONARY_DIR
    flagged = _compute_citation_flags(df)
    annotated = _annotate_high_citation(flagged, root)

    grouped = (
        flagged.groupby("document_chembl_id")["is_citation"]
        .agg(n_citation="sum", n_non_citation=lambda s: (~s).sum())
        .reset_index()
    )
    grouped["N"] = grouped["n_citation"] + grouped["n_non_citation"]
    valid = grouped[(grouped["n_citation"] > 0) & (grouped["n_non_citation"] > 0)]
    citation_fraction = _load_citation_fraction(root)
    valid = valid.merge(
        citation_fraction[["N", "K_min_significant"]],
        on="N",
        how="left",
    )
    stats_map = valid.set_index("document_chembl_id")
    grouped_map = grouped.set_index("document_chembl_id")

    explanations: list[str] = []
    for row in annotated.itertuples(index=False):
        doc_id = getattr(row, "document_chembl_id")
        current = getattr(row, "high_citation_rate")
        if doc_id in stats_map.index:
            stats = stats_map.loc[doc_id]
            threshold = stats.get("K_min_significant")
            n_cit = int(stats["n_citation"])
            n_non = int(stats["n_non_citation"])
            total = int(stats["N"])
            if pd.notna(threshold):
                condition = "met" if n_cit >= int(threshold) else "not met"
                explanations.append(
                    f"doc {doc_id}: citations={n_cit}, non_citations={n_non}, N={total}, threshold={int(threshold)} -> {condition} => high_citation_rate={_format_bool(current)}"
                )
            else:
                explanations.append(
                    f"doc {doc_id}: threshold missing for N={total} -> default False => high_citation_rate={_format_bool(current)}"
                )
        else:
            if doc_id in grouped_map.index:
                stats = grouped_map.loc[doc_id]
                n_cit = int(stats["n_citation"])
                n_non = int(stats["n_non_citation"])
                explanations.append(
                    f"doc {doc_id}: citations={n_cit}, non_citations={n_non} -> insufficient mix => high_citation_rate={_format_bool(current)}"
                )
            else:
                explanations.append(
                    f"doc {doc_id}: no citation statistics available -> default False => high_citation_rate={_format_bool(current)}"
                )

    return pd.DataFrame(
        {
            "high_citation_rate": annotated["high_citation_rate"].astype("boolean"),
            "__explain_high_citation_rate": pd.Series(explanations, index=annotated.index, dtype="string"),
        }
    )


def explain_original_activity_flags(df: pd.DataFrame) -> pd.DataFrame:
    index = df.index
    approx_series = (
        df["original_activity_approx"].astype("string")
        if "original_activity_approx" in df.columns
        else pd.Series(pd.NA, index=index, dtype="string")
    )
    exact_series = (
        df["original_activity_exact"].astype("string")
        if "original_activity_exact" in df.columns
        else pd.Series(pd.NA, index=index, dtype="string")
    )

    approx_explain = [
        (
            "original_activity_approx missing -> <NA>"
            if pd.isna(value)
            else f"original_activity_approx preserved literal {value!r}"
        )
        for value in approx_series
    ]
    exact_explain = [
        (
            "original_activity_exact missing -> <NA>"
            if pd.isna(value)
            else f"original_activity_exact preserved literal {value!r}"
        )
        for value in exact_series
    ]

    return pd.DataFrame(
        {
            "original_activity_approx": approx_series,
            "__explain_original_activity_approx": pd.Series(approx_explain, index=index, dtype="string"),
            "original_activity_exact": exact_series,
            "__explain_original_activity_exact": pd.Series(exact_explain, index=index, dtype="string"),
        }
    )


__all__ = [
    "explain_unknown_chirality",
    "explain_multmol_assay",
    "explain_exact_data_citation",
    "explain_higly_correlated_assay",
    "explain_shuffled_assay",
    "explain_review",
    "explain_rounded_data_citation",
    "explain_high_citation_rate",
    "explain_original_activity_flags",
]

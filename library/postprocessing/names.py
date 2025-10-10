"""Target names post-processing mirroring the legacy Power Query workbook.

The ``SSOT_CONTEXT`` constant holds the original Power Query script used as the
single source of truth for these transformations so that future refactors can
validate behaviour against the historical pipeline.
"""

from __future__ import annotations

import json
import re
import sys
from collections.abc import Iterable, Mapping, MutableMapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from library.common.csv_utils import write_csv_deterministic

from ..common.log import logger
from . import helpers
from .helpers import (
    ENCODING_FALLBACKS,
    normalise_export_basename,
    normalise_text,
    read_csv_with_fallbacks,
)

# ruff: noqa: E501
SSOT_CONTEXT = """
The implementation follows the Single Source of Truth (SSoT) captured in the
Power Query (M) script that shipped with the historical Excel workbook. The
script is reproduced below so future refactors can be validated against the
original transformation pipeline:

contrion=()=>let
    Source = #"Target Aggregated",
    #"Expanded target_components" = Table.ExpandListColumn(Source, "target_components"),
    #"Expanded target_components1" = Table.ExpandRecordColumn(#"Expanded target_components", "target_components", {"molecule_chembl_id", "component_id", "component_type", "component_name", "accession", "organism", "sequence", "description", "component_synonyms", "molecule", "molecule_pref_name", "molecule_synonyms", "molecule_structures", "molecule_properties", "active_component", "active_component_type", "salt_chembl_id", "relationships", "contrion"}, {"molecule_chembl_id", "component_id", "component_type", "component_name", "accession", "organism", "sequence", "description", "component_synonyms", "molecule", "molecule_pref_name", "molecule_synonyms", "molecule_structures", "molecule_properties", "active_component", "active_component_type", "salt_chembl_id", "relationships", "contrion"}),
    #"Expanded component_synonyms" = Table.ExpandTableColumn(#"Expanded target_components1", "component_synonyms", {"synonym"}, {"component_synonym"}),
    #"Expanded molecule_synonyms" = Table.ExpandTableColumn(#"Expanded component_synonyms", "molecule_synonyms", {"synonyms"}, {"molecule_synonym"}),
    #"Expanded relationships" = Table.ExpandTableColumn(#"Expanded molecule_synonyms", "relationships", {"relationship", "related_component"}, {"relationship", "related_component"}),
    #"Added Is Hydrate" = Table.AddColumn(#"Expanded relationships", "is_hydrate", each is_hidrate([component_synonym])),
    #"Removed Hydrate" = Table.TransformColumns(#"Added Is Hydrate", {{"component_synonym", remove_hidrate}}),
    #"Sorted Names" = Table.Sort(#"Removed Hydrate", {{"target_chembl_id", Order.Ascending}, {"molecule_chembl_id", Order.Ascending}, {"component_synonym", Order.Ascending}}),
    #"Grouped Rows" = Table.Group(#"Sorted Names", {"target_chembl_id", "molecule_chembl_id"}, {{"All Names", each sort_my_list([component_synonym]), type text}, {"Active Component", each Text.FirstN(Text.Combine(List.RemoveNulls([active_component]), "|"), 32766), type text}, {"Active Component Type", each Text.FirstN(Text.Combine(List.RemoveNulls([active_component_type]), "|"), 32766), type text}, {"Contrion", each sort_my_list([contrion]), type text}}),
    #"Merged Molecule" = Table.NestedJoin(#"Grouped Rows", {"molecule_chembl_id"}, Table.Distinct(#"Expanded relationships"[["molecule_chembl_id"], ["molecule"]]), {"molecule_chembl_id"}, "molecule", JoinKind.LeftOuter),
    #"Expanded molecule" = Table.ExpandRecordColumn(#"Merged Molecule", "molecule", {"molecule_structures", "molecule_properties", "unknown_chirality"}, {"molecule_structures", "molecule_properties", "unknown_chirality"}),
    #"Expanded molecule_structures" = Table.ExpandRecordColumn(#"Expanded molecule", "molecule_structures", {"canonical_smiles", "standard_inchi_key", "standard_inchi", "standard_inchi_stereo"}, {"canonical_smiles", "standard_inchi_key", "standard_inchi", "standard_inchi_stereo"}),
    #"Expanded molecule_properties" = Table.ExpandRecordColumn(#"Expanded molecule_structures", "molecule_properties", {"mw_freebase"}, {"mw_freebase"}),
    #"Added Reference SMILES" = Table.AddColumn(#"Expanded molecule_properties", "canonical_smiles_reference", each reference_SMILES([molecule_chembl_id])),
    #"Replaced Null SMILES" = Table.ReplaceValue(#"Added Reference SMILES", null, each [canonical_smiles_reference], Replacer.ReplaceValue, {"canonical_smiles"}),
    #"Removed Helper" = Table.RemoveColumns(#"Replaced Null SMILES", {"canonical_smiles_reference"}),
    #"Filled Down" = Table.FillDown(#"Removed Helper", {"mw_freebase", "standard_inchi_key", "standard_inchi", "standard_inchi_stereo", "unknown_chirality"}),
    #"Changed Type" = Table.TransformColumnTypes(#"Filled Down", {{"mw_freebase", type number}, {"unknown_chirality", type text}}),
    #"Sorted Rows" = Table.Sort(#"Changed Type", {{"target_chembl_id", Order.Ascending}, {"molecule_chembl_id", Order.Ascending}})
in
    #"Sorted Rows"

Helpers for generating target name exports. The helper mirrors the legacy Power
Query workbook that produced auxiliary name tables for reporting. It extracts
textual identifiers from the merged targets export, normalises them and emits a
long-form table where each row represents a distinct name attributed to a
target.
"""

__all__ = [
    "SSOT_CONTEXT",
    "process_target_names",
    "is_hidrate",
    "remove_hidrate",
    "sort_my_list",
    "reference_SMILES",
    "write_csv_deterministic",
]

_DEFAULT_SEARCH_DIR = Path("data/output")
_OUTPUT_COLUMNS: tuple[str, ...] = (
    "target_chembl_id",
    "molecule_chembl_id",
    "all_names",
    "canonical_smiles",
    "mw_freebase",
    "unknown_chirality",
    "salt_chembl_id",
    "standard_inchi_key",
    "standard_inchi_skeleton",
    "standard_inchi_stereo",
    "is_hydrate",
    "is_saltform",
    "active_component",
    "contrion",
    "active_component_type",
)


class TargetNamesError(RuntimeError):
    """Raised when the target names post-processing cannot proceed."""


def _current_default_search_dir() -> Path:
    package = sys.modules.get(__name__)
    if package is not None and hasattr(package, "_DEFAULT_SEARCH_DIR"):
        override = package._DEFAULT_SEARCH_DIR
        if override is not None:
            return Path(override)
    return _DEFAULT_SEARCH_DIR


def _matches_target_filename(name: str) -> bool:
    return bool(re.match(r"output\.target_\d{8}\.csv\Z", name))


def _latest_target_file(search_dir: Path) -> Path:
    candidates = sorted(
        (
            path
            for path in search_dir.iterdir()
            if path.is_file() and _matches_target_filename(path.name)
        ),
        key=lambda item: item.name,
    )
    if not candidates:
        raise TargetNamesError(
            f"No target exports matching 'output.target_YYYYMMDD.csv' found in {search_dir!s}"
        )
    return candidates[-1]


def is_hidrate(value: object | None) -> bool:
    """Return ``True`` when ``value`` contains a hydrate annotation."""

    text = normalise_text(value)
    if not text:
        return False
    lowered = text.lower()
    return any(
        token in lowered
        for token in (
            "hydrate",
            "hemihydrate",
            "sesquihydrate",
            "trihydrate",
            "dihydrate",
        )
    )


def remove_hidrate(value: object | None) -> str:
    """Strip hydrate-related qualifiers from ``value``."""

    text = normalise_text(value)
    if not text:
        return ""
    cleaned = re.sub(r"\b[\w\-]*hydrate\b", "", text, flags=re.IGNORECASE)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    return cleaned


def sort_my_list(values: Iterable[object] | object | None) -> str:
    """Return a deterministic pipe-delimited list derived from ``values``."""

    tokens: list[str] = []
    if values is None:
        return ""
    if isinstance(values, str):
        values = values.split("|")
    for item in values:
        text = normalise_text(item)
        if text and text not in tokens:
            tokens.append(text)
    decorated = [(token.lower(), idx, token) for idx, token in enumerate(tokens)]
    decorated.sort(key=lambda item: (item[0], item[1]))
    return "|".join(item[2] for item in decorated)


@dataclass
class MoleculeStructures:
    canonical_smiles: str | None = None
    standard_inchi_key: str | None = None
    standard_inchi: str | None = None
    standard_inchi_stereo: str | None = None


@dataclass
class MoleculeProperties:
    mw_freebase: str | None = None


def _coerce_dict(value: Any) -> dict[str, Any]:
    if isinstance(value, dict):
        return value
    return {}


def _coerce_list(value: Any) -> list[Any]:
    if isinstance(value, list):
        return value
    if value in (None, "", float("nan")):
        return []
    return [value]


def _parse_components(value: Any) -> list[dict[str, Any]]:
    text = normalise_text(value)
    if not text:
        return []
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        try:
            data = json.loads(f"[{text}]")
        except json.JSONDecodeError as exc:
            raise TargetNamesError(
                f"Unable to parse target_components payload: {text!r}"
            ) from exc
    if isinstance(data, dict):
        return [data]
    if isinstance(data, list):
        return [item for item in data if isinstance(item, dict)]
    raise TargetNamesError("Parsed target_components payload is not a list of records")


def _extract_synonyms(component: Mapping[str, Any]) -> list[str]:
    names: list[str] = []
    for key in ("component_synonyms", "molecule_synonyms", "synonyms"):
        entries = component.get(key)
        for entry in _coerce_list(entries):
            if isinstance(entry, str):
                text = normalise_text(entry)
                if text:
                    names.append(text)
            elif isinstance(entry, Mapping):
                text = normalise_text(
                    entry.get("synonym") or entry.get("name") or entry.get("synonyms")
                )
                if text:
                    names.append(text)
    return names


def _extract_contrion(component: Mapping[str, Any]) -> list[str]:
    contrion_value = component.get("contrion")
    if isinstance(contrion_value, str):
        return [contrion_value]
    if isinstance(contrion_value, Iterable) and not isinstance(
        contrion_value, str | bytes
    ):
        result: list[str] = []
        for item in contrion_value:
            text = normalise_text(item)
            if text:
                result.append(text)
        return result
    return []


def _extract_structures(component: Mapping[str, Any]) -> MoleculeStructures:
    structures = _coerce_dict(
        component.get("molecule_structures") or component.get("structures")
    )
    return MoleculeStructures(
        canonical_smiles=normalise_text(structures.get("canonical_smiles")),
        standard_inchi_key=normalise_text(structures.get("standard_inchi_key")),
        standard_inchi=normalise_text(structures.get("standard_inchi")),
        standard_inchi_stereo=normalise_text(structures.get("standard_inchi_stereo")),
    )


def _extract_properties(component: Mapping[str, Any]) -> MoleculeProperties:
    properties = _coerce_dict(
        component.get("molecule_properties") or component.get("properties")
    )
    value = properties.get("mw_freebase")
    if value is None:
        return MoleculeProperties()
    return MoleculeProperties(mw_freebase=str(value))


def _extract_unknown_chirality(component: Mapping[str, Any]) -> str:
    value = component.get("unknown_chirality")
    if isinstance(value, bool):
        return "Y" if value else "N"
    text = normalise_text(value)
    if not text:
        return ""
    if text.lower() in {"true", "yes", "1"}:
        return "Y"
    if text.lower() in {"false", "no", "0"}:
        return "N"
    return text


def _salt_identifier(component: Mapping[str, Any]) -> str:
    salt = component.get("salt_chembl_id")
    if salt:
        return normalise_text(salt)
    hierarchy = _coerce_dict(component.get("molecule") or {}).get("molecule_hierarchy")
    if isinstance(hierarchy, Mapping):
        salt_candidate = hierarchy.get("parent_molecule_chembl_id")
        return normalise_text(salt_candidate)
    return ""


def _component_identifier(component: Mapping[str, Any]) -> str:
    for key in ("molecule_chembl_id", "component_chembl_id", "molecule_id"):
        text = normalise_text(component.get(key))
        if text:
            return text
    molecule = _coerce_dict(component.get("molecule"))
    text = normalise_text(
        molecule.get("molecule_chembl_id") or molecule.get("chembl_id")
    )
    if text:
        return text
    raise TargetNamesError("Component record is missing molecule_chembl_id")


def _component_name(component: Mapping[str, Any]) -> str:
    text = normalise_text(
        component.get("component_name")
        or component.get("component_pref_name")
        or component.get("molecule_pref_name")
        or component.get("name")
    )
    return text


def _active_component(component: Mapping[str, Any]) -> str:
    text = normalise_text(component.get("active_component"))
    if text:
        return text
    name = _component_name(component)
    return remove_hidrate(name)


def _active_component_type(component: Mapping[str, Any]) -> str:
    text = normalise_text(
        component.get("active_component_type") or component.get("component_type")
    )
    return text


_REFERENCE_CACHE: MutableMapping[Path, pd.Series] = {}


def reference_SMILES(
    molecule_id: object | None,
    *,
    overrides: Mapping[str, str] | None = None,
    reference_path: str | Path | None = None,
) -> str | None:
    text = normalise_text(molecule_id)
    if not text:
        return None
    if overrides and text in overrides:
        return normalise_text(overrides[text]) or None
    path = Path(reference_path or Path("data/reference/Table6.csv"))
    if not path.exists():
        raise TargetNamesError(f"Reference SMILES table not found at {path!s}")
    if path not in _REFERENCE_CACHE:
        frame = read_csv_with_fallbacks(path, encodings=ENCODING_FALLBACKS)
        required = {"molecule_chembl_id", "canonical_smiles"}
        missing = required - set(frame.columns)
        if missing:
            raise TargetNamesError(
                f"Reference table {path!s} is missing required columns: {', '.join(sorted(missing))}"
            )
        series = (
            frame.set_index("molecule_chembl_id")["canonical_smiles"]
            .astype("string")
            .map(normalise_text)
        )
        _REFERENCE_CACHE[path] = series
    series = _REFERENCE_CACHE[path]
    value = series.get(text)
    if value is None or value == "":
        return None
    return str(value)


def _component_rows(frame: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in frame.to_dict("records"):
        target_id = normalise_text(record.get("target_chembl_id"))
        components = _parse_components(record.get("target_components"))
        for component in components:
            molecule_id = _component_identifier(component)
            structures = _extract_structures(component)
            properties = _extract_properties(component)
            canonical_smiles = structures.canonical_smiles
            if not canonical_smiles:
                canonical_smiles = reference_SMILES(molecule_id)
            salt_id = _salt_identifier(component)
            synonyms = _extract_synonyms(component)
            component_name = _component_name(component)
            if component_name:
                synonyms.insert(0, component_name)
            all_names = sort_my_list(synonyms)
            hydrate_flag = "Y" if is_hidrate(all_names) else "N"
            rows.append(
                {
                    "target_chembl_id": target_id,
                    "molecule_chembl_id": molecule_id,
                    "all_names": all_names,
                    "canonical_smiles": canonical_smiles or "",
                    "mw_freebase": properties.mw_freebase or "",
                    "unknown_chirality": _extract_unknown_chirality(component),
                    "salt_chembl_id": salt_id,
                    "standard_inchi_key": structures.standard_inchi_key or "",
                    "standard_inchi_skeleton": structures.standard_inchi or "",
                    "standard_inchi_stereo": structures.standard_inchi_stereo or "",
                    "is_hydrate": hydrate_flag,
                    "is_saltform": "Y" if salt_id and salt_id != molecule_id else "N",
                    "active_component": _active_component(component),
                    "contrion": sort_my_list(_extract_contrion(component)),
                    "active_component_type": _active_component_type(component),
                }
            )
    return rows


def _ensure_columns(frame: pd.DataFrame, columns: Sequence[str]) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise TargetNamesError(
            "Input target table is missing required columns: " + ", ".join(missing)
        )


# Columns in the targets table used to derive names.  Each entry maps the column
# name to a descriptive label stored alongside the extracted token so consumers
# can identify the provenance of a particular name.
NAME_COLUMN_SOURCES: tuple[tuple[str, str], ...] = (
    ("pref_name", "chembl_preferred"),
    ("protein_name_canonical", "uniprot_canonical"),
    ("protein_name_alt", "uniprot_alternative"),
    ("synonyms", "chembl_synonym"),
    ("protein_synonym_list", "uniprot_synonym"),
    ("gtop_synonyms", "gtop_synonym"),
    ("gene_symbol", "gene_symbol"),
    ("gene_symbol_list", "gene_symbol"),
    ("isoform_names", "isoform_name"),
    ("isoform_synonyms", "isoform_synonym"),
    ("recommendedName", "uniprot_recommended"),
    ("secondaryAccessionNames", "uniprot_secondary"),
)

# Columns that should be split on the pipe delimiter.  Other columns are treated
# as single textual values even if a pipe character is present.
PIPE_SPLIT_COLUMNS: frozenset[str] = frozenset(
    {
        "synonyms",
        "protein_synonym_list",
        "gtop_synonyms",
        "gene_symbol_list",
        "isoform_names",
        "isoform_synonyms",
        "secondaryAccessionNames",
    }
)

# Canonical column order for the emitted names table.
TARGET_NAMES_COLUMNS: tuple[str, ...] = (
    "target_chembl_id",
    "uniprot_id_primary",
    "name",
    "name_type",
    "source_column",
)

# Sentinel values considered as missing tokens when splitting pipe-delimited
# strings.  The list mirrors guards in the legacy workbook.
EMPTY_TOKEN_MARKERS: frozenset[str] = frozenset({"", "-", "n/a", "none"})


def _iter_name_tokens(value: object, *, split: bool) -> Iterable[str]:
    """Yield normalised name tokens extracted from ``value``."""

    text = helpers.normalise_text(value)
    if not text:
        return []

    if not split:
        return [text]

    tokens: list[str] = []
    for raw in text.split("|"):
        token = helpers.normalise_text(raw)
        if not token:
            continue
        lower = token.lower()
        if lower in EMPTY_TOKEN_MARKERS:
            continue
        tokens.append(token)
    return tokens


def _build_names_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Project ``frame`` onto the long-form target names table."""

    if frame.empty:
        return pd.DataFrame(columns=TARGET_NAMES_COLUMNS, dtype="string")

    working = frame.copy()
    working["target_chembl_id"] = working["target_chembl_id"].map(
        helpers.normalise_text
    )
    if "uniprot_id_primary" in working.columns:
        uniprot_ids = working["uniprot_id_primary"].map(helpers.normalise_text)
    else:
        uniprot_ids = pd.Series("", index=working.index, dtype="string")
    working["uniprot_id_primary"] = uniprot_ids

    working = working[working["target_chembl_id"] != ""].copy()
    if working.empty:
        return pd.DataFrame(columns=TARGET_NAMES_COLUMNS, dtype="string")

    column_set = set(working.columns)
    extracted_frames: list[pd.DataFrame] = []

    for column, label in NAME_COLUMN_SOURCES:
        if column not in column_set:
            continue

        series = working[column].map(helpers.normalise_text)
        if column in PIPE_SPLIT_COLUMNS:
            tokens = series.str.split("|").explode().map(helpers.normalise_text)
            if tokens.empty:
                continue
            tokens = tokens[tokens != ""]
            if tokens.empty:
                continue
            tokens = tokens[~tokens.str.lower().isin(EMPTY_TOKEN_MARKERS)]
        else:
            tokens = series[series != ""]

        if tokens.empty:
            continue

        column_frame = pd.DataFrame(
            {
                "target_chembl_id": working.loc[tokens.index, "target_chembl_id"],
                "uniprot_id_primary": working.loc[tokens.index, "uniprot_id_primary"],
                "name": tokens.astype("string"),
                "name_type": label,
                "source_column": column,
            }
        )
        extracted_frames.append(column_frame)

    if not extracted_frames:
        return pd.DataFrame(columns=TARGET_NAMES_COLUMNS, dtype="string")

    names_df = pd.concat(extracted_frames, ignore_index=True)
    names_df = helpers.ensure_string_columns(names_df, TARGET_NAMES_COLUMNS)
    names_df = names_df.drop_duplicates().sort_values(
        by=["target_chembl_id", "name", "name_type"], kind="mergesort"
    )
    names_df = names_df.reset_index(drop=True)
    return names_df


def _summarise_contrion(series: pd.Series | None) -> dict[str, int]:
    """Return counts for the ``contrion`` column if present."""

    if series is None:
        return {"contrion_unique": 0, "contrion_non_empty": 0, "contrion_total": 0}

    tokens: set[str] = set()
    non_empty = 0
    total = 0
    for value in series.dropna().astype(str):
        text = helpers.normalise_text(value)
        if not text or text in EMPTY_TOKEN_MARKERS:
            continue
        parts = [helpers.normalise_text(part) for part in text.split("|")]
        parts = [part for part in parts if part and part not in EMPTY_TOKEN_MARKERS]
        if not parts:
            continue
        non_empty += 1
        total += len(parts)
        tokens.update(parts)
    return {
        "contrion_unique": len(tokens),
        "contrion_non_empty": non_empty,
        "contrion_total": total,
    }


def _summarise_active_component_type(series: pd.Series | None) -> dict[str, int]:
    """Return counts per ``active_component_type`` value."""

    if series is None:
        return {}

    normalised = (
        series.fillna("<NA>")
        .astype(str)
        .map(lambda value: helpers.normalise_text(value) or "<EMPTY>")
    )
    counts = normalised.value_counts(dropna=False)
    return {str(key): int(value) for key, value in counts.items()}


def process_target_names(
    input_path: str | Path, *, verbose: bool = False
) -> dict[str, Any]:
    """Process ``input_path`` and emit the target names table."""

    source_path = Path(input_path)
    frame = helpers.read_csv_with_fallbacks(source_path)
    frame = helpers.ensure_string_columns(frame, frame.columns)

    names_df = _build_names_table(frame)
    base = normalise_export_basename(source_path)
    output_path = source_path.with_name(f"names.{base}").with_suffix(".csv")
    helpers.write_csv(names_df, output_path, columns=TARGET_NAMES_COLUMNS)

    summary: dict[str, Any] = {
        "rows_before": int(len(frame)),
        "rows_after": int(len(names_df)),
    }
    contrion_summary = _summarise_contrion(frame.get("contrion"))
    summary.update(contrion_summary)
    active_summary = _summarise_active_component_type(
        frame.get("active_component_type")
    )
    if active_summary:
        summary["active_component_type"] = active_summary
    if verbose:
        logger.info(
            "target_names_helper_done",
            path=str(output_path),
            rows=len(names_df),
        )
    return {"path": str(output_path), "summary": summary}

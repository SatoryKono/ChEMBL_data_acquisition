"""Normalization helpers for OpenAlex and CrossRef responses."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

from ..pubmed.parsing import combine

__all__ = [
    "normalize_openalex_response",
    "normalize_crossref_response",
]


_OPENALEX_TEMPLATE: dict[str, str] = {
    "OpenAlex.PublicationTypes": "",
    "OpenAlex.TypeCrossref": "",
    "OpenAlex.Genre": "",
    "OpenAlex.Id": "",
    "OpenAlex.Venue": "",
    "OpenAlex.MeshDescriptors": "",
    "OpenAlex.MeshQualifiers": "",
    "OpenAlex.Error": "",
}


def _extract_mesh(raw_entries: Any) -> tuple[list[str], list[str]]:
    """Return descriptor and qualifier names from ``raw_entries``."""

    descriptors: list[str] = []
    qualifiers: list[str] = []
    if not isinstance(raw_entries, Iterable):
        return descriptors, qualifiers

    for entry in raw_entries:
        if not isinstance(entry, Mapping):
            continue
        descriptor = entry.get("descriptor_name")
        if isinstance(descriptor, str) and descriptor:
            descriptors.append(descriptor)
        qualifier_entries = entry.get("qualifiers")
        if not isinstance(qualifier_entries, Iterable):
            continue
        for qualifier in qualifier_entries:
            if not isinstance(qualifier, Mapping):
                continue
            qualifier_name = qualifier.get("qualifier_name")
            if isinstance(qualifier_name, str) and qualifier_name:
                qualifiers.append(qualifier_name)
    return descriptors, qualifiers


def normalize_openalex_response(
    raw: Mapping[str, Any] | None,
    error: str | None,
) -> dict[str, str]:
    """Normalize OpenAlex API payload into a flat dictionary."""

    result = _OPENALEX_TEMPLATE.copy()
    if error or not isinstance(raw, Mapping):
        result["OpenAlex.Error"] = error or "Invalid response"
        return result

    descriptors, qualifiers = _extract_mesh(raw.get("mesh"))
    host_venue = raw.get("host_venue")
    venue_name = ""
    if isinstance(host_venue, Mapping):
        venue = host_venue.get("display_name")
        if isinstance(venue, str):
            venue_name = venue

    result.update(
        {
            "OpenAlex.PublicationTypes": str(raw.get("type", "")) if raw.get("type") is not None else "",
            "OpenAlex.TypeCrossref": str(raw.get("type_crossref", ""))
            if raw.get("type_crossref") is not None
            else "",
            "OpenAlex.Genre": str(raw.get("genre", "")) if raw.get("genre") is not None else "",
            "OpenAlex.Id": str(raw.get("id", "")) if raw.get("id") is not None else "",
            "OpenAlex.Venue": venue_name,
            "OpenAlex.MeshDescriptors": combine(descriptors),
            "OpenAlex.MeshQualifiers": combine(qualifiers),
            "OpenAlex.Error": "",
        }
    )
    return result


_CROSSREF_TEMPLATE: dict[str, str] = {
    "crossref.Type": "",
    "crossref.Subtype": "",
    "crossref.Title": "",
    "crossref.Subtitle": "",
    "crossref.Subject": "",
    "crossref.Error": "",
}


def _first_text(entry: Any) -> str:
    """Return first textual element from CrossRef list fields."""

    if isinstance(entry, str):
        return entry
    if isinstance(entry, Iterable):
        for item in entry:
            if isinstance(item, str) and item:
                return item
    return ""


def _join_subjects(subjects: Any) -> str:
    """Join subject list with ``; `` separator."""

    if not isinstance(subjects, Iterable) or isinstance(subjects, (str, bytes)):
        if isinstance(subjects, str):
            return subjects
        return ""
    values: list[str] = []
    for subject in subjects:
        if isinstance(subject, str) and subject:
            values.append(subject)
    return "; ".join(values)


def normalize_crossref_response(
    raw: Mapping[str, Any] | None,
    error: str | None,
) -> dict[str, str]:
    """Normalize CrossRef API payload into a flat dictionary."""

    result = _CROSSREF_TEMPLATE.copy()
    if error or not isinstance(raw, Mapping):
        result["crossref.Error"] = error or "Invalid response"
        return result

    message = raw.get("message", {})
    if not isinstance(message, Mapping):
        result["crossref.Error"] = "Invalid response"
        return result

    result.update(
        {
            "crossref.Type": str(message.get("type", "")) if message.get("type") is not None else "",
            "crossref.Subtype": str(message.get("subtype", ""))
            if message.get("subtype") is not None
            else "",
            "crossref.Title": _first_text(message.get("title")),
            "crossref.Subtitle": _first_text(message.get("subtitle")),
            "crossref.Subject": _join_subjects(message.get("subject")),
            "crossref.Error": "",
        }
    )
    return result

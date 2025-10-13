"""Shared helper utilities and constants for the test suite."""

from __future__ import annotations

import sys
import types
from collections.abc import Iterable, Mapping
from typing import Any

ASSAY_ENRICHMENT_MIN_RATIO = 0.8
"""Minimum acceptable non-null ratio for enriched assay attributes."""


def ensure_normalizer_stub() -> None:
    """Install a lightweight stub for ``library.integration._normalizers``.

    The production module depends on :mod:`library.pubmed.parsing`, which in turn
    imports the normalizer helpers creating a circular dependency during test
    collection. The stub provides minimal compatible implementations required by
    the CLI modules exercised in the integration tests.
    """

    if "library.integration._normalizers" in sys.modules:
        return

    stub = types.ModuleType("library.integration._normalizers")

    def _combine(values: Any) -> str:
        if isinstance(values, str):
            return values
        if not isinstance(values, Iterable) or isinstance(values, bytes):
            return ""
        tokens = [str(item) for item in values if isinstance(item, str) and item]
        return "; ".join(tokens)

    def _extract_mesh(raw_entries: Any) -> tuple[list[str], list[str]]:
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

    openalex_template: dict[str, str] = {
        "OpenAlex.PublicationTypes": "",
        "OpenAlex.TypeCrossref": "",
        "OpenAlex.Genre": "",
        "OpenAlex.Id": "",
        "OpenAlex.Venue": "",
        "OpenAlex.MeshDescriptors": "",
        "OpenAlex.MeshQualifiers": "",
        "OpenAlex.Error": "",
    }

    def normalize_openalex_response(
        raw: Mapping[str, Any] | None,
        error: str | None,
    ) -> dict[str, str]:
        result = dict(openalex_template)
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
                "OpenAlex.PublicationTypes": str(raw.get("type", ""))
                if raw.get("type") is not None
                else "",
                "OpenAlex.TypeCrossref": str(raw.get("type_crossref", ""))
                if raw.get("type_crossref") is not None
                else "",
                "OpenAlex.Genre": str(raw.get("genre", ""))
                if raw.get("genre") is not None
                else "",
                "OpenAlex.Id": str(raw.get("id", ""))
                if raw.get("id") is not None
                else "",
                "OpenAlex.Venue": venue_name,
                "OpenAlex.MeshDescriptors": _combine(descriptors),
                "OpenAlex.MeshQualifiers": _combine(qualifiers),
                "OpenAlex.Error": "",
            }
        )
        return result

    crossref_template: dict[str, str] = {
        "crossref.Type": "",
        "crossref.Subtype": "",
        "crossref.Title": "",
        "crossref.Subtitle": "",
        "crossref.Subject": "",
        "crossref.Error": "",
    }

    def _first_text(entry: Any) -> str:
        if isinstance(entry, str):
            return entry
        if isinstance(entry, Iterable):
            for item in entry:
                if isinstance(item, str) and item:
                    return item
        return ""

    def _join_subjects(subjects: Any) -> str:
        if isinstance(subjects, str):
            return subjects
        if not isinstance(subjects, Iterable) or isinstance(subjects, bytes):
            return ""
        tokens = [str(item) for item in subjects if isinstance(item, str) and item]
        return "; ".join(tokens)

    def normalize_crossref_response(
        raw: Mapping[str, Any] | None,
        error: str | None,
    ) -> dict[str, str]:
        result = dict(crossref_template)
        if error or not isinstance(raw, Mapping):
            result["crossref.Error"] = error or "Invalid response"
            return result

        message = raw.get("message", {})
        if not isinstance(message, Mapping):
            result["crossref.Error"] = "Invalid response"
            return result

        result.update(
            {
                "crossref.Type": str(message.get("type", ""))
                if message.get("type") is not None
                else "",
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

    stub.normalize_openalex_response = normalize_openalex_response
    stub.normalize_crossref_response = normalize_crossref_response
    stub.__all__ = ["normalize_openalex_response", "normalize_crossref_response"]

    sys.modules["library.integration._normalizers"] = stub

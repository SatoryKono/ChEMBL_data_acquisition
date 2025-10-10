"""Cellularity classification replicating the Power Query logic."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree

import numpy as np
import pandas as pd

from .helpers import first_element_text, normalize_text

ADD_CELLULARITY_SMART_FIELD_NAME = "AddCellularitySmart "

FetchLineageCallable = Callable[[Any, str | None], Sequence[str]]


@dataclass
class Cellularity:
    """Encapsulate the SSoT cellularity rules."""

    fetcher: FetchLineageCallable | None = None
    timeout: float = 10.0
    base_url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    default_email: str = "noreply@example.com"

    unicell_phyla: tuple[str, ...] = (
        "ciliophora",
        "apicomplexa",
        "euglenozoa",
        "dinoflagellata",
        "alveolata",
        "euglenozoa",
        "bacillariophyta",
        "parabasalia",
        "metamonada",
        "choanoflagellata",
        "chlorophyta",
    )
    multicell_animal_phyla: tuple[str, ...] = (
        "chordata",
        "arthropoda",
        "mollusca",
        "nematoda",
        "annelida",
        "ecdysozoa",
        "echinodermata",
        "cnidaria",
        "platyhelminthes",
        "porifera",
        "spiralia",
        "ctenophora",
        "tardigrada",
        "onychophora",
        "bryozoa",
        "brachiopoda",
        "echinodermata",
        "hemichordata",
        "rotifera",
        "nemertea",
        "sipuncula",
        "priapulida",
        "dikarya",
        "loricifera",
        "gastrotricha",
        "kinorhyncha",
    )
    multicell_plant_phyla: tuple[str, ...] = (
        "streptophyta",
        "tracheophyta",
        "bryophyta",
        "marchantiophyta",
    )
    mostly_multicell_phyla: tuple[str, ...] = ("rhodophyta",)
    fungi_phyla: tuple[str, ...] = (
        "ascomycota",
        "basidiomycota",
        "mucoromycota",
        "microsporidia",
        "chytridiomycota",
    )

    add_cellularity_smart_field: str = field(
        init=False, default=ADD_CELLULARITY_SMART_FIELD_NAME
    )

    def __post_init__(self) -> None:
        if self.fetcher is None:
            self.fetcher = self._fetch_lineage_default

    @staticmethod
    def normalize(value: Any) -> str:
        return normalize_text(value)

    def classify_by_lineage(self, superkingdom: Any, phylum: Any) -> str:
        sk = normalize_text(superkingdom)
        ph = normalize_text(phylum)
        if sk == "viruses":
            return "acellular (virus)"
        if sk in {"bacteria", "archaea"}:
            return "unicellular"
        if sk == "eukaryota":
            if ph in self.multicell_animal_phyla or ph in self.multicell_plant_phyla:
                return "multicellular"
            if ph in self.unicell_phyla:
                return "unicellular"
            if ph in self.fungi_phyla:
                return "multicellular"
            if ph in self.mostly_multicell_phyla:
                return "multicellular"
            return "ambiguous"
        if sk:
            return "ambiguous"
        return "ambiguous"

    def get_lineage_names(self, tax_id: Any, email: str | None = None) -> list[str]:
        missing_value: Any = False
        try:
            if bool(pd.isna(tax_id)):
                return []
            missing_value = pd.isna(tax_id)
        except (TypeError, ValueError):
            # ``pd.isna`` may return array-like objects for non-scalar inputs;
            # those should simply fall back to the fetcher.
            missing_value = pd.isna(tax_id)
        except Exception:
            missing_value = False

        if tax_id is None:
            return []

        if isinstance(missing_value, bool | np.bool_) and missing_value:
            return []

        if getattr(missing_value, "all", None) is not None:
            try:
                all_missing = missing_value.all()
            except TypeError:
                all_missing = False
            if isinstance(all_missing, bool | np.bool_) and all_missing:
                return []
            if all_missing is pd.NA:
                return []

        if tax_id is pd.NA:
            return []

        if isinstance(tax_id, float | np.floating) and pd.isna(tax_id):
            return []
        values = list(self.fetcher(tax_id, email))  # type: ignore[arg-type]
        return values

    def _fetch_lineage_default(self, tax_id: Any, email: str | None) -> list[str]:
        params = {
            "db": "taxonomy",
            "id": str(tax_id),
            "retmode": "xml",
            "tool": "cellularity-checker",
            "email": email or self.default_email,
        }
        query = urlencode(params)
        request = Request(
            f"{self.base_url}?{query}",
            headers={"User-Agent": "PowerQuery-M/1.0"},
        )
        try:
            with urlopen(request, timeout=self.timeout) as response:
                payload = response.read()
        except Exception:
            return []
        try:
            root = ElementTree.fromstring(payload)
        except ElementTree.ParseError:
            return []
        taxa = root.findall(".//Taxon")
        if not taxa:
            return []
        taxon = taxa[0]
        lineage_text = first_element_text(list(taxon), "Lineage")
        sci_name = first_element_text(list(taxon), "ScientificName")
        lineage_names: list[str] = []
        if lineage_text is not None:
            trimmed_lineage = lineage_text.strip()
            if len(trimmed_lineage) > 0:
                segments = lineage_text.split(";")
                lineage_names = [segment.strip() for segment in segments]
        if sci_name is not None and sci_name not in lineage_names:
            lineage_names.append(sci_name)
        return lineage_names

    def classify_by_fetch(self, tax_id: Any, email: str | None = None) -> str:
        names = [
            self._lower_token(name) for name in self.get_lineage_names(tax_id, email)
        ]
        if "viruses" in names:
            return "acellular (virus)"
        if "bacteria" in names or "archaea" in names:
            return "unicellular"
        if "eukaryota" in names:
            if (
                self._has_any(names, self.multicell_animal_phyla)
                or self._has_any(names, self.multicell_plant_phyla)
                or "metazoa" in names
            ):
                return "multicellular"
            if self._has_any(names, self.unicell_phyla):
                return "unicellular"
            if "fungi" in names:
                return "multicellular"
            if self._has_any(names, self.mostly_multicell_phyla):
                return "multicellular"
            return "ambiguous"
        if names:
            return "ambiguous"
        return "ambiguous"

    @staticmethod
    def _has_any(names: Sequence[str], candidates: Sequence[str]) -> bool:
        lookup = set(names)
        for candidate in candidates:
            if candidate in lookup:
                return True
        return False

    @staticmethod
    def _lower_token(value: Any) -> str:
        if value is None:
            return ""
        if isinstance(value, float) and pd.isna(value):
            return ""
        if value is pd.NA:
            return ""
        return str(value).lower()


def add_cellularity_smart(
    frame: pd.DataFrame,
    tax_id_column: str,
    superkingdom_column: str,
    phylum_column: str,
    *,
    email: str | None = None,
    fetcher: FetchLineageCallable | None = None,
) -> pd.DataFrame:
    classifier = Cellularity(fetcher=fetcher)
    result = frame.copy()

    def _classify(row: pd.Series) -> str:
        superkingdom = row.get(superkingdom_column)
        phylum = row.get(phylum_column)
        tax_id = row.get(tax_id_column)
        by_lineage = classifier.classify_by_lineage(superkingdom, phylum)
        if (
            by_lineage != "ambiguous"
            and superkingdom is not None
            and str(superkingdom).strip() != ""
        ):
            return by_lineage
        return classifier.classify_by_fetch(tax_id, email=email)

    result["cellularity"] = result.apply(_classify, axis=1)
    return result


__all__ = [
    "Cellularity",
    "add_cellularity_smart",
    "ADD_CELLULARITY_SMART_FIELD_NAME",
]

"""Cellularity classification translated from the Power Query workbook."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Iterable
from xml.etree import ElementTree as ET

import requests

from .helpers import first_element_text, normalize_text

__all__ = ["Cellularity"]


@dataclass(frozen=True)
class Cellularity:
    """Container for cellularity classification rules."""

    # Rule sets are kept as tuples to preserve the original order from the M script.
    UNICELL_PHYLA: tuple[str, ...] = (
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
    MULTICELL_ANIMAL_PHYLA: tuple[str, ...] = (
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
    MULTICELL_PLANT_PHYLA: tuple[str, ...] = (
        "streptophyta",
        "tracheophyta",
        "bryophyta",
        "marchantiophyta",
    )
    MOSTLY_MULTICELL_PHYLA: tuple[str, ...] = ("rhodophyta",)
    FUNGI_PHYLA: tuple[str, ...] = (
        "ascomycota",
        "basidiomycota",
        "mucoromycota",
        "microsporidia",
        "chytridiomycota",
    )

    BASE_URL: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    DEFAULT_EMAIL: str = "noreply@example.com"

    @staticmethod
    def _has_any(values: Iterable[str], targets: Iterable[str]) -> bool:
        values_set = set(values)
        for target in targets:
            if target in values_set:
                return True
        return False

    @classmethod
    def classify_by_lineage(cls, superkingdom: Any, phylum: Any) -> str:
        """Replicate ``ClassifyByLineage`` from the M record."""

        sk = normalize_text(superkingdom)
        ph = normalize_text(phylum)

        if sk == "viruses":
            return "acellular (virus)"
        if sk in {"bacteria", "archaea"}:
            return "unicellular"
        if sk == "eukaryota":
            if ph in cls.MULTICELL_ANIMAL_PHYLA or ph in cls.MULTICELL_PLANT_PHYLA:
                return "multicellular"
            if ph in cls.UNICELL_PHYLA:
                return "unicellular"
            if ph in cls.FUNGI_PHYLA:
                return "multicellular"
            if ph in cls.MOSTLY_MULTICELL_PHYLA:
                return "multicellular"
            return "ambiguous"
        return "ambiguous"

    @classmethod
    @lru_cache(maxsize=1024)
    def get_lineage_names(cls, tax_id: Any, email: str | None = None) -> tuple[str, ...]:
        """Fetch lineage names from NCBI Taxonomy using ``efetch``."""

        if tax_id is None:
            return ()
        text_id = str(tax_id).strip()
        if not text_id:
            return ()

        params = {
            "db": "taxonomy",
            "id": text_id,
            "retmode": "xml",
            "tool": "cellularity-checker",
            "email": email or cls.DEFAULT_EMAIL,
        }
        headers = {"User-Agent": "PowerQuery-M/1.0"}
        try:
            response = requests.get(cls.BASE_URL, params=params, headers=headers, timeout=30)
            response.raise_for_status()
        except requests.RequestException:
            return ()

        content = response.content
        if not content:
            return ()

        try:
            root = ET.fromstring(content)
        except ET.ParseError:
            return ()

        taxon = root.find(".//Taxon")
        if taxon is None:
            return ()

        children = list(taxon)
        lineage_text = first_element_text(children, "Lineage")
        scientific_name = first_element_text(children, "ScientificName")

        lineage_names: list[str] = []
        if lineage_text:
            pieces = [segment.strip() for segment in lineage_text.split(";")]
            lineage_names = [piece for piece in pieces if piece]

        if scientific_name and scientific_name not in lineage_names:
            lineage_names.append(scientific_name)

        return tuple(lineage_names)

    @classmethod
    def classify_by_fetch(cls, tax_id: Any, email: str | None = None) -> str:
        """Mirror ``ClassifyByFetch`` using :meth:`get_lineage_names`."""

        names = [normalize_text(name) for name in cls.get_lineage_names(tax_id, email=email)]
        if not names:
            return "ambiguous"

        if "viruses" in names:
            return "acellular (virus)"
        if "bacteria" in names or "archaea" in names:
            return "unicellular"
        if "eukaryota" in names:
            if (
                cls._has_any(names, cls.MULTICELL_ANIMAL_PHYLA)
                or cls._has_any(names, cls.MULTICELL_PLANT_PHYLA)
                or "metazoa" in names
            ):
                return "multicellular"
            if cls._has_any(names, cls.UNICELL_PHYLA):
                return "unicellular"
            if "fungi" in names:
                return "multicellular"
            if cls._has_any(names, cls.MOSTLY_MULTICELL_PHYLA):
                return "multicellular"
            return "ambiguous"
        return "ambiguous"

    @classmethod
    def add_cellularity_smart(
        cls,
        frame,
        tax_id_column: str,
        superkingdom_column: str,
        phylum_column: str,
        *,
        email: str | None = None,
    ):
        """Add the ``cellularity`` column replicating ``#"AddCellularitySmart "``."""

        import pandas as pd

        def classify_row(row: pd.Series) -> str:
            sk_value = row.get(superkingdom_column)
            ph_value = row.get(phylum_column)
            by_lineage = cls.classify_by_lineage(sk_value, ph_value)
            if (
                by_lineage != "ambiguous"
                and sk_value is not None
                and str(sk_value).strip() != ""
            ):
                return by_lineage
            return cls.classify_by_fetch(row.get(tax_id_column), email=email)

        result = frame.copy()
        result["cellularity"] = result.apply(classify_row, axis=1)
        return result

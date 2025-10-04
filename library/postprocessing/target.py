"""Post-processing helpers for target exports.

This module ports the Power Query post-processing logic used in
``target_all`` reports to Python. The implementation mirrors the
behaviour of the original M queries to ensure byte-for-byte identical
results for downstream analytics.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
import math
from pathlib import Path
from typing import Iterable
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

import pandas as pd

from ..common.log import logger

ENCODINGS_TO_TRY: tuple[str, ...] = ("utf-8", "utf-8-sig", "cp1252")
NCBI_EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_TOOL_NAME = "cellularity-checker"
DEFAULT_EMAIL = "noreply@example.com"
USER_AGENT = "Python-Postproc/1.0"


def normalize_text(value: str | None | object) -> str:
    """Return a strictly normalised textual value.

    Power Query's ``NormalizeText`` helper converts ``null`` to an empty
    string, trims whitespace, and lowercases the result. The Python
    implementation mirrors that behaviour exactly.
    """

    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value)
    return text.strip().lower()


def first_element_text(nodes: Iterable[ET.Element] | None, element_name: str) -> str | None:
    """Return the text value of the first element matching ``element_name``.

    The Power Query pipeline inspects a flattened XML table. In Python we
    work directly with ``Element`` objects, iterating through the provided
    nodes while preserving the semantics of returning ``null``/``None``
    when the node or its textual payload is absent.
    """

    if nodes is None:
        return None
    for node in nodes:
        if node.tag != element_name:
            continue
        text = node.text
        if text is None:
            return None
        return text
    return None


def _distinct(sequence: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for item in sequence:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result


def _http_get(url: str, params: dict[str, str], *, timeout: float) -> bytes | None:
    query = urlencode(params)
    request = Request(f"{url}?{query}", headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request, timeout=timeout) as response:  # nosec: B310 - controlled URL
            return response.read()
    except URLError:
        return None


@dataclass(frozen=True)
class _CellularityRules:
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


class CellularityClassifier:
    """Encapsulate the classification rules for cellularity labels."""

    def __init__(
        self,
        *,
        email: str | None = None,
        timeout: float = 15.0,
        offline: bool = False,
    ) -> None:
        self.email = email or DEFAULT_EMAIL
        self.timeout = timeout
        self.offline = offline
        self.rules = _CellularityRules()
        self.http_requests = 0
        self._fetch_cache: dict[str, list[str]] = {}

    def classify_by_lineage(self, superkingdom: str | None, phylum_or_class: str | None) -> str:
        sk = normalize_text(superkingdom)
        ph = normalize_text(phylum_or_class)
        if sk == "viruses":
            return "acellular (virus)"
        if sk in {"bacteria", "archaea"}:
            return "unicellular"
        if sk == "eukaryota":
            if ph in self.rules.multicell_animal_phyla or ph in self.rules.multicell_plant_phyla:
                return "multicellular"
            if ph in self.rules.unicell_phyla:
                return "unicellular"
            if ph in self.rules.fungi_phyla:
                return "multicellular"
            if ph in self.rules.mostly_multicell_phyla:
                return "multicellular"
            return "ambiguous"
        if sk:
            return "ambiguous"
        return "ambiguous"

    def _fetch_lineage_names(self, tax_id: str) -> list[str]:
        if tax_id in self._fetch_cache:
            return self._fetch_cache[tax_id]

        if self.offline:
            self._fetch_cache[tax_id] = []
            return []

        params = {
            "db": "taxonomy",
            "id": tax_id,
            "retmode": "xml",
            "tool": NCBI_TOOL_NAME,
            "email": self.email,
        }
        self.http_requests += 1
        payload = _http_get(NCBI_EUTILS_URL, params, timeout=self.timeout)
        if payload is None:
            self._fetch_cache[tax_id] = []
            return []

        try:
            tree = ET.fromstring(payload)
        except ET.ParseError:
            self._fetch_cache[tax_id] = []
            return []

        taxon_nodes = tree.findall(".//Taxon")
        if not taxon_nodes:
            self._fetch_cache[tax_id] = []
            return []
        taxon = taxon_nodes[0]
        children = list(taxon)
        lineage_text = first_element_text(children, "Lineage")
        sci_name = first_element_text(children, "ScientificName")
        names: list[str] = []
        if lineage_text is not None:
            stripped = lineage_text.strip()
            if stripped:
                names = [segment.strip() for segment in stripped.split(";") if segment.strip()]
        if sci_name and sci_name not in names:
            names.append(sci_name)

        self._fetch_cache[tax_id] = names
        return names

    def classify_by_fetch(self, tax_id: object) -> str:
        if tax_id is None:
            return "ambiguous"
        if isinstance(tax_id, float) and math.isnan(tax_id):
            return "ambiguous"

        tax_id_str = str(tax_id)
        if not tax_id_str.strip():
            return "ambiguous"

        names = [name.lower() for name in self._fetch_lineage_names(tax_id_str)]
        if not names:
            return "ambiguous"

        if "viruses" in names:
            return "acellular (virus)"
        if "bacteria" in names or "archaea" in names:
            return "unicellular"
        if "eukaryota" in names:
            if any(name in self.rules.multicell_animal_phyla for name in names) or any(
                name in self.rules.multicell_plant_phyla for name in names
            ) or "metazoa" in names:
                return "multicellular"
            if any(name in self.rules.unicell_phyla for name in names):
                return "unicellular"
            if "fungi" in names:
                return "multicellular"
            if any(name in self.rules.mostly_multicell_phyla for name in names):
                return "multicellular"
            return "ambiguous"
        return "ambiguous"

    def add_cellularity_smart(
        self,
        frame: pd.DataFrame,
        tax_id_column: str,
        superkingdom_column: str,
        phylum_or_class_column: str,
    ) -> pd.Series:
        def _classify(row: pd.Series) -> str:
            sk_value = row.get(superkingdom_column)
            by_lineage = self.classify_by_lineage(sk_value, row.get(phylum_or_class_column))
            sk_text = "" if sk_value is None else str(sk_value)
            if by_lineage != "ambiguous" and sk_text.strip():
                return by_lineage
            return self.classify_by_fetch(row.get(tax_id_column))

        return frame.apply(_classify, axis=1)


def _transform_reaction_ec_numbers(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, float) and math.isnan(value):
        return []
    text = str(value)
    parts = text.split("|")
    distinct_parts = _distinct(parts)
    first_segments = [segment.split(".")[0] if segment else "" for segment in distinct_parts]
    return _distinct(first_segments)


MULTIFUNCTIONAL_DROP_COLUMNS: tuple[str, ...] = (
    "isoform_ids",
    "isoform_names",
    "isoform_synonyms",
    "organism",
    "taxon_id",
    "lineage_superkingdom",
    "lineage_phylum",
    "lineage_class",
    "xref_chembl",
    "gtop_synonyms",
    "gtop_natural_ligands_n",
    "gtop_interactions_n",
    "gtop_function_text_short",
    "uniProtkbId",
    "hgnc_name",
    "secondaryAccessionNames",
    "transmembrane",
    "intramembrane",
    "hgnc_id",
    "features_signal_peptide",
    "tax_id",
    "species_group_flag",
    "family",
    "xref_uniprot",
    "xref_ensembl",
    "pfam",
    "interpro",
    "xref_pdb",
    "xref_alphafold",
    "SUPFAM",
    "PROSITE",
    "InterPro",
    "Pfam",
    "PRINTS",
    "TCDB",
    "target_type",
    "protein_classifications",
    "protein_name_alt",
    "recommendedName",
    "pref_name",
    "target_components",
    "cross_references",
    "gene_symbol_list",
    "protein_synonym_list",
    "ptm_glycosylation",
    "ptm_lipidation",
    "ptm_disulfide_bond",
    "ptm_modified_residue",
    "glycosylation",
    "lipidation",
    "disulfide_bond",
    "modified_residue",
    "phosphorylation",
    "acetylation",
    "ubiquitination",
    "signal_peptide",
    "propeptide",
    "gene_symbol",
    "protein_name_canonical",
    "sequence_length",
    "features_transmembrane",
    "features_topology",
    "xref_iuphar",
    "gtop_target_id",
    "uniprot_last_update",
    "uniprot_version",
    "pipeline_version",
    "timestamp_utc",
    "secondaryAccessions",
    "geneName",
    "molecular_function",
    "cellular_component",
    "subcellular_location",
    "topology",
    "GuidetoPHARMACOLOGY",
    "protein_class_pred_L1",
    "protein_class_pred_L2",
    "protein_class_pred_L3",
    "protein_class_pred_rule_id",
    "protein_class_pred_evidence",
    "protein_class_pred_confidence",
    "iuphar_target_id",
    "iuphar_family_id",
    "iuphar_type",
    "iuphar_class",
    "iuphar_subclass",
    "iuphar_chain",
    "iuphar_name",
    "iuphar_full_id_path",
    "iuphar_full_name_path",
)


def _prepare_multifunctional_frame(source: pd.DataFrame) -> pd.DataFrame:
    frame = source.drop(columns=[col for col in MULTIFUNCTIONAL_DROP_COLUMNS if col in source.columns])
    if "reaction_ec_numbers" in frame.columns:
        frame = frame.copy()
        frame["reaction_ec_numbers"] = frame["reaction_ec_numbers"].map(_transform_reaction_ec_numbers)
    else:
        frame = frame.copy()
        frame["reaction_ec_numbers"] = [[] for _ in range(len(frame))]
    frame["multifunctional_enzyme"] = frame["reaction_ec_numbers"].map(lambda values: len(values) > 1)
    return frame[["target_chembl_id", "multifunctional_enzyme"]]


TARGET_OUTPUT_COLUMNS: tuple[str, ...] = (
    "target_chembl_id",
    "uniprot_id_primary",
    "organism",
    "taxon_id",
    "lineage_superkingdom",
    "lineage_phylum",
    "lineage_class",
    "cellularity",
    "multifunctional_enzyme",
)


def _find_latest_output_csv(prefix: str = "output.target_", search_roots: Iterable[Path] | None = None) -> Path | None:
    if search_roots is None:
        cwd = Path.cwd()
        search_roots = (
            cwd,
            cwd / "data",
            cwd / "reports",
        )

    latest_path: Path | None = None
    for root in search_roots:
        if not root.exists():
            continue
        for path in root.rglob(f"{prefix}*.csv"):
            if latest_path is None or path.stat().st_mtime > latest_path.stat().st_mtime:
                latest_path = path
    return latest_path


def _read_csv_with_fallback(path: Path) -> pd.DataFrame:
    last_error: Exception | None = None
    for encoding in ENCODINGS_TO_TRY:
        try:
            return pd.read_csv(path, encoding=encoding)
        except UnicodeDecodeError as exc:  # pragma: no cover - depends on input encoding
            last_error = exc
    if last_error is not None:
        raise last_error
    return pd.read_csv(path)


def process_targets(
    input_csv: str | None = None,
    output_csv: str | None = None,
    output_prefix: str = "organism",
    ncbi_email: str | None = None,
    http_timeout_s: float = 15.0,
    offline: bool = False,
    verbose: bool = True,
) -> Path:
    """Post-process the target export according to the Power Query logic."""

    log = logger if verbose else logging.getLogger("null")
    if not verbose:
        logging.getLogger("null").addHandler(logging.NullHandler())

    input_path = Path(input_csv) if input_csv is not None else _find_latest_output_csv()
    if input_path is None:
        raise FileNotFoundError("No input CSV provided and no output.target_*.csv found")

    if verbose:
        log.info("Reading target CSV", extra={"path": str(input_path)})
    source = _read_csv_with_fallback(input_path)

    classifier = CellularityClassifier(email=ncbi_email, timeout=http_timeout_s, offline=offline)

    source_base = source[[
        "target_chembl_id",
        "uniprot_id_primary",
        "organism",
        "taxon_id",
        "lineage_superkingdom",
        "lineage_phylum",
        "lineage_class",
    ]].copy()

    for column in ("lineage_superkingdom", "lineage_phylum", "lineage_class"):
        source_base[column] = source_base[column].map(
            lambda value: value.lower() if isinstance(value, str) else value
        )

    source_base["cellularity"] = classifier.add_cellularity_smart(
        source_base,
        "taxon_id",
        "lineage_superkingdom",
        "lineage_class",
    )

    multifunctional_frame = _prepare_multifunctional_frame(source)
    merged = source_base.merge(multifunctional_frame, on="target_chembl_id", how="left")
    merged = merged.reindex(columns=TARGET_OUTPUT_COLUMNS)

    if output_csv is None:
        output_csv = f"{output_prefix}.{input_path.name}"
        output_path = input_path.with_name(output_csv)
    else:
        output_path = Path(output_csv)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, index=False, encoding="utf-8")

    ambiguous_count = int((merged["cellularity"] == "ambiguous").sum())
    if verbose:
        log.info(
            "Wrote organism-level target classification",
            extra={
                "rows": len(merged),
                "path": str(output_path),
                "http_requests": classifier.http_requests,
                "ambiguous": ambiguous_count,
            },
        )

    return output_path


__all__ = [
    "process_targets",
    "normalize_text",
    "first_element_text",
    "CellularityClassifier",
]


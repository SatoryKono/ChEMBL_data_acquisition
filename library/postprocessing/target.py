"""Target organism post-processing helpers.

This module provides a lightweight replica of the organism-centric workbook
used during the early iterations of the target pipeline.  The helpers here aim
to be deterministic, avoid interactive dependencies and offer an offline mode
that skips network lookups while still producing a stable CSV artefact.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
import re
from typing import Iterable, Mapping, MutableMapping, Sequence
from xml.etree import ElementTree

import pandas as pd
import requests

LOGGER = logging.getLogger(__name__)

# Power Query snapshots occasionally stored the aggregated target table in
# various encodings.  We mirror the legacy fallbacks so older exports can be
# replayed without manual re-encoding.
DEFAULT_ENCODINGS: Sequence[str] = (
    "utf-8-sig",
    "utf-8",
    "cp1251",
    "cp1252",
    "latin-1",
)

NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_TIMEOUT = 30
USER_AGENT = "Python-Postproc/1.0"

# Canonical column order for the generated artefact.  Additional columns (when
# present in the source CSV) are preserved after these core fields sorted
# alphabetically to guarantee stable diffs.
CORE_OUTPUT_COLUMNS: list[str] = [
    "target_chembl_id",
    "organism",
    "taxon_id",
    "organism_type",
    "type",
    "unicellular_organism",
    "multifunctional_enzyme",
    "gene_index",
    "taxon_index",
    "target_sort_order",
]

NULL_LITERALS = frozenset({"", "-", "na", "n/a", "none", "null", "nan"})
TOKEN_PATTERN = re.compile(r"[^|]+")

TYPE_MULTICELLULAR = "Multicellular organism"
TYPE_UNICELLULAR = "Unicellular organism"
TYPE_VIRAL = "Viruses"


def NormalizeText(
    value: object | None,
    *,
    lowercase: bool = False,
    uppercase: bool = False,
    default: str = "",
) -> str:
    """Return a trimmed representation of *value*.

    The helper mirrors the steps present in the Power Query workbook.  Values
    such as ``"NA"`` or ``"None"`` collapse to ``default`` and optional case
    normalisation can be requested via ``lowercase`` or ``uppercase`` (the two
    modes are mutually exclusive).
    """

    if lowercase and uppercase:
        raise ValueError("lowercase and uppercase modes are mutually exclusive")

    if value is None:
        return default

    if isinstance(value, float) and pd.isna(value):
        return default

    text = str(value).strip()
    if not text:
        return default

    collapsed = re.sub(r"\s+", " ", text)
    literal = collapsed.lower()
    if literal in NULL_LITERALS:
        return default

    if lowercase:
        return collapsed.lower()
    if uppercase:
        return collapsed.upper()
    return collapsed


def FirstElementText(
    value: object | None,
    *,
    delimiter: str = "|",
    uppercase: bool = False,
) -> str:
    """Return the first non-empty token extracted from *value*.

    Tokens are split using ``delimiter`` (default ``"|"``) and normalised using
    :func:`NormalizeText`.  The helper preserves the legacy behaviour where
    placeholder values collapse to the empty string.
    """

    if value is None:
        return ""

    normalized = NormalizeText(value)
    if not normalized:
        return ""

    for match in TOKEN_PATTERN.findall(normalized):
        token = NormalizeText(match, uppercase=uppercase)
        if token:
            return token
    return ""


@dataclass(slots=True)
class Cellularity:
    """Container tracking taxonomy fields and the inferred label."""

    genus: str
    superkingdom: str
    phylum: str
    klass: str

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        *,
        offline: bool,
        session: requests.Session | None,
    ) -> "Cellularity":
        genus = NormalizeText(record.get("genus"), lowercase=True)
        superkingdom = NormalizeText(
            record.get("lineage_superkingdom") or record.get("superkingdom"),
            lowercase=True,
        )
        phylum = NormalizeText(
            record.get("lineage_phylum") or record.get("phylum"),
            lowercase=True,
        )
        klass = NormalizeText(
            record.get("lineage_class") or record.get("class"),
            lowercase=True,
        )

        if offline:
            return cls(genus=genus, superkingdom=superkingdom, phylum=phylum, klass=klass)

        needs_lookup = not (genus and superkingdom and phylum and klass)
        tax_id = NormalizeText(record.get("taxon_id")) or NormalizeText(record.get("tax_id"))
        if needs_lookup and tax_id and session is not None:
            fetched = _fetch_taxonomy(session, tax_id)
            genus = genus or fetched.get("genus", "")
            superkingdom = superkingdom or fetched.get("superkingdom", "")
            phylum = phylum or fetched.get("phylum", "")
            klass = klass or fetched.get("class", "")

        return cls(genus=genus, superkingdom=superkingdom, phylum=phylum, klass=klass)

    def classify(self) -> str:
        if self.superkingdom == "viruses":
            return TYPE_VIRAL
        if self.superkingdom in {"bacteria", "archaea"}:
            return TYPE_UNICELLULAR

        unicellular_phyly = {
            "apicomplexa",
            "amoebozoa",
            "ciliophora",
            "chlorophyta",
            "euglenozoa",
            "myzozoa",
            "metamonada",
            "microsporidia",
        }
        unicellular_classes = {
            "aconoidasida",
            "conoidasida",
            "kinetoplastea",
            "saccharomycetes",
            "pneumocystidomycetes",
            "malasseziomycetes",
            "chlorophyceae",
        }
        unicellular_genera = {
            "candida",
            "malassezia",
            "pneumocystis",
            "chlamydomonas",
        }

        if self.phylum in unicellular_phyly:
            return TYPE_UNICELLULAR
        if self.klass in unicellular_classes:
            return TYPE_UNICELLULAR
        if self.genus in unicellular_genera:
            return TYPE_UNICELLULAR
        return TYPE_MULTICELLULAR

    def target_codes(self) -> tuple[str, str]:
        """Return taxonomy codes for ``target_sort_order`` construction."""

        superkingdom_code = _encode_taxonomy_token(self.superkingdom)
        phylum_code = _encode_taxonomy_token(self.phylum, suffix="-", width=4)
        return superkingdom_code, phylum_code

    @property
    def unicellular_flag(self) -> bool:
        label = self.classify()
        return label in {TYPE_UNICELLULAR, TYPE_VIRAL}


def multifunctional(record: Mapping[str, object]) -> bool:
    """Heuristic matching the Power Query ``multifunctional`` step."""

    candidates = [
        NormalizeText(record.get("IUPHAR_type"), lowercase=True),
        NormalizeText(record.get("target_type"), lowercase=True),
        NormalizeText(record.get("organism_type"), lowercase=True),
    ]
    return any("multifunctional" in candidate for candidate in candidates if candidate)


def _encode_taxonomy_token(value: str, *, width: int = 4, suffix: str = "--") -> str:
    if not value:
        return "0" * width + suffix

    base = value[: width].ljust(width, "0")
    return base.replace(" ", "0").upper() + suffix


def _fetch_taxonomy(session: requests.Session, tax_id: str) -> Mapping[str, str]:
    params = {"db": "taxonomy", "id": tax_id}
    headers = {"User-Agent": USER_AGENT}
    try:
        response = session.get(NCBI_URL, params=params, headers=headers, timeout=NCBI_TIMEOUT)
        response.raise_for_status()
    except requests.RequestException as exc:  # pragma: no cover - defensive
        LOGGER.warning("NCBI taxonomy request failed for %s: %s", tax_id, exc)
        return {}

    try:
        root = ElementTree.fromstring(response.text)
    except ElementTree.ParseError as exc:  # pragma: no cover - defensive
        LOGGER.warning("Unable to parse taxonomy response for %s: %s", tax_id, exc)
        return {}

    lineage: MutableMapping[str, str] = {}
    for node in root.findall(".//Taxon"):
        rank = NormalizeText(node.findtext("Rank"), lowercase=True)
        name = NormalizeText(node.findtext("ScientificName"), lowercase=True)
        if rank in {"superkingdom", "phylum", "class", "genus"} and name:
            lineage[rank] = name

    return lineage


def _read_csv(path: Path, *, encodings: Iterable[str], sep: str) -> pd.DataFrame:
    last_error: Exception | None = None
    for encoding in encodings:
        try:
            LOGGER.info("Reading target table from %s (encoding=%s)", path, encoding)
            return pd.read_csv(path, dtype=str, encoding=encoding, sep=sep)
        except UnicodeDecodeError as exc:  # pragma: no cover - defensive
            last_error = exc
            continue
    if last_error is not None:
        raise last_error
    raise UnicodeDecodeError("utf-8", b"", 0, 1, "unable to decode input")


def _resolve_paths(
    source: Path | str,
    *,
    output_dir: Path | str | None,
) -> tuple[Path, Path]:
    source_path = Path(source)
    if source_path.is_dir():
        matches = sorted(source_path.glob("output.target_*.csv"))
        if not matches:
            raise FileNotFoundError("No output.target_*.csv file found in directory")
        if len(matches) > 1:
            raise FileExistsError("Multiple output.target_*.csv files found; specify one explicitly")
        source_path = matches[0]

    if not source_path.name.startswith("output.target_"):
        raise ValueError("Input file must match pattern output.target_*.csv")

    if output_dir is None:
        output_path = source_path.with_name(f"organism.{source_path.name}")
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"organism.{source_path.name}"

    return source_path, output_path


def process_targets(
    source: Path | str,
    *,
    output_dir: Path | str | None = None,
    sep: str = ",",
    encodings: Sequence[str] | None = None,
    offline: bool = False,
    session: requests.Session | None = None,
) -> Path:
    """Process the aggregated target CSV and emit an organism-centric snapshot.

    Parameters
    ----------
    source:
        Directory containing ``output.target_*.csv`` or a path to the CSV file
        itself.
    output_dir:
        Optional directory receiving the generated artefact.  Defaults to the
        location of *source*.
    sep:
        Field separator used by the CSV file.
    encodings:
        Ordered list of encodings that are attempted when reading the input.
    offline:
        When ``True`` skips NCBI lookups and relies solely on embedded taxonomy
        columns.
    session:
        Optional :class:`requests.Session` reused for taxonomy lookups.
    """

    input_path, output_path = _resolve_paths(source, output_dir=output_dir)

    if encodings is None:
        encodings = DEFAULT_ENCODINGS

    frame = _read_csv(input_path, encodings=encodings, sep=sep)
    if frame.empty:
        LOGGER.warning("Input CSV %s is empty", input_path)

    working = frame.copy()

    should_close = False
    if not offline:
        if session is None:
            session = requests.Session()
            should_close = True

    taxon_col = "taxon_id" if "taxon_id" in working.columns else "tax_id"
    if taxon_col not in working.columns:
        working[taxon_col] = ""

    def _compute_row(row: pd.Series) -> pd.Series:
        cell = Cellularity.from_record(row, offline=offline, session=session if not offline else None)
        label = cell.classify()
        gene_symbol = FirstElementText(
            row.get("gene_symbol")
            or row.get("gene_name")
            or row.get("gene")
            or row.get("pref_name"),
            uppercase=True,
        )
        token = gene_symbol or "-"
        gene_index = token.ljust(8, "#")[:8]

        taxon_index = NormalizeText(row.get(taxon_col))
        taxon_index = taxon_index.zfill(6) if taxon_index.isdigit() else "000000"

        superkingdom_code, phylum_code = cell.target_codes()
        class_code = NormalizeText(row.get("IUPHAR_class")) or superkingdom_code
        subclass_code = NormalizeText(row.get("IUPHAR_subclass")) or phylum_code

        target_sort_order = f"{class_code}:{subclass_code}:{gene_index}:{taxon_index}"
        return pd.Series(
            {
                "organism_type": label,
                "type": label,
                "unicellular_organism": cell.unicellular_flag,
                "multifunctional_enzyme": multifunctional(row),
                "gene_index": gene_index,
                "taxon_index": taxon_index,
                "target_sort_order": target_sort_order,
            }
        )

    try:
        derived = working.apply(_compute_row, axis=1)
        for column in derived.columns:
            working[column] = derived[column]
    finally:
        if should_close and session is not None:
            session.close()

    for column in CORE_OUTPUT_COLUMNS:
        if column not in working.columns:
            working[column] = ""

    extra_columns = sorted(col for col in working.columns if col not in CORE_OUTPUT_COLUMNS)
    ordered = working[CORE_OUTPUT_COLUMNS + extra_columns]

    ordered = ordered.fillna("").astype(str)

    type_counts = ordered["type"].value_counts().to_dict()
    LOGGER.info("Organism classification summary: %s", type_counts)
    LOGGER.info(
        "Multifunctional enzymes: %d", (ordered["multifunctional_enzyme"] == "True").sum()
    )

    ordered.to_csv(output_path, index=False, encoding="utf-8", lineterminator="\n")
    LOGGER.info("Wrote organism snapshot to %s", output_path)
    return output_path


__all__ = [
    "NormalizeText",
    "FirstElementText",
    "Cellularity",
    "multifunctional",
    "process_targets",
]


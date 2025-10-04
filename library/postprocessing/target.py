"""Post-processing helpers for target exports.

The implementation mirrors the Single Source of Truth (SSoT) workbook that is
used by the analytics team to derive the *organism* enriched target tables.
The goal of this module is to provide a deterministic Python port that can be
executed in automated environments without depending on Excel.

The transformation performs the following high level steps:

* Locate the latest ``output.target_*.csv`` file within the provided directory
  (or process a specific file when the exact path is supplied).
* Load the CSV using a small catalogue of legacy encodings that were observed
  in historical exports.
* Normalise key textual columns using the :func:`NormalizeText` helper to
  ensure that whitespace, sentinels and casing are processed consistently.
* Derive additional helper columns used downstream by analytics including the
  organism cellularity classification and the multifunctional enzyme flag.
* Persist the transformed dataset as ``organism.output.target_*.csv`` encoded
  in UTF-8 with Unix line terminators, matching the legacy SSoT workbook.

Where possible the logic intentionally stays close to the semantics of the
original spreadsheet.  For example, the :class:`Cellularity` helper implements
the same lookup strategy by querying the NCBI taxonomy API (with graceful
offline fallbacks) and the :func:`FirstElementText` helper emulates the Excel
formulae that extract primary values from pipe separated lists.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pandas as pd
import requests

from ..schemas.targets import TARGETS_COLUMN_ORDER

LOGGER = logging.getLogger(__name__)


# ===== Constants ============================================================

DEFAULT_ENCODINGS: Sequence[str] = (
    "utf-8-sig",
    "utf-8",
    "cp1251",
    "cp1252",
    "latin-1",
)

INPUT_PATTERN = "output.target_*.csv"
OUTPUT_PREFIX = "organism.output."
NCBI_TAXON_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{tax_id}"
NCBI_HEADERS = {"User-Agent": "Python-Postproc/1.0"}
NCBI_TIMEOUT = 10

_NORMALISE_WS_RE = re.compile(r"\s+")
_SPLIT_VALUE_RE = re.compile(r"[|;,/]+")


# ===== Helper utilities =====================================================

def NormalizeText(value: object) -> str:
    """Return a cleaned textual representation of ``value``.

    Historically the target workbook stored a mixture of Excel ``NULL``
    sentinels, ``N/A`` values and values padded with whitespace.  The helper
    performs a conservative normalisation by collapsing whitespace, stripping
    the result and mapping a small set of sentinel strings to the empty string.
    All other values are preserved verbatim to avoid altering the analytical
    semantics of the table.
    """

    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):  # pragma: no cover - defensive
        return ""
    if value is pd.NA or value is pd.NaT:  # pragma: no cover - defensive
        return ""
    text = str(value).strip()
    if not text:
        return ""
    text = _NORMALISE_WS_RE.sub(" ", text)
    lowered = text.lower()
    if lowered in {"na", "n/a", "none", "null", "nan", "-"}:
        return ""
    return text


def FirstElementText(value: object) -> str:
    """Extract the first textual entry from a delimited string.

    ``value`` is expected to use the pipe separator (``|``) although historical
    exports also contained commas, semi-colons and forward slashes depending on
    the upstream loader.  The helper mirrors the Excel ``TEXTSPLIT`` based logic
    by considering each of these separators.
    """

    text = NormalizeText(value)
    if not text:
        return ""
    if "|" not in text and ";" not in text and "," not in text and "/" not in text:
        return text
    parts = [part for part in _SPLIT_VALUE_RE.split(text) if part]
    if not parts:
        return ""
    return NormalizeText(parts[0])


@dataclass(frozen=True)
class Cellularity:
    """Representation of the organism cellularity classification.

    Parameters
    ----------
    value:
        Human readable classification label such as ``"unicellular"``.
    method:
        Indicates whether the value was inferred locally or via the NCBI API.
    lineage:
        Tuple of lineage tokens used during inference.  This is primarily
        exposed for logging and testing purposes to guarantee determinism.
    """

    value: str
    method: str
    lineage: tuple[str, ...]

    @classmethod
    def from_row(
        cls,
        row: Mapping[str, object],
        *,
        offline: bool = False,
        session: requests.Session | None = None,
    ) -> "Cellularity":
        """Infer a :class:`Cellularity` instance from a row.

        The logic mirrors the Excel implementation which first relies on local
        lineage columns and only queries NCBI when a ``taxon_id`` is available
        and the caller allows network access.  Any network errors are treated as
        soft failures and the function falls back to the locally inferred value.
        """

        lineage_superkingdom = NormalizeText(row.get("lineage_superkingdom"))
        lineage_phylum = NormalizeText(row.get("lineage_phylum"))
        lineage_class = NormalizeText(row.get("lineage_class"))
        local_lineage = tuple(
            value
            for value in (lineage_superkingdom, lineage_phylum, lineage_class)
            if value
        )
        local_value = _classify_lineage(local_lineage)

        if offline:
            return cls(local_value, "local", local_lineage)

        taxon_id = NormalizeText(row.get("taxon_id") or row.get("tax_id"))
        if not taxon_id:
            return cls(local_value, "local", local_lineage)

        try:
            lineage = _fetch_taxonomy_lineage(taxon_id, session=session)
        except requests.RequestException as exc:  # pragma: no cover - defensive
            LOGGER.warning(
                "Unable to resolve taxonomy %s via NCBI (%s). Falling back to local classification.",
                taxon_id,
                exc,
            )
            return cls(local_value, "local", local_lineage)

        if not lineage:
            return cls(local_value, "local", local_lineage)

        remote_value = _classify_lineage(lineage)
        return cls(remote_value, "ncbi", tuple(lineage))


def multifunctional(row: Mapping[str, object]) -> bool:
    """Return ``True`` if the row represents a multifunctional enzyme.

    The SSoT workbook considers several columns when deriving the flag.  The
    Python port mirrors this behaviour by scanning a selection of classification
    columns for the literal ``"multifunctional enzyme"`` token.
    """

    CANDIDATE_COLUMNS = (
        "target_type",
        "protein_classifications",
        "protein_class_pred_L1",
        "protein_class_pred_L2",
        "protein_class_pred_L3",
    )
    for column in CANDIDATE_COLUMNS:
        value = NormalizeText(row.get(column))
        if value and "multifunctional enzyme" in value.lower():
            return True
    return False


# ===== Internal helpers =====================================================

def _classify_lineage(lineage: Iterable[str]) -> str:
    values = [value.lower() for value in lineage]
    if not values:
        return "unknown"

    superkingdom = values[0]
    if superkingdom in {"bacteria", "archaea"}:
        return "prokaryote"
    if superkingdom == "virus" or "viruses" in values:
        return "viral"
    if superkingdom == "eukaryota":
        if any(token in values for token in {"metazoa", "chordata", "mammalia"}):
            return "multicellular"
        if any(token in values for token in {"fungi", "ascomycota", "basidiomycota"}):
            return "fungal"
        return "unicellular"
    return "unknown"


def _fetch_taxonomy_lineage(
    taxon_id: str,
    *,
    session: requests.Session | None = None,
) -> tuple[str, ...] | tuple[()]:
    url = NCBI_TAXON_URL.format(taxon_id=taxon_id)
    sess = session or requests.Session()
    response = sess.get(url, headers=NCBI_HEADERS, timeout=NCBI_TIMEOUT)
    response.raise_for_status()

    payload = response.json()
    try:
        nodes = payload["taxonomy_nodes"]
    except (KeyError, TypeError) as exc:  # pragma: no cover - defensive
        LOGGER.warning("Unexpected taxonomy payload for %s: %s", taxon_id, exc)
        return tuple()

    if not nodes:
        return tuple()

    taxonomy = nodes[0].get("taxonomy", {})
    lineage = taxonomy.get("lineage", [])
    lineage_names: list[str] = []
    for entry in lineage:
        name = entry.get("taxon_name") or entry.get("name")
        if name:
            lineage_names.append(NormalizeText(name))

    if taxonomy.get("organism_name"):
        lineage_names.append(NormalizeText(taxonomy["organism_name"]))

    return tuple(value for value in lineage_names if value)


def _read_csv(
    path: Path,
    *,
    encodings: Sequence[str],
    sep: str,
) -> pd.DataFrame:
    last_error: Exception | None = None
    for encoding in encodings:
        try:
            LOGGER.info(
                "Reading target export %s with encoding %s", path, encoding
            )
            return pd.read_csv(path, dtype=str, encoding=encoding, sep=sep)
        except UnicodeDecodeError as exc:  # pragma: no cover - defensive
            last_error = exc
            continue
    if last_error is not None:
        raise last_error
    raise UnicodeDecodeError("utf-8", b"", 0, 1, "unable to decode input")


def _stabilise_columns(df: pd.DataFrame) -> pd.DataFrame:
    known_order = [col for col in TARGETS_COLUMN_ORDER if col in df.columns]
    remaining = [col for col in df.columns if col not in known_order]
    ordered_columns = known_order + sorted(remaining)
    return df.loc[:, ordered_columns]


# ===== Public API ===========================================================

def process_targets(
    source: Path | str,
    *,
    offline: bool = False,
    session: requests.Session | None = None,
    sep: str = ",",
    encodings: Sequence[str] | None = None,
) -> Path:
    """Process aggregated target data and persist the SSoT compliant export.

    Parameters
    ----------
    source:
        Either a directory containing the ``output.target_*.csv`` snapshot or
        the full path to the file that should be processed.
    offline:
        When ``True`` the function skips all network lookups and relies solely
        on the lineage information present in the CSV.
    session:
        Optional :class:`requests.Session` instance used for network requests.
        When omitted a new session is created for each invocation.
    sep:
        Column separator used by the CSV file.  Defaults to a comma.
    encodings:
        Optional sequence of encodings attempted while reading the CSV.  When
        ``None`` the :data:`DEFAULT_ENCODINGS` catalogue is used.

    Returns
    -------
    pathlib.Path
        Path to the generated ``organism.output.target_*.csv`` file.
    """

    source_path = Path(source)
    if source_path.is_dir():
        candidates = sorted(source_path.glob(INPUT_PATTERN))
        if not candidates:
            raise FileNotFoundError(
                f"No files matching '{INPUT_PATTERN}' found in {source_path}"
            )
        input_path = candidates[-1]
    else:
        input_path = source_path
        if not input_path.exists():
            raise FileNotFoundError(f"Input file does not exist: {input_path}")

    base_dir = input_path.parent
    input_name = input_path.name
    if input_name.startswith("output."):
        suffix = input_name[len("output.") :]
        output_name = f"{OUTPUT_PREFIX}{suffix}"
    else:
        output_name = f"{OUTPUT_PREFIX}{input_name}"
    output_path = base_dir / output_name

    encodings = encodings or DEFAULT_ENCODINGS
    df = _read_csv(input_path, encodings=encodings, sep=sep)

    LOGGER.info("Loaded %d target rows", len(df))

    # Normalise core textual columns used by analytics.
    for column in [
        "target_chembl_id",
        "organism",
        "target_names",
        "gene_symbol",
        "protein_name_canonical",
        "protein_name_alt",
    ]:
        if column in df.columns:
            df[column] = df[column].map(NormalizeText)

    if "target_names" in df.columns:
        df["target_name_primary"] = df["target_names"].map(FirstElementText)

    cellularity_series = df.apply(
        lambda row: Cellularity.from_row(row, offline=offline, session=session),
        axis=1,
    )
    df["organism_cellularity"] = [
        c.value if isinstance(c, Cellularity) else ""
        for c in cellularity_series
    ]
    df["organism_cellularity_method"] = [
        c.method if isinstance(c, Cellularity) else "local"
        for c in cellularity_series
    ]
    df["organism_cellularity_lineage"] = [
        "|".join(c.lineage) if isinstance(c, Cellularity) else ""
        for c in cellularity_series
    ]

    df["is_multifunctional_enzyme"] = df.apply(multifunctional, axis=1)

    cellularity_counts = (
        df["organism_cellularity"].value_counts().sort_index().to_dict()
    )
    LOGGER.info("Cellularity distribution: %s", json.dumps(cellularity_counts, sort_keys=True))
    LOGGER.info(
        "Multifunctional enzymes: %d of %d",
        int(df["is_multifunctional_enzyme"].sum()),
        len(df),
    )

    df = _stabilise_columns(df)

    LOGGER.info("Writing processed target export to %s", output_path)
    base_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, encoding="utf-8", line_terminator="\n")
    return output_path


__all__ = [
    "NormalizeText",
    "FirstElementText",
    "Cellularity",
    "multifunctional",
    "process_targets",
]


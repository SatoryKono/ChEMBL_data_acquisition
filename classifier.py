"""Deterministic classifier for scientific document types.

The module implements normalization, scoring and conflict resolution logic
according to project specification. It exposes :class:`DocumentClassifier`
which accepts :class:`io_schemas.DocumentRecord` and returns
:class:`io_schemas.ClassificationResult`.
"""

from __future__ import annotations

import logging
from typing import Iterable, List, Tuple

import io_schemas as io
from rules import (
    LIST_SPLIT_RE,
    MESH_A,
    MESH_B,
    MESH_C,
    MESH_D,
    MESH_E,
    PT_CLINICAL,
    PT_OTHER,
    PT_REVIEW,
    RE_BEHAV,
    RE_CELLLINES,
    RE_CLINICAL,
    RE_DOSE,
    RE_IN_VITRO,
    RE_PROTOCOL,
    RE_RESULTS,
    RE_ROUTE,
)

logger = logging.getLogger(__name__)


def _norm_str(value: str) -> str:
    """Normalize strings: lowercase, collapse whitespace, replace micro sign."""

    value = value.lower().strip()
    value = value.replace("μ", "u")
    value = " ".join(value.split())
    return value


def _norm_list(values: Iterable[str]) -> List[str]:
    """Normalize a list of strings using ``LIST_SPLIT_RE``."""

    items: List[str] = []
    for v in values:
        for part in LIST_SPLIT_RE.split(v):
            part = _norm_str(part)
            if part:
                items.append(part)
    dedup = list(dict.fromkeys(items))
    dedup.sort()
    return dedup


def _text(record: io.DocumentRecord) -> str:
    return f"{_norm_str(record.title)} {_norm_str(record.abstract)}".strip()


class DocumentClassifier:
    """Main classifier implementing deterministic rules."""

    def classify(self, record: io.DocumentRecord) -> io.ClassificationResult:
        """Classify a single record.

        Parameters
        ----------
        record: DocumentRecord
            Input data model.
        """

        text = _text(record)
        pubmed_pt = _norm_list(record.pubmed_publicationtype)
        scholar_pt = _norm_list(record.scholar_publicationtypes)
        openalex_pt = _norm_list(record.openalex_publicationtypes)
        mesh_desc = _norm_list(
            record.pubmed_mesh_descriptors + record.openalex_meshdescriptors
        )
        mesh_qual = _norm_list(
            record.pubmed_mesh_qualifiers + record.openalex_meshqualifiers
        )
        chem_list = _norm_list(record.pubmed_chemicallist)

        pt_class, source_priority = self._class_from_pt(
            pubmed_pt, scholar_pt, openalex_pt
        )
        final_class = pt_class

        if final_class is None and RE_CLINICAL.search(text) and "humans" in mesh_desc:
            final_class = "Clinical trials"
            source_priority = "-"

        protocol_veto = bool(RE_PROTOCOL.search(text) and not RE_RESULTS.search(text))
        if protocol_veto:
            final_class = "Other non experimental publication"

        vivo_score, vivo_hits = self._score_vivo(
            mesh_desc, mesh_qual, text, chem_list, record.optional_experiment_kind
        )
        vitro_score, vitro_hits = self._score_vitro(
            mesh_desc,
            mesh_qual,
            text,
            chem_list,
            pubmed_pt,
            record.optional_experiment_kind,
        )

        if not protocol_veto and (
            final_class is None or final_class == "Other non experimental publication"
        ):
            experimental_class = self._resolve_experimental(
                vivo_score, vitro_score, vivo_hits, vitro_hits
            )
            if experimental_class:
                final_class = experimental_class

        if final_class is None:
            final_class = "Other non experimental publication"

        pt_comparative = "comparative study" in pubmed_pt
        pt_binding = bool("binding" in text)

        xrf_type = _norm_str(record.crossref_type or "")
        xrf_is_generic = xrf_type == "journal-article"
        xrf_is_specific = bool(xrf_type and not xrf_is_generic)

        conflicts: List[str] = []
        if final_class == "Clinical trials" and xrf_is_generic:
            conflicts.append("xrf_generic_vs_clinical")
        if (
            final_class in {"Review", "Other non experimental publication"}
            and xrf_is_generic
            and pt_class is not None
        ):
            conflicts.append("xrf_generic_vs_nonexperimental")
        if vivo_score >= 5 and vitro_score >= 4:
            conflicts.append("vivo_vs_vitro_tie")

        confidence = 100 - 20 * (
            len([c for c in conflicts if c != "vivo_vs_vitro_tie"])
        )
        if "vivo_vs_vitro_tie" in conflicts:
            confidence -= 10
        confidence = max(0, confidence)

        explain_short = f"{final_class} via {source_priority or 'rules'}"

        return io.ClassificationResult(
            final_class=final_class,
            S_in_vivo=vivo_score,
            S_in_vitro=vitro_score,
            pt_comparative=pt_comparative,
            pt_binding=pt_binding,
            vivo_hits=vivo_hits,
            vitro_hits=vitro_hits,
            conflicts=conflicts,
            confidence=confidence,
            source_priority_used=source_priority or "-",
            xrf_is_generic=xrf_is_generic,
            xrf_is_specific=xrf_is_specific,
            explain_short=explain_short,
        )

    # ------------------------------------------------------------------
    @staticmethod
    def _class_from_pt(
        pubmed_pt: List[str], scholar_pt: List[str], openalex_pt: List[str]
    ) -> Tuple[str | None, str]:
        """Determine primary class from publication types with priority."""

        for pts, source in (
            (pubmed_pt, "pubmed"),
            (scholar_pt, "scholar"),
            (openalex_pt, "openalex"),
        ):
            veterinary = "veterinary" in pts
            for pt in pts:
                if pt in PT_CLINICAL and not veterinary:
                    return "Clinical trials", source
            for pt in pts:
                if pt in PT_REVIEW:
                    return "Review", source
            for pt in pts:
                if pt in PT_OTHER:
                    return "Other non experimental publication", source
        return None, "-"

    # ------------------------------------------------------------------
    @staticmethod
    def _score_vivo(
        mesh_desc: List[str],
        mesh_qual: List[str],
        text: str,
        chem_list: List[str],
        optional_kind: str | None,
    ) -> Tuple[int, List[str]]:
        """Compute in vivo score and collect hits."""

        hits: List[str] = []
        A = int(any(m in MESH_A for m in mesh_desc))
        B = int(any(m in MESH_B for m in mesh_qual))
        C = int(any(m in MESH_C for m in mesh_desc))
        TxtRouteDose = int(bool(RE_ROUTE.search(text) or RE_DOSE.search(text)))
        TxtBehavior = int(bool(RE_BEHAV.search(text)))
        D_or_E = int(any(m in MESH_D or m in MESH_E for m in mesh_desc + mesh_qual))
        Neg_vivo = int(not A and D_or_E and not TxtRouteDose)

        score = (
            2 * A + 2 * B + 1 * C + 2 * TxtRouteDose + 1 * TxtBehavior - 2 * Neg_vivo
        )
        if optional_kind == "F" and A:
            score += 1

        if A:
            hits.append("A")
        if B:
            hits.append("B")
        if C:
            hits.append("C")
        if TxtRouteDose:
            hits.append("route/dose")
        if TxtBehavior:
            hits.append("behavior")
        if Neg_vivo:
            hits.append("neg_vivo")

        return score, hits

    # ------------------------------------------------------------------
    @staticmethod
    def _score_vitro(
        mesh_desc: List[str],
        mesh_qual: List[str],
        text: str,
        chem_list: List[str],
        pubmed_pt: List[str],
        optional_kind: str | None,
    ) -> Tuple[int, List[str]]:
        """Compute in vitro score and collect hits."""

        hits: List[str] = []
        G = int(any(m in MESH_D for m in mesh_desc))
        H = int(any(m in MESH_E for m in mesh_qual))
        text_hits = min(
            len(RE_IN_VITRO.findall(text)) + len(RE_CELLLINES.findall(text)), 3
        )
        J = int(
            (
                any(m in MESH_A for m in mesh_desc)
                or bool(RE_ROUTE.search(text) or RE_BEHAV.search(text))
            )
            and not (G or H)
        )
        Chem = int(bool(chem_list))
        pt_comparative = int("comparative study" in pubmed_pt)
        pt_binding = int("binding" in text)

        score = 2 * G + 2 * H + text_hits + Chem - 2 * J + pt_comparative + pt_binding
        if optional_kind == "B":
            score += 1

        if G:
            hits.append("D")
        if H:
            hits.append("E")
        if text_hits:
            hits.append("assay")
        if Chem:
            hits.append("chem")
        if J:
            hits.append("neg_vitro")
        if pt_comparative:
            hits.append("comparative")
        if pt_binding:
            hits.append("binding")

        return score, hits

    # ------------------------------------------------------------------
    @staticmethod
    def _resolve_experimental(
        vivo_score: int, vitro_score: int, vivo_hits: List[str], vitro_hits: List[str]
    ) -> str | None:
        """Resolve final experimental class from scores."""

        if vivo_score >= 5 and vitro_score < 4:
            return "Experimental in vivo bioactivity"
        if vitro_score >= 4 and vivo_score < 5:
            return "Experimental in vitro bioactivity"
        if vivo_score >= 5 and vitro_score >= 4:
            if "behavior" in vivo_hits:
                return "Experimental in vivo bioactivity"
            return "Experimental in vitro bioactivity"
        return None

"""Utilities for deriving protein classification predictions."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from ...config import Config
from ...iuphar.integration.iuphar_library import (
    ClassificationRecord,
    IUPHARClassifier,
    IUPHARData,
)

PREDICTION_COLUMNS: tuple[str, ...] = (
    "protein_class_pred_L1",
    "protein_class_pred_L2",
    "protein_class_pred_L3",
    "protein_class_pred_rule_id",
    "protein_class_pred_evidence",
    "protein_class_pred_confidence",
)

_CONFIDENCE_MAP: Mapping[str, str] = {
    "target_id": "high",
    "family_id": "medium",
    "ec_number": "medium",
    "molecular_function": "low",
    "name": "low",
    "N/A": "low",
}


def classifier_from_config(cfg: Config) -> IUPHARClassifier:
    data = IUPHARData.from_files(
        target_path=cfg.target.iuphar.target_csv,
        family_path=cfg.target.iuphar.family_csv,
        encoding=cfg.io.csv_encoding,
    )
    return IUPHARClassifier(data)


def _normalise(value: object | None) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na"}:
        return ""
    return text


def _select_evidence(row: pd.Series, status: str) -> str:
    evidence_map: Mapping[str, str] = {
        "target_id": _normalise(
            row.get("target_id")
            or row.get("IUPHAR_target_id")
            or row.get("GuidetoPHARMACOLOGY")
        ),
        "family_id": _normalise(row.get("IUPHAR_family_id") or row.get("family_id")),
        "ec_number": _normalise(row.get("ec_number")),
        "name": _normalise(
            row.get("iuphar_name") or row.get("IUPHAR_name") or row.get("pref_name")
        ),
        "molecular_function": _normalise(row.get("molecular_function")),
    }
    return evidence_map.get(status, "") or "-"


@dataclass
class Prediction:
    L1: str
    L2: str
    L3: str
    rule_id: str
    evidence: str
    confidence: str

    @classmethod
    def from_record(cls, record: ClassificationRecord, *, evidence: str) -> Prediction:
        status = record.STATUS or "N/A"
        confidence = _CONFIDENCE_MAP.get(status, "low")
        rule_id = status if status not in {"", "N/A"} else "-"
        return cls(
            L1=record.IUPHAR_class or "-",
            L2=record.IUPHAR_subclass or "-",
            L3=record.IUPHAR_type or "-",
            rule_id=rule_id,
            evidence=evidence,
            confidence=confidence,
        )

    def as_dict(self) -> dict[str, str]:
        return {
            "protein_class_pred_L1": self.L1 or "-",
            "protein_class_pred_L2": self.L2 or "-",
            "protein_class_pred_L3": self.L3 or "-",
            "protein_class_pred_rule_id": self.rule_id or "-",
            "protein_class_pred_evidence": self.evidence or "-",
            "protein_class_pred_confidence": self.confidence or "-",
        }


def append_protein_class_predictions(
    df: pd.DataFrame, classifier: IUPHARClassifier
) -> pd.DataFrame:
    result = df.copy()
    if result.empty:
        for column in PREDICTION_COLUMNS:
            result[column] = pd.Series(dtype="string")
        return result

    def _predict(row: pd.Series) -> pd.Series:
        target_id = _normalise(row.get("target_id"))
        family_id = _normalise(row.get("IUPHAR_family_id") or row.get("family_id"))
        ec_number = _normalise(row.get("ec_number"))
        name = _normalise(
            row.get("iuphar_name") or row.get("IUPHAR_name") or row.get("pref_name")
        )

        record = classifier.get(target_id, family_id, ec_number, name)
        if (record.STATUS in {"", "N/A"}) or not record.IUPHAR_class:
            molecular_function = _normalise(row.get("molecular_function"))
            if molecular_function:
                record = classifier.by_molecular_function(molecular_function)

        evidence = _select_evidence(row, record.STATUS or "N/A")
        prediction = Prediction.from_record(record, evidence=evidence)
        return pd.Series(prediction.as_dict())

    predictions = result.apply(_predict, axis=1)
    for column in PREDICTION_COLUMNS:
        result[column] = (
            predictions[column].fillna("-").replace("", "-").astype("string")
        )
    return result

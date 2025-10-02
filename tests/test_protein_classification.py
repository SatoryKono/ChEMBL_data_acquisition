"""Tests for :mod:`library.protein_classification`."""

from __future__ import annotations

import pandas as pd

from library import protein_classification as pc
from library.integration.iuphar_library import ClassificationRecord


class DummyClassifier:
    """Simple stand-in for :class:`IUPHARClassifier` used in tests."""

    def __init__(
        self,
        record: ClassificationRecord,
        *,
        fallback: ClassificationRecord | None = None,
    ) -> None:
        self.record = record
        self.fallback = fallback or record

    def get(self, *_args: object, **_kwargs: object) -> ClassificationRecord:
        return self.record

    def by_molecular_function(
        self, *_args: object, **_kwargs: object
    ) -> ClassificationRecord:
        return self.fallback


def test_append_predictions_uses_classifier_record() -> None:
    """Classifier output is translated into ``protein_class_pred_*`` columns."""

    record = ClassificationRecord(
        IUPHAR_target_id="123",
        IUPHAR_family_id="456",
        IUPHAR_class="Enzyme",
        IUPHAR_subclass="Transferase",
        IUPHAR_type="Enzyme.Transferase",
        IUPHAR_name="Example",
        STATUS="target_id",
    )
    classifier = DummyClassifier(record)
    df = pd.DataFrame({"target_id": ["123"], "ec_number": ["1.1.1.1"]})

    result = pc.append_protein_class_predictions(df, classifier)

    row = result.iloc[0]
    assert row["protein_class_pred_L1"] == "Enzyme"
    assert row["protein_class_pred_L2"] == "Transferase"
    assert row["protein_class_pred_L3"] == "Enzyme.Transferase"
    assert row["protein_class_pred_rule_id"] == "target_id"
    assert row["protein_class_pred_confidence"] == "high"
    assert row["protein_class_pred_evidence"] == "123"


def test_append_predictions_falls_back_to_molecular_function() -> None:
    """Missing classifications fall back to molecular function keywords."""

    empty = ClassificationRecord()
    fallback = ClassificationRecord(
        IUPHAR_class="Transporter",
        IUPHAR_subclass="SLC superfamily of solute carrier",
        IUPHAR_type="Transporter.SLC superfamily of solute carrier",
        STATUS="molecular_function",
    )
    classifier = DummyClassifier(empty, fallback=fallback)
    df = pd.DataFrame({"molecular_function": ["solute carrier activity"]})

    result = pc.append_protein_class_predictions(df, classifier)

    row = result.iloc[0]
    assert row["protein_class_pred_rule_id"] == "molecular_function"
    assert row["protein_class_pred_confidence"] == "low"
    assert row["protein_class_pred_L1"] == "Transporter"
    assert row["protein_class_pred_L2"] == "SLC superfamily of solute carrier"

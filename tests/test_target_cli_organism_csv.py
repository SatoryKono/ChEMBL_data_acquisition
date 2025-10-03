"""CLI integration tests for organism classification in the target pipeline."""

from __future__ import annotations

import pandas as pd
import pytest

from library.pipelines.target import protein_classification as pc
from library.config import Config
from scripts import get_target_data as gtd


class DummyRecord:
    """Lightweight stand-in for :class:`IUPHARClassificationRecord`."""

    def __init__(self, status: str = "target_id") -> None:
        self.STATUS = status
        self.IUPHAR_class = "ClassA"
        self.IUPHAR_subclass = "SubclassA"
        self.IUPHAR_type = "TypeA"


class DummyClassifier:
    """Minimal classifier returning deterministic predictions."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, tuple[str, ...]]] = []

    def get(
        self,
        target_id: str,
        family_id: str,
        ec_number: str,
        name: str,
    ) -> DummyRecord:
        self.calls.append(("get", (target_id, family_id, ec_number, name)))
        return DummyRecord()

    def by_molecular_function(self, molecular_function: str) -> DummyRecord:
        self.calls.append(("by_molecular_function", (molecular_function,)))
        return DummyRecord(status="molecular_function")


def test_all_subcommand_rejects_organism_csv() -> None:
    """The ``all`` sub-command no longer accepts ``--organism-csv``."""

    parser, _log_cfg = gtd.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["all", "--organism-csv", "ignored.csv"])


def test_merge_results_uses_classifier(
    monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """``merge_results`` loads the classifier when none is provided."""

    combined_df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprot_id": ["P12345"],
            "uniProtkbId": ["P12345-1"],
            "pref_name": ["Alpha"],
            "component_description": ["Alpha component"],
            "gene": ["GENEA|ALPHA"],
            "chembl_alternative_name": ["Alt"],
            "names": ["Name1|Name2"],
            "secondaryAccessionNames": ["Sec1"],
            "ec_numbers": ["1.1.1.1"],
            "reaction_ec_numbers": ["2.2.2.2"],
            "genus": ["Homo"],
            "lineage_superkingdom": ["Eukaryota"],
            "lineage_phylum": ["Chordata"],
            "lineage_class": ["Mammalia"],
        }
    )
    iuphar_df = pd.DataFrame(
        {
            "uniprot_id": ["P12345"],
            "IUPHAR_class": ["Kinase"],
            "target_id": ["T1"],
        }
    )

    sentinel = DummyClassifier()
    recorded: dict[str, object] = {}

    def fake_classifier(config: Config) -> DummyClassifier:
        recorded["config"] = config
        return sentinel

    def fake_append(df: pd.DataFrame, classifier: DummyClassifier) -> pd.DataFrame:
        recorded["classifier"] = classifier
        df = df.copy()
        for column in pc.PREDICTION_COLUMNS:
            df[column] = "ClassA"
        return df

    monkeypatch.setattr(pc, "classifier_from_config", fake_classifier)
    monkeypatch.setattr(pc, "append_protein_class_predictions", fake_append)

    orig_finalise = gtd.tp.finalise_targets

    def patched_finalise(df: pd.DataFrame, **kwargs: object) -> pd.DataFrame:
        df = df.drop(
            columns=[col for col in ("type", "target_type") if col in df.columns]
        )
        return orig_finalise(df, **kwargs)

    monkeypatch.setattr(gtd.tp, "finalise_targets", patched_finalise)

    result = gtd.merge_results(combined_df, iuphar_df, cfg)

    assert recorded["config"] is cfg
    assert recorded["classifier"] is sentinel
    for column in pc.PREDICTION_COLUMNS:
        assert column in result.columns
    assert result.loc[0, "target_type"] == "Multicellular organism"
    assert result.loc[0, "protein_class_pred_L1"] == "ClassA"

"""End-to-end tests for :func:`scripts.get_target_data.run_all`.

The test executes the combined target pipeline using small, static input
files for ChEMBL, UniProt and IUPHAR data. External HTTP requests are
mocked by replacing the individual source runners with functions that copy
pre-generated outputs. The final result is validated for column ordering and
expected content.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd
from pytest import MonkeyPatch

from library import protein_classification as pc
from library.config import Config
from library.constants import TargetsSchema
from library.constants import TARGETS_COLUMN_ORDER
from scripts import get_target_data as gtd


class DummyRecord:
    """Lightweight classification record for deterministic predictions."""

    def __init__(self, status: str = "target_id") -> None:
        self.STATUS = status
        self.IUPHAR_class = "ClassA"
        self.IUPHAR_subclass = "SubclassA"
        self.IUPHAR_type = "TypeA"


class DummyClassifier:
    """Minimal classifier returning deterministic values for tests."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, tuple[str, ...]]] = []

    def get(
        self, target_id: str, family_id: str, ec_number: str, name: str
    ) -> DummyRecord:
        self.calls.append(("get", (target_id, family_id, ec_number, name)))
        return DummyRecord()

    def by_molecular_function(self, molecular_function: str) -> DummyRecord:
        self.calls.append(("by_molecular_function", (molecular_function,)))
        return DummyRecord(status="molecular_function")


def test_prepare_targets_for_schema_adds_missing_columns() -> None:
    """Helper fills optional columns and preserves schema order."""

    df = pd.DataFrame({"target_chembl_id": ["CHEMBL1"]})

    prepared, missing_required, missing_optional = gtd._prepare_targets_for_schema(df)

    assert missing_required == set()
    assert "uniprot_id_primary" in prepared.columns
    assert prepared.columns.tolist() == TARGETS_COLUMN_ORDER
    assert prepared.loc[0, "uniprot_id_primary"] == "-"
    assert "uniprot_id_primary" in missing_optional
    assert list(df.columns) == ["target_chembl_id"]


def test_run_all_uses_local_inputs(
    tmp_path: Path, monkeypatch: MonkeyPatch, cfg: Config
) -> None:
    """Run ``run_all`` on local mini inputs and verify outputs.

    Parameters
    ----------
    tmp_path:
        Temporary directory provided by ``pytest``.
    monkeypatch:
        Fixture for dynamically modifying objects during the test.
    cfg:
        Base :class:`~library.config.Config` instance.
    """

    # ------------------------------------------------------------------
    # Configure paths for mini input data
    # ------------------------------------------------------------------
    chembl_data = Path("tests/data/chembl_targets_min.csv")
    uniprot_data = Path("tests/data/uniprot_targets_min.csv")
    iuphar_data = Path("tests/data/iuphar_targets_min.csv")
    original_chunk = cfg.target.chembl.chunk_size
    cfg.target.all.chunk_size = original_chunk + 2
    cfg.target.all.chembl_out = tmp_path / "chembl_out.csv"
    cfg.target.all.uniprot_out = tmp_path / "uniprot_out.csv"
    cfg.target.all.iuphar_out = tmp_path / "iuphar_out.csv"

    # ------------------------------------------------------------------
    # Patch network-dependent functions to use local files
    # ------------------------------------------------------------------
    recorded: dict[str, int] = {}

    def fake_run_chembl(cfg: Config, args: argparse.Namespace) -> int:
        recorded["chunk_size"] = cfg.target.chembl.chunk_size
        df = pd.read_csv(chembl_data)
        df["type"] = "Legacy"
        df.to_csv(args.output_csv, index=False)
        return 0

    def fake_run_uniprot(cfg: Config, args: argparse.Namespace) -> int:
        df = pd.read_csv(uniprot_data)
        df["lineage_superkingdom"] = "Eukaryota"
        df["lineage_phylum"] = "Chordata"
        df["lineage_class"] = "Mammalia"
        df.to_csv(args.output_csv, index=False)
        return 0

    def fake_run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
        shutil.copy(iuphar_data, args.output_csv)
        return 0

    monkeypatch.setattr(gtd, "run_chembl", fake_run_chembl)
    monkeypatch.setattr(gtd, "run_uniprot", fake_run_uniprot)
    monkeypatch.setattr(gtd, "run_iuphar", fake_run_iuphar)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *a, **k: None)
    classifier = DummyClassifier()
    monkeypatch.setattr(pc, "classifier_from_config", lambda cfg: classifier)

    orig_finalise = gtd.tp.finalise_targets

    def patched_finalise(df: pd.DataFrame, **kwargs: object) -> pd.DataFrame:
        df = df.drop(
            columns=[col for col in ("type", "target_type") if col in df.columns]
        )
        return orig_finalise(df, **kwargs)

    monkeypatch.setattr(gtd.tp, "finalise_targets", patched_finalise)

    # Skip strict schema validation; the goal is to exercise orchestration.
    monkeypatch.setattr(
        TargetsSchema, "validate", staticmethod(lambda df, lazy=True: df)
    )

    # ------------------------------------------------------------------
    # Prepare input identifiers and run the combined pipeline
    # ------------------------------------------------------------------
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding)
    output_csv = tmp_path / "out.csv"

    args = argparse.Namespace(input_csv=input_csv, output_csv=output_csv)
    exit_code = gtd.run_all(cfg, args)
    assert exit_code == 0
    assert recorded["chunk_size"] == cfg.target.all.chunk_size
    assert cfg.target.chembl.chunk_size == original_chunk

    # ------------------------------------------------------------------
    # Intermediate files are persisted and match expectations
    # ------------------------------------------------------------------
    assert cfg.target.all.chembl_out.exists()
    assert cfg.target.all.uniprot_out.exists()
    assert cfg.target.all.iuphar_out.exists()

    # ------------------------------------------------------------------
    # Validate final output
    # ------------------------------------------------------------------
    result = pd.read_csv(output_csv, dtype=str)
    assert result.columns.tolist() == TARGETS_COLUMN_ORDER
    row = result.loc[0]
    assert row["target_chembl_id"] == "CHEMBL1"
    assert row["uniprot_id_primary"] == "P12345"
    assert row["gene_symbol"] == "GENEA"
    assert row["target_type"] == "Multicellular organism"
    assert row["protein_class_pred_L1"] == "ClassA"
    assert row["protein_class_pred_rule_id"] == "target_id"
    assert classifier.calls
    assert (
        row["protein_synonym_list"]
        == "gene1|genea|sec1|alpha component|alt|recommended|name1|name2|alpha"
    )

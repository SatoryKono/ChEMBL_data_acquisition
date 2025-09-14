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

from library import target_postprocessing as tp
from library.config import Config
from schemas import TargetsSchema
from scripts import get_target_data as gtd


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
    organism_csv = Path("tests/data/organism_min.csv")

    cfg.target.all.organism_csv = organism_csv
    cfg.target.all.chembl_out = tmp_path / "chembl_out.csv"
    cfg.target.all.uniprot_out = tmp_path / "uniprot_out.csv"
    cfg.target.all.iuphar_out = tmp_path / "iuphar_out.csv"

    # ------------------------------------------------------------------
    # Patch network-dependent functions to use local files
    # ------------------------------------------------------------------
    def fake_run_chembl(cfg: Config, args: argparse.Namespace) -> int:
        shutil.copy(chembl_data, args.output_csv)
        return 0

    def fake_run_uniprot(cfg: Config, args: argparse.Namespace) -> int:
        shutil.copy(uniprot_data, args.output_csv)
        return 0

    def fake_run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
        shutil.copy(iuphar_data, args.output_csv)
        return 0

    monkeypatch.setattr(gtd, "run_chembl", fake_run_chembl)
    monkeypatch.setattr(gtd, "run_uniprot", fake_run_uniprot)
    monkeypatch.setattr(gtd, "run_iuphar", fake_run_iuphar)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *a, **k: None)

    # ``finalise_targets`` expects the ``type`` column to be missing prior to
    # merging with organism data. The wrapper removes it to mimic production
    # behaviour and keeps the original logic otherwise.
    orig_finalise = tp.finalise_targets

    def patched_finalise(
        df: pd.DataFrame, organism: pd.DataFrame, **kw: object
    ) -> pd.DataFrame:
        df = df.drop(columns=["type"], errors="ignore")
        return orig_finalise(df, organism, **kw)

    monkeypatch.setattr(tp, "finalise_targets", patched_finalise)

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
    assert result.columns.tolist() == sorted(result.columns)
    row = result.loc[0]
    assert row["target_chembl_id"] == "CHEMBL1"
    assert row["uniprotkb_Id"] == "P12345"
    assert row["gene_name"] == "GENEA"
    assert row["type"] == "Type1"
    assert (
        row["synonyms"]
        == "gene1|genea|sec1|alpha component|alt|recommended|name1|name2|alpha"
    )

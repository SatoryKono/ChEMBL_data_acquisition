"""Tests for data fetching utilities in get_target_data."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd
import yaml
from pytest import MonkeyPatch

from library import protein_classification as pc
from library.config import Config
from schemas import TargetsSchema
from schemas.targets import TARGETS_COLUMN_ORDER
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


def test_fetch_chembl(monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config) -> None:
    out = tmp_path / "chembl.csv"
    inp = tmp_path / "in.csv"

    def fake_run_chembl(cfg: Config, args: argparse.Namespace) -> int:
        pd.DataFrame({"target_chembl_id": ["C1"], "uniprot_id": ["P1"]}).to_csv(
            args.output_csv, index=False
        )
        return 0

    monkeypatch.setattr(gtd, "run_chembl", fake_run_chembl)
    df = gtd.fetch_chembl(cfg, inp, out)
    assert df.loc[0, "uniprot_id"] == "P1"


def test_fetch_chembl_custom_chunk_size(
    monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    out = tmp_path / "chembl.csv"
    inp = tmp_path / "in.csv"

    recorded: dict[str, int] = {}
    original_chunk = cfg.target.chembl.chunk_size

    def fake_run_chembl(cfg: Config, args: argparse.Namespace) -> int:
        recorded["chunk_size"] = cfg.target.chembl.chunk_size
        pd.DataFrame({"target_chembl_id": ["C1"]}).to_csv(args.output_csv, index=False)
        return 0

    monkeypatch.setattr(gtd, "run_chembl", fake_run_chembl)
    chunk_size = original_chunk + 3
    df = gtd.fetch_chembl(cfg, inp, out, chunk_size=chunk_size)

    assert df.loc[0, "target_chembl_id"] == "C1"
    assert recorded["chunk_size"] == chunk_size
    assert cfg.target.chembl.chunk_size == original_chunk


def test_run_chembl_streams_chunks(
    monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    cfg.target.chembl.chunk_size = 2
    cfg.target.chembl.limit = None
    cfg.target.chembl.column = "target_chembl_id"

    input_csv = tmp_path / "targets.csv"
    input_csv.write_text(
        "target_chembl_id\nCHEMBL1\nCHEMBL2\nCHEMBL3\nCHEMBL4\nCHEMBL5\n",
        encoding=cfg.io.csv_encoding,
    )
    output_csv = tmp_path / "chembl.csv"
    args = argparse.Namespace(input_csv=input_csv, output_csv=output_csv, offset=0)

    chunk_calls: list[list[str]] = []

    def fake_get_targets(
        chunk_ids: list[str],
        *,
        cfg: object,
        client: object,
        mapping_cfg: object,
        chunk_size: int,
        timeout: float | None,
    ) -> pd.DataFrame:
        ids = list(chunk_ids)
        chunk_calls.append(ids)
        return pd.DataFrame({"target_chembl_id": ids})

    class DummyClient:
        def __init__(self, *args: object, **kwargs: object) -> None:
            pass

        def __enter__(self) -> "DummyClient":
            return self

        def __exit__(
            self,
            exc_type: object,
            exc: object,
            tb: object,
        ) -> bool:
            return False

    frames_to_write: list[pd.DataFrame] = []

    def fake_write_csv(
        df: pd.DataFrame | list[pd.DataFrame],
        path: Path,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        chunksize: int | None = None,
    ) -> Path:
        assert not isinstance(df, pd.DataFrame)
        chunks = list(df)
        frames_to_write.extend(chunks)
        destination = Path(path)
        destination.write_text("dummy", encoding=encoding or cfg.io.csv_encoding)
        return destination

    monkeypatch.setattr(gtd.cl, "get_targets", fake_get_targets)
    monkeypatch.setattr(gtd, "ChemblClient", DummyClient)
    monkeypatch.setattr(gtd.io, "write_csv", fake_write_csv)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(
        TargetsSchema, "validate", staticmethod(lambda df, lazy=True: df)
    )

    exit_code = gtd.run_chembl(cfg, args)

    assert exit_code == 0
    assert chunk_calls == [
        ["CHEMBL1", "CHEMBL2"],
        ["CHEMBL3", "CHEMBL4"],
        ["CHEMBL5"],
    ]
    assert [df["target_chembl_id"].tolist() for df in frames_to_write] == chunk_calls


def test_fetch_uniprot(monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config) -> None:
    chembl_df = pd.DataFrame({"uniprot_id": ["P12345"]})
    out = tmp_path / "uniprot.csv"

    def fake_run_uniprot(cfg: Config, args: argparse.Namespace) -> int:
        pd.DataFrame({"uniprot_id": ["P12345"], "names": ["Foo"]}).to_csv(
            args.output_csv, index=False
        )
        return 0

    monkeypatch.setattr(gtd, "run_uniprot", fake_run_uniprot)
    df = gtd.fetch_uniprot(cfg, chembl_df, out)
    assert list(df["original_id"]) == ["P12345"]


def test_run_uniprot_initialises_session(
    monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("uniprot_id\nP12345\n", encoding=cfg.io.csv_encoding)
    output_csv = tmp_path / "uniprot.csv"
    cfg.target.uniprot.data_dir = tmp_path

    called: dict[str, object] = {}

    def fake_init_session(api: object, retry: object) -> None:
        called["init"] = (api, retry)

    def fake_process(
        input_csv: str,
        output_csv: str,
        data_dir: Path | str | None = None,
        *,
        cfg: object,
        gtop_cfg: object | None = None,
        sep: str = ",",
        encoding: str = "utf-8",
    ) -> None:
        called["cfg"] = cfg
        called["gtop_cfg"] = gtop_cfg
        pd.DataFrame({"uniprot_id": ["P12345"], "names": ["Foo"]}).to_csv(
            output_csv, index=False
        )

    monkeypatch.setattr(gtd.uu, "init_session", fake_init_session)
    monkeypatch.setattr(gtd.uu, "process", fake_process)

    def fake_write_csv(
        df: pd.DataFrame,
        path: Path,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
        **__: object,
    ) -> Path:
        return path

    monkeypatch.setattr(gtd.io, "write_csv", fake_write_csv)

    args = argparse.Namespace(input_csv=input_csv, output_csv=output_csv)
    rc = gtd.run_uniprot(cfg, args)

    assert rc == 0
    assert called["init"] == (cfg.api, cfg.retry)
    assert called["cfg"] is cfg.uniprot
    assert called["gtop_cfg"] is cfg.iuphar


def test_run_uniprot_writes_sidecar(
    monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("uniprot_id\nP12345\nP23456\n", encoding=cfg.io.csv_encoding)
    output_csv = tmp_path / "uniprot.csv"
    cfg.target.uniprot.data_dir = tmp_path

    def fake_init_session(api: object, retry: object) -> None:
        pass

    def fake_process(
        input_csv: str,
        output_csv: str,
        data_dir: Path | str | None = None,
        *,
        cfg: object,
        gtop_cfg: object | None = None,
        sep: str = ",",
        encoding: str = "utf-8",
    ) -> None:
        pd.DataFrame(
            {
                "uniprot_id": ["P23456", "P12345"],
                "names": ["Beta", "Alpha"],
            }
        ).to_csv(output_csv, index=False)

    monkeypatch.setattr(gtd.uu, "init_session", fake_init_session)
    monkeypatch.setattr(gtd.uu, "process", fake_process)

    args = argparse.Namespace(input_csv=input_csv, output_csv=output_csv)
    rc = gtd.run_uniprot(cfg, args)

    assert rc == 0

    df = pd.read_csv(output_csv)
    assert list(df["uniprot_id"]) == ["P12345", "P23456"]

    meta_path = output_csv.with_name(output_csv.name + ".meta.yaml")
    assert meta_path.exists()
    metadata = yaml.safe_load(meta_path.read_text(encoding="utf-8"))

    assert "uniprot_id" in metadata.get("columns", [])
    assert metadata["schema"] == "UniProtExport"

    stats = metadata["stats"]
    assert stats["rows_total"] == 2
    assert stats["rows_kept"] == 2
    assert stats["rows_dropped"] == 0
    assert len(stats["output_sha256"]) == 64

    inputs = metadata["inputs"]
    assert inputs["input_csv"] == str(input_csv)
    assert inputs["data_dir"] == str(cfg.target.uniprot.data_dir)


def test_fetch_iuphar(monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config) -> None:
    chembl_df = pd.DataFrame(
        {
            "target_chembl_id": ["C1"],
            "uniprot_id": ["P1"],
            "mapping_uniprot_id": ["P1|ALT1"],
            "pref_name": ["pref"],
            "component_description": ["desc"],
            "gene": ["GENE"],
            "chembl_alternative_name": ["alt"],
            "ec_numbers": ["1.1.1.1"],
            "reaction_ec_numbers": ["2.2.2.2"],
        }
    )
    uniprot_df = pd.DataFrame(
        {
            "uniprot_id": ["P1"],
            "original_id": ["P1"],
            "names": ["name"],
            "secondaryAccessionNames": ["sec"],
            "ec_numbers": ["3.3.3.3"],
            "reaction_ec_numbers": ["4.4.4.4"],
        }
    )
    out = tmp_path / "iuphar.csv"

    def fake_run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
        input_df = pd.read_csv(
            args.input_csv,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            dtype=str,
        )
        assert "mapping_uniprot_id" in input_df.columns
        assert input_df.loc[0, "mapping_uniprot_id"] == "P1|ALT1"
        pd.DataFrame({"uniprot_id": ["P1"], "IUPHAR_class": ["Enzyme"]}).to_csv(
            args.output_csv, index=False
        )
        return 0

    monkeypatch.setattr(gtd, "run_iuphar", fake_run_iuphar)
    combined_df, iuphar_df = gtd.fetch_iuphar(cfg, chembl_df, uniprot_df, out)
    assert "synonyms" in combined_df.columns
    assert iuphar_df.loc[0, "IUPHAR_class"] == "Enzyme"


def test_fetch_iuphar_missing_uniprot_column(
    monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    chembl_df = pd.DataFrame({"target_chembl_id": ["C1"]})
    uniprot_df = pd.DataFrame(columns=["uniprot_id", "original_id"])
    out = tmp_path / "iuphar.csv"

    def fake_run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
        pd.DataFrame(columns=["uniprot_id"]).to_csv(args.output_csv, index=False)
        return 0

    monkeypatch.setattr(gtd, "run_iuphar", fake_run_iuphar)
    combined_df, iuphar_df = gtd.fetch_iuphar(cfg, chembl_df, uniprot_df, out)

    merge_column = cfg.target.all.uniprot_column
    assert merge_column in combined_df.columns
    assert combined_df.loc[0, merge_column] == gtd.UNIPROT_MISSING_VALUE
    assert iuphar_df.empty


def test_run_all_preserves_reaction_ec_numbers(
    monkeypatch: MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    chembl_data = Path("tests/data/chembl_targets_min.csv")
    uniprot_data = Path("tests/data/uniprot_targets_min.csv")
    iuphar_data = Path("tests/data/iuphar_targets_min.csv")
    cfg.target.all.chembl_out = tmp_path / "chembl_out.csv"
    cfg.target.all.uniprot_out = tmp_path / "uniprot_out.csv"
    cfg.target.all.iuphar_out = tmp_path / "iuphar_out.csv"
    cfg.target.all.uniprot_column = "uniprot_id"

    def fake_run_chembl(cfg: Config, args: argparse.Namespace) -> int:
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
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)
    classifier = DummyClassifier()
    monkeypatch.setattr(pc, "classifier_from_config", lambda cfg: classifier)

    orig_finalise = gtd.tp.finalise_targets

    def patched_finalise(df: pd.DataFrame, **kwargs: object) -> pd.DataFrame:
        df = df.drop(
            columns=[col for col in ("type", "target_type") if col in df.columns]
        )
        return orig_finalise(df, **kwargs)

    monkeypatch.setattr(gtd.tp, "finalise_targets", patched_finalise)
    monkeypatch.setattr(
        TargetsSchema, "validate", staticmethod(lambda df, lazy=True: df)
    )

    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding)
    output_csv = tmp_path / "out.csv"

    args = argparse.Namespace(input_csv=input_csv, output_csv=output_csv)
    exit_code = gtd.run_all(cfg, args)
    assert exit_code == 0

    result = pd.read_csv(output_csv, dtype=str)
    assert result.loc[0, "reaction_ec_numbers"] == "2.2.2.2|4.4.4.4"
    assert result.loc[0, "target_type"] == "Multicellular organism"
    assert result.loc[0, "protein_class_pred_L1"] == "ClassA"
    assert result.loc[0, "protein_class_pred_rule_id"] == "target_id"
    assert classifier.calls


def test_validate_and_write_skips_quality_for_empty(
    tmp_path: Path, cfg: Config
) -> None:
    df = pd.DataFrame(columns=TARGETS_COLUMN_ORDER)
    output = tmp_path / "targets.csv"

    exit_code = gtd.validate_and_write(df, output, cfg)

    assert exit_code == 0
    assert output.exists()

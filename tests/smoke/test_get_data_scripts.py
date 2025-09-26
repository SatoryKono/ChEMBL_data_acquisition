from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Callable, Iterable

import pandas as pd
import pytest
from pandas.api import types as ptypes

from library import chembl_library as cl
from library import pubchem_library as pl
from scripts import (
    get_activity_data,
    get_assay_data,
    get_document_data,
    get_target_data,
    get_testitem_data,
)

TypeCheck = Callable[[pd.Series], bool]


def _cleanup_output(path: Path) -> None:
    """Remove previous outputs and sidecars for deterministic assertions."""

    candidates = [
        path,
        path.with_name(path.name + ".meta.yaml"),
        path.with_suffix(".quality.json"),
        path.with_name(f"{path.stem}_failure_cases.csv"),
    ]
    for candidate in candidates:
        if candidate.exists() and candidate.is_file():
            candidate.unlink()


def _assert_columns(df: pd.DataFrame, expected: Iterable[str]) -> None:
    missing = set(expected) - set(df.columns)
    assert not missing, f"missing columns: {sorted(missing)}"


def _assert_column_types(
    df: pd.DataFrame, expectations: dict[str, Iterable[TypeCheck] | TypeCheck]
) -> None:
    for column, validators in expectations.items():
        assert column in df.columns, f"column {column} not present"
        series = df[column]
        checks: Iterable[TypeCheck]
        if callable(validators):
            checks = (validators,)
        else:
            checks = tuple(validators)
        if not any(check(series) for check in checks):
            raise AssertionError(
                f"unexpected dtype for {column}: {series.dtype!s}"
            )


def test_get_activity_data_smoke(smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Smoke-test ``get_activity_data`` with bundled identifiers."""

    input_csv = Path("data/input-smoke/activity.csv")
    output_csv = smoke_output_dir / "activity.csv"
    _cleanup_output(output_csv)

    def fake_get_activities(ids, cfg, client, chunk_size, timeout):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for idx, raw_id in enumerate(ids, start=1):
            activity_id = int(str(raw_id))
            rows.append(
                {
                    "activity_id": activity_id,
                    "molecule_chembl_id": f"CHEMBL_MOL_{idx}",
                    "assay_chembl_id": f"CHEMBL_ASSAY_{idx}",
                    "standard_type": "IC50",
                    "standard_value": float(idx),
                    "document_chembl_id": f"CHEMBL_DOC_{idx}",
                }
            )
        return pd.DataFrame(rows)

    monkeypatch.setattr(cl, "get_activities", fake_get_activities)
    monkeypatch.setattr(get_activity_data, "analyze_table_quality", lambda *_, **__: None)

    exit_code = get_activity_data.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            str(Path("config.yaml")),
        ]
    )
    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(output_csv)
    assert not df.empty
    _assert_columns(
        df,
        [
            "activity_id",
            "molecule_chembl_id",
            "assay_chembl_id",
            "standard_type",
            "standard_value",
            "pipeline_version",
            "timestamp_utc",
        ],
    )
    _assert_column_types(
        df,
        {
            "activity_id": (ptypes.is_integer_dtype, ptypes.is_numeric_dtype),
            "standard_value": (ptypes.is_float_dtype, ptypes.is_numeric_dtype),
            "molecule_chembl_id": ptypes.is_object_dtype,
            "assay_chembl_id": ptypes.is_object_dtype,
            "standard_type": ptypes.is_object_dtype,
            "pipeline_version": ptypes.is_object_dtype,
            "timestamp_utc": ptypes.is_object_dtype,
        },
    )


def test_get_assay_data_smoke(smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Smoke-test ``get_assay_data`` using canned assay identifiers."""

    input_csv = Path("data/input-smoke/assay.csv")
    output_csv = smoke_output_dir / "assay.csv"
    _cleanup_output(output_csv)

    def fake_get_assays(ids, cfg, client, chunk_size, timeout):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for idx, assay_id in enumerate(ids, start=1):
            rows.append(
                {
                    "assay_chembl_id": str(assay_id),
                    "document_chembl_id": f"CHEMBL_DOC_{idx}",
                    "target_chembl_id": f"CHEMBL_TGT_{idx}",
                    "description": f"Assay description {idx}",
                }
            )
        return pd.DataFrame(rows)

    monkeypatch.setattr(cl, "get_assays", fake_get_assays)
    monkeypatch.setattr(get_assay_data.ap, "postprocess_assays", lambda df: df)
    monkeypatch.setattr(get_assay_data, "analyze_table_quality", lambda *_, **__: None)

    exit_code = get_assay_data.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            str(Path("config.yaml")),
        ]
    )
    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(output_csv)
    assert not df.empty
    _assert_columns(
        df,
        [
            "assay_chembl_id",
            "document_chembl_id",
            "target_chembl_id",
            "description",
            "pipeline_version",
            "timestamp_utc",
        ],
    )
    _assert_column_types(
        df,
        {
            "assay_chembl_id": ptypes.is_object_dtype,
            "document_chembl_id": ptypes.is_object_dtype,
            "target_chembl_id": ptypes.is_object_dtype,
            "description": ptypes.is_object_dtype,
            "pipeline_version": ptypes.is_object_dtype,
            "timestamp_utc": ptypes.is_object_dtype,
        },
    )


def test_get_document_data_smoke(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Smoke-test the ``chembl`` sub-command of ``get_document_data``."""

    input_csv = Path("data/input-smoke/documents.csv")
    output_csv = smoke_output_dir / "documents.csv"
    _cleanup_output(output_csv)

    def fake_get_documents(ids, cfg, client, chunk_size, timeout):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for idx, doc_id in enumerate(ids, start=1):
            rows.append(
                {
                    "document_chembl_id": str(doc_id),
                    "title": f"Document title {idx}",
                    "abstract": f"Summary {idx}",
                    "doi": f"10.1000/sample{idx}",
                    "year": 2020 + idx,
                    "journal": "Journal of Testing",
                    "publication_class": "experimental",
                }
            )
        return pd.DataFrame(rows)

    monkeypatch.setattr(cl, "get_documents", fake_get_documents)
    monkeypatch.setattr(get_document_data, "analyze_table_quality", lambda *_, **__: None)

    exit_code = get_document_data.main(
        [
            "chembl",
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            str(Path("config.yaml")),
        ]
    )
    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(output_csv)
    assert not df.empty
    _assert_columns(
        df,
        [
            "ChEMBL.document_chembl_id",
            "ChEMBL.title",
            "ChEMBL.doi",
            "ChEMBL.journal",
            "publication_class",
        ],
    )
    _assert_column_types(
        df,
        {
            "ChEMBL.document_chembl_id": ptypes.is_object_dtype,
            "ChEMBL.title": ptypes.is_object_dtype,
            "ChEMBL.doi": ptypes.is_object_dtype,
            "ChEMBL.journal": ptypes.is_object_dtype,
            "publication_class": ptypes.is_object_dtype,
        },
    )


def test_get_target_data_smoke(smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Smoke-test ``get_target_data`` using the ``chembl`` pipeline."""

    input_csv = Path("data/input-smoke/targets.csv")
    output_csv = smoke_output_dir / "targets.csv"
    _cleanup_output(output_csv)

    def fake_get_targets(ids, cfg, client, mapping_cfg, chunk_size, timeout):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for idx, target_id in enumerate(ids, start=1):
            rows.append(
                {
                    "target_chembl_id": str(target_id),
                    "gene_symbol": f"GENE{idx}",
                    "organism": "Homo sapiens",
                }
            )
        return pd.DataFrame(rows)

    monkeypatch.setattr(cl, "get_targets", fake_get_targets)
    monkeypatch.setattr(get_target_data, "analyze_table_quality", lambda *_, **__: None)

    exit_code = get_target_data.main(
        [
            "chembl",
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            str(Path("config.yaml")),
        ]
    )
    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(output_csv)
    assert not df.empty
    _assert_columns(
        df,
        [
            "target_chembl_id",
            "gene_symbol",
            "organism",
            "pipeline_version",
            "timestamp_utc",
        ],
    )
    _assert_column_types(
        df,
        {
            "target_chembl_id": ptypes.is_object_dtype,
            "gene_symbol": ptypes.is_object_dtype,
            "organism": ptypes.is_object_dtype,
            "pipeline_version": ptypes.is_object_dtype,
            "timestamp_utc": ptypes.is_object_dtype,
        },
    )


def test_get_testitem_data_smoke(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Smoke-test ``get_testitem_data`` with PubChem augmentation patched."""

    input_csv = Path("data/input-smoke/testitem.csv")
    output_csv = smoke_output_dir / "testitem.csv"
    _cleanup_output(output_csv)

    def fake_get_testitem(ids, cfg, client, chunk_size, timeout):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for idx, mol_id in enumerate(ids, start=1):
            rows.append(
                {
                    "molecule_chembl_id": str(mol_id),
                    "canonical_smiles": f"C{idx}H{2 * idx + 2}",
                    "max_phase": "4",
                    "standard_inchi": f"InChI=1S/C{idx}H{2 * idx + 2}",
                }
            )
        return pd.DataFrame(rows)

    def fake_get_cid(smiles: str, cfg):  # type: ignore[no-untyped-def]
        return "12345"

    def fake_get_properties(cid: str, cfg):  # type: ignore[no-untyped-def]
        return SimpleNamespace(
            IUPACName="Mock acid",
            MolecularFormula="C2H6O",
            iSMILES="CCO",
            cSMILES="CCO",
            InChI="InChI=1S/C2H6O",
            InChIKey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        )

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(pl, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pl, "get_cid_from_smiles", fake_get_cid)
    monkeypatch.setattr(pl, "get_properties", fake_get_properties)
    monkeypatch.setattr(
        get_testitem_data,
        "load_parent_catalog",
        lambda **__: {"CHEMBL1": "CHEMBL1_PARENT"},
    )
    monkeypatch.setattr(get_testitem_data, "analyze_table_quality", lambda *_, **__: None)

    exit_code = get_testitem_data.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            str(Path("config.yaml")),
        ]
    )
    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(output_csv)
    assert not df.empty
    _assert_columns(
        df,
        [
            "molecule_chembl_id",
            "parent_molecule_chembl_id",
            "canonical_smiles",
            "pubchem_cid",
            "pubchem_inchi",
            "pipeline_version",
            "timestamp_utc",
        ],
    )
    _assert_column_types(
        df,
        {
            "molecule_chembl_id": ptypes.is_object_dtype,
            "parent_molecule_chembl_id": (
                ptypes.is_object_dtype,
                ptypes.is_float_dtype,
            ),
            "canonical_smiles": ptypes.is_object_dtype,
            "pubchem_cid": (ptypes.is_object_dtype, ptypes.is_numeric_dtype),
            "pubchem_inchi": ptypes.is_object_dtype,
            "pipeline_version": ptypes.is_object_dtype,
            "timestamp_utc": ptypes.is_object_dtype,
        },
    )

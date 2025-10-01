from __future__ import annotations

import argparse
import io
import json
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from pandas.api import types as ptypes

from library import chembl_library as cl
from library import pubchem_library as pl
from library.clients import pubchem as pc
from library.logging_setup import LoggerConfig, configure_logger
from library.utils.config import DEFAULT_CONFIG_RELATIVE
from scripts import (
    get_activity_data,
    get_assay_data,
    get_document_data,
    get_data,
    get_target_data,
    get_testitem_data,
)

TypeCheck = Callable[[pd.Series], bool]

CONFIG_CLI_PATH = str(DEFAULT_CONFIG_RELATIVE)

DEFAULT_DATE = "20240101"


def _cli_args(
    *extra: str, date: str = DEFAULT_DATE, log_level: str = "ERROR"
) -> list[str]:
    """Return base CLI arguments extended with ``extra`` tokens."""

    return [
        "--base-path",
        ".",
        "--input-dir",
        "data/input-smoke",
        "--output-dir",
        "data/output-smoke",
        "--date",
        date,
        "--log-level",
        log_level,
        "--config",
        CONFIG_CLI_PATH,
        *extra,
    ]


def _expected_output(output_dir: Path, stem: str, *, date: str = DEFAULT_DATE) -> Path:
    """Return the expected output path for ``stem`` and ``date``."""

    return output_dir / f"output.{stem}_{date}.csv"


def test_get_data_main_smoke(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure ``scripts.get_data`` orchestrates all pipelines with mocked steps."""

    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()

    config_path = base_path / "config.yaml"
    config_path.write_text("test: true\n")

    invocations: list[str] = []

    def make_stub(
        name: str, subcommand: str | None
    ) -> Callable[[Sequence[str] | None], int]:
        def _main(argv: Sequence[str] | None) -> int:
            args = list(argv or [])
            if subcommand is not None:
                assert args, f"missing subcommand for {name}"
                assert args[0] == subcommand
                args = args[1:]
            parser = argparse.ArgumentParser(add_help=False)
            parser.add_argument("--config")
            parser.add_argument("--input")
            parser.add_argument("--output")
            parser.add_argument("--log-level")
            ns = parser.parse_args(args)
            input_path = Path(ns.input)
            output_path = Path(ns.output)
            assert input_path.exists(), f"missing input for {name}"
            output_path.write_text(f"{name} output\n")
            invocations.append(name)
            return 0

        return _main

    for step in get_data._PIPELINE_STEPS:
        input_path = input_dir / get_data._DEFAULT_INPUT_FILES[step.name]
        input_path.write_text("id\n1\n")

    stub_steps = tuple(
        get_data.PipelineStep(
            step.name, make_stub(step.name, step.subcommand), step.subcommand
        )
        for step in get_data._PIPELINE_STEPS
    )
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", stub_steps)

    class DummyLogger:
        def info(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def debug(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def warning(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def error(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

    dummy_logger = DummyLogger()

    def fake_configure(level_name: str, *, run_id: str | None = None) -> DummyLogger:
        return dummy_logger

    monkeypatch.setattr(get_data, "_configure_logging", fake_configure)

    exit_code = get_data.main(
        [
            "--base-path",
            str(base_path),
            "--input-dir",
            "input",
            "--output-dir",
            "output",
            "--config",
            str(config_path),
            "--date",
            "20240101",
            "--log-level",
            "ERROR",
        ]
    )

    assert exit_code == 0
    assert invocations == [step.name for step in stub_steps]
    for step in stub_steps:
        path = output_dir / (
            f"output.{get_data._DEFAULT_OUTPUT_STEMS[step.name]}_20240101.csv"
        )
        assert path.exists()
    assert path.read_text() == f"{step.name} output\n"


def test_get_data_limit_zero_skips_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``--limit 0`` should forward the skip to every pipeline without writes."""

    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()

    config_path = base_path / "config.yaml"
    config_path.write_text("test: true\n")

    for step in get_data._PIPELINE_STEPS:
        input_path = input_dir / get_data._DEFAULT_INPUT_FILES[step.name]
        input_path.write_text("id\nCHEMBL1\n", encoding="utf8")

    invocations: list[str] = []

    def make_stub(
        name: str, subcommand: str | None
    ) -> Callable[[Sequence[str] | None], int]:
        def _main(argv: Sequence[str] | None) -> int:
            args = list(argv or [])
            if subcommand is not None:
                assert args, f"missing subcommand for {name}"
                assert args[0] == subcommand
                args = args[1:]
            parser = argparse.ArgumentParser(add_help=False)
            parser.add_argument("--config")
            parser.add_argument("--input")
            parser.add_argument("--output")
            parser.add_argument("--log-level")
            parser.add_argument("--limit")
            ns = parser.parse_args(args)
            assert Path(ns.input).exists(), f"missing input for {name}"
            assert ns.limit == "0"
            invocations.append(name)
            return 0

        return _main

    stub_steps = tuple(
        get_data.PipelineStep(
            step.name, make_stub(step.name, step.subcommand), step.subcommand
        )
        for step in get_data._PIPELINE_STEPS
    )
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", stub_steps)

    class DummyLogger:
        def info(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def debug(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def warning(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def error(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

    dummy_logger = DummyLogger()

    def fake_configure(level_name: str, *, run_id: str | None = None) -> DummyLogger:
        return dummy_logger

    monkeypatch.setattr(get_data, "_configure_logging", fake_configure)

    exit_code = get_data.main(
        [
            "--base-path",
            str(base_path),
            "--input-dir",
            "input",
            "--output-dir",
            "output",
            "--config",
            str(config_path),
            "--date",
            "20240101",
            "--log-level",
            "ERROR",
            "--limit",
            "0",
        ]
    )

    assert exit_code == 0
    assert invocations == [step.name for step in stub_steps]

    for step in stub_steps:
        final_output = output_dir / (
            f"output.{get_data._DEFAULT_OUTPUT_STEMS[step.name]}_20240101.csv"
        )
        working_output = final_output.with_name(f".{final_output.name}.tmp")
        sentinel = final_output.with_name(f"{final_output.name}.failed")
        stemless = final_output.with_suffix("")
        assert not final_output.exists()
        assert not Path(f"{final_output}.meta.yaml").exists()
        assert not final_output.with_name(f"{final_output.stem}_failure_cases.csv").exists()
        assert not final_output.with_suffix(".quality.json").exists()
        assert not Path(f"{stemless}_quality_report_table.csv").exists()
        assert not Path(f"{stemless}_data_correlation_report_table.csv").exists()
        assert not working_output.exists()
        assert not sentinel.exists()


def test_get_data_forwards_skip_existing_flag(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify orchestrator forwards ``--skip-existing`` to every pipeline."""

    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()

    config_path = base_path / "config.yaml"
    config_path.write_text("test: true\n")

    recorded_args: dict[str, list[str]] = {}

    def make_stub(
        name: str, subcommand: str | None
    ) -> Callable[[Sequence[str] | None], int]:
        def _main(argv: Sequence[str] | None) -> int:
            args = list(argv or [])
            if subcommand is not None:
                assert args, f"missing subcommand for {name}"
                assert args[0] == subcommand
                args = args[1:]
            recorded_args[name] = args
            parser = argparse.ArgumentParser(add_help=False)
            parser.add_argument("--config")
            parser.add_argument("--input")
            parser.add_argument("--output")
            parser.add_argument("--log-level")
            parser.add_argument("--limit")
            parser.add_argument("--force", action="store_true")
            parser.add_argument("--skip-existing", action="store_true")
            ns = parser.parse_args(args)
            Path(ns.output).write_text(f"{name} output\n")
            return 0

        return _main

    for step in get_data._PIPELINE_STEPS:
        input_path = input_dir / get_data._DEFAULT_INPUT_FILES[step.name]
        input_path.write_text("id\n1\n")

    stub_steps = tuple(
        get_data.PipelineStep(
            step.name, make_stub(step.name, step.subcommand), step.subcommand
        )
        for step in get_data._PIPELINE_STEPS
    )
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", stub_steps)

    class DummyLogger:
        def info(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def debug(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def warning(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def error(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

    dummy_logger = DummyLogger()

    def fake_configure(level_name: str, *, run_id: str | None = None) -> DummyLogger:
        return dummy_logger

    monkeypatch.setattr(get_data, "_configure_logging", fake_configure)

    exit_code = get_data.main(
        [
            "--base-path",
            str(base_path),
            "--input-dir",
            "input",
            "--output-dir",
            "output",
            "--config",
            str(config_path),
            "--date",
            "20240101",
            "--log-level",
            "ERROR",
            "--skip-existing",
        ]
    )

    assert exit_code == 0
    assert set(recorded_args) == {step.name for step in stub_steps}
    for args in recorded_args.values():
        assert "--skip-existing" in args
        assert "--force" not in args


def test_get_data_pipeline_events_include_run_id(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Structured log events include the configured ``run_id``."""

    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()

    config_path = base_path / "config.yaml"
    config_path.write_text("test: true\n")

    for step in get_data._PIPELINE_STEPS:
        input_path = input_dir / get_data._DEFAULT_INPUT_FILES[step.name]
        input_path.write_text("id\n1\n")

    def make_stub(
        name: str, subcommand: str | None
    ) -> Callable[[Sequence[str] | None], int]:
        def _main(argv: Sequence[str] | None) -> int:
            args = list(argv or [])
            if subcommand is not None:
                assert args
                assert args[0] == subcommand
                args = args[1:]
            parser = argparse.ArgumentParser(add_help=False)
            parser.add_argument("--config")
            parser.add_argument("--input")
            parser.add_argument("--output")
            parser.add_argument("--log-level")
            ns = parser.parse_args(args)
            Path(ns.output).write_text(f"{name} output\n")
            return 0

        return _main

    stub_steps = tuple(
        get_data.PipelineStep(
            step.name, make_stub(step.name, step.subcommand), step.subcommand
        )
        for step in get_data._PIPELINE_STEPS
    )
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", stub_steps)

    stream = io.StringIO()
    expected_run_id = "test-run-id"

    def fake_configure(level_name: str, *, run_id: str | None = None):
        return configure_logger(
            LoggerConfig(
                level=level_name.upper(),
                run_id=run_id or expected_run_id,
                stream=stream,
            )
        )

    monkeypatch.setattr(get_data, "_configure_logging", fake_configure)

    exit_code = get_data.main(
        [
            "--base-path",
            str(base_path),
            "--input-dir",
            "input",
            "--output-dir",
            "output",
            "--config",
            str(config_path),
            "--date",
            "20240101",
            "--log-level",
            "INFO",
        ]
    )

    assert exit_code == 0

    records = [
        json.loads(line) for line in stream.getvalue().splitlines() if line.strip()
    ]

    events = {record["event"] for record in records}
    assert {
        "pipeline_start",
        "pipeline_done",
        "workflow_complete",
        "workflow_succeeded",
    } <= events

    interesting = [
        record
        for record in records
        if record["event"].startswith("pipeline_")
        or record["event"].startswith("step_")
        or record["event"].startswith("workflow_")
    ]
    assert interesting
    assert all(record["run_id"] == expected_run_id for record in interesting)


def _cleanup_output(path: Path) -> None:
    """Remove previous outputs and sidecars for deterministic assertions."""

    candidates = [
        path,
        path.with_name(path.name + ".meta.yaml"),
        path.with_suffix(".quality.json"),
        path.with_name(f"{path.stem}_failure_cases.csv"),
        Path(f"{path.with_suffix('')}_quality_report_table.csv"),
        Path(f"{path.with_suffix('')}_data_correlation_report_table.csv"),
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
            raise AssertionError(f"unexpected dtype for {column}: {series.dtype!s}")


def test_get_activity_data_smoke(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Smoke-test ``get_activity_data`` with bundled identifiers."""

    output_csv = _expected_output(smoke_output_dir, "activities")
    _cleanup_output(output_csv)

    def fake_get_activities(ids, cfg, client, chunk_size, timeout, **kwargs):  # type: ignore[no-untyped-def]
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
    monkeypatch.setattr(
        get_activity_data, "analyze_table_quality", lambda *_, **__: None
    )

    exit_code = get_activity_data.main(_cli_args())
    assert exit_code == 0
    assert output_csv.exists()
    assert output_csv.name == f"output.activities_{DEFAULT_DATE}.csv"
    assert output_csv.parent == smoke_output_dir

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


def test_get_assay_data_smoke(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Smoke-test ``get_assay_data`` using canned assay identifiers."""

    output_csv = _expected_output(smoke_output_dir, "assays")
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

    exit_code = get_assay_data.main(_cli_args())
    assert exit_code == 0
    assert output_csv.exists()
    assert output_csv.name == f"output.assays_{DEFAULT_DATE}.csv"
    assert output_csv.parent == smoke_output_dir

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

    output_csv = _expected_output(smoke_output_dir, "documents")
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
    monkeypatch.setattr(
        get_document_data, "analyze_table_quality", lambda *_, **__: None
    )

    exit_code = get_document_data.main(["chembl", *_cli_args()])
    assert exit_code == 0
    assert output_csv.exists()
    assert output_csv.name == f"output.documents_{DEFAULT_DATE}.csv"
    assert output_csv.parent == smoke_output_dir

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


def test_get_target_data_smoke(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Smoke-test ``get_target_data`` using the ``chembl`` pipeline."""

    output_csv = _expected_output(smoke_output_dir, "targets")
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

    exit_code = get_target_data.main(["chembl", *_cli_args()])
    assert exit_code == 0
    assert output_csv.exists()
    assert output_csv.name == f"output.targets_{DEFAULT_DATE}.csv"
    assert output_csv.parent == smoke_output_dir

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
    smoke_output_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Smoke-test ``get_testitem_data`` with PubChem augmentation patched."""

    output_csv = _expected_output(smoke_output_dir, "testitems")
    _cleanup_output(output_csv)

    polymer_smiles = "POLY-SMILES"
    mixture_smiles = "MIX-SMILES"

    def fake_get_testitem(ids, cfg, client, chunk_size, timeout, **kwargs):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for idx, mol_id in enumerate(ids, start=1):
            base: dict[str, object] = {
                "molecule_chembl_id": str(mol_id),
                "max_phase": "4",
                "standard_inchi": f"InChI=1S/C{idx}H{2 * idx + 2}",
            }
            if idx == 1:
                base.update(
                    {
                        "molecule_type": "Polymer",
                        "canonical_smiles": polymer_smiles,
                        "pubchem_cid": "POLY-CID",
                        "pubchem_iupac_name": "polymer existing",
                    }
                )
            elif idx == 2:
                base.update(
                    {
                        "molecule_type": "Mixture",
                        "canonical_smiles": mixture_smiles,
                        "pubchem_cid": "MIX-CID",
                        "pubchem_iupac_name": "mixture existing",
                    }
                )
            else:
                base.update(
                    {
                        "molecule_type": "Small molecule",
                        "canonical_smiles": f"C{idx}H{2 * idx + 2}",
                    }
                )
            rows.append(base)
        return pd.DataFrame(rows)

    def fake_get_cid(smiles: str, cfg):  # type: ignore[no-untyped-def]
        smiles_calls.append(smiles)
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
    monkeypatch.setattr(pc, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pl, "get_cid_from_smiles", fake_get_cid)
    monkeypatch.setattr(pl, "get_properties", fake_get_properties)
    monkeypatch.setattr(
        get_testitem_data,
        "load_parent_catalog",
        lambda **__: {"CHEMBL1": "CHEMBL1_PARENT"},
    )
    monkeypatch.setattr(
        get_testitem_data, "analyze_table_quality", lambda *_, **__: None
    )

    original_apply = get_testitem_data.cli.apply_config_overrides

    def patched_apply(*args, **kwargs):  # type: ignore[no-untyped-def]
        cfg = original_apply(*args, **kwargs)
        cfg.pubchem.allow_polymer = False
        return cfg

    monkeypatch.setattr(get_testitem_data.cli, "apply_config_overrides", patched_apply)

    smiles_calls: list[str] = []
    warning_events: list[tuple[str, dict[str, object]]] = []
    monkeypatch.setattr(
        get_testitem_data.logger,
        "warning",
        lambda event, **kwargs: warning_events.append((event, kwargs)),
    )

    exit_code = get_testitem_data.main(_cli_args())
    assert exit_code == 0
    assert output_csv.exists()
    assert output_csv.name == f"output.testitems_{DEFAULT_DATE}.csv"
    assert output_csv.parent == smoke_output_dir
    assert polymer_smiles not in smiles_calls
    assert mixture_smiles not in smiles_calls
    assert any(event == "pubchem_skip_polymers" for event, _ in warning_events)

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

    polymer_row = df.loc[df["molecule_type"] == "Polymer"].iloc[0]
    assert polymer_row["pubchem_cid"] == "POLY-CID"
    mixture_row = df.loc[df["molecule_type"] == "Mixture"].iloc[0]
    assert mixture_row["pubchem_cid"] == "MIX-CID"


def test_get_activity_data_skip_existing(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure ``--skip-existing`` prevents pipeline execution."""

    output_csv = _expected_output(smoke_output_dir, "activities")
    _cleanup_output(output_csv)
    output_csv.write_text("sentinel")

    called = False

    def _unexpected_call(*_args, **_kwargs):  # type: ignore[no-untyped-def]
        nonlocal called
        called = True
        raise AssertionError("pipeline should have been skipped")

    monkeypatch.setattr(cl, "get_activities", _unexpected_call)

    exit_code = get_activity_data.main(_cli_args("--skip-existing"))
    assert exit_code == 0
    assert not called
    assert output_csv.read_text() == "sentinel"


def test_get_activity_data_force_overrides_skip(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify ``--force`` re-runs the pipeline despite ``--skip-existing``."""

    output_csv = _expected_output(smoke_output_dir, "activities")
    _cleanup_output(output_csv)
    output_csv.write_text("sentinel")

    def fake_get_activities(ids, cfg, client, chunk_size, timeout, **kwargs):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for raw_id in ids:
            rows.append(
                {
                    "activity_id": int(str(raw_id)),
                    "molecule_chembl_id": "CHEMBL_FORCE",
                    "assay_chembl_id": "CHEMBL_ASSAY_FORCE",
                    "standard_type": "IC50",
                    "standard_value": 1.0,
                    "document_chembl_id": "CHEMBL_DOC_FORCE",
                }
            )
        return pd.DataFrame(rows)

    monkeypatch.setattr(cl, "get_activities", fake_get_activities)
    monkeypatch.setattr(
        get_activity_data, "analyze_table_quality", lambda *_, **__: None
    )

    exit_code = get_activity_data.main(_cli_args("--skip-existing", "--force"))
    assert exit_code == 0
    assert output_csv.exists()
    assert output_csv.read_text() != "sentinel"


def test_run_pipeline_failure_removes_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure failed orchestrator steps clean artefacts and emit a sentinel."""

    base_dir = tmp_path
    input_dir = base_dir / "input"
    output_dir = base_dir / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path = base_dir / "config.yaml"
    config_path.write_text("pipeline: test\n")

    input_csv = input_dir / "activity.csv"
    input_csv.write_text("activity_id\n1\n")

    cfg = get_data.PipelineRunConfig(
        base_path=base_dir,
        input_dir=input_dir,
        output_dir=output_dir,
        config_path=config_path,
        date_prefix="20240101",
        log_level="ERROR",
        limit=None,
        force=False,
        skip_existing=False,
    )

    final_output = cfg.output_path("activity")
    working_output = final_output.with_name(f".{final_output.name}.tmp")
    sentinel_path = final_output.with_name(f"{final_output.name}.failed")

    def failing_main(argv: list[str] | None) -> int:
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument("--output")
        parser.add_argument("--input")
        parser.add_argument("--config")
        parser.add_argument("--log-level")
        ns = parser.parse_args(argv)
        Path(ns.output).write_text("temporary output\n")
        final_output.write_text("unexpected final output\n")
        failure_tmp = Path(ns.output).with_name(
            f"{Path(ns.output).stem}_failure_cases.csv"
        )
        failure_tmp.write_text("errors\n", encoding="utf-8")
        return 2

    class DummyLogger:
        def info(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def debug(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def error(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def warning(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def exception(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

    monkeypatch.setattr(
        get_data,
        "_PIPELINE_STEPS",
        (get_data.PipelineStep("activity", failing_main, None),),
    )
    monkeypatch.setattr(get_data, "_LOGGER", DummyLogger())

    status = get_data.run_pipeline(cfg)

    assert status == 2
    assert not final_output.exists()
    assert not working_output.exists()
    assert sentinel_path.exists()

    failure_final = final_output.with_name(f"{final_output.stem}_failure_cases.csv")
    failure_working = working_output.with_name(
        f"{working_output.stem}_failure_cases.csv"
    )
    assert failure_final.exists()
    assert not failure_working.exists()


def test_run_pipeline_system_exit_cleans_up(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """SystemExit raised by a step should trigger cleanup and reporting."""

    base_dir = tmp_path
    input_dir = base_dir / "input"
    output_dir = base_dir / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path = base_dir / "config.yaml"
    config_path.write_text("pipeline: test\n")

    input_csv = input_dir / "activity.csv"
    input_csv.write_text("activity_id\n1\n")

    cfg = get_data.PipelineRunConfig(
        base_path=base_dir,
        input_dir=input_dir,
        output_dir=output_dir,
        config_path=config_path,
        date_prefix="20240101",
        log_level="ERROR",
        limit=None,
        force=False,
        skip_existing=False,
    )

    final_output = cfg.output_path("activity")
    working_output = final_output.with_name(f".{final_output.name}.tmp")
    sentinel_path = final_output.with_name(f"{final_output.name}.failed")

    def aborting_main(argv: list[str] | None) -> int:
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument("--output")
        parser.add_argument("--input")
        parser.add_argument("--config")
        parser.add_argument("--log-level")
        ns = parser.parse_args(argv)
        Path(ns.output).write_text("temporary output\n")
        raise SystemExit(3)

    class DummyLogger:
        def info(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def debug(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def error(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def warning(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def exception(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

    monkeypatch.setattr(
        get_data,
        "_PIPELINE_STEPS",
        (get_data.PipelineStep("activity", aborting_main, None),),
    )
    monkeypatch.setattr(get_data, "_LOGGER", DummyLogger())

    status = get_data.run_pipeline(cfg)

    assert status == 3
    assert not final_output.exists()
    assert not working_output.exists()
    assert sentinel_path.exists()


def test_run_pipeline_success_promotes_sidecars(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure successful steps promote sidecars from the working prefix."""

    base_dir = tmp_path
    input_dir = base_dir / "input"
    output_dir = base_dir / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    config_path = base_dir / "config.yaml"
    config_path.write_text("pipeline: test\n")

    input_csv = input_dir / "activity.csv"
    input_csv.write_text("activity_id\n1\n")

    cfg = get_data.PipelineRunConfig(
        base_path=base_dir,
        input_dir=input_dir,
        output_dir=output_dir,
        config_path=config_path,
        date_prefix="20240101",
        log_level="ERROR",
        limit=None,
        force=False,
        skip_existing=False,
    )

    final_output = cfg.output_path("activity")
    working_output = final_output.with_name(f".{final_output.name}.tmp")

    def successful_main(argv: list[str] | None) -> int:
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument("--output")
        parser.add_argument("--input")
        parser.add_argument("--config")
        parser.add_argument("--log-level")
        ns = parser.parse_args(argv)
        out_path = Path(ns.output)
        out_path.write_text("activity_id\n1\n")
        meta_path = out_path.with_name(out_path.name + ".meta.yaml")
        meta_path.write_text("meta: 1\n", encoding="utf-8")
        quality_json = out_path.with_suffix(".quality.json")
        quality_json.write_text("{}\n", encoding="utf-8")
        quality_table = Path(f"{out_path.with_suffix('')}_quality_report_table.csv")
        quality_table.write_text("col\n", encoding="utf-8")
        corr_table = Path(
            f"{out_path.with_suffix('')}_data_correlation_report_table.csv"
        )
        corr_table.write_text("corr\n", encoding="utf-8")
        return 0

    class DummyLogger:
        def info(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def debug(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def error(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def warning(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

        def exception(self, *args, **kwargs) -> None:  # pragma: no cover - simple stub
            pass

    monkeypatch.setattr(
        get_data,
        "_PIPELINE_STEPS",
        (get_data.PipelineStep("activity", successful_main, None),),
    )
    monkeypatch.setattr(get_data, "_LOGGER", DummyLogger())

    status = get_data.run_pipeline(cfg)

    assert status == 0
    assert final_output.exists()
    assert not working_output.exists()

    meta_final = final_output.with_name(final_output.name + ".meta.yaml")
    quality_json_final = final_output.with_suffix(".quality.json")
    quality_table_final = Path(
        f"{final_output.with_suffix('')}_quality_report_table.csv"
    )
    corr_table_final = Path(
        f"{final_output.with_suffix('')}_data_correlation_report_table.csv"
    )

    assert meta_final.exists()
    assert quality_json_final.exists()
    assert quality_table_final.exists()
    assert corr_table_final.exists()

    meta_working = working_output.with_name(working_output.name + ".meta.yaml")
    quality_json_working = working_output.with_suffix(".quality.json")
    quality_table_working = Path(
        f"{working_output.with_suffix('')}_quality_report_table.csv"
    )
    corr_table_working = Path(
        f"{working_output.with_suffix('')}_data_correlation_report_table.csv"
    )

    assert not meta_working.exists()
    assert not quality_json_working.exists()
    assert not quality_table_working.exists()
    assert not corr_table_working.exists()

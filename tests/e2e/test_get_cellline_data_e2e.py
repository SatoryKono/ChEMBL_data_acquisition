from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from library.config import Config
from library.pipelines.cellline import CellLinePipelineOptions, run_cellline_pipeline
from library.pipelines.cellline.chembl import CELL_LINE_COLUMN_ORDER
from library.pipelines.common.metadata import get_pipeline_version, get_timestamp_utc
from scripts import get_cellline_data


class _MemoryLogger:
    """Capture structured log events emitted by the CLI."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def info(self, event: str, **payload: object) -> None:
        self.events.append(("info", event, dict(payload)))

    def warning(self, event: str, **payload: object) -> None:
        self.events.append(("warning", event, dict(payload)))

    def error(self, event: str, **payload: object) -> None:
        self.events.append(("error", event, dict(payload)))


class _StubChemblClient:
    def __init__(self, payloads: dict[str, dict[str, object]]) -> None:
        self.payloads = payloads

    def request_json(
        self, url: str, *, cfg, timeout: float | None = None
    ) -> dict[str, object]:
        records: list[dict[str, object]] = []
        for key, payload in self.payloads.items():
            if key in url:
                records.extend(payload.get("cell_lines", []))
        return {"cell_lines": records, "page_meta": {}}

    def close(self) -> None:  # pragma: no cover - compatibility
        return None


@pytest.fixture()
def cellline_payloads() -> dict[str, dict[str, object]]:
    return {
        "CHEMBL3307636": {
            "cell_lines": [
                {
                    "cell_chembl_id": "CHEMBL3307636",
                    "cell_name": "22Rv1",
                    "cell_description": "Example cell line",
                    "cell_id": 101,
                    "cell_source_organism": "Homo sapiens",
                    "cell_source_tax_id": 9606,
                    "cell_source_tissue": "Prostate",
                    "cellosaurus_id": "CVCL_1045",
                    "cl_lincs_id": "LCL-1234",
                    "clo_id": "CLO_0001001",
                    "efo_id": "EFO_0001663",
                }
            ],
            "page_meta": {},
        },
        "CHEMBL3307790": {
            "cell_lines": [
                {
                    "cell_chembl_id": "CHEMBL3307790",
                    "cell_name": "Granta",
                    "cell_description": "Granta cell line",
                    "cell_id": 820,
                    "cell_source_organism": "Homo sapiens",
                    "cell_source_tax_id": 9606,
                    "cell_source_tissue": None,
                    "cellosaurus_id": None,
                    "cl_lincs_id": None,
                    "clo_id": None,
                    "efo_id": None,
                }
            ],
            "page_meta": {},
        },
    }


def test_get_cellline_data_main_success(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cellline_payloads: dict[str, dict[str, object]],
) -> None:
    input_csv = tmp_path / "cellline.csv"
    input_csv.write_text(
        "cell_chembl_id\nCHEMBL3307636\nCHEMBL3307790\n",
        encoding="utf-8",
    )
    output_csv = tmp_path / "output.csv"

    logger_stub = _MemoryLogger()
    monkeypatch.setattr(get_cellline_data, "logger", logger_stub)

    original_pipeline = run_cellline_pipeline

    def fake_run_pipeline(cfg: Config, options: CellLinePipelineOptions):
        client = _StubChemblClient(cellline_payloads)
        return original_pipeline(cfg, options, client=client)

    monkeypatch.setattr(get_cellline_data, "run_cellline_pipeline", fake_run_pipeline)

    def fake_run_cli_command(
        *, args, parser, log_cfg, mapping, run, logger, base_parser=None
    ):
        args._config_metadata = None
        if isinstance(args.input_csv, (str, Path)):
            args.input_csv = Path(args.input_csv)
        if isinstance(args.final_out, (str, Path)):
            args.final_out = Path(args.final_out)
        cfg = Config()
        return run(cfg, args)

    monkeypatch.setattr(get_cellline_data, "run_cli_command", fake_run_cli_command)

    @contextmanager
    def fake_setup_cli_logging(*args, **kwargs):
        yield SimpleNamespace(log_cfg=SimpleNamespace(level="INFO", run_id="test"))

    monkeypatch.setattr(get_cellline_data, "setup_cli_logging", fake_setup_cli_logging)

    exit_code = get_cellline_data.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
            "--limit",
            "2",
            "--offset",
            "0",
        ]
    )

    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(
        output_csv,
        dtype=str,
        keep_default_na=False,
    ).replace({"": pd.NA})
    df = df.reindex(columns=CELL_LINE_COLUMN_ORDER)

    expected = pd.DataFrame(
        [
            {
                "cell_chembl_id": "CHEMBL3307636",
                "cell_name": "22Rv1",
                "cell_description": "Example cell line",
                "cell_id": "101",
                "cell_source_organism": "Homo sapiens",
                "cell_source_tax_id": "9606",
                "cell_source_tissue": "Prostate",
                "cellosaurus_id": "CVCL_1045",
                "cl_lincs_id": "LCL-1234",
                "clo_id": "CLO_0001001",
                "efo_id": "EFO_0001663",
                "pipeline_version": get_pipeline_version(),
                "timestamp_utc": get_timestamp_utc(),
            },
            {
                "cell_chembl_id": "CHEMBL3307790",
                "cell_name": "Granta",
                "cell_description": "Granta cell line",
                "cell_id": "820",
                "cell_source_organism": "Homo sapiens",
                "cell_source_tax_id": "9606",
                "cell_source_tissue": pd.NA,
                "cellosaurus_id": pd.NA,
                "cl_lincs_id": pd.NA,
                "clo_id": pd.NA,
                "efo_id": pd.NA,
                "pipeline_version": get_pipeline_version(),
                "timestamp_utc": get_timestamp_utc(),
            },
        ]
    )
    expected = expected.astype("string")
    pd.testing.assert_frame_equal(df.astype("string"), expected)

    summary_events = [
        event for event in logger_stub.events if event[1] == "cellline_pipeline_summary"
    ]
    assert summary_events
    assert summary_events[0][2]["status"] == "OK"

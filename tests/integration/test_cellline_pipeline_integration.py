from pathlib import Path

import pandas as pd
import pytest
import requests

from library.config import Config
from library.pipelines.cellline import (
    CELL_LINE_COLUMN_ORDER,
    CellLinePipelineOptions,
    run_cellline_pipeline,
)
from library.pipelines.common.metadata import (
    get_pipeline_version,
    get_timestamp_utc,
)


class _StubChemblClient:
    """Return canned responses for the cell line endpoint."""

    def __init__(self, payloads: dict[str, dict[str, object]]) -> None:
        self.payloads = payloads
        self.calls: list[str] = []

    def request_json(
        self, url: str, *, cfg, timeout: float | None = None
    ) -> dict[str, object]:
        self.calls.append(url)
        records: list[dict[str, object]] = []
        for key, payload in self.payloads.items():
            if key in url:
                records.extend(payload.get("cell_lines", []))
        return {"cell_lines": records, "page_meta": {}}

    def close(self) -> None:  # pragma: no cover - interface compatibility
        return None


@pytest.fixture()
def sample_payloads() -> dict[str, dict[str, object]]:
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


def test_run_cellline_pipeline__writes_expected_csv(
    tmp_path: Path,
    cfg: Config,
    sample_payloads: dict[str, dict[str, object]],
) -> None:
    input_csv = tmp_path / "cellline.csv"
    input_csv.write_text(
        "cell_chembl_id\nCHEMBL3307636\nCHEMBL3307790\n",
        encoding="utf-8",
    )
    output_csv = tmp_path / "output.csv"
    options = CellLinePipelineOptions(
        input_csv=input_csv,
        output_csv=output_csv,
        column=cfg.cellline.column,
        batch_size=1,
    )

    client = _StubChemblClient(sample_payloads)

    result = run_cellline_pipeline(cfg, options, client=client)

    assert result.exit_code == 0
    assert result.records == 2
    assert result.written is True
    assert result.missing_ids == ()
    assert client.calls  # ensure endpoint was queried

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


def test_run_cellline_pipeline__api_error(tmp_path: Path, cfg: Config) -> None:
    input_csv = tmp_path / "cellline.csv"
    input_csv.write_text("cell_chembl_id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"
    options = CellLinePipelineOptions(
        input_csv=input_csv,
        output_csv=output_csv,
        column=cfg.cellline.column,
        batch_size=1,
    )

    class _ErrorClient:
        def request_json(self, *args, **kwargs):
            raise requests.RequestException("boom")

        def close(self) -> None:  # pragma: no cover - compatibility
            return None

    client = _ErrorClient()

    result = run_cellline_pipeline(cfg, options, client=client)

    assert result.exit_code == 1
    assert result.written is False
    assert result.records == 0

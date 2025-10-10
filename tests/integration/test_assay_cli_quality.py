from __future__ import annotations

import argparse
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pandas as pd
import pytest

from library.config import Config
from library.pipelines.assay.chembl_assay import MAX_ASSAY_CHUNK_SIZE
from scripts import get_assay_data
from tests.helpers import ASSAY_ENRICHMENT_MIN_RATIO


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


@pytest.mark.integration
def test_get_assay_cli__enrichment_quality(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    data_dir = Path(__file__).resolve().parents[1] / "data"
    input_csv = tmp_path / "assay.csv"
    input_csv.write_text(
        (data_dir / "assay.csv").read_text(encoding="utf-8"), encoding="utf-8"
    )
    output_csv = tmp_path / "out" / "assays.csv"
    dictionary_path = data_dir / "assay_dictionary.csv"

    def _stub_run_chembl(config: Config, args: argparse.Namespace) -> int:
        frame = pd.read_csv(args.input_csv)
        dictionary = pd.read_csv(dictionary_path)
        dictionary["assay_chembl_id"] = dictionary["assay_chembl_id"].astype("string")
        enriched = frame.merge(dictionary, on="assay_chembl_id", how="left")
        enriched["description"] = enriched["description"].astype("string").str.strip()
        enriched["description_length"] = (
            enriched["description"].str.len().astype("Int64")
        )
        enriched["year"] = pd.to_numeric(enriched["year"], errors="coerce").astype(
            "Int64"
        )
        quality_columns = ["assay_strain", "assay_group", "year", "accession"]
        completeness = 1.0 - enriched[quality_columns].isna().mean()
        if float(completeness.min()) < ASSAY_ENRICHMENT_MIN_RATIO:
            raise AssertionError(
                "assay enrichment below threshold "
                f"(threshold={ASSAY_ENRICHMENT_MIN_RATIO}, completeness={completeness.to_dict()})"
            )
        _ensure_parent(args.final_out)
        columns = [
            "assay_chembl_id",
            "target_chembl_id",
            "document_chembl_id",
            "description",
            "description_length",
            "assay_strain",
            "assay_group",
            "year",
            "accession",
        ]
        enriched.to_csv(args.final_out, index=False, columns=columns)
        return 0

    monkeypatch.setattr(get_assay_data, "run_chembl", _stub_run_chembl)

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
    )

    exit_code = get_assay_data.run(cfg, args)

    assert exit_code == 0
    result = pd.read_csv(output_csv)
    expected_columns = [
        "assay_chembl_id",
        "target_chembl_id",
        "document_chembl_id",
        "description",
        "description_length",
        "assay_strain",
        "assay_group",
        "year",
        "accession",
    ]
    assert list(result.columns) == expected_columns
    assert (result["description_length"] == result["description"].str.len()).all()
    quality_columns = ["assay_strain", "assay_group", "year", "accession"]
    completeness = 1.0 - result[quality_columns].isna().mean()
    assert (completeness >= ASSAY_ENRICHMENT_MIN_RATIO).all(), completeness.to_dict()


@pytest.mark.integration
def test_get_assay_cli__clamps_batch_size(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    cfg.assay.batch_size = MAX_ASSAY_CHUNK_SIZE * 2

    identifiers = [f"CHEMBL{i}" for i in range(1, MAX_ASSAY_CHUNK_SIZE + 4)]
    request_chunk_sizes: list[int] = []

    class _StubChemblClient:
        def __init__(self, *args, **kwargs) -> None:
            del args, kwargs

        def __enter__(self) -> _StubChemblClient:
            return self

        def __exit__(self, exc_type, exc, tb) -> bool:  # noqa: D401 - context protocol
            return False

        def close(self) -> None:
            return None

        def request_json(
            self,
            url: str,
            *,
            cfg: Config,
            timeout: float | None,
        ) -> dict[str, object]:
            del cfg, timeout
            parsed = urlparse(url)
            params = parse_qs(parsed.query)
            ids_param = params.get("assay_chembl_id__in")
            if ids_param:
                identifiers_chunk = [
                    value for value in ids_param[0].split(",") if value
                ]
                request_chunk_sizes.append(len(identifiers_chunk))
                return {
                    "assays": [
                        {"assay_chembl_id": identifier}
                        for identifier in identifiers_chunk
                    ],
                    "page_meta": {},
                }
            identifier_param = params.get("assay_chembl_id")
            if identifier_param:
                identifier = identifier_param[0]
                return {
                    "assays": [{"assay_chembl_id": identifier}],
                    "page_meta": {},
                }
            return {"assays": [], "page_meta": {}}

    class _StubChunkFailureTracker:
        def __init__(self) -> None:
            self.stats: dict[str, object] = {}

        def add_failure(self, chunk_ids: Sequence[str], error_message: str) -> None:
            del chunk_ids, error_message

        def save(self, path: Path, cfg: Config) -> None:
            del path, cfg

    written_frames: list[pd.DataFrame] = []

    def _fake_prepare_chunked_pipeline(
        *,
        fetch_config,
        fetch_chunk,
        csv_writer,
    ):
        del csv_writer
        ids_list = list(fetch_config.ids)

        def _fetcher():
            frames: list[pd.DataFrame] = []
            for start in range(0, len(ids_list), fetch_config.chunk_size):
                chunk = ids_list[start : start + fetch_config.chunk_size]
                frames.append(fetch_chunk(chunk))
            return frames

        def _writer(
            frames: Iterable[pd.DataFrame],
            destination: Path,
            col_order: Sequence[str] | None,
            key_cols: Sequence[str],
        ) -> Path:
            del col_order, key_cols
            chunk_list = list(frames)
            written_frames.extend(chunk_list)
            return Path(destination)

        return _fetcher, _writer

    def _fake_run_pipeline(
        *,
        definition,
        fetcher,
        output_path,
        failure_path,
        **kwargs,
    ) -> int:
        del failure_path, kwargs
        frames = fetcher()
        definition.writer(frames, output_path, None, ())
        return 0

    def _fake_read_ids(path: Path, *, column: str, cfg) -> Iterator[str]:
        del path, column, cfg
        return iter(identifiers)

    monkeypatch.setattr(
        "library.orchestration.context.ChemblClient",
        _StubChemblClient,
    )
    monkeypatch.setattr(get_assay_data, "ChunkFailureTracker", _StubChunkFailureTracker)
    monkeypatch.setattr(
        get_assay_data, "prepare_chunked_pipeline", _fake_prepare_chunked_pipeline
    )
    monkeypatch.setattr(get_assay_data, "run_pipeline", _fake_run_pipeline)
    monkeypatch.setattr(get_assay_data.io, "read_ids", _fake_read_ids)

    input_csv = tmp_path / "assay_ids.csv"
    input_csv.write_text("assay_chembl_id\n" + "\n".join(identifiers), encoding="utf-8")
    output_csv = tmp_path / "assays.csv"

    args = argparse.Namespace(
        input_csv=input_csv,
        final_out=output_csv,
        skip_existing=False,
        force=False,
        offset=0,
        limit=None,
        batch_size=cfg.assay.batch_size,
        timeout=cfg.assay.timeout,
    )

    exit_code = get_assay_data.run(cfg, args)

    assert exit_code == 0
    assert request_chunk_sizes, "Expected get_assays to perform batched requests"
    assert max(request_chunk_sizes) == MAX_ASSAY_CHUNK_SIZE
    remainder = len(identifiers) - MAX_ASSAY_CHUNK_SIZE
    assert sorted(request_chunk_sizes) == sorted([MAX_ASSAY_CHUNK_SIZE, remainder])
    assert written_frames, "Expected writer to receive frames"

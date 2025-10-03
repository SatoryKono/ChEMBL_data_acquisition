from __future__ import annotations

import argparse
import logging
import shutil
from collections.abc import Callable, Sequence
from pathlib import Path

import pandas as pd
import pandas.api.types as ptypes
import pytest

from library.schemas.normalize import (
    normalize_activities,
    normalize_assays,
    normalize_documents,
    normalize_targets,
    normalize_testitems,
)
from scripts import get_data
from scripts.get_data import PipelineStep


PostProcess = Callable[[pd.DataFrame, logging.Logger], pd.DataFrame]


def _make_stub_pipeline(
    name: str,
    *,
    required_columns: set[str],
    normalizer: Callable[[pd.DataFrame], pd.DataFrame],
    key_columns: Sequence[str],
    column_order: Sequence[str],
    postprocess: PostProcess | None = None,
    warning_sink: list[str] | None = None,
) -> Callable[[Sequence[str] | None], int]:
    """Return a deterministic stub ``main`` callable for a pipeline step."""

    logger = logging.getLogger(f"chembl_da.tests.stub.{name}")

    def _main(argv: Sequence[str] | None = None) -> int:
        args = list(argv or [])
        if args and not args[0].startswith("-"):
            command = args.pop(0)
            if command != "all":
                raise SystemExit(f"unsupported subcommand: {command}")
        parser = argparse.ArgumentParser(prog=f"stub_{name}")
        parser.add_argument("--config", dest="config", required=False)
        parser.add_argument("--input", dest="input_csv", required=True)
        parser.add_argument("--final-out", dest="final_out", required=False)
        parser.add_argument("--output", dest="final_out", required=False)
        parser.add_argument("--log-level", dest="log_level", required=False)
        parser.add_argument("--limit", dest="limit", type=int, required=False)
        parser.add_argument("--force", dest="force", action="store_true")
        parser.add_argument(
            "--skip-existing", dest="skip_existing", action="store_true"
        )
        parser.add_argument("--dry-run", dest="dry_run", action="store_true")
        parsed = parser.parse_args(args)
        if parsed.dry_run:
            logger.info("dry_run", extra={"step": name})
            return 0
        if parsed.limit == 0:
            logger.info("limit_skip", extra={"step": name})
            return 0

        input_path = Path(parsed.input_csv)
        output_path = Path(parsed.final_out)
        frame = pd.read_csv(input_path)
        missing = required_columns - set(frame.columns)
        if missing:
            logger.error(
                "schema_missing_columns",
                extra={"step": name, "missing": sorted(missing)},
            )
            return 1

        normalised = normalizer(frame)
        processed = normalised if postprocess is None else postprocess(normalised, logger)

        duplicates = processed.duplicated(subset=key_columns)
        if duplicates.any():
            if warning_sink is not None:
                warning_sink.append("deduplicated_rows")
            logger.warning(
                "deduplicated_rows",
                extra={"step": name, "removed": int(duplicates.sum())},
            )
            processed = processed.loc[~duplicates].copy()

        processed = processed.loc[:, list(column_order)]
        processed = processed.sort_values(by=list(key_columns)).reset_index(drop=True)

        for column in processed.columns:
            series = processed[column]
            if ptypes.is_string_dtype(series.dtype) or ptypes.is_object_dtype(series.dtype):
                processed[column] = series.fillna("")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        processed.to_csv(output_path, index=False)
        return 0

    return _main



@pytest.mark.e2e
def test_get_data_e2e__miniature_pipeline(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    base_dir = tmp_path / "chembl"
    input_dir = base_dir / "input"
    output_dir = base_dir / "output"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)

    fixture_dir = Path(__file__).resolve().parents[1] / "data"
    expected_dir = Path(__file__).resolve().parents[1] / "resources" / "expected_get_data"
    resource_dir = Path(__file__).resolve().parents[1] / "resources"

    for name in ("document", "target", "assay", "testitem", "activity"):
        shutil.copy(fixture_dir / f"{name}.csv", input_dir / f"{name}.csv")

    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_dir))
    monkeypatch.setattr(get_data, "_warm_parent_catalog", lambda cfg: None)

    warning_events: list[str] = []

    def _documents_postprocess(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
        result = df.copy()
        result["pmid"] = result["pmid"].astype("string")
        result["title"] = result["title"].astype("string").str.strip()
        missing = result["title"].isna() | (result["title"] == "")
        if missing.any():
            warning_events.append("documents_missing_title")
            logger.warning(
                "documents_missing_title",
                extra={"pmids": sorted(result.loc[missing, "pmid"].tolist())},
            )
            result.loc[missing, "title"] = "Unknown Title"
        result["relation"] = result["relation"].astype("string")
        result["units"] = result["units"].astype("string")
        return result

    def _targets_postprocess(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
        result = df.copy()
        result["target_chembl_id"] = result["target_chembl_id"].astype("string")
        result["target_pref_name"] = (
            result["target_pref_name"].astype("string").str.strip()
        )
        missing_name = result["target_pref_name"].isna() | (result["target_pref_name"] == "")
        if missing_name.any():
            warning_events.append("targets_missing_name")
            logger.warning(
                "targets_missing_name",
                extra={
                    "target_ids": sorted(
                        result.loc[missing_name, "target_chembl_id"].tolist()
                    )
                },
            )
            result.loc[missing_name, "target_pref_name"] = "Unknown Target"
        tax_map = {"Human": 9606}
        mapped = result["organism"].map(tax_map)
        result["tax_id"] = pd.Series(mapped, dtype="Int64")
        missing_tax = result["tax_id"].isna()
        if missing_tax.any():
            warning_events.append("targets_missing_taxonomy")
            logger.warning(
                "targets_missing_taxonomy",
                extra={
                    "target_ids": sorted(
                        result.loc[missing_tax, "target_chembl_id"].tolist()
                    )
                },
            )
        return result

    def _assays_postprocess(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
        result = df.copy()
        result["assay_chembl_id"] = result["assay_chembl_id"].astype("string")
        result["description"] = result["description"].astype("string").str.strip()
        missing = result["description"].isna() | (result["description"] == "")
        if missing.any():
            warning_events.append("assays_missing_description")
            logger.warning(
                "assays_missing_description",
                extra={
                    "assay_ids": sorted(
                        result.loc[missing, "assay_chembl_id"].tolist()
                    )
                },
            )
            result.loc[missing, "description"] = "N/A"
        result["relation"] = result["relation"].astype("string")
        result["units"] = result["units"].astype("string")
        return result

    hierarchy_path = resource_dir / "molecule_hierarchy.csv"
    catalog_path = resource_dir / "molecule_catalog.csv"

    def _testitems_postprocess(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
        result = df.copy()
        result["molecule_chembl_id"] = result["molecule_chembl_id"].astype("string")
        hierarchy = pd.read_csv(hierarchy_path)
        catalog = pd.read_csv(catalog_path)
        result = result.merge(hierarchy, on="molecule_chembl_id", how="left")
        result = result.merge(
            catalog[["molecule_chembl_id", "natural_product"]],
            on="molecule_chembl_id",
            how="left",
        )
        missing_parent = result["parent_molecule_chembl_id"].isna() | (
            result["parent_molecule_chembl_id"] == ""
        )
        if missing_parent.any():
            warning_events.append("testitems_missing_parent")
            logger.warning(
                "testitems_missing_parent",
                extra={
                    "molecule_ids": sorted(
                        result.loc[missing_parent, "molecule_chembl_id"].tolist()
                    )
                },
            )
        result["relation"] = result["relation"].astype("string")
        result["units"] = result["units"].astype("string")
        result["natural_product"] = result["natural_product"].astype("string")
        result["parent_molecule_chembl_id"] = (
            result["parent_molecule_chembl_id"].astype("string").fillna("")
        )
        return result

    def _activities_postprocess(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
        result = df.copy()
        result["activity_id"] = result["activity_id"].astype("string")
        result["value"] = pd.to_numeric(result["value"], errors="coerce")
        missing = result["value"].isna()
        if missing.any():
            warning_events.append("activities_missing_value")
            logger.warning(
                "activities_missing_value",
                extra={
                    "activity_ids": sorted(
                        result.loc[missing, "activity_id"].tolist()
                    )
                },
            )
            result.loc[missing, "value"] = -1
        result["value"] = result["value"].astype("Int64")
        result["relation"] = result["relation"].astype("string")
        result["standard_units"] = result["standard_units"].astype("string")
        return result

    monkeypatch.setattr(
        get_data,
        "_PIPELINE_STEPS",
        (
            PipelineStep(
                "document",
                _make_stub_pipeline(
                    "document",
                    required_columns={"pmid", "title", "relation", "units"},
                    normalizer=normalize_documents,
                    key_columns=["pmid"],
                    column_order=["pmid", "title", "relation", "units"],
                    postprocess=_documents_postprocess,
                    warning_sink=warning_events,
                ),
                "all",
            ),
            PipelineStep(
                "target",
                _make_stub_pipeline(
                    "target",
                    required_columns={
                        "target_chembl_id",
                        "target_pref_name",
                        "organism",
                    },
                    normalizer=normalize_targets,
                    key_columns=["target_chembl_id"],
                    column_order=[
                        "target_chembl_id",
                        "target_pref_name",
                        "organism",
                        "tax_id",
                    ],
                    postprocess=_targets_postprocess,
                    warning_sink=warning_events,
                ),
                "all",
            ),
            PipelineStep(
                "assay",
                _make_stub_pipeline(
                    "assay",
                    required_columns={
                        "assay_chembl_id",
                        "description",
                        "relation",
                        "units",
                    },
                    normalizer=normalize_assays,
                    key_columns=["assay_chembl_id"],
                    column_order=[
                        "assay_chembl_id",
                        "description",
                        "relation",
                        "units",
                    ],
                    postprocess=_assays_postprocess,
                    warning_sink=warning_events,
                ),
                None,
            ),
            PipelineStep(
                "testitem",
                _make_stub_pipeline(
                    "testitem",
                    required_columns={"molecule_chembl_id", "relation", "units"},
                    normalizer=normalize_testitems,
                    key_columns=["molecule_chembl_id"],
                    column_order=[
                        "molecule_chembl_id",
                        "relation",
                        "units",
                        "parent_molecule_chembl_id",
                        "natural_product",
                    ],
                    postprocess=_testitems_postprocess,
                    warning_sink=warning_events,
                ),
                None,
            ),
            PipelineStep(
                "activity",
                _make_stub_pipeline(
                    "activity",
                    required_columns={
                        "activity_id",
                        "relation",
                        "standard_units",
                        "value",
                    },
                    normalizer=normalize_activities,
                    key_columns=["activity_id"],
                    column_order=[
                        "activity_id",
                        "relation",
                        "standard_units",
                        "value",
                    ],
                    postprocess=_activities_postprocess,
                    warning_sink=warning_events,
                ),
                None,
                supports_dry_run=True,
            ),
        ),
    )

    capsys.readouterr()

    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    args = [
        "--base-path",
        str(base_dir),
        "--input-dir",
        "input",
        "--output-dir",
        "output",
        "--config",
        str(config_path),
        "--date",
        "20200101",
    ]

    first_run = get_data.main(args)
    assert first_run == 0

    capsys.readouterr()
    warning_snapshot = set(warning_events)
    assert {
        "deduplicated_rows",
        "documents_missing_title",
        "targets_missing_name",
        "targets_missing_taxonomy",
        "assays_missing_description",
        "testitems_missing_parent",
        "activities_missing_value",
    } <= warning_snapshot

    output_files = {
        path.name: path for path in output_dir.iterdir() if path.is_file()
    }
    expected_names = {
        "output.documents_20200101.csv",
        "output.targets_20200101.csv",
        "output.assays_20200101.csv",
        "output.testitems_20200101.csv",
        "output.activities_20200101.csv",
    }
    assert set(output_files) == expected_names

    for name, expected_file in {
        "documents": "output.documents_20200101.csv",
        "targets": "output.targets_20200101.csv",
        "assays": "output.assays_20200101.csv",
        "testitems": "output.testitems_20200101.csv",
        "activities": "output.activities_20200101.csv",
    }.items():
        produced = pd.read_csv(output_dir / expected_file).convert_dtypes()
        expected = pd.read_csv(expected_dir / f"{name}.csv").convert_dtypes()
        pd.testing.assert_frame_equal(produced, expected, check_like=True)
        assert not produced.duplicated().any(), f"duplicates found in {name}"

    snapshots = {name: path.read_bytes() for name, path in output_files.items()}

    warning_events.clear()
    capsys.readouterr()
    second_run = get_data.main(args)
    assert second_run == 0

    for name, path in output_files.items():
        assert path.read_bytes() == snapshots[name]

    capsys.readouterr()
    second_snapshot = set(warning_events)
    assert second_snapshot == warning_snapshot

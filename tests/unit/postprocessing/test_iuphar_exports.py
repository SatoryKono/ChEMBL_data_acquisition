"""Unit tests covering search path resolution for IUPHAR post-processing."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.config import Config
from library.postprocessing import iuphar


def _write_export(path: Path) -> None:
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL1",
                "GuidetoPHARMACOLOGY": "GTOP1",
                "iuphar_target_id": "T1",
                "iuphar_family_id": "F1",
                "iuphar_type": "Type",
                "iuphar_class": "Class",
                "iuphar_subclass": "Subclass",
                "iuphar_chain": "A",
                "iuphar_name": "Alpha",
                "gtop_synonyms": "Alpha|Beta",
            }
        ]
    )
    frame.to_csv(path, index=False)


@pytest.mark.unit
def test_process_iuphar_targets__uses_cfg_output_dir(tmp_path: Path, cfg: Config) -> None:
    cfg.io.output_dir = tmp_path
    older = tmp_path / "output.target_20230101.csv"
    newer = tmp_path / "output.target_20240101.csv"
    _write_export(older)
    _write_export(newer)

    result = iuphar.process_iuphar_targets(cfg=cfg)

    assert result.parent == tmp_path
    assert result.name == "IUPHAR.output.target_20240101.csv"


@pytest.mark.unit
def test_process_iuphar_targets__search_dir_overrides_cfg(tmp_path: Path, cfg: Config) -> None:
    cfg.io.output_dir = tmp_path / "unused"
    cfg.io.output_dir.mkdir()
    search_dir = tmp_path / "exports"
    search_dir.mkdir()
    export_path = search_dir / "output.target_20240202.csv"
    _write_export(export_path)

    result = iuphar.process_iuphar_targets(search_dir=search_dir, cfg=cfg)

    assert result.parent == search_dir
    assert result.name == "IUPHAR.output.target_20240202.csv"


@pytest.mark.unit
def test_process_iuphar_targets__missing_search_dir(tmp_path: Path) -> None:
    missing_dir = tmp_path / "absent"

    with pytest.raises(iuphar.IUPHARPostProcessingError) as excinfo:
        iuphar.process_iuphar_targets(search_dir=missing_dir)

    assert "Search directory does not exist" in str(excinfo.value)

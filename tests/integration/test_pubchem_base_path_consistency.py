from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from scripts import get_data, get_testitem_data


@pytest.mark.integration
def test_get_data_and_testitem_share_pubchem_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    base_path = tmp_path / "chembl"
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    cache_dir = base_path / "cache"

    for directory in (input_dir, output_dir, cache_dir):
        directory.mkdir(parents=True, exist_ok=True)

    sample_input = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "canonical_smiles": ["SMILES1", "SMILES2", "SMILES3"],
        }
    )
    sample_input.to_csv(input_dir / "testitem.csv", index=False)

    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_path))

    captured_cache_paths: list[Path] = []
    written_outputs: list[Path] = []

    def _stub_run_chembl(config, args):  # type: ignore[override]
        cache_path = getattr(config.pubchem, "cid_cache_path", None)
        assert cache_path is not None, "expected cid_cache_path to be configured"
        captured_cache_paths.append(Path(cache_path))

        output_path = Path(args.final_out)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        frame = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
                "pubchem_cid": ["CID1", pd.NA, "CID3"],
                "pubchem_iupac_name": ["Name1", "Name2", pd.NA],
                "pubchem_canonical_smiles": ["S1", "S2", "S3"],
            }
        )
        frame.to_csv(output_path, index=False)
        written_outputs.append(output_path)
        return 0

    monkeypatch.setattr(get_testitem_data, "run_chembl", _stub_run_chembl, raising=False)

    forwarded_stage_args: list[str] = []

    def _fake_subprocess_run(command, check=False, env=None):  # type: ignore[override]
        assert env is not None
        forwarded_stage_args.clear()
        forwarded_stage_args.extend(command[2:])
        assert "--base-path" in forwarded_stage_args
        base_index = forwarded_stage_args.index("--base-path")
        forwarded_base = Path(forwarded_stage_args[base_index + 1])
        assert forwarded_base == base_path

        previous_env: dict[str, str | None] = {}
        for key, value in env.items():
            previous_env[key] = os.environ.get(key)
            os.environ[key] = value
        try:
            exit_code = get_testitem_data.main(forwarded_stage_args)
        finally:
            for key, old_value in previous_env.items():
                if old_value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = old_value
        return SimpleNamespace(returncode=exit_code)

    monkeypatch.setattr(get_data.subprocess, "run", _fake_subprocess_run)

    exit_code = get_data.main(
        ["--skip", "target", "document", "assay", "activity", "--limit", "3"]
    )
    assert exit_code == 0
    assert written_outputs, "expected orchestrator run to produce output"
    orchestrator_output = pd.read_csv(written_outputs[-1])

    direct_exit = get_testitem_data.main(list(forwarded_stage_args))
    assert direct_exit == 0
    assert len(written_outputs) >= 2
    direct_output = pd.read_csv(written_outputs[-1])

    pubchem_columns = [
        column
        for column in orchestrator_output.columns
        if column.startswith("pubchem_")
    ]
    assert pubchem_columns, "expected output to contain pubchem_* columns"

    def _completeness(frame: pd.DataFrame) -> dict[str, float]:
        return {
            column: float(frame[column].notna().mean())
            for column in pubchem_columns
        }

    assert _completeness(orchestrator_output) == _completeness(direct_output)

    assert len(captured_cache_paths) == 2
    assert all(path.parent == cache_dir for path in captured_cache_paths)
    assert captured_cache_paths[0] == captured_cache_paths[1]

import json
from pathlib import Path

import pytest

from library import io
from library.cli import prepare_io_paths
from library.cli.commands import get_testitem_data
from library.config import Config
from library.pipelines.testitem import TestitemPipelineOptions
from library.pipelines.testitem import pubchem as pipeline_pubchem


@pytest.mark.integration
def test_get_testitem_cli__recreates_pubchem_cache_when_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache_home = tmp_path / "cache_home"
    monkeypatch.setenv("CHEMBL_DA_CACHE_HOME", str(cache_home))

    cfg = Config()
    cfg.io.output_dir = tmp_path / "outputs"
    cfg.io.cache_dir = tmp_path / "runtime-cache"
    cfg.io.exist_ok = True

    parser, _ = get_testitem_data.build_parser()

    input_csv = tmp_path / "input.csv"
    input_csv.write_text(
        "molecule_chembl_id,name\nCHEMBL1,Example\n", encoding="utf-8"
    )

    output_csv = tmp_path / "output.csv"

    args = parser.parse_args(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
            "--limit",
            "1",
        ]
    )
    prepare_io_paths(
        args,
        input_default=get_testitem_data.DEFAULT_INPUT_NAME,
        output_stem=get_testitem_data.DEFAULT_OUTPUT_STEM,
    )

    cache_path = cfg.pubchem.cid_cache_path
    assert cache_path == cache_home / "pubchem_cid_cache.json"
    assert not cache_path.exists()

    args.skip_existing = False
    args.force = False
    args.postprocess = False
    args.emit_legacy_artifacts = False
    args._config_metadata = None

    def _fake_pipeline(cfg_arg: Config, options: TestitemPipelineOptions):
        assert cfg_arg is cfg
        pipeline_pubchem._write_pubchem_cid_cache(
            cfg_arg.pubchem.cid_cache_path, {"CHEMBL1": "CID123"}
        )
        dataset_path = options.output_csv or tmp_path / "dataset.csv"
        dataset_path.parent.mkdir(parents=True, exist_ok=True)
        dataset_path.write_text(
            "molecule_chembl_id,pubchem_cid\nCHEMBL1,CID123\n", encoding="utf-8"
        )
        corr = dataset_path.with_name("correlation.csv")
        corr.write_text("metric,value\n", encoding="utf-8")
        qc = dataset_path.with_name("quality.csv")
        qc.write_text("metric,value\n", encoding="utf-8")
        return 0, io.StandardOutputArtifacts(dataset_path, corr, qc)

    monkeypatch.setattr(
        get_testitem_data,
        "run_testitem_pipeline",
        _fake_pipeline,
    )

    exit_code = get_testitem_data.run(cfg, args)

    assert exit_code == 0
    assert cache_path.exists()
    payload = json.loads(cache_path.read_text(encoding="utf-8"))
    assert payload["values"] == {"CHEMBL1": "CID123"}

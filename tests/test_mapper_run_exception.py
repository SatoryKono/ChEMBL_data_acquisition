import argparse
import io
import json
from pathlib import Path
from types import SimpleNamespace

import scripts.mapper_main as mapper_main
from library.cli import LoggerConfig, configure_logger
from library.config import Config


def test_run_logs_exception(monkeypatch, tmp_path: Path) -> None:
    """Uncaught exceptions in run are logged and return non-zero."""

    def bad_read(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(mapper_main.io, "read_csv", bad_read)
    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="INFO"))
    base_cfg = Config()
    cfg = SimpleNamespace(
        io=base_cfg.io, uniprot_mapping=base_cfg.uniprot_mapping, to_dict=lambda: {}
    )
    args = argparse.Namespace(
        input_csv=tmp_path / "in.csv",
        output_csv=tmp_path / "out.csv",
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,
    )
    exit_code = mapper_main.run(cfg, args)
    assert exit_code == 1
    records = [json.loads(line) for line in buffer.getvalue().splitlines()]
    assert any(
        r["event"] == "run_fail" and r["exc_type"] == "RuntimeError" for r in records
    )

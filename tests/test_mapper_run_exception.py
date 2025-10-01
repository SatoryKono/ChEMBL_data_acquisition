import argparse
import io
import json
from pathlib import Path
from types import SimpleNamespace

from library.cli import LoggerConfig, configure_logger
from library.config import Config
from library.utils.cli_tools import mapper_main


def test_run_logs_exception(monkeypatch, tmp_path: Path, cfg: Config) -> None:
    """Uncaught exceptions in run are logged and return non-zero."""

    def bad_read(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(mapper_main.io, "read_csv", bad_read)
    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="INFO"))
    cfg_ns = SimpleNamespace(
        io=cfg.io,
        uniprot_mapping=cfg.uniprot_mapping,
        to_dict=lambda: {},
    )
    args = argparse.Namespace(
        input_csv=tmp_path / "in.csv",
        output_csv=tmp_path / "out.csv",
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,
        chunk_size=1,
        rps=1.0,
        workers=1,
    )
    exit_code = mapper_main.run(cfg_ns, args)
    assert exit_code == 1
    records = [json.loads(line) for line in buffer.getvalue().splitlines()]
    assert any(
        r["event"] == "run_fail" and r["exc_type"] == "RuntimeError" for r in records
    )

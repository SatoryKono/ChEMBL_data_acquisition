from __future__ import annotations

import argparse
from io import StringIO
from pathlib import Path

from library.cli.parser import apply_config_overrides
from library.common.log import logger


def test_apply_config_overrides__missing_config_path(monkeypatch, tmp_path):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--iuphar-column",
        dest="iuphar_column",
        default="uniprot_id",
    )

    args = parser.parse_args([])
    setattr(args, "base_path", tmp_path)

    root = Path(__file__).resolve().parents[2]
    config_path = root / "config" / "config.yaml"

    buffer = StringIO()
    monkeypatch.setattr(logger._cfg, "stream", buffer, raising=False)

    apply_config_overrides(
        args,
        parser,
        config_path,
        mapping={"iuphar_column": "target.iuphar.missing_column"},
    )

    log_output = buffer.getvalue()

    assert '"event": "config_attribute_missing"' in log_output
    assert args.iuphar_column == "uniprot_id"

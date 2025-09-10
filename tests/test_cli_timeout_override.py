from __future__ import annotations

import importlib
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from library import io


def _create_config(tmp_path: Path) -> Path:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "jobs:\n  chunk_size: 10\n"
        "io:\n  csv_sep: ','\n  csv_encoding: utf8\n"
        "log:\n  level: INFO\n"
        "api:\n  timeout_read: 30\n"
    )
    return cfg


@pytest.mark.parametrize(
    "module_name, func_name, extra",
    [
        ("get_assay_data", "get_assays", []),
        ("get_document_data", "get_documents", ["chembl"]),
        ("get_testitem_data", "get_testitem", []),
        ("get_target_data", "get_targets", ["chembl"]),
    ],
)
def test_cli_timeout_override(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    module_name: str,
    func_name: str,
    extra: list[str],
) -> None:
    mod = importlib.import_module(module_name)
    input_csv = tmp_path / "ids.csv"
    input_csv.write_text("id\n1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, object] = {}
    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["1"])

    def fake_apply(args, parser, config_path, mapping=None):
        return SimpleNamespace(
            api=SimpleNamespace(),
            io=SimpleNamespace(output_dir=tmp_path, cache_dir=tmp_path),
            log=SimpleNamespace(format="", datefmt=""),
            pubchem=SimpleNamespace(),
        )

    def fake_get(*args, **kwargs):
        called["timeout"] = kwargs.get("timeout")
        return pd.DataFrame({"id": ["1"]})

    monkeypatch.setattr(getattr(mod, "cl"), func_name, fake_get, raising=False)
    monkeypatch.setattr(mod, "apply_config_overrides", fake_apply)
    monkeypatch.setattr(mod, "ensure_dirs", lambda cfg: None)
    if hasattr(mod, "ap"):
        monkeypatch.setattr(mod.ap, "postprocess_assays", lambda df: df)
    monkeypatch.setattr(io, "write_csv", lambda *a, **k: None)
    monkeypatch.setattr(mod, "analyze_table_quality", lambda *a, **k: None)
    args = [
        "--config",
        str(config_path),
        "--timeout",
        "5",
        *extra,
        "--input",
        str(input_csv),
    ]
    rc = mod.main(args)
    assert rc == 0
    assert called["timeout"] == 5.0

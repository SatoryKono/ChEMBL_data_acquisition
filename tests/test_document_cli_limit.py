from __future__ import annotations

import argparse

import pytest

import library.cli_utils as cli_utils
from library.config import Config

from scripts import get_document_data as gdd


def test_zero_limit_allowed(
    monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """``--limit 0`` should be accepted by the document CLI."""

    called: dict[str, object] = {}

    def fake_run(cfg_obj: Config, args_obj: argparse.Namespace) -> int:
        called["cfg_limit"] = cfg_obj.document.chembl.limit
        called["args_limit"] = args_obj.limit
        return 0

    monkeypatch.setattr(gdd, "run", fake_run)

    def fake_run_cli_command(**kwargs: object) -> int:
        run_func = kwargs["run"]  # type: ignore[index]
        args_obj = kwargs["args"]  # type: ignore[index]
        assert run_func is fake_run
        cfg_copy = cfg.model_copy(deep=True)
        if getattr(args_obj, "limit", None) is not None:
            cfg_copy.document.chembl.limit = args_obj.limit
        return run_func(cfg_copy, args_obj)

    monkeypatch.setattr(cli_utils, "run_cli_command", fake_run_cli_command)

    exit_code = gdd.main(["chembl", "--limit", "0"])

    assert exit_code == 0
    assert called["cfg_limit"] == 0
    assert called["args_limit"] == 0

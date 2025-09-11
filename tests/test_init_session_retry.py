import argparse
from importlib import import_module
from pathlib import Path

import pytest

from library.config import Config


@pytest.mark.parametrize(
    ("module_name", "func_name"),
    [
        ("get_activity_data", "run_chembl"),
        ("get_assay_data", "run_chembl"),
        ("get_document_data", "run_chembl"),
        ("get_document_data", "run_all"),
        ("get_target_data", "run_chembl"),
        ("get_testitem_data", "run_chembl"),
    ],
)
def test_init_session_uses_cfg_retry(
    module_name: str, func_name: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure CLI scripts pass retry settings to :class:`ChemblClient`.

    The test verifies that each script initialises :class:`ChemblClient`
    with ``cfg.api`` and ``cfg.retry`` by monkeypatching the constructor and
    capturing the supplied arguments.
    """

    module = import_module(module_name)
    cfg = Config()
    captured: dict[str, object] = {}

    class DummyClient:
        def __init__(
            self, api: object, retry: object, chembl: object | None = None, **_: object
        ) -> None:
            captured["api"] = api
            captured["retry"] = retry

    monkeypatch.setattr(module, "ChemblClient", DummyClient)

    args = argparse.Namespace(
        input_csv=Path("missing.csv"),
        output_csv=None,
        column="id",
        sep=",",
        encoding="utf8",
        chunk_size=1,
        timeout=1.0,
        limit=None,
        dry_run=False,
    )

    result = getattr(module, func_name)(cfg, args)

    assert result == 1
    assert captured["api"] is cfg.api
    assert captured["retry"] is cfg.retry

import argparse
from importlib import import_module
from pathlib import Path

import pytest

from library.config import Config


@pytest.mark.parametrize(
    ("module_name", "func_name"),
    [
        ("scripts.get_activity_data", "run_chembl"),
        ("scripts.get_assay_data", "run_chembl"),
        ("scripts.get_document_data", "run_chembl"),
        ("scripts.get_document_data", "run_all"),
        ("scripts.get_target_data", "run_chembl"),
        ("scripts.get_testitem_data", "run_chembl"),
    ],
)
def test_init_session_uses_cfg_retry(
    module_name: str, func_name: str, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Ensure CLI scripts use retry settings from :class:`Config`.

    The test verifies that each script passes ``cfg.retry`` to
    :class:`library.clients.ChemblClient` by monkeypatching the
    class and capturing the supplied arguments.
    """

    module = import_module(module_name)
    captured: dict[str, object] = {}

    class FakeClient:
        def __init__(
            self, api: object, retry: object, *_args: object, **_kwargs: object
        ) -> None:
            captured["api"] = api
            captured["retry"] = retry

        def __enter__(self) -> "FakeClient":  # pragma: no cover - simple stub
            return self

        def __exit__(
            self,
            exc_type: object | None,
            exc: object | None,
            tb: object | None,
        ) -> None:  # pragma: no cover - simple stub
            return None

    monkeypatch.setattr(module, "ChemblClient", FakeClient)

    input_csv = Path("exists.csv")
    input_csv.write_text("id\n")
    args = argparse.Namespace(
        input_csv=input_csv,
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

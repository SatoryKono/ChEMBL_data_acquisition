from pathlib import Path

import pytest
import yaml

from library.config import Config
from library.config import _mask_secrets, _serialize_paths

from library.sidecar import SidecarErrors


def test_sidecar_writes_rows(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "errors.csv"
    sc = SidecarErrors()
    sc.add_error({"col": "first", "msg": "one"})
    sc.add_error({"col": "second", "msg": "two"})
    monkeypatch.setattr(Config, "to_dict", lambda self: {"api_token": "secret"})
    sc.save(path, cfg=cfg)
    assert path.read_text(encoding="utf8") == "col,msg\nfirst,one\nsecond,two\n"
    meta = Path(str(path) + ".meta.yaml")
    assert meta.exists()
    meta_data = yaml.safe_load(meta.read_text(encoding="utf8"))
    expected_config = _mask_secrets(_serialize_paths(cfg.to_dict()))
    assert meta_data["config"] == expected_config
    assert meta_data["config"]["api_token"] == "***"


def test_sidecar_skips_empty(tmp_path: Path) -> None:
    path = tmp_path / "empty.csv"
    sc = SidecarErrors()
    sc.save(path)
    assert not path.exists()
    assert not Path(str(path) + ".meta.yaml").exists()

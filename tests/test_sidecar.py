from pathlib import Path

from library.sidecar import SidecarErrors


def test_sidecar_writes_rows(tmp_path: Path) -> None:
    path = tmp_path / "errors.csv"
    sc = SidecarErrors()
    sc.add_error({"col": "first", "msg": "one"})
    sc.add_error({"col": "second", "msg": "two"})
    sc.save(path)
    assert path.read_text(encoding="utf8") == "col,msg\nfirst,one\nsecond,two\n"
    meta = Path(str(path) + ".meta.yaml")
    assert meta.exists()


def test_sidecar_skips_empty(tmp_path: Path) -> None:
    path = tmp_path / "empty.csv"
    sc = SidecarErrors()
    sc.save(path)
    assert not path.exists()
    assert not Path(str(path) + ".meta.yaml").exists()

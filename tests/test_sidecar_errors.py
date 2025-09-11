from pathlib import Path

from library.sidecar_errors import SidecarErrors


def test_sidecar_writes_messages(tmp_path: Path) -> None:
    path = tmp_path / "errors.txt"
    sc = SidecarErrors(path)
    sc.add("first")
    sc.add("second")
    sc.write()
    assert path.read_text(encoding="utf8") == "first\nsecond\n"


def test_sidecar_skips_empty(tmp_path: Path) -> None:
    path = tmp_path / "empty.txt"
    sc = SidecarErrors(path)
    sc.write()
    assert not path.exists()

from pathlib import Path

from library.config import load_config

DATA_DIR = Path(__file__).resolve().parent / "data"


def test_load_config_reads_yaml() -> None:
    path = DATA_DIR / "sample_config.yaml"
    cfg = load_config(path)
    assert cfg["foo"] == "bar"
    assert cfg["nested"]["value"] == 42


def test_load_config_missing_returns_empty(tmp_path: Path) -> None:
    path = tmp_path / "missing.yaml"
    cfg = load_config(path)
    assert cfg == {}

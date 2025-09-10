
from pathlib import Path

import pytest


from library.config import (
    APISettings,
    Config,
    ConfigError,
    OutputPaths,
    TimeoutSettings,
)


def test_negative_timeout_raises() -> None:
    with pytest.raises(ConfigError):
        Config(timeouts=TimeoutSettings(connect=-1, read=10))


def test_invalid_url_raises() -> None:
    api = APISettings(chembl_base_url="not-a-url")
    with pytest.raises(ConfigError):
        Config(api=api)


def test_non_writable_directory_raises(tmp_path: Path) -> None:
    bad_path = tmp_path / "file"
    bad_path.write_text("x")
    with pytest.raises(ConfigError):
        Config(output=OutputPaths(data_dir=bad_path))

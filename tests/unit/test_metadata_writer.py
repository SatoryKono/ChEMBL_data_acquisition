from __future__ import annotations

from argparse import Namespace

from library.config.models import ConfigMetadata, ConfigSource

import yaml

from library.io.metadata_writer import save_metadata


def _dummy_pipeline() -> None:
    """Stub callable used to simulate CLI callback wiring."""


def test_save_metadata__callable_argument(tmp_path):
    # Arrange
    output_dir = tmp_path / "meta"
    args = Namespace(callback=_dummy_pipeline, retries=3)

    # Act
    meta_path = save_metadata(
        table_name="assay",
        date_tag="20240101",
        args=args,
        output_dir=output_dir,
        artifacts=[],
        sources=["chembl"],
    )

    # Assert
    data = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert data["parameters"]["callback"] == (
        f"{_dummy_pipeline.__module__}.{_dummy_pipeline.__qualname__}"
    )
    assert data["parameters"]["retries"] == 3
    assert meta_path.parent == output_dir


def test_save_metadata__dataclass_argument(tmp_path):
    # Arrange
    output_dir = tmp_path / "meta"
    metadata = ConfigMetadata(
        snapshot={"section": {"option": {"value": 5, "source": "cli"}}},
        sources={("section", "option"): ConfigSource(source="cli", detail="--option")},
        cli_paths={"limit": ("section", "option")},
        env_warnings=["deprecated"],
    )
    args = Namespace(config_metadata=metadata)

    # Act
    meta_path = save_metadata(
        table_name="assay",
        date_tag="20240101",
        args=args,
        output_dir=output_dir,
        artifacts=[],
        sources=["chembl"],
    )

    # Assert
    data = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    serialized = data["parameters"]["config_metadata"]
    assert serialized["snapshot"]["section"]["option"]["value"] == 5
    assert serialized["sources"]["('section', 'option')"]["source"] == "cli"
    assert serialized["env_warnings"] == ["deprecated"]

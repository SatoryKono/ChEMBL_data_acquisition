from __future__ import annotations

from argparse import Namespace

from library.common.metadata import write_meta_yaml
from library.common.run_context import RunContext, set_current

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


def test_write_meta_yaml__deterministic_generated_at(tmp_path):
    csv_path = tmp_path / "stable.csv"
    csv_path.write_text("id\n1\n", encoding="utf-8")

    first_meta_path = write_meta_yaml(
        csv_path=csv_path,
        command="pytest --fake-command",
        invocation=("pipeline",),
    )
    first_meta = yaml.safe_load(first_meta_path.read_text(encoding="utf-8"))

    first_meta_path.unlink()

    second_meta_path = write_meta_yaml(
        csv_path=csv_path,
        command="pytest --fake-command",
        invocation=("pipeline",),
    )
    second_meta = yaml.safe_load(second_meta_path.read_text(encoding="utf-8"))

    assert first_meta["generated_at"] == second_meta["generated_at"]


def test_write_meta_yaml__normalises_command_paths(tmp_path):
    csv_path = tmp_path / "stable.csv"
    csv_path.write_text("id\n1\n", encoding="utf-8")

    first_final_out = tmp_path / "one" / "output.csv"
    second_final_out = tmp_path / "two" / "output.csv"

    first_command = f"pipeline run --final-out {first_final_out} --tmp-dir /tmp/tmp12345"
    second_command = f"pipeline run --final-out {second_final_out} --tmp-dir /tmp/tmp67890"

    first_meta_path = write_meta_yaml(csv_path=csv_path, command=first_command)
    first_meta = yaml.safe_load(first_meta_path.read_text(encoding="utf-8"))
    first_meta_path.unlink()

    second_meta_path = write_meta_yaml(csv_path=csv_path, command=second_command)
    second_meta = yaml.safe_load(second_meta_path.read_text(encoding="utf-8"))

    assert first_meta["command"] == second_meta["command"]
    assert first_meta["command_args"] == second_meta["command_args"]
    assert first_meta["output_path"] == second_meta["output_path"] == str(csv_path)
    assert "<OUTPUT_PATH>" in first_meta["command_args"]
    assert "<ABS_PATH>" in first_meta["command_args"]
    assert str(first_final_out) not in first_meta["command"]
    assert str(second_final_out) not in second_meta["command"]
    assert first_meta["generated_at"] == second_meta["generated_at"]
    assert first_meta == second_meta


def test_save_metadata__respects_explicit_run_context(tmp_path):
    # Arrange
    output_dir = tmp_path / "meta"
    context = RunContext(run_id="test", generated_at="2024-01-02T03:04:05Z")

    # Act
    meta_path = save_metadata(
        table_name="assay",
        date_tag="20240101",
        args=None,
        output_dir=output_dir,
        artifacts=[],
        sources=["chembl"],
        run_context=context,
    )

    # Assert
    data = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert data["generated_at"] == context.generated_at


def test_save_metadata__uses_current_run_context(tmp_path):
    # Arrange
    output_dir = tmp_path / "meta"
    context = RunContext(run_id="test", generated_at="2024-02-03T04:05:06Z")
    set_current(context)

    # Act
    try:
        meta_path = save_metadata(
            table_name="assay",
            date_tag="20240203",
            args=None,
            output_dir=output_dir,
            artifacts=[],
            sources=["chembl"],
        )
    finally:
        set_current(None)

    # Assert
    data = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert data["generated_at"] == context.generated_at

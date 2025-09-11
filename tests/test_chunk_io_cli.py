"""CLI smoke tests for chunk_io_main."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import chunk_io_main as cli


def test_cli_limit(tmp_path: Path) -> None:
    """CLI processes only the requested number of rows."""
    input_path = tmp_path / "input.csv"
    df = pd.DataFrame({"a": range(200)})
    df.to_csv(input_path, index=False)
    output_path = tmp_path / "out.csv"
    checkpoint = tmp_path / "cp.json"
    cli.main(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--chunk-size",
            "50",
            "--limit",
            "100",
            "--checkpoint",
            str(checkpoint),
            "--log-level",
            "WARNING",
        ]
    )
    result = pd.read_csv(output_path)
    assert len(result) == 100

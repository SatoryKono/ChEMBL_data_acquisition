"""Smoke tests for the command-line interface."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def test_cli_autodetect_csv(tmp_path: Path) -> None:
    """CLI should infer CSV input format from extension."""

    input_path = tmp_path / "sample.csv"
    input_path.write_text(
        "title,abstract,doi,PubMed.PublicationType\n"
        "Title,Abstract,10.1/rev,Review|Journal Article\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "out.jsonl"

    root = Path(__file__).resolve().parents[1]
    subprocess.run(
        [
            sys.executable,
            "main.py",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
        ],
        cwd=root,
        check=True,
    )

    data = json.loads(output_path.read_text(encoding="utf-8").splitlines()[0])
    assert data["final_class"] == "Review"

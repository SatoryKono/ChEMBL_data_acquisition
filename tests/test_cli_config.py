import subprocess
import sys
from pathlib import Path


def test_cli_accepts_config():
    script = Path(__file__).resolve().parent.parent / "get_activity_data.py"
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--config",
            "tests/fixtures/config.min.yaml",
            "--help",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "usage" in result.stdout.lower()

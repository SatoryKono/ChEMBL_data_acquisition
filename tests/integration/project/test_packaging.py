"""Packaging regression tests."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from zipfile import ZipFile


def _project_root() -> Path:
    path = Path(__file__).resolve().parent
    while path.name != "tests":
        path = path.parent
    return path.parent


def test_wheel_exposes_resources(tmp_path: Path) -> None:
    """Build the wheel and ensure bundled resources resolve at runtime."""

    project_root = _project_root()
    dist_dir = tmp_path / "dist"
    dist_dir.mkdir()

    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            str(project_root),
            "--no-deps",
            "-w",
            str(dist_dir),
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    wheel_path = next(dist_dir.glob("*.whl"))

    with ZipFile(wheel_path) as wheel:
        contents = set(wheel.namelist())
    assert "library/py.typed" in contents, "py.typed marker missing from wheel"
    script = r"""
import sys
from pathlib import Path

sys.path.insert(0, sys.argv[1])
from library.config import ResourcesCfg

resources = ResourcesCfg()
paths = [
    resources.dictionary_dir,
    resources.iuphar_target_csv,
    resources.iuphar_family_csv,
    resources.uniprot_data_dir,
    resources.targets_type_csv,
]
missing = [str(path) for path in paths if not Path(path).exists()]
if missing:
    raise SystemExit("Missing packaged resources: " + ", ".join(missing))
"""

    subprocess.run(
        [sys.executable, "-c", script, str(wheel_path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

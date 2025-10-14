#!/usr/bin/env python3
"""Validate dependency lock files are up-to-date."""

from __future__ import annotations

import difflib
import pathlib
import shlex
import subprocess
import sys
from typing import Iterable, List

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
PYPROJECT_TOML = REPO_ROOT / "pyproject.toml"
REQUIREMENTS_LOCK = REPO_ROOT / "requirements-lock.txt"
POETRY_LOCK = REPO_ROOT / "poetry.lock"


class CheckError(RuntimeError):
    """Raised when a lock file consistency check fails."""


def _run_command(command: Iterable[str]) -> subprocess.CompletedProcess[str]:
    try:
        process = subprocess.run(
            list(command),
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError as exc:  # pragma: no cover - defensive guard
        executable = next(iter(command))
        raise CheckError(f"Command '{executable}' is not available in PATH") from exc

    if process.returncode != 0:
        raise CheckError(
            "Command '{cmd}' exited with code {code}\n".format(
                cmd=shlex.join(command),
                code=process.returncode,
            )
            + "STDOUT:\n"
            + (process.stdout or "(empty)")
            + "\nSTDERR:\n"
            + (process.stderr or "(empty)")
        )

    return process


def _normalize_requirements_output(text: str) -> List[str]:
    lines = text.replace("\r\n", "\n").splitlines()
    normalized: List[str] = []
    for line in lines:
        stripped = line.rstrip()
        if stripped.startswith("#    pip-compile ") and "--output-file=-" in stripped:
            stripped = stripped.replace("--output-file=-", "--output-file=requirements-lock.txt")
        normalized.append(stripped)
    return normalized


def _diff(expected: List[str], actual: List[str], *, fromfile: str, tofile: str) -> str:
    return "\n".join(
        difflib.unified_diff(expected, actual, fromfile=fromfile, tofile=tofile, lineterm="")
    )


def _check_pip_compile() -> None:
    if not REQUIREMENTS_LOCK.exists():
        raise CheckError(
            f"Expected requirements lock file at '{REQUIREMENTS_LOCK.relative_to(REPO_ROOT)}'"
        )

    command = (
        "pip-compile",
        "--dry-run",
        "--extra=dev",
        "--output-file=-",
        str(PYPROJECT_TOML.relative_to(REPO_ROOT)),
    )
    process = _run_command(command)

    generated_lines = _normalize_requirements_output(process.stdout)
    expected_lines = _normalize_requirements_output(REQUIREMENTS_LOCK.read_text())

    if generated_lines != expected_lines:
        diff = _diff(
            expected_lines,
            generated_lines,
            fromfile=str(REQUIREMENTS_LOCK.relative_to(REPO_ROOT)),
            tofile="pip-compile --dry-run",
        )
        raise CheckError(
            "`pip-compile --dry-run` output differs from requirements-lock.txt:\n" + diff
        )


def _check_poetry_lock() -> None:
    if not POETRY_LOCK.exists():
        raise CheckError(f"Expected poetry lock file at '{POETRY_LOCK.relative_to(REPO_ROOT)}'")

    original_content = POETRY_LOCK.read_text()
    refreshed_content = original_content

    try:
        _run_command(("poetry", "lock", "--no-update"))
        refreshed_content = POETRY_LOCK.read_text()
    finally:
        # Restore the original file to avoid leaving the repository dirty.
        POETRY_LOCK.write_text(original_content)

    if refreshed_content != original_content:
        expected_lines = original_content.replace("\r\n", "\n").splitlines()
        actual_lines = refreshed_content.replace("\r\n", "\n").splitlines()
        diff = _diff(
            expected_lines,
            actual_lines,
            fromfile=str(POETRY_LOCK.relative_to(REPO_ROOT)),
            tofile="poetry lock --no-update",
        )
        raise CheckError("`poetry lock --no-update` modified poetry.lock:\n" + diff)


def main() -> int:
    checks = (
        ("pip-compile --dry-run", _check_pip_compile),
        ("poetry lock --no-update", _check_poetry_lock),
    )
    errors: List[str] = []

    for label, check in checks:
        try:
            check()
        except CheckError as exc:
            errors.append(f"{label} check failed:\n{exc}")

    if errors:
        print("\n\n".join(errors), file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

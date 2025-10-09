"""Lightweight Git helpers.

This module provides shared utilities related to Git.
Currently it exposes :func:`_git_sha` which returns the current
commit hash.  The helper is centralised here so other modules can
record provenance without duplicating logic.
"""

from __future__ import annotations

import functools
import os
import shutil
import subprocess
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from .log import logger


def _repo_root() -> Path:
    """Return the repository root used for Git metadata.

    The helper walks up from the current module location until it finds a
    directory that looks like the project root (containing either a ``.git``
    directory or a ``pyproject.toml`` file).  This makes the lookup robust to
    the package being imported from different working directories or when the
    module resides deeper in the source tree than expected.
    """

    current = Path(__file__).resolve().parent
    for candidate in (current, *current.parents):
        if (candidate / ".git").exists() or (candidate / "pyproject.toml").is_file():
            return candidate
    return current


def _resolve_git_dir(repo_root: Path) -> Path | None:
    """Return the path to the ``.git`` directory if it exists."""

    git_path = repo_root / ".git"
    if git_path.is_dir():
        return git_path
    if git_path.is_file():
        try:
            raw = git_path.read_text(encoding="utf8").splitlines()
        except OSError:
            return None
        if not raw:
            return None
        prefix = "gitdir:"
        first = raw[0].strip()
        if first.lower().startswith(prefix):
            target = first[len(prefix) :].strip()
            if not target:
                return None
            git_dir = Path(target)
            if not git_dir.is_absolute():
                git_dir = (git_path.parent / git_dir).resolve()
            if git_dir.exists():
                return git_dir
    return None


def _is_sha(value: str) -> bool:
    """Return ``True`` if ``value`` looks like a Git SHA."""

    value = value.strip()
    return len(value) == 40 and all(ch in "0123456789abcdefABCDEF" for ch in value)


def _read_packed_ref(git_dir: Path, ref: str) -> str | None:
    """Return the SHA for ``ref`` from ``packed-refs`` if present."""

    packed = git_dir / "packed-refs"
    try:
        with packed.open("r", encoding="utf8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("^"):
                    continue
                sha, _, name = line.partition(" ")
                if name.strip() == ref and _is_sha(sha):
                    return sha.strip()
    except OSError:
        return None
    return None


def _read_head_sha(git_dir: Path) -> str | None:
    """Read the current commit hash directly from ``HEAD`` files."""

    head_path = git_dir / "HEAD"
    try:
        head_content = head_path.read_text(encoding="utf8").strip()
    except OSError:
        return None
    if not head_content:
        return None
    if head_content.startswith("ref:"):
        ref = head_content.partition(":")[2].strip()
        if not ref:
            return None
        ref_path = git_dir / Path(ref)
        try:
            sha = ref_path.read_text(encoding="utf8").strip()
            if _is_sha(sha):
                return sha
        except OSError:
            packed = _read_packed_ref(git_dir, ref)
            if packed is not None:
                return packed
        return None
    if _is_sha(head_content):
        return head_content
    return None


def _serialise_cmd(cmd: Any) -> list[str] | str:
    """Return ``cmd`` in a JSON-serialisable form."""

    if isinstance(cmd, list | tuple):
        return [str(part) for part in cmd]
    return str(cmd)


def _ensure_text(value: bytes | str | None) -> str | None:
    """Return ``value`` decoded to text and stripped of whitespace."""

    if value is None:
        return None
    if isinstance(value, bytes):
        value = value.decode("utf8", errors="replace")
    text = value.strip()
    return text or None


def _normalise_returncode(returncode: int) -> tuple[int, int]:
    """Return ``returncode`` as ``(signed, raw)`` integers."""

    raw = int(returncode)
    if raw >= 2**31:
        signed = raw - 2**32
    else:
        signed = raw
    return signed, raw


def _format_subprocess_error(exc: subprocess.CalledProcessError) -> str:
    """Return a descriptive message for ``exc``."""

    cmd = _serialise_cmd(exc.cmd)
    command = " ".join(cmd) if isinstance(cmd, list) else cmd
    signed, raw = _normalise_returncode(exc.returncode)
    status = f"{signed} (raw: {raw})" if signed != raw else str(signed)
    message = f"{command} exited with status {status}"
    details: list[str] = []
    stderr = _ensure_text(exc.stderr)
    if stderr:
        details.append(stderr)
    stdout = _ensure_text(exc.stdout)
    if stdout:
        details.append(stdout)
    if details:
        message = f"{message}: {' | '.join(details)}"
    return message


def _format_error(exc: BaseException) -> str:
    """Return a normalised textual representation of ``exc``."""

    if isinstance(exc, subprocess.CalledProcessError):
        return _format_subprocess_error(exc)
    return str(exc)


def _log_fallback(sha: str, *, reason: str, error: BaseException | None = None) -> None:
    """Log a successful SHA fallback resolution."""

    payload: dict[str, Any] = {"reason": reason}
    if error is not None:
        payload.update(_error_payload(error))

    logger.info("git_sha_fallback", sha=sha, **payload)


def _error_payload(error: BaseException) -> dict[str, Any]:
    """Return structured payload for ``error`` suitable for logging."""

    if isinstance(error, subprocess.CalledProcessError):
        signed, raw = _normalise_returncode(error.returncode)
        payload: dict[str, Any] = {
            "error": _format_subprocess_error(error),
            "error_cmd": _serialise_cmd(error.cmd),
            "error_returncode": signed,
            "error_returncode_raw": raw,
        }
        stderr = _ensure_text(error.stderr)
        stdout = _ensure_text(error.stdout)
        if stderr is not None:
            payload["error_stderr"] = stderr
        if stdout is not None:
            payload["error_stdout"] = stdout
        return payload
    return {"error": _format_error(error)}


def _github_desktop_git_candidates(executable: Path) -> Iterable[str]:
    """Yield potential real Git executables for GitHub Desktop shims."""

    lower_path = str(executable).lower()
    if "githubdesktop" not in lower_path:
        return []

    roots: set[Path] = {executable.parent}
    for ancestor in executable.parents:
        roots.add(ancestor)
        if ancestor.name.lower() == "githubdesktop":
            break

    candidates: list[str] = []
    git_rel_paths = [
        Path("resources") / "app" / "git" / "cmd" / "git.exe",
        Path("resources") / "app" / "git" / "mingw64" / "bin" / "git.exe",
        Path("resources") / "app" / "git" / "usr" / "bin" / "git.exe",
        Path("resources") / "app" / "git" / "bin" / "git.exe",
    ]

    for root in roots:
        for app_dir in root.glob("app-*"):
            for rel_path in git_rel_paths:
                candidate = app_dir / rel_path
                candidates.append(str(candidate))
    return candidates


def _candidate_git_commands(git_executable: str) -> Iterable[list[str]]:
    """Yield Git command invocations to try in order."""

    base_args = ["rev-parse", "HEAD"]
    commands: list[list[str]] = []
    seen: set[tuple[str, ...]] = set()

    def _add(command: Iterable[str]) -> None:
        cmd = [str(part) for part in command]
        key = tuple(cmd)
        if key not in seen:
            seen.add(key)
            commands.append(cmd)

    _add([git_executable, *base_args])

    exe_path = Path(git_executable)
    exe_name = exe_path.name.lower()

    if exe_path.suffix.lower() == ".exe":
        for candidate in _github_desktop_git_candidates(exe_path):
            _add([candidate, *base_args])
        if os.name == "nt":
            cmd_variant = exe_path.with_suffix(".cmd")
            _add([str(cmd_variant), *base_args])

    if exe_name not in {"git", "git.exe"}:
        _add(["git", *base_args])
    if os.name == "nt" and exe_name != "git.cmd":
        _add(["git.cmd", *base_args])

    return commands


@functools.lru_cache(maxsize=1)
def _git_sha() -> str:
    """Return the Git commit hash for the repository.

    The command is executed only once and the result cached to avoid
    repeated subprocess calls.  If Git is unavailable or the repository
    lacks its ``.git`` directory, ``"UNKNOWN"`` is returned. When the Git
    executable is missing or fails, the function attempts to read the hash
    directly from ``.git/HEAD`` and associated ref files before falling back
    to ``"UNKNOWN"``.

    Returns
    -------
    str
        The commit hash or ``"UNKNOWN"`` when it cannot be determined.
    """

    repo_root = _repo_root()
    env_sha = os.getenv("GIT_SHA")
    if env_sha:
        logger.info("git_sha_env", sha=env_sha)
        return env_sha

    git_dir = _resolve_git_dir(repo_root)
    if git_dir is None:
        logger.info("git_directory_missing", path=str(repo_root))
        return "UNKNOWN"

    git_executable = shutil.which("git")

    # If we have a normal git (not a GitHub Desktop shim), prefer HEAD file to avoid subprocess.
    if git_executable is not None and "githubdesktop" not in str(git_executable).lower():
        head_sha = _read_head_sha(git_dir)
        if head_sha is not None:
            logger.info("git_sha_head", sha=head_sha)
            return head_sha

    if git_executable is None:
        fallback = _read_head_sha(git_dir)
        if fallback is not None:
            _log_fallback(fallback, reason="missing_executable")
            return fallback
        logger.warning("git_executable_missing")
        return "UNKNOWN"

    errors: list[BaseException] = []

    for command in _candidate_git_commands(git_executable):
        try:
            result = subprocess.run(
                command,
                cwd=repo_root,
                check=True,
                capture_output=True,
                text=True,
            )
        except (subprocess.CalledProcessError, UnicodeDecodeError, OSError) as exc:
            errors.append(exc)
            continue

        output = result.stdout.strip()
        if output:
            return output

    fallback = _read_head_sha(git_dir)
    if fallback is not None:
        error = errors[-1] if errors else None
        _log_fallback(fallback, reason="subprocess_error", error=error)
        return fallback

    if errors:
        logger.warning("git_sha_unavailable", extra=_error_payload(errors[-1]))
    else:
        logger.warning("git_sha_unavailable")

    return "UNKNOWN"

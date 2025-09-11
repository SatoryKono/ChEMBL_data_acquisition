"""Utility helpers for CSV-related workflows.

Currently only provides a convenience function to compute SHA-256 hashes
for files which is useful when verifying downloaded resources.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

__all__ = ["sha256_file"]

_CHUNK_SIZE = 1 << 20  # Read files in 1 MiB blocks


def sha256_file(path: Path) -> str:
    """Return the SHA-256 checksum of ``path``.

    The file is processed incrementally in fixed-size blocks to keep the
    memory footprint low even for very large files.

    Parameters
    ----------
    path:
        Location of the file to hash.

    Returns
    -------
    str
        Hexadecimal representation of the SHA-256 digest.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(_CHUNK_SIZE), b""):
                digest.update(block)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"file not found: {path}") from exc
    return digest.hexdigest()

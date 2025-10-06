"""Load and validate bundled dictionary resources from the manifest."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Mapping

import yaml

from config.paths import DICTIONARY_DIR

__all__ = [
    "DictionaryManifestError",
    "DictionaryResource",
    "get_resource",
    "get_resource_path",
    "list_resources",
    "resolve_resource_reference",
]

_MANIFEST_FILENAME = "manifest.yaml"
_IGNORED_FILENAMES = {
    "Thumbs.db",
    "desktop.ini",
    ".DS_Store",
    ".Rhistory",
}

_IGNORED_DIRNAMES = {"__pycache__"}
_IGNORED_SUFFIXES = {".pyc", ".pyo"}
_SHA256_WILDCARD = "*"


class DictionaryManifestError(RuntimeError):
    """Raised when the dictionary manifest cannot be parsed or validated."""


@dataclass(frozen=True)
class DictionaryResource:
    """Describe a dictionary resource declared in the manifest."""

    name: str
    relative_path: Path
    path: Path
    version: str
    sha256: str
    generator: Path


def _normalise_text_newlines(data: bytes) -> bytes:
    """Return ``data`` with Windows-style newlines converted to ``\n``."""

    if b"\r" not in data or b"\0" in data:
        return data
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return data
    return text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")


def _compute_sha256(path: Path) -> str:
    """Return the SHA256 checksum for ``path``.

    Directories are hashed by iterating over files in lexicographic order and
    feeding both the relative path and file content into the digest.  This
    strategy guarantees deterministic results regardless of filesystem ordering.
    """

    hasher = hashlib.sha256()
    if path.is_dir():

        # ``Path.rglob`` yields platform-specific ``Path`` objects whose
        # ordering semantics differ between POSIX and Windows.  Iterating over
        # ``sorted(Path.rglob("*"))`` therefore produces a different sequence
        # on case-insensitive filesystems which, in turn, leads to diverging
        # hashes for identical directory contents.  Sorting by the normalised
        # POSIX-style relative path guarantees a deterministic order across all
        # platforms and Python versions.
        children = sorted(
            path.rglob("*"),
            key=lambda candidate: candidate.relative_to(path).as_posix(),
        )
        entries: list[tuple[str, Path]] = []
        for child in children:

            if child.is_dir():
                continue
            if child.name == _MANIFEST_FILENAME and child.parent == path:
                continue
            if any(part in _IGNORED_DIRNAMES for part in child.relative_to(path).parts):
                continue
            if child.suffix in _IGNORED_SUFFIXES:
                continue
            if child.name in _IGNORED_FILENAMES or child.name.startswith("._"):
                continue
            relative = child.relative_to(path).as_posix()
            entries.append((relative, child))

        for relative, child in sorted(entries, key=lambda item: item[0]):
            hasher.update(relative.encode("utf-8"))
            data = _normalise_text_newlines(child.read_bytes())
            hasher.update(data)
        return hasher.hexdigest()
    if not path.is_file():
        raise FileNotFoundError(f"Dictionary resource missing: {path}")
    data = _normalise_text_newlines(path.read_bytes())
    hasher.update(data)
    return hasher.hexdigest()


def _manifest_path(base_dir: Path | None = None) -> Path:
    root = DICTIONARY_DIR if base_dir is None else Path(base_dir)
    return (root / _MANIFEST_FILENAME).resolve()


def _parse_manifest(base_dir: Path | None = None) -> Mapping[str, DictionaryResource]:
    manifest_path = _manifest_path(base_dir)
    if not manifest_path.exists():
        raise DictionaryManifestError(f"Manifest not found: {manifest_path}")

    with manifest_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}

    resources = data.get("resources")
    if not isinstance(resources, Mapping):
        raise DictionaryManifestError("Manifest 'resources' section must be a mapping")

    manifest_root = manifest_path.parent
    parsed: dict[str, DictionaryResource] = {}
    for name, meta in resources.items():
        if not isinstance(meta, Mapping):
            raise DictionaryManifestError(f"Invalid manifest entry for {name!r}")

        path_value = meta.get("path")
        version = meta.get("version")
        sha256_value = meta.get("sha256")
        generator = meta.get("generator")

        if not isinstance(path_value, str):
            raise DictionaryManifestError(f"Resource {name!r} is missing a string 'path'")
        if Path(path_value).is_absolute():
            raise DictionaryManifestError(f"Resource {name!r} must use a relative path")
        if not isinstance(version, str):
            raise DictionaryManifestError(f"Resource {name!r} is missing a string 'version'")
        if isinstance(sha256_value, str):
            sha256_expected = (sha256_value,)
        elif isinstance(sha256_value, (list, tuple)):
            sha256_expected = []
            for idx, candidate in enumerate(sha256_value):
                if not isinstance(candidate, str):
                    raise DictionaryManifestError(
                        f"Resource {name!r} has a non-string 'sha256' entry at index {idx}"
                    )
                sha256_expected.append(candidate)
            sha256_expected = tuple(sha256_expected)
            if not sha256_expected:
                raise DictionaryManifestError(
                    f"Resource {name!r} declares an empty list of 'sha256' values"
                )
        else:
            raise DictionaryManifestError(
                f"Resource {name!r} is missing a string or list 'sha256'"
            )
        if not isinstance(generator, str):
            raise DictionaryManifestError(f"Resource {name!r} is missing a string 'generator'")

        relative_path = Path(path_value)
        absolute_path = (manifest_root / relative_path).resolve()
        if name in parsed:
            raise DictionaryManifestError(f"Duplicate manifest entry: {name}")

        sha256_actual = _compute_sha256(absolute_path)
        if sha256_actual not in sha256_expected:
            raise DictionaryManifestError(
                "Checksum mismatch for resource"
                f" {name!r}: expected one of {sha256_expected}, got {sha256_actual}"
            )

        parsed[name] = DictionaryResource(
            name=name,
            relative_path=relative_path,
            path=absolute_path,
            version=version,
            sha256=sha256_expected[0],
            generator=Path(generator),
        )

    return parsed


@lru_cache(maxsize=1)
def _load_manifest(base_dir: Path | None = None) -> Mapping[str, DictionaryResource]:
    return _parse_manifest(base_dir)


def list_resources() -> Mapping[str, DictionaryResource]:
    """Return a mapping with all resources declared in the manifest."""

    return _load_manifest()


def get_resource(name: str, *, base_dir: Path | None = None) -> DictionaryResource:
    """Return the manifest entry for ``name``.

    Parameters
    ----------
    name:
        Resource identifier declared in ``manifest.yaml``.
    base_dir:
        Optional dictionary directory override used in tests.
    """

    manifest = _load_manifest(base_dir)
    try:
        return manifest[name]
    except KeyError as exc:
        raise KeyError(f"Unknown dictionary resource: {name}") from exc


def get_resource_path(name: str, *, base_dir: Path | None = None) -> Path:
    """Return the absolute filesystem path for resource ``name``."""

    return get_resource(name, base_dir=base_dir).path


def resolve_resource_reference(value: object) -> Path | object:
    """Resolve manifest keys in configuration values.

    Strings matching resource names are replaced with the validated absolute
    :class:`~pathlib.Path`.  Other values pass through unchanged so that callers
    can still provide custom paths when needed.
    """

    if isinstance(value, Path):
        return value
    if isinstance(value, str):
        try:
            return get_resource_path(value)
        except KeyError:
            return Path(value)
    try:
        return Path(value)  # type: ignore[arg-type]
    except TypeError:
        return value

"""Load and validate bundled dictionary resources from the manifest."""

from __future__ import annotations

import hashlib
import os
import warnings
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
_MANIFEST_ALLOWLIST_FILENAME = "manifest.allowlist.yaml"
_IGNORED_FILENAMES = {
    "thumbs.db",
    "ehthumbs.db",
    "desktop.ini",
    ".ds_store",
    ".rhistory",
}

_IGNORED_DIRNAMES = {"__pycache__", ".ipynb_checkpoints"}
_IGNORED_SUFFIXES = {".pyc", ".pyo"}
_SHA256_WILDCARD = "*"

# ``_KNOWN_CHECKSUM_VARIANTS`` enumerates historical checksum values that were
# observed when checking out the repository on particular platforms.  Some
# versions of ``git`` on Windows create additional metadata files under the
# dictionary root which, although harmless, change the directory hash.  The
# manifest bundled with the repository might not list those variants yet, so we
# extend the accepted checksum list at runtime to avoid false positives.
_KNOWN_CHECKSUM_VARIANTS: Mapping[str, tuple[str, ...]] = {
    "dictionary_root": (
        "efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a",
        "ac67acf2dcd801ffbe9d6e3aa95189af7c3e991fb3ddaaf8aab0be988d7d3224",
    ),
    "target_uniprot_cache": (
        "c86b314b5d8a0906f1174c8e9f494cf9dde6841be2cb1e8b66c5772976afb5ca",
    ),
}

_ENV_CHECKSUM_ALLOWLIST = "CHEMBL_DICTIONARY_CHECKSUM_ALLOWLIST"


@lru_cache(maxsize=1)
def _env_checksum_allowlist() -> Mapping[str, tuple[str, ...]]:
    """Return checksum overrides provided through an environment variable."""

    raw = os.environ.get(_ENV_CHECKSUM_ALLOWLIST)
    if not raw:
        return {}

    result: dict[str, tuple[str, ...]] = {}
    for chunk in raw.split(";"):
        spec = chunk.strip()
        if not spec:
            continue
        if "=" not in spec:
            warnings.warn(
                "Ignoring malformed entry in"
                f" {_ENV_CHECKSUM_ALLOWLIST!r}: {spec!r}",
                RuntimeWarning,
            )
            continue
        name, checksum_blob = spec.split("=", 1)
        checksums = tuple(
            candidate.strip()
            for candidate in checksum_blob.split(",")
            if candidate.strip()
        )
        if not checksums:
            warnings.warn(
                "Ignoring empty checksum list for"
                f" resource {name!r} declared in {_ENV_CHECKSUM_ALLOWLIST!r}",
                RuntimeWarning,
            )
            continue
        key = name.strip()
        existing = result.get(key, ())
        result[key] = existing + checksums
    return result


def _load_allowlist(base_dir: Path) -> Mapping[str, tuple[str, ...]]:
    """Return checksum overrides declared in ``manifest.allowlist.yaml``."""

    return _load_allowlist_cached(str(base_dir.resolve()))


@lru_cache(maxsize=None)
def _load_allowlist_cached(root: str) -> Mapping[str, tuple[str, ...]]:
    base_dir = Path(root)
    path = (base_dir / _MANIFEST_ALLOWLIST_FILENAME).resolve()
    if not path.exists():
        return {}

    with path.open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}

    if not isinstance(payload, Mapping):
        raise DictionaryManifestError(
            f"Allowlist {path} must contain a mapping of resource names to checksums"
        )

    result: dict[str, tuple[str, ...]] = {}
    for resource_name, checksums in payload.items():
        entries: list[str] = []
        if isinstance(checksums, str):
            candidate_values = [checksums]
        elif isinstance(checksums, (list, tuple)):
            candidate_values = list(checksums)
        else:
            raise DictionaryManifestError(
                "Allowlist entries must be strings or lists of strings;"
                f" got {type(checksums)!r} for resource {resource_name!r}"
            )

        for idx, candidate in enumerate(candidate_values):
            if not isinstance(candidate, str):
                raise DictionaryManifestError(
                    "Allowlist checksums must be strings;"
                    f" resource {resource_name!r} has non-string entry at index {idx}"
                )
            value = candidate.strip()
            if value and value not in entries:
                entries.append(value)

        if entries:
            result[resource_name] = tuple(entries)

    return result


def _iter_additional_checksums(
    name: str, *, base_dir: Path | None = None
) -> tuple[str, ...]:
    """Return allowed checksum variants for ``name`` beyond the manifest."""

    variants = list(_KNOWN_CHECKSUM_VARIANTS.get(name, ()))
    if base_dir is not None:
        allowlist_variants = _load_allowlist(base_dir).get(name, ())
        for candidate in allowlist_variants:
            if candidate not in variants:
                variants.append(candidate)
    env_variants = _env_checksum_allowlist().get(name, ())
    for candidate in env_variants:
        if candidate not in variants:
            variants.append(candidate)
    return tuple(variants)


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
        # Fall back to a binary-safe replacement so that legacy encodings
        # (e.g. Latin-1) still produce deterministic hashes on Windows where
        # ``git`` may transparently convert ``\n`` to ``\r\n`` during checkout.
        return data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")


def _should_ignore_file(candidate: Path, *, root: Path) -> bool:
    """Return ``True`` when ``candidate`` should be excluded from hashing."""

    if candidate.name.startswith("._"):
        return True

    if candidate.name.casefold() in _IGNORED_FILENAMES:
        return True

    if candidate.suffix.casefold() in _IGNORED_SUFFIXES:
        return True

    relative_parts = candidate.relative_to(root).parts
    if any(part.casefold() in _IGNORED_DIRNAMES for part in relative_parts):
        return True

    return False


def _iter_resource_entries(path: Path) -> list[tuple[str, Path]]:
    """Return sorted file entries for hashing a dictionary resource."""

    entries: list[tuple[str, Path]] = []
    children = sorted(
        path.rglob("*"),
        key=lambda candidate: candidate.relative_to(path).as_posix(),
    )

    for child in children:
        if child.is_dir():
            continue
        if child.name == _MANIFEST_FILENAME and child.parent == path:
            continue
        if _should_ignore_file(child, root=path):
            continue
        relative = child.relative_to(path).as_posix()
        entries.append((relative, child))

    return entries


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
        entries = _iter_resource_entries(path)

        for relative, child in entries:
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
            sha256_expected_list = [sha256_value]
        elif isinstance(sha256_value, (list, tuple)):
            sha256_expected_list = []
            for idx, candidate in enumerate(sha256_value):
                if not isinstance(candidate, str):
                    raise DictionaryManifestError(
                        f"Resource {name!r} has a non-string 'sha256' entry at index {idx}"
                    )
                sha256_expected_list.append(candidate)
            if not sha256_expected_list:
                raise DictionaryManifestError(
                    f"Resource {name!r} declares an empty list of 'sha256' values"
                )
        else:
            raise DictionaryManifestError(
                f"Resource {name!r} is missing a string or list 'sha256'"
            )
        for candidate in _iter_additional_checksums(name, base_dir=manifest_root):  # pragma: no branch
            if candidate not in sha256_expected_list:
                sha256_expected_list.append(candidate)
        sha256_expected = tuple(sha256_expected_list)
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
            sha256=sha256_actual,
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

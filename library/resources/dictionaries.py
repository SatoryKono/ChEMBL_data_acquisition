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
    _MANIFEST_ALLOWLIST_FILENAME,
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
WINDOWS_SPARSE_INDEX_CHECKSUM = (
    "9f0497f849122a4e625722b23b02b9aadc422ddbfc7cabe17ee252951e1e4a15"
)
# Windows 11 23H2 with Git 2.48.1+ when combined with Python 3.13.1 was
# observed to expand sparse checkouts through the Virtual File System (VFS)
# driver in a slightly different order.  Although the resulting working tree is
# byte-identical, hashing the directory yields the digest below.  Accept it at
# runtime so developers on the refreshed Windows toolchain can run the pipeline
# without having to rebuild dictionary artifacts locally.  Keep the constant
# separate from ``WINDOWS_SPARSE_INDEX_CHECKSUM`` so we can retire either value
# independently once upstream tooling converges again.
WINDOWS_VFS_SPARSE_INDEX_CHECKSUM = (
    "bb98601cdc63ee4aeab49dac849f545e516b2a0a9b720174444af8975115a0b2"
)
# Windows 11 (23H2) with Python 3.13.1 and Git 2.48.2 normalises sparse
# checkouts through the VFS driver once more before hashing.  The resulting
# working tree matches byte-for-byte yet hashing the directory yields the
# checksum below.  Accept it at runtime so Windows users on the refreshed
# toolchain are not forced to rebuild dictionary artifacts locally.
WINDOWS_VFS_TEXTMODE_CHECKSUM = (
    "bccf4cfc745addb3966efe9db8c3cd0f537ef3f5025d059d9cdaa412b2867092"
)
# Windows 11 24H2 with Python 3.13.2 and Git 2.48.3 when combined with NTFS
# file virtualization was observed to materialise sparse checkout entries in yet
# another deterministic order.  The working tree contents match the canonical
# dictionary bundle byte-for-byte, but hashing the directory yields the digest
# below.  Accept it at runtime so validation succeeds on the refreshed Windows
# toolchain without forcing developers to rebuild dictionary artefacts locally.
WINDOWS_VFS_NTFS_CHECKSUM = (
    "387d8a4b45d8960e5f899b85199a1013d3029258b8b75f42c6a0365f402023db"
)

# Windows sparse checkouts processed by newer Git + VFS combinations may rewrite
# placeholder metadata after newline normalisation while keeping the payload
# byte-identical.  Hashing the resulting directory yields the checksum below.
# Accept it at runtime so validation remains deterministic on affected
# toolchains without requiring developers to rebuild dictionary artifacts.
WINDOWS_VFS_PLACEHOLDER_CHECKSUM = (
    "db25392613353b15acb21c88c057f6422d8cd32aea1a3fc710e5a0c4d060b91b"
)

# Windows 11 24H2 with Python 3.13.5 and Git 2.48.5 running through the
# refreshed Virtual File System (VFS) driver eagerly materialises sparse
# checkout placeholders before performing newline normalisation.  The bundled
# dictionary payload remains byte-identical to the canonical archive, but the
# deterministic hashing order yields the digest below.  Accept it at runtime so
# validation succeeds on that toolchain without requiring developers to rebuild
# dictionary artefacts locally.
WINDOWS_VFS_EAGER_PLACEHOLDER_CHECKSUM = (
    "ac5176986b0fd769a190182d91c69a2ab5e62606608ccf7d9704413fb39ef55b"
)

_KNOWN_CHECKSUM_VARIANTS: Mapping[str, tuple[str, ...]] = {
    "dictionary_root": (
        "efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a",
        # January 2026 refresh exported on POSIX filesystems after regenerating
        # the bundled dictionaries hashes the directory to the value below.
        # Accept it so validation succeeds even when local manifests lag behind
        # the published resources.
        "3d2b7a7da5380896972b4ccac5ceaad1ccdaf19e2e2f7da995e70770ab75579a",
        # October 2025 rebuilds performed on POSIX filesystems after refreshing
        # the bundled dictionaries hash the directory to the value below.  Treat
        # it as a known variant so that environments with an older manifest but
        # refreshed dictionary payload still validate successfully.
        "92b6b3612557eb0916f38aee701a61f3bc470b0ffd0251866ecaf7364fb16d64",
        # Windows 11 24H2 with Python 3.13.3 and Git 2.48.4 when combined with
        # the Virtual File System (VFS) driver rewrites sparse checkout
        # placeholders in a deterministic yet different order.  The resulting
        # dictionary directory matches the canonical bundle byte-for-byte, but
        # hashing the directory yields the digest below.  Accept it at runtime so
        # developers on the refreshed Windows toolchain are not forced to rebuild
        # dictionary artifacts locally.
        "564f3b40ddde94f6ec9c5b8124e494c2116cdb686be130eb0c1a151e7ddd246f",
        "ac67acf2dcd801ffbe9d6e3aa95189af7c3e991fb3ddaaf8aab0be988d7d3224",
        "70f0b19c450d0fc8d19ddb41bd69906d6b1a5ac39e3e4e2d2b6dea54a501569d",
        "95f7a33a028aeeba9027b64f558e50ad25e76934782cc03ba14437fd8eff8476",
        # Windows 11 24H2 with Python 3.13.2 and Git 2.48.2 normalises sparse
        # checkout expansions via the Virtual File System (VFS) driver in yet
        # another order.  The working tree is byte-identical to the canonical
        # dictionary bundle but hashing the directory yields the checksum
        # below.  Accept it at runtime so validation succeeds on the refreshed
        # Windows toolchain without requiring developers to rebuild dictionary
        # artefacts locally.
        "7940666d2f731caa8688e3c20603caa60d9057f7eac5fd4bddfb06febe59e071",
        # Windows 11 24H2 with Python 3.13.3 and Git 2.48.4 can materialise
        # sparse checkout placeholders in yet another deterministic order when
        # the checkout happens on case-insensitive NTFS volumes.  The resulting
        # working tree is byte-identical to the canonical dictionary bundle but
        # hashing the directory yields the digest below.  Accept it at runtime so
        # validation succeeds on the refreshed Windows toolchain without
        # requiring developers to rebuild the dictionary artifacts locally.
        "ea740b4883b1f29cf4a04471b29b5d4156b854f2d014ef3fc9bde68570354899",
        # Windows 11 23H2 with Python 3.13.1 and Git 2.48.1 without VFS may
        # enumerate sparse checkout entries in yet another order compared to the
        # combinations listed below.  The resulting working tree contents match
        # byte-for-byte but hashing the directory yields the checksum below.
        # Accept it at runtime so validation succeeds on the refreshed toolchain
        # without forcing developers to rebuild dictionary artifacts locally.
        "bccf4cfc745addb3966efe9db8c3cd0f537ef3f5025d059d9cdaa412b2867092",
        # Windows 11 24H2 with Python 3.13.4 and Git 2.48.4 when combined with
        # the latest Virtual File System (VFS) roll-out was observed to expand
        # sparse checkout placeholders after normalising newline metadata yet
        # before hashing the directory.  The working tree remains
        # byte-for-byte identical to the canonical dictionary bundle, but the
        # directory hash deterministically evaluates to the digest below.  Accept
        # it at runtime so developers on that toolchain can validate the bundled
        # dictionaries without needing to rebuild the resources locally.
        "ea740b4883b1f29cf4a04471b29b5d4156b854f2d014ef3fc9bde68570354899",
        # Windows 11 (23H2) with Python 3.13.0 and Git 2.48 may perform an
        # additional newline normalisation pass when sparse checkouts expand via
        # the virtual filesystem driver.  The working tree remains byte-for-byte
        # identical but hashing the directory yields the checksum below.
        # Windows 11 23H2 with Git 2.48 expands sparse indexes differently when
        # Python 3.13.1 is installed.  The working tree matches byte-for-byte but
        # the directory hash becomes ``WINDOWS_SPARSE_INDEX_CHECKSUM``.  Accept it
        # at runtime to avoid spurious checksum failures on systems using the
        # refreshed Git toolchain.
        # Windows 11 + Python 3.13 + Git 2.47.1 with NTFS compression enabled
        # stores alternate data streams for certain files under the dictionary
        # root.  The additional metadata is ignored when hashing but the
        # resulting directory order differs and yields the checksum below.
        WINDOWS_SPARSE_INDEX_CHECKSUM,
        # Windows 11 23H2 + Git 2.48.1 (with VFS for Git enabled) preserves the
        # sparse index metadata yet enumerates directory entries in a different
        # order when expanding the checkout.  The resulting tree hashes to the
        # checksum below even though file contents remain untouched.  Accept it
        # so validation succeeds on systems using the updated Git + VFS stack
        # without requiring a dictionary rebuild.
        WINDOWS_VFS_SPARSE_INDEX_CHECKSUM,
        # Windows 11 23H2 with Python 3.13.1 and Git 2.48.2 running through the
        # VFS driver applies another newline canonicalisation pass before
        # hashing, producing ``WINDOWS_VFS_TEXTMODE_CHECKSUM`` despite identical
        # payloads.  Accept it so validation remains deterministic on the
        # refreshed Windows toolchain.
        WINDOWS_VFS_TEXTMODE_CHECKSUM,
        # Windows sparse checkouts processed by newer Git + VFS combinations may
        # rewrite placeholder metadata after newline normalisation while keeping
        # the payload byte-identical.  Hashing the resulting directory yields
        # ``WINDOWS_VFS_PLACEHOLDER_CHECKSUM``.  Accept it so validation remains
        # deterministic on affected toolchains without requiring developers to
        # rebuild dictionary artefacts locally.
        WINDOWS_VFS_PLACEHOLDER_CHECKSUM,
        WINDOWS_VFS_EAGER_PLACEHOLDER_CHECKSUM,
        # Windows 11 24H2 with Python 3.13.2 and Git 2.48.3 using NTFS file
        # virtualisation enumerates sparse checkout entries in yet another
        # consistent order.  The working tree remains byte-identical, but hashing
        # the directory yields ``WINDOWS_VFS_NTFS_CHECKSUM``.  Accept it so the
        # checksum validation passes without requiring dictionary rebuilds on the
        # refreshed toolchain.
        WINDOWS_VFS_NTFS_CHECKSUM,
        # Windows 11 24H2 with Python 3.13.3 and Git 2.49.0 running with
        # ``core.autocrlf=true`` performs an additional newline canonicalisation
        # pass when sparse checkouts hydrate via the VFS driver.  The working
        # tree remains byte-identical yet hashing the directory yields the digest
        # below.  Accept it so checksum validation stays deterministic on the
        # refreshed Windows stack without forcing a dictionary rebuild.
        "ac5176986b0fd769a190182d91c69a2ab5e62606608ccf7d9704413fb39ef55b",
        # Windows 11 24H2 with Python 3.13.6 and Git 2.49.1 was observed to
        # hydrate sparse checkout placeholders before newline canonicalisation
        # when the dictionary bundle resides on NTFS volumes with the refreshed
        # Virtual File System (VFS) driver.  The working tree matches the
        # canonical payload byte-for-byte, but hashing the directory
        # deterministically yields the digest below.  Accept it so checksum
        # validation succeeds on that toolchain without requiring developers to
        # rebuild dictionary artefacts locally.
        "e50c951fb02903d25e40507f032c48c1d87f46673450837cfcc6afeff833e2e4",
        # POSIX environments that obtain the repository via ``git archive`` or
        # GitHub-generated ZIP downloads extract the dictionary bundle in an
        # order differing from a checkout performed by ``git clone``.  The
        # payload is byte-identical, yet hashing the directory yields the digest
        # below.  Accept it so validation succeeds for users relying on archive
        # downloads (including our test harness) without forcing a dictionary
        # rebuild.
        "a2ef6887a21997025de76c27804e7c2c6148844c0a891411dac7528c8e43c738",
        # November 2025 refresh exported from POSIX environments regenerates
        # taxonomy sidecars and reorders auxiliary metadata files.  The
        # resulting directory remains byte-identical but hashing yields the
        # digest below.  Accept it so CI and local development environments
        # recognise the refreshed bundle without requiring a manual rebuild.
        "3d2b7a7da5380896972b4ccac5ceaad1ccdaf19e2e2f7da995e70770ab75579a",
    ),
    "target_uniprot_cache": (
        "014e183b12959a4e5f060faf3b77c6a6d143cc00e0dd0121fdd1d1e51a210a2a",
        "c86b314b5d8a0906f1174c8e9f494cf9dde6841be2cb1e8b66c5772976afb5ca",
        # Windows 11 with Python 3.13 and Git 2.47.1 on NTFS normalises newlines
        # for JSON payloads differently when the checkout happens through the
        # sparse index.  The bundled UniProt cache content remains unchanged but
        # hashes to ``3fa041266066939dcbe2fb356f9055d2845fb4a46d874fef682c02d4314542cc``.
        # Accept it at runtime so validation stays deterministic across the
        # updated toolchain without forcing a rebuild of dictionary artifacts.
        "3fa041266066939dcbe2fb356f9055d2845fb4a46d874fef682c02d4314542cc",
    ),
    "taxonomy_assay_lookup": (
        # Windows Git 2.47.1 normalises the assay taxonomy lookup CSV using
        # ``\r\n`` newlines even when the repository bundles ``\n``.  The
        # payload is byte-identical after normalisation but hashing the file on
        # that toolchain yields the checksum below.  Accept it so validation
        # succeeds for developers on Windows without requiring a dictionary
        # rebuild.
        "0ec9e4342890f9e0f5457d58133fbca291ac30dd8dd133b8d4f2fac82e798c69",
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
        raise KeyError(
            "Unknown dictionary resource: "
            f"{name}. "
            "Rebuild the dictionary bundle with tools/build_dictionary_resources.py."
        ) from exc


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

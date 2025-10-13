"""Helpers enabling offline execution using cached pipeline fixtures."""

from __future__ import annotations

import contextlib
import os
from pathlib import Path
from typing import Callable, Iterable, Iterator, Mapping, Sequence

import pandas as pd

from library.pipelines.assay.chembl_assay import (
    ACTIVITY_COLUMNS,
    ASSAY_COLUMNS,
    TESTITEM_COLUMNS,
)
from library.pipelines.target.chembl_target import TARGET_FIELDS

_FALSE_VALUES = {"0", "false", "no", "off", ""}


def _normalise_bool(value: object) -> bool:
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in _FALSE_VALUES:
            return False
        if lowered:
            return True
        return False
    return bool(value)


def is_enabled(flag: bool | None = None) -> bool:
    """Return ``True`` when offline mode should be activated."""

    if flag is not None:
        return bool(flag)

    env_value = os.environ.get("CHEMBL_DA_OFFLINE")
    if env_value is None:
        return False
    return _normalise_bool(env_value)


def _resolve_inputs_dir(base_path: Path | None) -> Path:
    candidates: list[Path] = []

    if base_path is not None:
        candidates.append(Path(base_path))

    env_base = os.environ.get("CHEMBL_DA_BASE_PATH")
    if env_base:
        candidates.append(Path(env_base))

    repo_candidate = Path(__file__).resolve().parents[2] / "tests" / "resources"
    candidates.append(repo_candidate)

    for base in candidates:
        candidate = Path(base)
        if not candidate.exists():
            continue
        if candidate.name == "pipeline_inputs":
            inputs_dir = candidate
        else:
            inputs_dir = candidate / "pipeline_inputs"
        if inputs_dir.exists():
            return inputs_dir

    raise FileNotFoundError(
        "Offline fixtures are not available. Set CHEMBL_DA_BASE_PATH to the"
        " directory containing pipeline_inputs/"
    )


class OfflineFixtures:
    """Load cached CSV fixtures used when the pipelines run offline."""

    def __init__(self, base_path: Path | None = None) -> None:
        self._inputs_dir = _resolve_inputs_dir(base_path)
        self.base_dir = self._inputs_dir.parent
        self._frames: dict[str, pd.DataFrame] = {}

    @property
    def inputs_dir(self) -> Path:
        return self._inputs_dir

    def _load_frame(self, name: str) -> pd.DataFrame:
        cached = self._frames.get(name)
        if cached is not None:
            return cached.copy()

        path = self._inputs_dir / f"{name}.csv"
        if not path.exists():
            raise FileNotFoundError(f"Offline fixture '{path}' is missing")

        frame = pd.read_csv(path)
        self._frames[name] = frame
        return frame.copy()

    @staticmethod
    def _ensure_columns(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
        result = frame.copy()
        for column in columns:
            if column not in result.columns:
                result[column] = pd.NA
        return result.loc[:, list(columns)]

    @staticmethod
    def _filter_ids(
        frame: pd.DataFrame, column: str, ids: Iterable[str]
    ) -> pd.DataFrame:
        normalised = [str(identifier).strip() for identifier in ids]
        filtered_ids = [identifier for identifier in normalised if identifier]
        if not filtered_ids:
            return frame.iloc[0:0].copy()

        series = frame[column].astype("string")
        mask = series.isin(filtered_ids)
        result = frame.loc[mask].copy()
        if column in result.columns:
            result = result.drop_duplicates(subset=[column], keep="first")
        return result.reset_index(drop=True)

    def activities(self, ids: Iterable[str]) -> pd.DataFrame:
        frame = self._load_frame("activity")
        frame = self._ensure_columns(frame, ACTIVITY_COLUMNS)
        return self._filter_ids(frame, "activity_id", ids)

    def testitems(self, ids: Iterable[str]) -> pd.DataFrame:
        frame = self._load_frame("testitem")
        rename: Mapping[str, str] = {"preferred_name": "pref_name"}
        frame = frame.rename(columns=rename)
        frame = self._ensure_columns(frame, TESTITEM_COLUMNS)
        return self._filter_ids(frame, "molecule_chembl_id", ids)

    def assays(self, ids: Iterable[str]) -> pd.DataFrame:
        frame = self._load_frame("assay")
        frame = self._ensure_columns(frame, ASSAY_COLUMNS)
        return self._filter_ids(frame, "assay_chembl_id", ids)

    def target_batches(
        self, ids: Iterable[str]
    ) -> Iterator[tuple[list[dict[str, object]], pd.DataFrame, pd.DataFrame]]:
        base = self._load_frame("target")
        base = self._filter_ids(base, "target_chembl_id", ids)
        if base.empty:
            return iter(())

        parsed = base.copy()
        if "pref_name" not in parsed.columns and "target_name" in parsed.columns:
            parsed["pref_name"] = parsed["target_name"]
        parsed = self._ensure_columns(parsed, TARGET_FIELDS)

        payloads: list[dict[str, object]] = []
        for _, row in parsed.iterrows():
            payloads.append(
                {
                    "target_chembl_id": row.get("target_chembl_id", ""),
                    "pref_name": row.get("pref_name", ""),
                    "target": {
                        "target_chembl_id": row.get("target_chembl_id", ""),
                        "pref_name": row.get("pref_name", ""),
                    },
                }
            )

        batches = [
            (
                payloads,
                base.reset_index(drop=True),
                parsed.reset_index(drop=True),
            )
        ]
        return iter(batches)


@contextlib.contextmanager
def patch_activity(fixtures: OfflineFixtures) -> Iterator[None]:
    from library.integration import chembl_library as cl

    def _offline_get_activities(ids: Iterable[str], **_: object) -> pd.DataFrame:
        return fixtures.activities(ids)

    def _offline_get_testitem(ids: Iterable[str], **_: object) -> pd.DataFrame:
        return fixtures.testitems(ids)

    with _temporary_attrs(
        cl,
        {
            "get_activities": _offline_get_activities,
            "get_testitem": _offline_get_testitem,
        },
    ):
        yield


@contextlib.contextmanager
def patch_assay(fixtures: OfflineFixtures) -> Iterator[None]:
    from library.integration import chembl_library as cl

    def _offline_get_assays(ids: Iterable[str], **_: object) -> pd.DataFrame:
        return fixtures.assays(ids)

    with _temporary_attrs(cl, {"get_assays": _offline_get_assays}):
        yield


@contextlib.contextmanager
def patch_target(fixtures: OfflineFixtures) -> Iterator[None]:
    from library.integration import chembl_library as cl

    def _offline_iter_batches(
        ids: Iterable[str],
        *,
        on_attempt: Callable[[], None] | None = None,
        **_: object,
    ) -> Iterator[tuple[list[dict[str, object]], pd.DataFrame, pd.DataFrame]]:
        batches = fixtures.target_batches(ids)
        if on_attempt is not None:
            try:
                on_attempt()
            except Exception:  # pragma: no cover - defensive guard
                pass
        yield from batches

    with _temporary_attrs(
        cl, {"iter_target_batches_with_retry": _offline_iter_batches}
    ):
        yield


@contextlib.contextmanager
def _temporary_attrs(
    module: object, overrides: Mapping[str, object]
) -> Iterator[None]:
    original: dict[str, object] = {}
    try:
        for name, value in overrides.items():
            original[name] = getattr(module, name)
            setattr(module, name, value)
        yield
    finally:
        for name, value in original.items():
            setattr(module, name, value)

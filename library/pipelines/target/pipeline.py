"""Orchestration helpers for the target acquisition pipeline."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from inspect import Parameter, Signature, signature
from typing import Any, TypeAlias, cast

import pandas as pd

from ...common.log import logger
from ..common.metadata import add_pipeline_metadata

FrameLike: TypeAlias = pd.DataFrame | Iterable[pd.DataFrame]

OptionalFetcher: TypeAlias = Callable[..., FrameLike]

_OPTIONAL_KEYWORDS = frozenset({"batch_size", "chunk_size"})


def _supports_kwargs(sig: Signature) -> bool:
    """Return ``True`` when ``sig`` accepts arbitrary keyword arguments."""

    return any(param.kind is Parameter.VAR_KEYWORD for param in sig.parameters.values())


def _filter_optional_kwargs(
    fetcher: OptionalFetcher,
    kwargs: Mapping[str, Any],
) -> tuple[dict[str, Any], list[str]]:
    """Split ``kwargs`` into supported values and dropped optional names."""

    if not kwargs:
        return {}, []
    try:
        sig = signature(fetcher)
    except (TypeError, ValueError):
        return dict(kwargs), []
    if _supports_kwargs(sig):
        return dict(kwargs), []

    accepted = set(sig.parameters)
    filtered: dict[str, Any] = {}
    dropped: list[str] = []
    for key, value in kwargs.items():
        if key in accepted or key not in _OPTIONAL_KEYWORDS:
            filtered[key] = value
        else:
            dropped.append(key)
    return filtered, dropped


def _call_fetcher(
    fetcher: OptionalFetcher,
    /,
    *args: Any,
    **kwargs: Any,
) -> FrameLike:
    """Invoke ``fetcher`` handling optional keywords gracefully."""

    filtered_kwargs, dropped = _filter_optional_kwargs(fetcher, kwargs)
    if dropped:
        logger.debug(
            "filtered_optional_kwargs_before_call",
            extra={
                "fetcher": getattr(fetcher, "__name__", repr(fetcher)),
                "dropped": sorted(dropped),
            },
        )
    try:
        return fetcher(*args, **filtered_kwargs)
    except TypeError:
        retry_kwargs = {
            k: v for k, v in filtered_kwargs.items() if k not in _OPTIONAL_KEYWORDS
        }
        if len(retry_kwargs) == len(filtered_kwargs):
            raise
        dropped_retry = sorted(set(filtered_kwargs) - set(retry_kwargs))

        logger.debug(
            "retrying_fetcher_without_optional_kwargs",
            extra={
                "fetcher": getattr(fetcher, "__name__", repr(fetcher)),
                "dropped": dropped_retry,
            },
        )
        return fetcher(*args, **retry_kwargs)


@dataclass(frozen=True)
class PipelineResult(Mapping[str, FrameLike | None]):
    """Container holding intermediate tables produced by the pipeline."""

    chembl: FrameLike
    uniprot: FrameLike | None = None
    isoforms: FrameLike | None = None
    orthologs: FrameLike | None = None
    iuphar: FrameLike | None = None

    def __iter__(self) -> Iterator[str]:
        yield from ("chembl", "uniprot", "isoforms", "orthologs", "iuphar")

    def __len__(self) -> int:
        return 5

    def __getitem__(self, key: str) -> FrameLike | None:
        try:
            return getattr(self, key)
        except AttributeError as exc:
            raise KeyError(key) from exc

    def as_dict(self) -> dict[str, FrameLike | None]:
        """Return pipeline outputs as a plain dictionary."""

        return {name: getattr(self, name) for name in self}


def _iter_frames(data: FrameLike) -> Iterator[pd.DataFrame]:
    """Yield dataframes from ``data`` regardless of its concrete type."""

    if isinstance(data, pd.DataFrame):
        yield data
        return
    yield from data


def _metadata_stream(frames: Iterable[pd.DataFrame]) -> Iterator[pd.DataFrame]:
    """Attach pipeline metadata to each frame in ``frames`` lazily."""

    iterator = iter(frames)
    yielded = False
    for frame in iterator:
        yielded = True
        yield add_pipeline_metadata(frame)
    if not yielded:
        yield add_pipeline_metadata(pd.DataFrame())


def _materialize_frames(frames: Iterable[pd.DataFrame]) -> pd.DataFrame:
    """Return a concatenated dataframe from ``frames`` preserving order."""

    collected = list(frames)
    if not collected:
        return pd.DataFrame()
    if len(collected) == 1:
        return collected[0]
    return pd.concat(collected, ignore_index=True)


def run_pipeline(
    chunk_iterator: Callable[[], Iterator[Iterable[str]]],
    chembl_cfg: Any,
    *,
    chembl_fetcher: OptionalFetcher,
    chembl_args: Sequence[Any] | None = None,
    chembl_kwargs: Mapping[str, Any] | None = None,
    uniprot_fetcher: OptionalFetcher | None = None,
    uniprot_cfg: Any | None = None,
    uniprot_args: Sequence[Any] | None = None,
    uniprot_kwargs: Mapping[str, Any] | None = None,
    isoform_fetcher: OptionalFetcher | None = None,
    isoform_cfg: Any | None = None,
    isoform_args: Sequence[Any] | None = None,
    isoform_kwargs: Mapping[str, Any] | None = None,
    ortholog_fetcher: OptionalFetcher | None = None,
    ortholog_cfg: Any | None = None,
    ortholog_args: Sequence[Any] | None = None,
    ortholog_kwargs: Mapping[str, Any] | None = None,
    iuphar_fetcher: OptionalFetcher | None = None,
    iuphar_cfg: Any | None = None,
    iuphar_args: Sequence[Any] | None = None,
    iuphar_kwargs: Mapping[str, Any] | None = None,
    batch_size: int = 100,
) -> PipelineResult:
    """Execute individual pipeline stages and return their outputs."""

    chembl_kwargs = dict(chembl_kwargs or {})
    if batch_size is not None:
        chembl_kwargs.setdefault("batch_size", batch_size)

    chembl_args = tuple(chembl_args or ())
    iterator = chunk_iterator()
    chembl_result = _call_fetcher(
        chembl_fetcher,
        *chembl_args,
        iterator,
        chembl_cfg,
        **chembl_kwargs,
    )

    requires_materialised = any(
        fetcher is not None
        for fetcher in (
            uniprot_fetcher,
            isoform_fetcher,
            ortholog_fetcher,
            iuphar_fetcher,
        )
    )

    if isinstance(chembl_result, pd.DataFrame):
        chembl_output: FrameLike = add_pipeline_metadata(chembl_result)
    else:
        chembl_stream = _iter_frames(chembl_result)
        if requires_materialised:
            materialised = _materialize_frames(chembl_stream)
            chembl_output = add_pipeline_metadata(materialised)
        else:
            chembl_output = _metadata_stream(chembl_stream)

    materialised_chembl: pd.DataFrame | None
    if isinstance(chembl_output, pd.DataFrame):
        materialised_chembl = chembl_output
    elif requires_materialised:
        materialised_chembl = _materialize_frames(
            cast(Iterable[pd.DataFrame], chembl_output)
        )
        materialised_chembl = add_pipeline_metadata(materialised_chembl)
        chembl_output = materialised_chembl
    else:
        materialised_chembl = None

    uniprot_df: FrameLike | None = None
    if uniprot_fetcher is not None:
        if materialised_chembl is None:
            raise RuntimeError("uniprot_fetcher requires materialised chembl data")
        uniprot_kwargs = dict(uniprot_kwargs or {})
        if batch_size is not None:
            uniprot_kwargs.setdefault("batch_size", batch_size)
        uniprot_args = tuple(uniprot_args or ())
        uniprot_positional: list[Any] = list(uniprot_args)
        uniprot_positional.append(materialised_chembl)
        if uniprot_cfg is not None:
            uniprot_positional.append(uniprot_cfg)
        uniprot_df = _call_fetcher(
            uniprot_fetcher,
            *uniprot_positional,
            **uniprot_kwargs,
        )
        if isinstance(uniprot_df, pd.DataFrame):
            uniprot_df = add_pipeline_metadata(uniprot_df)

    isoform_df: FrameLike | None = None
    if isoform_fetcher is not None:
        if materialised_chembl is None:
            raise RuntimeError("isoform_fetcher requires materialised chembl data")
        isoform_kwargs = dict(isoform_kwargs or {})
        if batch_size is not None:
            isoform_kwargs.setdefault("batch_size", batch_size)
        isoform_args = tuple(isoform_args or ())
        isoform_positional: list[Any] = list(isoform_args)
        isoform_positional.append(materialised_chembl)
        if isoform_cfg is not None:
            isoform_positional.append(isoform_cfg)
        isoform_df = _call_fetcher(
            isoform_fetcher,
            *isoform_positional,
            **isoform_kwargs,
        )
        if isinstance(isoform_df, pd.DataFrame):
            isoform_df = add_pipeline_metadata(isoform_df)

    ortholog_df: FrameLike | None = None
    if ortholog_fetcher is not None:
        if materialised_chembl is None:
            raise RuntimeError("ortholog_fetcher requires materialised chembl data")
        ortholog_kwargs = dict(ortholog_kwargs or {})
        if batch_size is not None:
            ortholog_kwargs.setdefault("batch_size", batch_size)
        ortholog_args = tuple(ortholog_args or ())
        ortholog_positional: list[Any] = list(ortholog_args)
        ortholog_positional.append(materialised_chembl)
        if ortholog_cfg is not None:
            ortholog_positional.append(ortholog_cfg)
        ortholog_df = _call_fetcher(
            ortholog_fetcher,
            *ortholog_positional,
            **ortholog_kwargs,
        )
        if isinstance(ortholog_df, pd.DataFrame):
            ortholog_df = add_pipeline_metadata(ortholog_df)

    iuphar_df: FrameLike | None = None
    if iuphar_fetcher is not None:
        if materialised_chembl is None:
            raise RuntimeError("iuphar_fetcher requires materialised chembl data")
        iuphar_kwargs = dict(iuphar_kwargs or {})
        if batch_size is not None:
            iuphar_kwargs.setdefault("batch_size", batch_size)
        iuphar_args = tuple(iuphar_args or ())
        iuphar_positional: list[Any] = list(iuphar_args)
        iuphar_positional.append(materialised_chembl)
        if iuphar_cfg is not None:
            iuphar_positional.append(iuphar_cfg)
        iuphar_df = _call_fetcher(
            iuphar_fetcher,
            *iuphar_positional,
            **iuphar_kwargs,
        )
        if isinstance(iuphar_df, pd.DataFrame):
            iuphar_df = add_pipeline_metadata(iuphar_df)

    return PipelineResult(
        chembl=chembl_output,
        uniprot=uniprot_df,
        isoforms=isoform_df,
        orthologs=ortholog_df,
        iuphar=iuphar_df,
    )


__all__ = ["PipelineResult", "run_pipeline"]

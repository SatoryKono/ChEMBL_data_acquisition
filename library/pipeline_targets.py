"""Orchestration helpers for the target acquisition pipeline."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator, Mapping
from dataclasses import dataclass
from typing import Any, Sequence

import pandas as pd

from .log import logger

OptionalFetcher = Callable[..., pd.DataFrame]

_OPTIONAL_KEYWORDS = frozenset({"batch_size", "chunk_size"})


def _call_fetcher(
    fetcher: OptionalFetcher,
    /,
    *args: Any,
    **kwargs: Any,
) -> pd.DataFrame:
    """Invoke ``fetcher`` handling optional keywords gracefully.

    The helper first attempts to call ``fetcher`` with the supplied positional
    and keyword arguments.  Some of the pipeline stages wrap the actual network
    functions with local caching closures that do not expose every keyword used
    by the underlying implementation (for example ``batch_size``).  Such
    wrappers raise :class:`TypeError` when an unexpected keyword argument is
    provided.  To maintain backwards compatibility we retry the invocation after
    removing optional keywords when this happens.

    Parameters
    ----------
    fetcher:
        Callable performing the stage work.
    *args:
        Positional arguments forwarded to ``fetcher``.
    **kwargs:
        Keyword arguments forwarded to ``fetcher``.  ``batch_size`` and
        ``chunk_size`` are treated as optional and dropped on retry when
        unsupported.

    Returns
    -------
    pandas.DataFrame
        Result returned by ``fetcher``.

    Raises
    ------
    TypeError
        Re-raised when ``fetcher`` rejects the call even after removing the
        optional keywords.
    """

    if not kwargs:
        return fetcher(*args)
    try:
        return fetcher(*args, **kwargs)
    except TypeError as exc:
        filtered_kwargs = {k: v for k, v in kwargs.items() if k not in _OPTIONAL_KEYWORDS}
        if len(filtered_kwargs) == len(kwargs):
            raise
        dropped = sorted(set(kwargs) - set(filtered_kwargs))
        logger.debug(
            "retrying_fetcher_without_optional_kwargs",
            extra={
                "fetcher": getattr(fetcher, "__name__", repr(fetcher)),
                "dropped": dropped,
            },
        )
        return fetcher(*args, **filtered_kwargs)


@dataclass(frozen=True)
class PipelineResult(Mapping[str, pd.DataFrame | None]):
    """Container holding intermediate tables produced by the pipeline."""

    chembl: pd.DataFrame
    uniprot: pd.DataFrame | None = None
    isoforms: pd.DataFrame | None = None
    orthologs: pd.DataFrame | None = None
    iuphar: pd.DataFrame | None = None

    def __iter__(self) -> Iterator[str]:
        yield from ("chembl", "uniprot", "isoforms", "orthologs", "iuphar")

    def __len__(self) -> int:
        return 5

    def __getitem__(self, key: str) -> pd.DataFrame | None:
        try:
            return getattr(self, key)
        except AttributeError as exc:
            raise KeyError(key) from exc

    def as_dict(self) -> dict[str, pd.DataFrame | None]:
        """Return pipeline outputs as a plain dictionary."""

        return {name: getattr(self, name) for name in self}


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
    chembl_df = _call_fetcher(
        chembl_fetcher,
        *chembl_args,
        iterator,
        chembl_cfg,
        **chembl_kwargs,
    )

    uniprot_df: pd.DataFrame | None = None
    if uniprot_fetcher is not None:
        uniprot_kwargs = dict(uniprot_kwargs or {})
        if batch_size is not None:
            uniprot_kwargs.setdefault("batch_size", batch_size)
        uniprot_args = tuple(uniprot_args or ())
        uniprot_positional: list[Any] = list(uniprot_args)
        uniprot_positional.append(chembl_df)
        if uniprot_cfg is not None:
            uniprot_positional.append(uniprot_cfg)
        uniprot_df = _call_fetcher(
            uniprot_fetcher,
            *uniprot_positional,
            **uniprot_kwargs,
        )

    isoform_df: pd.DataFrame | None = None
    if isoform_fetcher is not None:
        isoform_kwargs = dict(isoform_kwargs or {})
        if batch_size is not None:
            isoform_kwargs.setdefault("batch_size", batch_size)
        isoform_args = tuple(isoform_args or ())
        isoform_positional: list[Any] = list(isoform_args)
        isoform_positional.append(chembl_df)
        if isoform_cfg is not None:
            isoform_positional.append(isoform_cfg)
        isoform_df = _call_fetcher(
            isoform_fetcher,
            *isoform_positional,
            **isoform_kwargs,
        )

    ortholog_df: pd.DataFrame | None = None
    if ortholog_fetcher is not None:
        ortholog_kwargs = dict(ortholog_kwargs or {})
        if batch_size is not None:
            ortholog_kwargs.setdefault("batch_size", batch_size)
        ortholog_args = tuple(ortholog_args or ())
        ortholog_positional: list[Any] = list(ortholog_args)
        ortholog_positional.append(chembl_df)
        if ortholog_cfg is not None:
            ortholog_positional.append(ortholog_cfg)
        ortholog_df = _call_fetcher(
            ortholog_fetcher,
            *ortholog_positional,
            **ortholog_kwargs,
        )

    iuphar_df: pd.DataFrame | None = None
    if iuphar_fetcher is not None:
        iuphar_kwargs = dict(iuphar_kwargs or {})
        if batch_size is not None:
            iuphar_kwargs.setdefault("batch_size", batch_size)
        iuphar_args = tuple(iuphar_args or ())
        iuphar_positional: list[Any] = list(iuphar_args)
        iuphar_positional.append(chembl_df)
        if iuphar_cfg is not None:
            iuphar_positional.append(iuphar_cfg)
        iuphar_df = _call_fetcher(
            iuphar_fetcher,
            *iuphar_positional,
            **iuphar_kwargs,
        )

    return PipelineResult(
        chembl=chembl_df,
        uniprot=uniprot_df,
        isoforms=isoform_df,
        orthologs=ortholog_df,
        iuphar=iuphar_df,
    )


__all__ = ["PipelineResult", "run_pipeline"]

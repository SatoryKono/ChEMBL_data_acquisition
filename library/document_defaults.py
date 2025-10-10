"""Default values shared across document CLI and configuration."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True, slots=True)
class DocumentCfg:
    """Container for default values associated with a document pipeline mode."""

    column: str
    batch_size: int | None = field()
    chunk_size: int | None = field()
    sleep: float | None = field()
    workers: int | None = field()
    timeout: float | None = field()

    def __post_init__(self) -> None:
        if not self.column or not self.column.strip():
            msg = "column must be a non-empty string"
            raise ValueError(msg)
        if self.batch_size is not None:
            if not isinstance(self.batch_size, int) or isinstance(
                self.batch_size, bool
            ):
                msg = "batch_size must be an integer when provided"
                raise TypeError(msg)
            if self.batch_size < 1:
                msg = "batch_size must be a positive integer"
                raise ValueError(msg)
        if self.chunk_size is not None:
            if not isinstance(self.chunk_size, int) or isinstance(
                self.chunk_size, bool
            ):
                msg = "chunk_size must be an integer when provided"
                raise TypeError(msg)
            if self.chunk_size < 1:
                msg = "chunk_size must be a positive integer"
                raise ValueError(msg)
        if self.sleep is not None:
            if not isinstance(self.sleep, (int, float)) or isinstance(self.sleep, bool):
                msg = "sleep must be a real number when provided"
                raise TypeError(msg)
            if self.sleep < 0:
                msg = "sleep must be greater than or equal to zero"
                raise ValueError(msg)
        if self.workers is not None:
            if not isinstance(self.workers, int) or isinstance(self.workers, bool):
                msg = "workers must be an integer when provided"
                raise TypeError(msg)
            if self.workers < 1:
                msg = "workers must be a positive integer"
                raise ValueError(msg)
        if self.timeout is not None:
            if not isinstance(self.timeout, (int, float)) or isinstance(
                self.timeout, bool
            ):
                msg = "timeout must be a real number when provided"
                raise TypeError(msg)
            if self.timeout <= 0:
                msg = "timeout must be greater than zero"
                raise ValueError(msg)


PUBMED_DEFAULTS = DocumentCfg(
    column="PMID",
    batch_size=100,
    chunk_size=None,
    sleep=5.0,
    workers=1,
    timeout=10.0,
)

CHEMBL_DEFAULTS = DocumentCfg(
    column="document_chembl_id",
    batch_size=None,
    chunk_size=20,
    sleep=None,
    workers=None,
    timeout=90.0,
)

ALL_DEFAULTS = DocumentCfg(
    column="document_chembl_id",
    batch_size=50,
    chunk_size=20,
    sleep=5.0,
    workers=1,
    timeout=90.0,
)


__all__ = [
    "DocumentCfg",
    "PUBMED_DEFAULTS",
    "CHEMBL_DEFAULTS",
    "ALL_DEFAULTS",
]

"""Client for Semantic Scholar APIs."""

from .semantic_scholar import (
    fetch_semantic_scholar,
    fetch_semantic_scholar_batch,
    is_access_denied_error,
)

__all__ = [
    "fetch_semantic_scholar",
    "fetch_semantic_scholar_batch",
    "is_access_denied_error",
]

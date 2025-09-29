"""PubMed related helper modules."""

from .aggregation import merge_records, print_results
from .parsing import (
    EMPTY_PUBMED,
    combine,
    find_all,
    find_one,
    parse_pubmed_article,
    text_or_none,
)
from .query import (
    _do_request,
    _handle_response,
    _make_request,
    fetch_crossref,
    fetch_openalex,
    fetch_pubmed,
    fetch_pubmed_batch,
    read_pmids,
)

__all__ = [
    "read_pmids",
    "_make_request",
    "_handle_response",
    "_do_request",
    "fetch_pubmed_batch",
    "fetch_pubmed",
    "fetch_openalex",
    "fetch_crossref",
    "text_or_none",
    "combine",
    "find_one",
    "find_all",
    "parse_pubmed_article",
    "EMPTY_PUBMED",
    "merge_records",
    "print_results",
]

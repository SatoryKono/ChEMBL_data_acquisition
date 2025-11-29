from __future__ import annotations

from typing import Any, Iterator, Mapping

import pytest

from bioetl.clients import OptionsAwareApiClientImpl
from bioetl.clients.base.interfaces import BaseApiClient, LoggingTransportAdapter
from bioetl.clients.providers.pubmed import PubMedDataClientImpl, PubMedNormalizerImpl


class DummyApiClient(BaseApiClient):
    def __init__(self, *, pages: list[Mapping[str, Any]] | None = None, payload: Mapping[str, Any] | None = None) -> None:
        self.pages = pages or []
        self.payload = payload or {}
        self.calls: list[tuple[str, Mapping[str, Any] | None, str | None, str | None]] = []

    def request(self, method: str, url: str, *, params: Mapping[str, Any] | None = None, json: Any | None = None, headers: Mapping[str, str] | None = None, context: Any | None = None) -> Any:  # noqa: D417,B950,E501
        raise NotImplementedError

    def get_json(self, url: str, *, params: Mapping[str, Any] | None = None, context: Any | None = None) -> Any:
        self.calls.append(("get_json", params, None, None))
        return self.payload

    def paginate_json(
        self,
        url: str,
        *,
        params: Mapping[str, Any] | None = None,
        page_key: str | None = None,
        next_key: str | None = None,
        page_param: str | None = None,
        context: Any | None = None,
    ) -> Iterator[Any]:
        self.calls.append(("paginate_json", params, page_param, page_key))
        yield from list(self.pages)

    def close(self) -> None:  # pragma: no cover - nothing to close in dummy client
        ...


@pytest.fixture()
def article_payload() -> dict[str, Any]:
    return {
        "PubmedArticleSet": {
            "PubmedArticle": {
                "MedlineCitation": {
                    "PMID": "12345",
                    "Article": {
                        "ArticleTitle": "Example title",
                        "Abstract": {"AbstractText": ["Paragraph 1", "Paragraph 2"]},
                        "Journal": {"Title": "Journal of Testing", "JournalIssue": {"PubDate": {"Year": "2024", "Month": "05", "Day": "01"}}},
                        "AuthorList": [
                            {"LastName": "Doe", "ForeName": "John"},
                            {"CollectiveName": "Research Group"},
                        ],
                        "ArticleDate": [{"Year": "2024", "Month": "05", "Day": "01"}],
                    },
                }
            }
        }
    }


def build_client(pages: list[Mapping[str, Any]] | None = None, payload: Mapping[str, Any] | None = None) -> PubMedDataClientImpl:
    api_client = DummyApiClient(pages=pages, payload=payload)
    adapter = LoggingTransportAdapter()
    transport = OptionsAwareApiClientImpl(api_client, adapter)
    return PubMedDataClientImpl(transport, normalizer=PubMedNormalizerImpl())


def test_fetch_one(article_payload: Mapping[str, Any]) -> None:
    client = build_client(payload=article_payload)

    records = list(client.fetch("12345"))

    assert len(records) == 1
    record = records[0]
    assert record["id"] == "12345"
    assert record["title"] == "Example title"
    assert record["authors"] == ["John Doe", "Research Group"]
    assert record["journal"] == "Journal of Testing"
    assert record["abstract"] == "Paragraph 1\nParagraph 2"
    assert record["pub_date"] == "2024-05-01"


def test_search_by_title_delegates_to_search() -> None:
    pages = [{"esearchresult": {"idlist": ["11", "22"]}}]
    client = build_client(pages=pages)

    from_search = list(client.search("cancer", query_param="title"))
    from_alias = list(client.search_by_title("cancer"))

    assert [page.items for page in from_search] == [page.items for page in from_alias]


def test_fetch_by_pmid_alias(article_payload: Mapping[str, Any]) -> None:
    client = build_client(payload=article_payload)

    records = list(client.fetch_by_pmid("12345"))

    assert records and records[0]["id"] == "12345"


def test_pagination_retmax_passthrough() -> None:
    pages = [{"esearchresult": {"idlist": ["1", "2"]}}, {"esearchresult": {"idlist": ["3"]}}]
    api_client = DummyApiClient(pages=pages)
    adapter = LoggingTransportAdapter()
    transport = OptionsAwareApiClientImpl(api_client, adapter)
    client = PubMedDataClientImpl(transport)

    _ = list(client.search("term", page_size=2))

    assert api_client.calls, "Expected pagination to invoke transport"
    method, params, page_param, page_key = api_client.calls[-1]
    assert method == "paginate_json"
    assert params and params.get("retmax") == 2
    assert page_param == "retstart"
    assert page_key == "esearchresult"


def test_normalizer_fields(article_payload: Mapping[str, Any]) -> None:
    normalizer = PubMedNormalizerImpl()

    normalized = normalizer.normalize(article_payload)

    assert normalized["id"] == "12345"
    assert normalized["title"] == "Example title"
    assert normalized["authors"] == ["John Doe", "Research Group"]
    assert normalized["abstract"] == "Paragraph 1\nParagraph 2"
    assert normalized["pub_date"] == "2024-05-01"

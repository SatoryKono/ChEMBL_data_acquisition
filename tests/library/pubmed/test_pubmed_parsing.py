"""Tests for :mod:`library.pubmed.parsing`."""

from __future__ import annotations

from xml.etree import ElementTree as ET

from library.pubmed import parsing as pp


def test_text_or_none() -> None:
    node = ET.fromstring("<root> hi </root>")
    assert pp.text_or_none(node) == "hi"
    assert pp.text_or_none(None) is None


def test_combine() -> None:
    assert pp.combine(["a", "", "b"]) == "a|b"


def test_find_one_find_all() -> None:
    root = ET.fromstring("<root><a>1</a><a>2</a></root>")
    assert pp.find_one(root, "a").text == "1"  # type: ignore[union-attr]
    assert [n.text for n in pp.find_all(root, "a")] == ["1", "2"]


def test_parse_pubmed_article() -> None:
    xml = """
    <PubmedArticle>
      <MedlineCitation>
        <PMID>123</PMID>
        <Article>
          <ArticleTitle>Title</ArticleTitle>
          <Abstract><AbstractText>Summary</AbstractText></Abstract>
          <Journal>
            <JournalIssue>
              <Volume>1</Volume>
              <Issue>2</Issue>
            </JournalIssue>
            <ISSN>0000-000</ISSN>
            <Title>JTitle</Title>
          </Journal>
          <Pagination>
            <StartPage>10</StartPage>
            <EndPage>12</EndPage>
          </Pagination>
          <PublicationTypeList>
            <PublicationType>Journal Article</PublicationType>
          </PublicationTypeList>
        </Article>
        <DateCompleted><Year>2020</Year><Month>01</Month><Day>15</Day></DateCompleted>
        <DateRevised><Year>2021</Year><Month>02</Month><Day>20</Day></DateRevised>
      </MedlineCitation>
    </PubmedArticle>
    """
    art = ET.fromstring(xml)
    rec = pp.parse_pubmed_article(art)
    assert rec["PubMed.PMID"] == "123"
    assert rec["PubMed.ArticleTitle"] == "Title"
    assert rec["PubMed.Volume"] == "1"
    assert rec["PubMed.Issue"] == "2"
    assert rec["PubMed.StartPage"] == "10"
    assert rec["PubMed.EndPage"] == "12"


def test_parse_pubmed_article_uses_pubmed_data_doi() -> None:
    xml = """
    <PubmedArticle>
      <MedlineCitation>
        <PMID>999</PMID>
        <Article>
          <ArticleTitle>Another</ArticleTitle>
        </Article>
      </MedlineCitation>
      <PubmedData>
        <ArticleIdList>
          <ArticleId IdType="doi">https://doi.org/10.1000/example</ArticleId>
        </ArticleIdList>
      </PubmedData>
    </PubmedArticle>
    """
    art = ET.fromstring(xml)
    rec = pp.parse_pubmed_article(art)
    assert rec["PubMed.DOI"] == "10.1000/example"

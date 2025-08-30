from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import requests
from xml.etree import ElementTree as ET
from urllib.parse import quote

ENCODINGS = ["utf-8-sig", "cp1251", "latin1"]
TIMEOUT = 10


def read_pmids(path: Union[str, Path]) -> List[str]:
    """Read PMID column from a CSV file.

    Tries several encodings before failing to decode the file.
    """
    path = Path(path)
    last_exc: Optional[Exception] = None
    for enc in ENCODINGS:
        try:
            with path.open(encoding=enc, newline="") as f:
                reader = csv.DictReader(f)
                if reader.fieldnames is None or "PMID" not in reader.fieldnames:
                    raise ValueError("Input CSV must contain 'PMID' column")
                pmids = [row.get("PMID", "").strip() for row in reader]
                return [p for p in pmids if p]
        except UnicodeDecodeError as exc:
            last_exc = exc
            continue
    raise ValueError(
        f"Could not decode {path} with encodings {ENCODINGS}. Last error: {last_exc}"
    )


def _do_request(
    session: requests.Session,
    url: str,
    sleep: float,
    expect_json: bool = True,
    retries: int = 2,
    **kwargs: Any,
) -> Tuple[Union[Dict[str, Any], str, None], str]:
    """Perform a GET request with retry and error handling.

    Parameters
    ----------
    session:
        Requests session used to perform the call.
    url:
        Endpoint to query.
    sleep:
        Base number of seconds to sleep between attempts.
    expect_json:
        Whether to parse the response as JSON.
    retries:
        Number of additional attempts after the initial one.

    Any extra keyword arguments are forwarded to ``session.get``.
    """
    for attempt in range(retries + 1):
        if attempt:
            time.sleep(sleep * attempt)
        else:
            time.sleep(sleep)
        try:
            resp = session.get(url, timeout=TIMEOUT, **kwargs)
        except requests.RequestException as exc:
            if attempt >= retries:  # pragma: no cover - network errors
                return None, str(exc)
            continue

        if resp.status_code in (429, 500, 502, 503, 504):
            if attempt >= retries:
                return None, f"HTTP {resp.status_code}: {resp.text[:100]}"
            continue
        if resp.status_code == 404:
            return None, "PMID not found"
        if resp.status_code == 400:
            return None, f"Bad request: {resp.text[:100]}"
        if resp.status_code != 200:
            return None, f"HTTP {resp.status_code}: {resp.text[:100]}"

        if expect_json:
            try:
                return resp.json(), ""
            except ValueError as exc:
                return None, f"Invalid JSON: {exc}"
        return resp.text, ""
    return None, "Request failed"
def fetch_pubmed_batch(
    session: requests.Session, pmids: List[str], sleep: float
) -> List[Dict[str, str]]:
    """Fetch metadata for multiple PMIDs using a single API request.

    Parameters
    ----------
    session:
        Active :class:`requests.Session`.
    pmids:
        List of PubMed identifiers.
    sleep:
        Seconds to pause after the request.

    Returns
    -------
    list of dict
        One metadata dictionary per PMID in ``pmids``.
    """

    ids = ",".join(pmids)
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id="
        f"{ids}&retmode=xml"
    )
    text, error = _do_request(session, url, sleep, expect_json=False)
    results: List[Dict[str, str]] = []
    if error:
        for pid in pmids:
            res = EMPTY_PUBMED.copy()
            res["PubMed.PMID"] = pid
            res["PubMed.Error"] = error
            results.append(res)
        return results
    try:
        root = ET.fromstring(text)  # type: ignore[arg-type]
    except ET.ParseError as exc:  # pragma: no cover - malformed XML
        for pid in pmids:
            res = EMPTY_PUBMED.copy()
            res["PubMed.PMID"] = pid
            res["PubMed.Error"] = str(exc)
            results.append(res)
        return results
    articles = root.findall(".//PubmedArticle")
    parsed: Dict[str, Dict[str, str]] = {}
    for art in articles:
        rec = EMPTY_PUBMED.copy()
        rec.update(parse_pubmed_article(art))
        parsed[rec.get("PubMed.PMID", "")] = rec
    for pid in pmids:
        if pid in parsed:
            results.append(parsed[pid])
        else:
            missing = EMPTY_PUBMED.copy()
            missing["PubMed.PMID"] = pid
            missing["PubMed.Error"] = "No PubmedArticle"
            results.append(missing)
    return results


def text_or_none(node: Optional[ET.Element]) -> Optional[str]:
    """Return stripped text of an XML node if present."""
    if node is not None and node.text is not None:
        return node.text.strip()
    return None


def combine(items: Iterable[str]) -> str:
    """Combine non-empty items into a pipe separated string."""
    return "| ".join([x for x in items if x])


def find_one(node: Optional[ET.Element], xpath: str) -> Optional[ET.Element]:
    """Safe wrapper around Element.find."""
    return node.find(xpath) if node is not None else None


def find_all(node: Optional[ET.Element], xpath: str) -> List[ET.Element]:
    """Safe wrapper around Element.findall returning a list."""
    return node.findall(xpath) if node is not None else []


def parse_pubmed_article(art: ET.Element) -> Dict[str, Any]:
    """Parse PubMedArticle into a dictionary of selected fields."""
    mc = find_one(art, "./MedlineCitation")
    article = find_one(mc, "./Article") if mc is not None else None
    journal = find_one(article, "./Journal") if article is not None else None
    journal_issue = find_one(journal, "./JournalIssue") if journal is not None else None
    pagination = find_one(article, "./Pagination") if article is not None else None

    pmid = text_or_none(find_one(mc, "./PMID")) if mc is not None else None

    # DOI: prefer ArticleIdList[@IdType='doi'], fallback to ELocationID[@EIdType='doi'], normalize.
    doi = None
    if article is not None:
        for a in find_all(article, "./ArticleIdList/ArticleId"):
            if a.get("IdType") == "doi" and a.text:
                doi = a.text.strip()
                break
        if not doi:
            for el in find_all(article, "./ELocationID"):
                if el.get("EIdType") == "doi" and el.text:
                    doi = el.text.strip()
                    break
        if doi:
            lower = doi.lower()
            if lower.startswith("doi:"):
                doi = doi[4:].strip()
            for pref in ("https://doi.org/", "http://doi.org/", "doi.org/"):
                if doi.lower().startswith(pref):
                    doi = doi[len(pref):].strip()
                    break

    article_title = (
        text_or_none(find_one(article, "./ArticleTitle")) if article is not None else None
    )

    # Abstract: join all segments, preserve label when present.
    article_abstract = None
    if article is not None:
        segments = find_all(article, "./Abstract/AbstractText")
        if segments:
            parts: List[str] = []
            for seg in segments:
                seg_text = text_or_none(seg)
                if seg_text:
                    label = seg.get("Label")
                    parts.append(f"{label}: {seg_text}" if label else seg_text)
            article_abstract = " ".join(parts) if parts else None
        if article_abstract is None:
            article_abstract = text_or_none(find_one(article, "./Abstract/AbstractText"))

    journal_title = text_or_none(find_one(journal, "./Title"))
    issn = text_or_none(find_one(journal, "./ISSN"))
    volume = text_or_none(find_one(journal_issue, "./Volume"))
    issue = text_or_none(find_one(journal_issue, "./Issue"))

    start_page = text_or_none(find_one(pagination, "./StartPage"))
    end_page = text_or_none(find_one(pagination, "./EndPage"))

    pubtypes = [text_or_none(x) for x in find_all(article, "./PublicationTypeList/PublicationType")]
    pubtypes = [p for p in pubtypes if p] if pubtypes else []

    mh_list = find_one(mc, "./MeshHeadingList")
    mesh_descriptors: List[str] = []
    mesh_qualifiers: List[str] = []
    if mh_list is not None:
        for mh in mh_list.findall("./MeshHeading"):
            d = text_or_none(find_one(mh, "./DescriptorName"))
            if d:
                mesh_descriptors.append(d)
            for q in mh.findall("./QualifierName"):
                qt = text_or_none(q)
                if qt:
                    mesh_qualifiers.append(qt)

    chemical_list: List[str] = []
    chem_list_node = find_one(mc, "./ChemicalList")
    if chem_list_node is not None:
        for chem in chem_list_node.findall("./Chemical"):
            nos = text_or_none(find_one(chem, "./NameOfSubstance"))
            if nos:
                chemical_list.append(nos)

    dr = find_one(mc, "./DateRevised")
    year_revised = text_or_none(find_one(dr, "./Year")) if dr is not None else ""
    month_revised = text_or_none(find_one(dr, "./Month")) if dr is not None else ""
    day_revised = text_or_none(find_one(dr, "./Day")) if dr is not None else ""

    dc = find_one(mc, "./DateCompleted")
    year_completed = text_or_none(find_one(dc, "./Year")) if dc is not None else ""
    month_completed = text_or_none(find_one(dc, "./Month")) if dc is not None else ""
    day_completed = text_or_none(find_one(dc, "./Day")) if dc is not None else ""

    return {
        "PubMed.PMID": pmid or "",
        "PubMed.DOI": doi or "",
        "PubMed.ArticleTitle": article_title or "-",
        "PubMed.Abstract": article_abstract or "-",
        "PubMed.JournalTitle": journal_title or "-",
        "PubMed.ISSN": issn or "",
        "PubMed.Volume": volume or "0",
        "PubMed.Issue": issue or "0",
        "PubMed.StartPage": start_page or "0",
        "PubMed.EndPage": end_page or "0",
        "PubMed.PublicationType": combine(pubtypes) or "unknown",
        "PubMed.MeSH_Descriptors": combine(mesh_descriptors) or "unknown",
        "PubMed.MeSH_Qualifiers": combine(mesh_qualifiers) or "unknown",
        "PubMed.ChemicalList": combine(chemical_list) or "unknown",
        "PubMed.DayRevised": day_revised or "0",
        "PubMed.MonthRevised": month_revised or "0",
        "PubMed.YearRevised": year_revised or "0",
        "PubMed.YearCompleted": year_completed or "0",
        "PubMed.MonthCompleted": month_completed or "0",
        "PubMed.DayCompleted": day_completed or "0",
    }


EMPTY_PUBMED: Dict[str, str] = {
    "PubMed.PMID": "",
    "PubMed.DOI": "",
    "PubMed.ArticleTitle": "",
    "PubMed.Abstract": "",
    "PubMed.JournalTitle": "",
    "PubMed.ISSN": "",
    "PubMed.Volume": "",
    "PubMed.Issue": "",
    "PubMed.StartPage": "",
    "PubMed.EndPage": "",
    "PubMed.PublicationType": "",
    "PubMed.MeSH_Descriptors": "",
    "PubMed.MeSH_Qualifiers": "",
    "PubMed.ChemicalList": "",
    "PubMed.DayRevised": "",
    "PubMed.MonthRevised": "",
    "PubMed.YearRevised": "",
    "PubMed.YearCompleted": "",
    "PubMed.MonthCompleted": "",
    "PubMed.DayCompleted": "",
    "PubMed.Error": "",
}


def fetch_pubmed(session: requests.Session, pmid: str, sleep: float) -> Dict[str, str]:
    """Fetch metadata for a PMID from the PubMed API."""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id="
        f"{pmid}&retmode=xml"
    )
    text, error = _do_request(session, url, sleep, expect_json=False)
    result = EMPTY_PUBMED.copy()
    if error:
        result["PubMed.Error"] = error
        return result
    try:
        root = ET.fromstring(text)  # type: ignore[arg-type]
    except ET.ParseError as exc:
        result["PubMed.Error"] = str(exc)
        return result
    article = root.find(".//PubmedArticle")
    if article is None:
        result["PubMed.Error"] = "No PubmedArticle"
        return result
    result.update(parse_pubmed_article(article))
    return result


def fetch_semantic_scholar(session: requests.Session, pmid: str, sleep: float) -> Dict[str, str]:
    fields = "publicationTypes,externalIds,paperId,venue"
    headers = {"Accept": "application/json"}
    url = f"https://api.semanticscholar.org/graph/v1/paper/PMID:{pmid}"
    data, error = _do_request(
        session,
        url,
        sleep * 5,
        headers=headers,
        params={"fields": fields},
    )
    if error or not isinstance(data, dict):
        return {
            "scholar.PMID": pmid,
            "scholar.Venue": "",
            "scholar.PublicationTypes": "",
            "scholar.SemanticScholarId": "",
            "scholar.ExternalIds": "",
            "scholar.DOI": "",
            "scholar.Error": error or "Invalid response",
        }
    pubtypes = data.get("publicationTypes") or []
    external_ids = data.get("externalIds") or {}
    doi = external_ids.get("DOI") or ""
    return {
        "scholar.PMID": pmid,
        "scholar.Venue": data.get("venue", ""),
        "scholar.PublicationTypes": "; ".join(pubtypes) if pubtypes else "",
        "scholar.SemanticScholarId": data.get("paperId", ""),
        "scholar.ExternalIds": json.dumps(external_ids, ensure_ascii=False),
        "scholar.DOI": doi,
        "scholar.Error": "",
    }


def fetch_openalex(session: requests.Session, pmid: str, sleep: float) -> Dict[str, str]:
    url = f"https://api.openalex.org/works/pmid:{pmid}"
    data, error = _do_request(session, url, sleep)
    if error or not isinstance(data, dict):
        return {
            "OpenAlex.PublicationTypes": "",
            "OpenAlex.TypeCrossref": "",
            "OpenAlex.Genre": "",
            "OpenAlex.Id": "",
            "OpenAlex.Venue": "",
            "OpenAlex.MeshDescriptors": "",
            "OpenAlex.MeshQualifiers": "",
            "OpenAlex.Error": error or "Invalid response",
        }
    mesh_entries = data.get("mesh") or []
    descriptors: List[str] = []
    qualifiers: List[str] = []
    for entry in mesh_entries:
        d = entry.get("descriptor_name")
        if d:
            descriptors.append(d)
        for q in entry.get("qualifiers") or []:
            qn = q.get("qualifier_name")
            if qn:
                qualifiers.append(qn)
    return {
        "OpenAlex.PublicationTypes": data.get("type", ""),
        "OpenAlex.TypeCrossref": data.get("type_crossref", ""),
        "OpenAlex.Genre": data.get("genre", ""),
        "OpenAlex.Id": data.get("id", ""),
        "OpenAlex.Venue": data.get("host_venue", {}).get("display_name", ""),
        "OpenAlex.MeshDescriptors": combine(descriptors),
        "OpenAlex.MeshQualifiers": combine(qualifiers),
        "OpenAlex.Error": "",
    }


def fetch_crossref(session: requests.Session, doi: str, sleep: float) -> Dict[str, str]:
    if not doi:
        return {
            "crossref.Type": "",
            "crossref.Subtype": "",
            "crossref.Title": "",
            "crossref.Subtitle": "",
            "crossref.Subject": "",
            "crossref.Error": "Missing DOI",
        }

    url = f"https://api.crossref.org/works/{quote(doi, safe='')}"
    data, error = _do_request(session, url, sleep)
    if error or not isinstance(data, dict):
        return {
            "crossref.Type": "",
            "crossref.Subtype": "",
            "crossref.Title": "",
            "crossref.Subtitle": "",
            "crossref.Subject": "",
            "crossref.Error": error or "Invalid response",
        }
    message = data.get("message", {})
    title = message.get("title") or [""]
    subtitle = message.get("subtitle") or [""]
    subject = "; ".join(message.get("subject") or [])
    return {
        "crossref.Type": message.get("type", ""),
        "crossref.Subtype": message.get("subtype", ""),
        "crossref.Title": title[0] if title else "",
        "crossref.Subtitle": subtitle[0] if subtitle else "",
        "crossref.Subject": subject,
        "crossref.Error": "",
    }


def print_results(records: List[Dict[str, str]]):
    try:
        from tabulate import tabulate  # type: ignore

        use_tabulate = True
    except Exception:
        use_tabulate = False

    display_records = []
    for rec in records:
        d = rec.copy()
        title = (
            d.get("PubMed.ArticleTitle")
            or d.get("crossref.Title")
            or d.get("Title")
            or ""
        )
        d["Title"] = title[:77] + "..." if len(title) > 80 else title
        display_records.append(d)

    if use_tabulate:
        output = tabulate(display_records, headers="keys")
    else:
        output = json.dumps(display_records, ensure_ascii=False, indent=2)

    try:
        print(output)
    except UnicodeEncodeError:
        encoded = output.encode(sys.stdout.encoding or "utf-8", errors="replace")
        sys.stdout.buffer.write(encoded + b"\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch publication metadata by PMID"
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV path with PMID column")
    parser.add_argument("-o", "--output", required=True, help="Output CSV path")
    parser.add_argument("--sleep", type=float, default=5.0, help="Sleep between requests")
    args = parser.parse_args()

    pmids = read_pmids(args.input)
    records: List[Dict[str, str]] = []
    with requests.Session() as session:
        for pmid in pmids:
            pubmed = fetch_pubmed(session, pmid, args.sleep)
            semsch = fetch_semantic_scholar(session, pmid, args.sleep)
            openalex = fetch_openalex(session, pmid, args.sleep)
            doi = pubmed.get("PubMed.DOI") or semsch.get("scholar.DOI") or ""
            crossref = fetch_crossref(session, doi, args.sleep)
            combined: Dict[str, str] = {}
            combined.update(pubmed)
            combined.update(semsch)
            combined.update(openalex)
            combined.update(crossref)
            # print_results expects a list of records, so wrap the single record
            print_results([combined])
            records.append(combined)

    all_keys = set()
    for rec in records:
        all_keys.update(rec.keys())
    fieldnames = sorted(all_keys)

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    main()
import io
import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent))

from mylib import transforms


def _classify(row_dict):
    row = pd.Series(row_dict)
    buffer = io.StringIO()
    return transforms.classify_row(row, buffer)


def test_meta_analysis_review():
    row = {
        "record_id": "1",
        "pubmed_publication_types": "Meta-Analysis",
        "crossref_type": "",
        "openalex_type": "",
        "scholar_type": "",
        "pubmed_mesh_descriptors": "review literature as topic",
        "pubmed_mesh_qualifiers": "",
        "openalex_mesh_descriptors": "",
        "openalex_mesh_qualifiers": "",
        "title": "Example",
    }
    res = _classify(row)
    assert res.label == "review"


def test_rct_non_review():
    row = {
        "record_id": "2",
        "pubmed_publication_types": "Randomized Controlled Trial",
        "crossref_type": "",
        "openalex_type": "",
        "scholar_type": "",
        "pubmed_mesh_descriptors": "animals",
        "pubmed_mesh_qualifiers": "",
        "openalex_mesh_descriptors": "",
        "openalex_mesh_qualifiers": "",
        "title": "Trial",
    }
    res = _classify(row)
    assert res.label == "non_review"


def test_no_signals_uncertain():
    row = {
        "record_id": "3",
        "pubmed_publication_types": "",
        "crossref_type": "",
        "openalex_type": "",
        "scholar_type": "",
        "pubmed_mesh_descriptors": "",
        "pubmed_mesh_qualifiers": "",
        "openalex_mesh_descriptors": "",
        "openalex_mesh_qualifiers": "",
        "title": "No signals",
    }
    res = _classify(row)
    assert res.label == "uncertain"

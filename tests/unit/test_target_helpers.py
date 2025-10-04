from __future__ import annotations

import xml.etree.ElementTree as ET

import pandas as pd
import pytest

from library.postprocessing.target.helpers import first_element_text, normalize_text


@pytest.mark.unit
def test_normalize_text__trims_and_lowercases() -> None:
    assert normalize_text("  HeLLo ") == "hello"
    assert normalize_text(None) == ""
    assert normalize_text(pd.NA) == ""


@pytest.mark.unit
def test_first_element_text__returns_first_match() -> None:
    root = ET.fromstring(
        """
        <Taxon>
            <Lineage>Viruses; Negarnaviricota</Lineage>
            <ScientificName>Influenza A virus</ScientificName>
        </Taxon>
        """
    )
    nodes = list(root)
    assert first_element_text(nodes, "ScientificName") == "Influenza A virus"
    assert first_element_text(nodes, "NonExisting") is None

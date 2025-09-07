"""Rule definitions for deterministic document classifier.

This module centralizes regex patterns and keyword sets so that the
classifier logic can remain readable. All terms are stored in lowercase
for case-insensitive comparisons.
"""

from __future__ import annotations

import re
from typing import Final, Pattern, Set

# Regular expressions used by text scoring
RE_IN_VITRO: Final[Pattern[str]] = re.compile(
    r"\b(in\s*vitro|cell[-\s]?based|cell[-\s]?free|enzyme assay|biochemical assay|"
    r"binding assay|spr|itc|fluorescence polarization|reporter assay|luciferase|"
    r"western blot|flow cytometry|ic50|ec50|ki|kd|patch clamp|whole[-\s]?cell)\b",
    re.I,
)
RE_CELLLINES: Final[Pattern[str]] = re.compile(
    r"\b(hek\s*293|hepg2|hela|cho|u2os|mcf[-\s]?7|thp[-\s]?1|293t|nhdf|huh[-\s]?7|"
    r"sf9|primary\s+(?:cells|neurons)|organoid[s]?)\b",
    re.I,
)
RE_ROUTE: Final[Pattern[str]] = re.compile(
    r"\b(ip|iv|po|sc|im|intrathecal|gavage)\b", re.I
)
RE_DOSE: Final[Pattern[str]] = re.compile(
    r"\b\d+(?:\.\d+)?\s*(?:mg|u?g)\s*/\s*kg\b", re.I
)
RE_BEHAV: Final[Pattern[str]] = re.compile(
    r"\b(von\s+frey|hot\s*plate|tail\s*flick|formalin|rotarod|paw\s*edema|"
    r"allodynia|hyperalgesia|carrageenan|cfa)\b",
    re.I,
)
RE_PROTOCOL: Final[Pattern[str]] = re.compile(
    r"\b(protocol|study protocol|statistical analysis plan|sap|trial registration)\b",
    re.I,
)
RE_RESULTS: Final[Pattern[str]] = re.compile(
    r"\b(result|observed|measured|significant|p\s*[<≤]\s*0?\.\d+|dose[-\s]?response|"
    r"endpoint[s]?)\b",
    re.I,
)
RE_CLINICAL: Final[Pattern[str]] = re.compile(
    r"\b(randomized|placebo|double[-\s]?blind|single[-\s]?blind|phase\s*(?:i|ii|iii|iv)\b|"
    r"multicenter|participants|subjects|cohort|arm[s]?)\b",
    re.I,
)

# Publication type mappings
PT_CLINICAL: Final[Set[str]] = {
    "randomized controlled trial",
    "controlled clinical trial",
    "clinical trial",
    "clinical trial phase i",
    "clinical trial phase ii",
    "clinical trial phase iii",
    "clinical trial phase iv",
    "pragmatic clinical trial",
    "multicenter study",
    "double-blind method",
    "single-blind method",
}
PT_REVIEW: Final[Set[str]] = {
    "review",
    "systematic review",
    "scoping review",
    "umbrella review",
    "literature review",
}
PT_OTHER: Final[Set[str]] = {
    "meta-analysis",
    "guideline",
    "practice guideline",
    "consensus",
    "editorial",
    "comment",
    "letter",
    "news",
    "erratum",
    "study protocol",
    "protocol",
    "statistical analysis plan",
    "case report",
    "case series",
}

# MeSH descriptor/qualifier keyword sets
MESH_A: Final[Set[str]] = {
    "animals",
    "mouse",
    "mice",
    "rat",
    "rabbit",
    "dog",
    "swine",
    "zebrafish",
    "hamster",
    "primate",
    "disease models",
    "animal",
}
MESH_B: Final[Set[str]] = {
    "pharmacology",
    "drug effects",
    "administration & dosage",
    "therapeutic use",
    "toxicity",
    "adverse effects",
    "pharmacokinetics",
    "metabolism",
}
MESH_C: Final[Set[str]] = {
    "pain measurement",
    "nociceptors",
    "hyperalgesia",
    "allodynia",
    "edema",
    "inflammation",
}
MESH_D: Final[Set[str]] = {
    "cells, cultured",
    "cell line",
    "recombinant proteins",
    "enzymes",
    "hek293",
    "hepg2",
    "hela",
    "cho",
    "u2os",
    "mcf-7",
    "thp-1",
    "293t",
    "nhdf",
    "huh7",
    "sf9",
    "organoid",
    "primary cell",
}
MESH_E: Final[Set[str]] = {
    "in vitro techniques",
    "pharmacology",
    "drug effects",
    "enzyme activation",
    "enzyme inhibition",
    "binding sites",
    "affinity",
}
MESH_F: Final[Set[str]] = {
    "kinetics",
    "dose-response relationship, drug",
    "signal transduction",
}

# Utility regex for splitting list fields
LIST_SPLIT_RE: Final[Pattern[str]] = re.compile(r"[;,|/]+")

"""Tests for helpers defined in ``library.processing.target``."""

from library.processing.target import _collect_reaction_ec_numbers


def test_collect_reaction_ec_numbers_excludes_reactome_like_xrefs() -> None:
    """Ensure excluded xref sources are ignored while valid ones remain."""

    components = [
        {
            "target_component_synonyms": {
                "target_component_synonym": [
                    {
                        "syn_type": "REACTION",
                        "component_synonym": "1.1.1.1",
                    },
                    {
                        "syn_type": "REACTION_NUMBER",
                        "component_synonym": "2.2.2.2",
                    },
                ]
            },
            "target_component_xrefs": {
                "target": [
                    {"xref_src_db": "Reactome", "xref_id": "R-HSA-111111"},
                    {"xref_src_db": "RHEA", "xref_id": "RHEA:12345"},
                    {"xref_src_db": "MetaCyc", "xref_id": "META:999"},
                    {"xref_src_db": "EC_REACTION", "xref_id": "EC:9.9.9.9"},
                    {"xref_src_db": "SABIO-RK", "xref_id": "SABIO:555"},
                ]
            },
        }
    ]

    assert _collect_reaction_ec_numbers(components) == "1.1.1.1|2.2.2.2"


def test_collect_reaction_ec_numbers_filters_non_ec_tokens() -> None:
    """Ensure only strict EC numbers survive from mixed inputs."""

    components = [
        {
            "target_component_synonyms": {
                "target_component_synonym": [
                    {
                        "syn_type": "reaction",
                        "component_synonym": "EC 3.1.1.1",
                    },
                    {
                        "syn_type": "reaction_number",
                        "component_synonym": "Not an EC",
                    },
                ]
            },
            "target_component_xrefs": {
                "target": [
                    {"xref_src_db": "OTHER", "xref_id": "1.1.-.-"},
                    {"xref_src_db": "OTHER", "xref_id": "SABIO:555"},
                    {"xref_src_db": "OTHER", "xref_id": "4.2.1.1"},
                ]
            },
        }
    ]

    assert _collect_reaction_ec_numbers(components) == "1.1.-.-|3.1.1.1|4.2.1.1"

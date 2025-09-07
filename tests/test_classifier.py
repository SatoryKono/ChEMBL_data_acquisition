from classifier import DocumentClassifier
from io_schemas import DocumentRecord


def make_record(**kwargs):
    base = {
        "title": "",
        "abstract": "",
        "doi": "",
        "pubmed_publicationtype": [],
        "scholar_publicationtypes": [],
        "openalex_publicationtypes": [],
        "crossref_type": None,
        "openalex_typecrossref": None,
        "pubmed_mesh_descriptors": [],
        "pubmed_mesh_qualifiers": [],
        "openalex_meshdescriptors": [],
        "openalex_meshqualifiers": [],
        "pubmed_chemicallist": [],
        "optional_experiment_kind": None,
    }
    base.update(kwargs)
    return DocumentRecord(**base)


classifier = DocumentClassifier()


def test_protocol_veto_overrides_clinical():
    record = make_record(
        title="Randomized Controlled Trial",
        abstract="This study protocol evaluates...",
        pubmed_publicationtype=["Randomized Controlled Trial"],
        pubmed_mesh_descriptors=["humans"],
    )
    result = classifier.classify(record)
    assert result.final_class == "Other non experimental publication"


def test_review_crossref_conflict():
    record = make_record(
        title="A review of therapy",
        pubmed_publicationtype=["Review"],
        crossref_type="journal-article",
    )
    result = classifier.classify(record)
    assert result.final_class == "Review"
    assert "xrf_generic_vs_nonexperimental" in result.conflicts


def test_veterinary_trial_in_vivo():
    record = make_record(
        title="Dose response in mice",
        abstract="IP 10 mg/kg reduced pain in rats",
        pubmed_publicationtype=["Clinical Trial, Veterinary"],
        pubmed_mesh_descriptors=["animals", "mice", "pain measurement"],
        pubmed_mesh_qualifiers=["pharmacology"],
        optional_experiment_kind="F",
    )
    result = classifier.classify(record)
    assert result.final_class == "Experimental in vivo bioactivity"


def test_clean_in_vitro():
    record = make_record(
        title="Binding assay in cell line",
        abstract="In vitro binding assay using HEK293 cells showed IC50",
        pubmed_publicationtype=["Comparative Study"],
        pubmed_mesh_descriptors=["cell line"],
        pubmed_mesh_qualifiers=["in vitro techniques"],
        pubmed_chemicallist=["aspirin"],
        optional_experiment_kind="B",
    )
    result = classifier.classify(record)
    assert result.final_class == "Experimental in vitro bioactivity"
    assert result.pt_comparative


def test_review_with_mouse_not_in_vivo():
    record = make_record(
        title="Review of mouse cell lines",
        abstract="We discuss mouse models and HEK293 cells in vitro",
        pubmed_publicationtype=["Review"],
    )
    result = classifier.classify(record)
    assert result.final_class == "Review"


def test_tie_with_behavior_goes_vivo():
    record = make_record(
        title="Effects in mice and cells",
        abstract="IP 5 mg/kg reduced paw edema and in vitro assay on HEK293 cells",
        pubmed_mesh_descriptors=["animals", "pain measurement", "cell line"],
        pubmed_mesh_qualifiers=["pharmacology", "in vitro techniques"],
    )
    result = classifier.classify(record)
    assert result.final_class == "Experimental in vivo bioactivity"
    assert "vivo_vs_vitro_tie" in result.conflicts


def test_tie_without_behavior_goes_vitro():
    record = make_record(
        title="Animal cell line kinetics",
        abstract="Study on mice cell line kinetics in vitro",
        pubmed_mesh_descriptors=["animals", "pain measurement", "cell line"],
        pubmed_mesh_qualifiers=["pharmacology", "in vitro techniques"],
    )
    result = classifier.classify(record)
    assert result.final_class == "Experimental in vitro bioactivity"
    assert "vivo_vs_vitro_tie" in result.conflicts

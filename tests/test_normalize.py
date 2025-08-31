from review_classifier.normalize import canonicalize_mesh, is_review, split_and_normalize


def test_split_and_normalize():
    assert split_and_normalize("Review|Journal Article", "|;,/") == ["review", "journal article"]


def test_is_review():
    assert is_review({"review"})
    assert is_review({"systematic review"})
    assert not is_review({"journal article", "letter"})


def test_canonicalize_mesh():
    assert canonicalize_mesh("Rats.") == "rat"
    assert canonicalize_mesh("cells (drug effects)") == "cell"
    assert canonicalize_mesh("DNA-binding proteins") == "dna-binding protein"

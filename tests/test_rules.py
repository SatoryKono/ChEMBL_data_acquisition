import pandas as pd

from review_classifier.mesh_probs import load_mesh_probabilities
from review_classifier.pipeline import Config, process_dataframe


def make_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "PubMed.PublicationType": ["Review", "Article", "Review", "Article"],
            "OpenAlex.PublicationTypes": ["Review", "Article", "Article", "Article"],
            "scholar.PublicationTypes": ["Article", "Article", "Article", "Review"],
            "PubMed.MeSH_Descriptors": ["", "", "rats|mice|cells", "humans"],
            "PubMed.MeSH_Qualifiers": ["", "", "", ""],
            "OpenAlex.MeshDescriptors": ["", "", "", ""],
            "OpenAlex.MeshQualifiers": ["", "", "", ""],
        }
    )


def test_vote_and_mesh(tmp_path):
    mesh_map = load_mesh_probabilities("tests/data/experimental_mesh.csv")
    df = make_df()
    cfg = Config()
    result = process_dataframe(df, mesh_map, cfg)
    labels = result["label"].tolist()
    assert labels[0] == "review"  # votes: 2
    assert labels[1] == "non-review"  # votes: 0
    assert labels[2] == "non-review"  # mesh strong experimental
    assert labels[3] == "review"  # mesh strong review

"""Property-based tests for document publication type parsing."""

from __future__ import annotations

import pytest

from library.pipelines.document import type_terms

hypothesis = pytest.importorskip("hypothesis")
st = hypothesis.strategies
SearchStrategy = st.SearchStrategy
given = hypothesis.given
settings = hypothesis.settings

_TOKEN_ALPHABET = list("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ -_")


def _token_strategy() -> SearchStrategy[str]:
    random_token = st.text(alphabet=_TOKEN_ALPHABET, min_size=0, max_size=12)
    synonym_token = st.sampled_from(tuple(type_terms._SYNONYMS.keys()))
    return st.one_of(random_token, synonym_token)


@settings(max_examples=50, derandomize=True, deadline=None)
@given(
    tokens=st.lists(_token_strategy(), min_size=0, max_size=6),
    separator=st.sampled_from(["|", ";", ",", "/"]),
)
def test_parse_terms__normalises_and_deduplicates(
    tokens: list[str], separator: str
) -> None:
    """Parsed tokens are normalised, unique and sorted."""

    value = separator.join(tokens)

    expected_unique: list[str] = []
    seen: set[str] = set()
    for token in tokens:
        normalised = type_terms._normalise_token(token)
        if normalised and normalised not in seen:
            seen.add(normalised)
            expected_unique.append(normalised)

    expected = sorted(expected_unique)
    assert type_terms.parse_terms(value) == expected


@settings(max_examples=50, derandomize=True, deadline=None)
@given(
    st.one_of(
        st.none(),
        st.integers(),
        st.floats(allow_nan=False),
        st.lists(st.integers()),
        st.dictionaries(st.text(min_size=1, max_size=3), st.integers()),
    )
)
def test_parse_terms__non_string_inputs_return_empty(value: object) -> None:
    """Non-string values yield an empty token list."""

    assert type_terms.parse_terms(value) == []

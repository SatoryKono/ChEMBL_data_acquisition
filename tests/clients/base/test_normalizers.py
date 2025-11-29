from bioetl.clients.base.normalizers import IdentityNormalizerImpl


def test_identity_normalizer_returns_copy():
    normalizer = IdentityNormalizerImpl()
    record = {"id": 1, "name": "foo"}

    normalized = normalizer.normalize(record)

    assert normalized == record
    assert normalized is not record


def test_identity_normalizer_batch_returns_copies():
    normalizer = IdentityNormalizerImpl()
    records = [{"id": 1}, {"id": 2}]

    normalized_batch = list(normalizer.normalize_batch(records))

    assert normalized_batch == records
    for original, normalized in zip(records, normalized_batch):
        assert normalized is not original

from __future__ import annotations

import textwrap

import pytest

from library.resources import dictionaries


@pytest.mark.unit
def test_parse_manifest__ignores_manifest_in_checksum(tmp_path):
    dictionary_dir = tmp_path / "dictionary"
    dictionary_dir.mkdir()

    data_file = dictionary_dir / "data.csv"
    data_file.write_text("id\n1\n", encoding="utf-8")

    manifest_path = dictionary_dir / "manifest.yaml"

    data_sha = dictionaries._compute_sha256(data_file)
    root_sha = dictionaries._compute_sha256(
        dictionary_dir, exclude=(manifest_path,)
    )

    manifest_template = textwrap.dedent(
        """
        version: 1
        resources:
          dictionary_root:
            path: .
            version: "1.0"
            sha256: "{root_sha}"
            generator: tools/build_dictionary_resources.py
          sample_table:
            path: data.csv
            version: "1.0"
            sha256: "{data_sha}"
            generator: tools/build_dictionary_resources.py
        """
    )
    manifest_path.write_text(
        manifest_template.format(root_sha=root_sha, data_sha=data_sha),
        encoding="utf-8",
    )

    resources = dictionaries._parse_manifest(base_dir=dictionary_dir)

    assert resources["dictionary_root"].sha256 == root_sha
    assert resources["sample_table"].sha256 == data_sha

    mutated_manifest = manifest_template.format(
        root_sha="deadbeef" * 8, data_sha=data_sha
    )
    manifest_path.write_text(mutated_manifest, encoding="utf-8")

    with pytest.raises(dictionaries.DictionaryManifestError):
        dictionaries._parse_manifest(base_dir=dictionary_dir)

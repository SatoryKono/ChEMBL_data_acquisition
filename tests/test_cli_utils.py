from __future__ import annotations

from library.cli_utils import build_parser


def test_cli_utils_flags_and_help() -> None:
    parser = build_parser()
    actions = parser._option_string_actions
    expected = {
        "-h",
        "--help",
        "--log-level",
        "--input",
        "--output",
        "--sep",
        "--encoding",
        "--col-order",
        "--key-cols",
        "--chunk-size",
        "--merge-chunk-size",
    }
    assert set(actions) == expected
    assert actions["--log-level"].help == "Logging level"
    assert actions["--input"].help == "Input CSV file"
    assert (
        actions["--output"].help
        == "Destination CSV file (default: output_<stem>_<YYYYMMDD>.csv)"
    )
    assert actions["--sep"].help == "CSV delimiter"
    assert actions["--encoding"].help == "File encoding"
    assert actions["--col-order"].help == "Preferred column order"
    assert actions["--key-cols"].help == "Columns used for sorting"
    assert actions["--chunk-size"].help == "Number of rows read per chunk"
    assert (
        actions["--merge-chunk-size"].help
        == "Rows loaded per temporary file during merge"
    )
    assert parser.description is not None
    assert parser.description.startswith(
        "CLI wrapper for :func:`write_csv_deterministic`"
    )

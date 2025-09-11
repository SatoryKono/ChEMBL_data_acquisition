import pytest

import get_activity_data as gad


def test_negative_limit_rejected(capsys: pytest.CaptureFixture[str]) -> None:
    """Ensure ``--limit`` rejects non-positive integers."""
    with pytest.raises(SystemExit) as excinfo:
        gad.main(["--limit", "-1"])
    assert excinfo.value.code == 2
    err = capsys.readouterr().err
    assert "--limit must be a positive integer" in err

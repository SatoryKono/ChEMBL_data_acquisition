from __future__ import annotations

import pytest
from pydantic import ValidationError

from library.parser_schema import CSVExportArgs


def test_extra_fields_forbidden() -> None:
    """Extra fields should raise a validation error."""
    with pytest.raises(ValidationError):
        CSVExportArgs.model_validate({"input_csv": "in.csv", "extra": 1})

# Activity output column pruning

## Summary
- Updated `scripts/get_activity_data.py` to remove deprecated activity citation columns from the final export using a whitelist approach with safe dropping.
- Added a unit smoke test ensuring the restricted columns are excluded and the removal is logged.

## Verification
- `pytest tests/test_activity_output_columns.py`

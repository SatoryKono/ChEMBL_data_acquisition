# Assay output column pruning

## Summary
- Added a final metadata hook in `scripts/get_assay_data.py` to drop the deprecated assay columns before export while preserving schema order.
- Logged the removed columns once per pipeline run to aid auditing.

## Testing
- `pytest tests/test_assay_output_schema.py` *(fails: file or directory not found in repository)*
